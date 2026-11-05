"""Append-only launch reservations, separate from mathematical checkpoints.

No solver, reader retry, orphan adoption or inferred exit status. Unknown
attempts retain their full reservation. Only a complete wait4 witness may
settle one explicitly born process. External plan and expected HEAD are trust
anchors; authenticating the supervisor/witness belongs to the caller.
"""
import copy
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import stat

MAX_RECORDS = 98304 + 1024
MAX_BYTES = 128 * 2**10


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def digest(value):
    return hashlib.sha256(value).hexdigest()


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def number(value, minimum=0):
    require(type(value) is int and value >= minimum, 'integer resource domain')
    return value


def sha(value):
    require(type(value) is str and re.fullmatch('[0-9a-f]{64}', value), 'SHA256 domain')
    return value


def keys(value, expected):
    require(type(value) is dict and set(value) == set(expected.split()), 'resource member inventory')


def identity(value):
    keys(value, 'pid ppid start_ticks boot_id argv')
    number(value['pid'], 1); number(value['ppid'], 1); number(value['start_ticks'], 1)
    require(type(value['boot_id']) is str and re.fullmatch('[0-9a-f-]{36}', value['boot_id']), 'process boot')
    require(type(value['argv']) is list and value['argv'] and all(type(x) is str and x and '\x00' not in x for x in value['argv']), 'complete argv')
    return value


def regular(path):
    info = path.lstat()
    require(stat.S_ISREG(info.st_mode) and info.st_nlink == 1 and info.st_size <= MAX_BYTES, 'regular resource file required')
    return path.read_bytes()


def directory(path):
    require(stat.S_ISDIR(path.lstat().st_mode), 'regular resource directory required')


def sync(path):
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def write_new(path, data):
    with path.open('xb') as stream:
        stream.write(data); stream.flush(); os.fsync(stream.fileno())


def decode(data):
    value = json.loads(data)
    require(canonical(value) == data, 'canonical resource record required')
    return value


def validate_contract(contract):
    keys(contract, 'schema plan_sha256 source_sha256 CPU_cap_ns wall_cap_ns prior_charges measurement')
    require(contract['schema'] == 1 and type(contract['schema']) is int, 'resource schema')
    sha(contract['plan_sha256']); sha(contract['source_sha256'])
    number(contract['CPU_cap_ns'], 1); number(contract['wall_cap_ns'], 1)
    require(contract['measurement'] == 'wait4_rusage_ceil_nanoseconds_v1', 'resource measurement policy')
    require(type(contract['prior_charges']) is list, 'explicit prior charges required')
    ids = set()
    for entry in contract['prior_charges']:
        keys(entry, 'id CPU_ns wall_ns evidence_sha256 kind')
        require(type(entry['id']) is str and entry['id'] and entry['id'] not in ids, 'prior charge identity')
        ids.add(entry['id']); number(entry['CPU_ns']); number(entry['wall_ns']); sha(entry['evidence_sha256'])
        require(entry['kind'] in ('complete_measured', 'unknown_full_reservation'), 'prior cost cannot be inferred from progress')
    for kind, cap in [('CPU_ns', 'CPU_cap_ns'), ('wall_ns', 'wall_cap_ns')]:
        require(sum(e[kind] for e in contract['prior_charges']) <= contract[cap], 'prior costs already exceed frozen cap')


class Ledger:
    """Every mutation CASes current HEAD; uncommitted files fail closed.

    reserve -> authorize -> born -> settle. Authorize must be committed BEFORE
    Popen/fork. A launch with missing birth or interrupted settlement remains
    fully charged. Failed/over-budget complete evidence permanently halts new
    dispatch. Mathematical state/branch selection is never an input.
    """
    def __init__(self, root, contract):
        self.root = Path(root); directory(self.root)
        validate_contract(contract)
        self.contract = copy.deepcopy(contract)
        self.binding = digest(canonical(contract))

    def head(self):
        return digest(regular(self.root / 'HEAD')) if (self.root / 'HEAD').exists() or (self.root / 'HEAD').is_symlink() else None

    def _chain(self, expected):
        require(expected is None or sha(expected), 'expected resource HEAD')
        require(self.head() == expected, 'stale resource HEAD')
        data = regular(self.root / 'HEAD') if expected else None
        chain = []
        while data:
            q = decode(data)
            keys(q, 'schema contract_sha256 sequence parent_sha256 event')
            require(q['schema'] == 1 and q['contract_sha256'] == self.binding, 'resource contract binding')
            seq = number(q['sequence'], 1)
            require(seq <= MAX_RECORDS, 'resource history budget')
            require(regular(self.root / f'{seq:08d}.json') == data, 'resource archive/HEAD mismatch')
            chain.append(q)
            if seq == 1:
                require(q['parent_sha256'] is None, 'initial resource ancestry')
                break
            data = regular(self.root / f'{seq-1:08d}.json')
            require(digest(data) == sha(q['parent_sha256']), 'resource ancestry SHA')
            require(decode(data)['sequence'] == seq - 1, 'resource ancestry sequence')
        allowed = {'SESSION.lock'} | ({'HEAD'} if chain else set()) | {f'{q["sequence"]:08d}.json' for q in chain}
        require({p.name for p in self.root.iterdir()} <= allowed, 'uncommitted/orphan resource record; seal and stop,never adopt/refund')
        return list(reversed(chain))

    def _state(self, chain):
        state = dict(attempts={}, CPU_ns=sum(e['CPU_ns'] for e in self.contract['prior_charges']),
                     wall_ns=sum(e['wall_ns'] for e in self.contract['prior_charges']), halted=False)
        for q in chain:
            e = q['event']; require(type(e) is dict and 'kind' in e, 'resource event')
            if e['kind'] == 'reserve':
                keys(e, 'kind id CPU_ns wall_ns plan_sha256 command supervisor')
                require(not state['halted'], 'dispatch after terminal resource/execution failure')
                require(type(e['id']) is str and re.fullmatch('[A-Za-z0-9_-]{1,80}', e['id']) and e['id'] not in state['attempts'], 'unique attempt ID')
                number(e['CPU_ns'], 1); number(e['wall_ns'], 1)
                require(e['plan_sha256'] == self.contract['plan_sha256'], 'reservation plan binding')
                identity(e['supervisor']); require(type(e['command']) is list and e['command'] and all(type(x) is str and x and '\x00' not in x for x in e['command']), 'reservation argv')
                state['CPU_ns'] += e['CPU_ns']; state['wall_ns'] += e['wall_ns']
                require(state['CPU_ns'] <= self.contract['CPU_cap_ns'] and state['wall_ns'] <= self.contract['wall_cap_ns'], 'cumulative resource admission')
                state['attempts'][e['id']] = dict(reservation=e, phase='reserved')
                continue
            require(e.get('id') in state['attempts'], 'resource event without reservation')
            attempt = state['attempts'][e['id']]
            if e['kind'] == 'authorize':
                keys(e, 'kind id'); require(attempt['phase'] == 'reserved', 'launch authorization is single use')
                attempt.update(phase='authorized', intent_sha256=digest(canonical(q)))
            elif e['kind'] == 'born':
                keys(e, 'kind id process')
                require(attempt['phase'] == 'authorized', 'birth without committed intent')
                identity(e['process']); r = attempt['reservation']
                require(e['process']['argv'] == r['command'] and e['process']['ppid'] == r['supervisor']['pid'] and e['process']['boot_id'] == r['supervisor']['boot_id'], 'child/supervisor/argv binding')
                attempt.update(phase='born', process=e['process'], birth_sha256=digest(canonical(q)))
            elif e['kind'] == 'settle':
                keys(e, 'kind id evidence evidence_sha256')
                require(attempt['phase'] == 'born', 'settlement requires explicit process birth')
                proof = e['evidence']; require(digest(canonical(proof)) == sha(e['evidence_sha256']), 'complete resource witness SHA')
                keys(proof, 'schema complete attempt_id plan_sha256 intent_sha256 birth_sha256 process supervisor measurement CPU_ns wall_ns start_monotonic_ns end_monotonic_ns return_code')
                require(proof['schema'] == 1 and proof['complete'] is True and type(proof['return_code']) is int, 'incomplete/unknown exit resource witness')
                r = attempt['reservation']
                require(proof['attempt_id'] == e['id'] and proof['plan_sha256'] == r['plan_sha256'] and proof['intent_sha256'] == attempt['intent_sha256'] and proof['birth_sha256'] == attempt['birth_sha256'], 'resource attempt/plan/commit binding')
                require(proof['process'] == attempt['process'] and proof['supervisor'] == r['supervisor'], 'resource identity reuse rejected')
                require(proof['measurement'] == self.contract['measurement'], 'resource measurement witness')
                number(proof['CPU_ns']); number(proof['wall_ns'], 1); number(proof['start_monotonic_ns']); number(proof['end_monotonic_ns'], 1)
                require(proof['wall_ns'] == proof['end_monotonic_ns'] - proof['start_monotonic_ns'], 'resource wall interval')
                state['CPU_ns'] += proof['CPU_ns'] - r['CPU_ns']; state['wall_ns'] += proof['wall_ns'] - r['wall_ns']
                failure = proof['return_code'] != 0 or proof['CPU_ns'] > r['CPU_ns'] or proof['wall_ns'] > r['wall_ns']
                attempt.update(phase='settled', evidence_sha256=e['evidence_sha256'], execution_or_resource_failure=failure)
                state['halted'] |= failure
            else:
                raise ValueError('unsupported resource event')
        return state

    def read(self, expected_head):
        return self._state(self._chain(expected_head))

    def append(self, event, expected_head, hook=lambda stage: None):
        """No child may launch until authorize append has returned successfully."""
        lock = self.root / 'SESSION.lock'
        if lock.exists() or lock.is_symlink():
            regular(lock)
        with lock.open('a+b') as owner:
            fcntl.flock(owner.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            chain = self._chain(expected_head)
            q = dict(schema=1, contract_sha256=self.binding, sequence=len(chain)+1,
                     parent_sha256=expected_head, event=event)
            require(q['sequence'] <= MAX_RECORDS, 'resource history budget')
            state = self._state(chain + [q])  # Semantic rejection before any append.
            data = canonical(q); require(len(data) <= MAX_BYTES, 'resource record byte budget')
            write_new(self.root / f'{q["sequence"]:08d}.json', data); sync(self.root)
            hook('before_HEAD')
            temp = self.root / 'HEAD.pending'; write_new(temp, data)
            os.replace(temp, self.root / 'HEAD'); sync(self.root)
            hook('after_HEAD')
            return dict(HEAD_sha256=digest(data), sequence=q['sequence'], state=state)

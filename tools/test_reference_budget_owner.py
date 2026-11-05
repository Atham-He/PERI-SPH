"""New cumulative-ledger and real-class memory-owner controls,zero physics.

Synthetic receipts test the ledger's semantics; one real wait4 process tests
resource-witness integration. Real Inspector classes use synthetic full-shaped
memory,NOT a constitutive trajectory or a qualified review restart.
"""
import copy
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time
import unittest
from fractions import Fraction
import decimal as dec

from peri_reference_budget import Ledger, canonical, digest, write_new


def process_identity(pid=None):
    pid = os.getpid() if pid is None else pid
    p = Path('/proc') / str(pid)
    fields = (p/'stat').read_text().rsplit(')', 1)[1].split()
    return dict(pid=pid, ppid=int(fields[1]), start_ticks=int(fields[19]),
                boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
                argv=(p/'cmdline').read_bytes().decode().rstrip('\x00').split('\x00'))


def ns(seconds):
    q = Fraction(seconds) * 10**9
    return -(-q.numerator // q.denominator)


class BudgetControls(unittest.TestCase):
    def setUp(self):
        self.root = Path(os.environ['REFERENCE_BUDGET_WORK']) / self._testMethodName
        self.root.mkdir(); (self.root/'ledger').mkdir()
        self.contract = dict(schema=1, plan_sha256=digest(b'prospectively_frozen_SYNTHETIC_plan'),
            source_sha256=digest(b'prospectively_frozen_SYNTHETIC_source'), CPU_cap_ns=100*10**9,
            wall_cap_ns=200*10**9, prior_charges=[dict(id='prior_unknown', CPU_ns=7*10**9, wall_ns=9*10**9,
            evidence_sha256=digest(b'SYNTHETIC_unknown_retained_full_reservation'), kind='unknown_full_reservation')],
            measurement='wait4_rusage_ceil_nanoseconds_v1')
        self.ledger = Ledger(self.root/'ledger', self.contract)
        self.parent = process_identity(); self.expected = None

    def event(self, kind, **fields):
        value = dict(kind=kind, id='attempt_1', **fields)
        result = self.ledger.append(value, self.expected); self.expected = result['HEAD_sha256']
        return result

    def reserve(self, cpu=10*10**9, wall=20*10**9, command=None):
        return self.event('reserve', CPU_ns=cpu, wall_ns=wall, plan_sha256=self.contract['plan_sha256'],
                          command=command or ['/synthetic/reader', '--fixed95'], supervisor=self.parent)

    def born(self):
        self.reserve(); intent = self.event('authorize')
        process = dict(pid=99999, ppid=self.parent['pid'], start_ticks=88888,
                       boot_id=self.parent['boot_id'], argv=['/synthetic/reader', '--fixed95'])
        birth = self.event('born', process=process)
        return dict(schema=1, complete=True, attempt_id='attempt_1', plan_sha256=self.contract['plan_sha256'],
            intent_sha256=intent['HEAD_sha256'], birth_sha256=birth['HEAD_sha256'], process=process,
            supervisor=self.parent, measurement=self.contract['measurement'], CPU_ns=2*10**9,
            wall_ns=3*10**9, start_monotonic_ns=100, end_monotonic_ns=3*10**9+100, return_code=0)

    def settle(self, proof):
        return self.event('settle', evidence=proof, evidence_sha256=digest(canonical(proof)))

    def files(self):
        return {p.name: p.read_bytes() for p in (self.root/'ledger').iterdir() if p.is_file()}

    def test_exact_unknown_charge_is_not_last_progress(self):
        self.reserve()
        state = self.ledger.read(self.expected)
        self.assertEqual(state['CPU_ns'], 17*10**9); self.assertEqual(state['wall_ns'], 29*10**9)
        self.assertEqual(state['attempts']['attempt_1']['phase'], 'reserved')

    def test_reserve_before_authorize_and_birth(self):
        with self.assertRaisesRegex(ValueError, 'without reservation'):
            self.event('authorize')
        self.reserve(); self.event('authorize')
        before = self.files()
        with self.assertRaisesRegex(ValueError, 'single use'):
            self.event('authorize')
        self.assertEqual(before, self.files())

    def test_complete_settlement_only_releases_measured_unused_capacity(self):
        proof = self.born(); state = self.settle(proof)['state']
        self.assertEqual(state['CPU_ns'], 9*10**9); self.assertEqual(state['wall_ns'], 12*10**9)
        self.assertFalse(state['halted'])

    def test_unknown_exit_and_partial_progress_cannot_settle(self):
        proof = self.born(); before = self.files()
        for field, value in [('complete', False), ('return_code', None)]:
            q = copy.deepcopy(proof); q[field] = value
            with self.assertRaisesRegex(ValueError, 'incomplete/unknown'):
                self.settle(q)
            self.assertEqual(before, self.files())
        self.assertEqual(self.ledger.read(self.expected)['CPU_ns'], 17*10**9)

    def test_correct_SHA_semantic_corruptions_leave_ledger_unchanged(self):
        proof = self.born(); before = self.files()
        mutations = [('wrong_attempt', lambda p: p.update(attempt_id='other')),
            ('wrong_plan', lambda p: p.update(plan_sha256=digest(b'wrong'))),
            ('wrong_intent', lambda p: p.update(intent_sha256=digest(b'wrong'))),
            ('wrong_birth', lambda p: p.update(birth_sha256=digest(b'wrong'))),
            ('reused_PID_ticks', lambda p: p['process'].update(start_ticks=1)),
            ('wrong_boot', lambda p: p['process'].update(boot_id='0'*36)),
            ('wrong_full_argv', lambda p: p['process'].update(argv=['/other'])),
            ('wrong_supervisor', lambda p: p['supervisor'].update(start_ticks=1)),
            ('negative_CPU', lambda p: p.update(CPU_ns=-1)),
            ('wrong_wall_interval', lambda p: p.update(wall_ns=1)),
            ('unknown_measurement', lambda p: p.update(measurement='progress_estimate')),
            ('incomplete_inventory', lambda p: p.pop('complete'))]
        evidence = []
        for label, mutate in mutations:
            p = copy.deepcopy(proof); mutate(p); raw = canonical(p)
            write_new(self.root/(label+'.correctSHA.json'), raw)
            with self.assertRaises(ValueError) as caught:
                self.settle(p)
            self.assertEqual(before, self.files())
            evidence.append(dict(name=label, sha256=digest(raw), rejected=True, reason=str(caught.exception)))
        write_new(self.root/'semantic_controls.json', canonical(evidence))

    def test_actual_cost_above_reservation_is_recorded_and_halts(self):
        proof = self.born(); proof['CPU_ns'] = 11*10**9
        state = self.settle(proof)['state']
        self.assertEqual(state['CPU_ns'], 18*10**9); self.assertTrue(state['halted'])
        with self.assertRaisesRegex(ValueError, 'terminal'):
            self.event('reserve', CPU_ns=1, wall_ns=1, plan_sha256=self.contract['plan_sha256'], command=['x'], supervisor=self.parent)

    def test_nonzero_return_keeps_real_charge_and_halts(self):
        proof = self.born(); proof['return_code'] = -9
        self.assertTrue(self.settle(proof)['state']['halted'])

    def test_cumulative_cap_rejects_before_append(self):
        before = self.files()
        with self.assertRaisesRegex(ValueError, 'admission'):
            self.reserve(cpu=94*10**9)
        # The empty lock is the only possible setup side effect,no ledger event.
        self.assertEqual({k:v for k,v in self.files().items() if k!='SESSION.lock'}, before)
        self.assertIsNone(self.ledger.head())

    def test_old_mathematical_selection_cannot_rewind_resource_HEAD(self):
        proof = self.born(); self.settle(proof); before = self.files()
        math = self.root/'mathematical_selection_only.json'
        write_new(math, canonical(dict(selected_step=8)))
        math.write_bytes(canonical(dict(selected_step=3)))
        self.assertEqual(self.files(), before)
        self.assertEqual(self.ledger.read(self.expected)['CPU_ns'], 9*10**9)

    def test_stale_HEAD_and_replayed_reservation_reject(self):
        self.reserve(); first = self.expected; self.event('authorize'); before = self.files()
        with self.assertRaisesRegex(ValueError, 'stale'):
            self.ledger.read(first)
        with self.assertRaisesRegex(ValueError, 'unique'):
            self.reserve()
        self.assertEqual(before, self.files())

    def test_rolled_back_HEAD_or_orphan_never_automatically_adopted(self):
        self.reserve(); first_data = (self.root/'ledger/HEAD').read_bytes(); self.event('authorize')
        (self.root/'ledger/HEAD').write_bytes(first_data)
        with self.assertRaisesRegex(ValueError, 'orphan'):
            self.ledger.read(digest(first_data))

    def test_wrong_contract_or_alias_rejected(self):
        self.reserve(); wrong = copy.deepcopy(self.contract); wrong['CPU_cap_ns'] += 1
        with self.assertRaisesRegex(ValueError, 'binding'):
            Ledger(self.root/'ledger', wrong).read(self.expected)
        (self.root/'ledger/foreign').symlink_to(self.root/'ledger/HEAD')
        with self.assertRaisesRegex(ValueError, 'orphan'):
            self.ledger.read(self.expected)

    def killed(self, stage):
        self.reserve(); command = [sys.executable, '-B', str(Path(__file__).with_name('reference_budget_child.py')),
                                  str(self.root), stage, self.expected]
        write_new(self.root/'contract.json', canonical(self.contract))
        with (self.root/'child.stdout').open('xb') as out, (self.root/'child.stderr').open('xb') as err:
            child = subprocess.Popen(command, stdout=out, stderr=err); code = child.wait(timeout=20)
        self.assertEqual(code, -9)
        write_new(self.root/'child.attempt.json', canonical(dict(command=command, pid=child.pid, return_code=code,
                    stage=stage, witness=(self.root/(stage+'.witness')).read_text())))

    def test_real_child_kill_before_HEAD_fails_closed_with_full_charge(self):
        self.killed('before_HEAD')
        with self.assertRaisesRegex(ValueError, 'orphan'):
            self.ledger.read(self.expected)
        self.assertEqual(self.ledger.head(), self.expected)

    def test_real_child_kill_after_HEAD_keeps_unsettled_charge(self):
        self.killed('after_HEAD'); head = self.ledger.head()
        state = self.ledger.read(head)
        self.assertEqual(state['CPU_ns'], 17*10**9)
        self.assertEqual(state['attempts']['attempt_1']['phase'], 'authorized')

    def test_real_wait4_witness_complete_process_accounting(self):
        command = [sys.executable, '-B', str(Path(__file__).with_name('reference_budget_child.py')), str(self.root), 'usage']
        self.reserve(command=command); intent = self.event('authorize'); start = time.monotonic_ns()
        with (self.root/'usage.stdout').open('xb') as out, (self.root/'usage.stderr').open('xb') as err:
            child = subprocess.Popen(command, stdout=out, stderr=err)
            for _ in range(200):
                if (self.root/'usage.identity.json').exists(): break
                time.sleep(.005)
            self.assertTrue((self.root/'usage.identity.json').exists())
            ident = json.loads((self.root/'usage.identity.json').read_text())
            self.assertEqual(ident, process_identity(child.pid))
            birth = self.event('born', process=ident); write_new(self.root/'usage.release', b'go')
            pid, status, usage = os.wait4(child.pid, 0); end = time.monotonic_ns()
            child.returncode = os.waitstatus_to_exitcode(status)
        self.assertEqual(pid, child.pid); self.assertEqual(child.returncode, 0)
        proof = dict(schema=1, complete=True, attempt_id='attempt_1', plan_sha256=self.contract['plan_sha256'],
            intent_sha256=intent['HEAD_sha256'], birth_sha256=birth['HEAD_sha256'], process=ident, supervisor=self.parent,
            measurement=self.contract['measurement'], CPU_ns=ns(Fraction(usage.ru_utime)+Fraction(usage.ru_stime)),
            wall_ns=end-start, start_monotonic_ns=start, end_monotonic_ns=end, return_code=child.returncode)
        write_new(self.root/'measured.wait4.json', canonical(proof))
        state = self.settle(proof)['state']; self.assertFalse(state['halted'])
        self.assertEqual(state['CPU_ns'], 7*10**9+proof['CPU_ns'])
        self.assertEqual(state['wall_ns'], 9*10**9+proof['wall_ns'])


class RealOwnerMemoryControls(unittest.TestCase):
    def setUp(self):
        from reference_checkpoint_fixture import initial, advance, inspector
        from peri_physical_time_friction_audit import PhysicalTimeFrictionInspector
        from peri_physical_friction_negative import OneStepCapture
        self.contract, self.state, self.case, self.header = initial()
        for _ in range(3): self.state = advance(self.state, self.contract, self.case, self.header)
        donor, cap = inspector(self.state, self.case, self.header)
        self.owner = PhysicalTimeFrictionInspector.__new__(PhysicalTimeFrictionInspector)
        self.owner.__dict__.update(copy.deepcopy(donor.__dict__))
        self.capture = OneStepCapture(); self.capture.__dict__.update(copy.deepcopy(cap.__dict__))
        self.owner.capture = self.capture

    def snapshot(self):
        from peri_reference_checkpoint import snapshot
        s = self.state
        return snapshot(self.owner, s['old'], s['reference_viscous'], s['reference_loss'], s['friction_history'],
                        self.capture, s['raw_cursor'], self.contract, s['cpu_completed'])

    def test_actual_class_19_member_snapshot_install_no_aliases(self):
        from peri_reference_checkpoint import encode, install, AUDIT_FIELDS
        self.assertEqual(set(self.owner.__dict__), AUDIT_FIELDS)
        image = self.snapshot(); before = encode(image, self.contract)
        # Static operators are reconstructed by the future reader,not replaced
        # by the serialized history; verify this object's identities survive.
        statics = {k:id(getattr(self.owner,k)) for k in AUDIT_FIELDS-set(image['audit'])-{'capture'}}
        self.owner.maximum = 9e-11; self.owner.statistics['sliding_rows'] += 10
        self.capture.selected = None
        old,z,loss,history = install(image, self.owner, self.capture, self.contract)
        after = self.snapshot(); self.assertEqual(encode(after, self.contract), before)
        self.assertEqual(statics, {k:id(getattr(self.owner,k)) for k in statics})
        self.owner.statistics['sliding_rows'] += 1
        self.assertNotEqual(image['audit']['statistics'], self.owner.statistics)
        self.assertIs(self.owner.capture, self.capture)
        self.assertEqual(old, image['old']); self.assertEqual(z, image['reference_viscous'])
        self.assertEqual(loss, image['reference_loss']); self.assertEqual(history, image['friction_history'])

    def test_real_owner_inherited_close_method_preserves_error_highwater(self):
        from peri_reference_checkpoint import install
        image = self.snapshot(); install(image,self.owner,self.capture,self.contract)
        previous = self.owner.maximum
        self.owner.close(dec.Decimal(1),dec.Decimal(1))
        self.assertEqual(self.owner.maximum,previous)
        self.assertEqual(type(self.owner).__name__,'PhysicalTimeFrictionInspector')
        self.assertEqual(type(self.capture).__name__,'OneStepCapture')

    def test_real_owner_rejects_extra_member_or_active_trial_cache(self):
        self.owner.unrecorded_history = []
        with self.assertRaisesRegex(ValueError,'inventory'):
            self.snapshot()
        del self.owner.unrecorded_history
        self.owner.cache[('uncommitted_trial',)] = dec.Decimal(1)
        with self.assertRaisesRegex(ValueError,'boundary cache'):
            self.snapshot()

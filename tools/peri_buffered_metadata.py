"""Same byte/SHA checks with one check per identical absolute record pin.

No path alias normalization;two differently spelled paths remain separate
pins. Conflicting records reject before IO. Whole-file identity is rechecked
after reading. These functions never grant a numerical qualification.
"""
import hashlib,json,os,re
from pathlib import Path

BUFFER=2*2**20


def unique_records(entries):
    unique={}
    for item in entries:
        if type(item) is not dict or set(item)!={'path','size','sha256'}:raise ValueError('record inventory')
        path=item['path']
        if type(path) is not str or not Path(path).is_absolute() or type(item['size']) is not int or item['size']<0:
            raise ValueError('absolute record path/size')
        if type(item['sha256']) is not str or not re.fullmatch('[0-9a-f]{64}',item['sha256']):raise ValueError('record SHA256')
        if path in unique and unique[path]!=item:raise ValueError('conflicting expected record for identical path')
        unique[path]=dict(item)
    return list(unique.values())


def fingerprint(st):return st.st_dev,st.st_ino,st.st_size,st.st_mtime_ns


def buffered_record(path):
    p=Path(path);before=p.stat();h=hashlib.sha256()
    with p.open('rb',buffering=BUFFER) as stream:
        opened=os.fstat(stream.fileno())
        if fingerprint(opened)!=fingerprint(before):raise ValueError('record identity changed before read')
        for block in iter(lambda:stream.read(BUFFER),b''):h.update(block)
        ended=os.fstat(stream.fileno())
    if fingerprint(before)!=fingerprint(ended) or fingerprint(before)!=fingerprint(p.stat()):raise ValueError('record identity changed during read')
    return dict(path=str(p),size=before.st_size,sha256=h.hexdigest())


def verify_unique(entries,recorder=buffered_record):
    pins=unique_records(entries)
    for expected in pins:
        if recorder(expected['path'])!=expected:raise ValueError('frozen record size/SHA mismatch:'+expected['path'])
    return dict(declared_records=len(entries),checked_unique_paths=len(pins),identical_redundancies=len(entries)-len(pins),
                unique_declared_bytes=sum(p['size'] for p in pins))


class BinaryReader:
    def __init__(self,path):
        self.path=Path(path)
        if self.path.is_symlink() or not self.path.is_file():raise ValueError('regular raw file required')
        self.before=self.path.stat();self.stream=self.path.open('rb',buffering=BUFFER);self.hash=hashlib.sha256()
    def close(self):self.stream.close()
    def row(self):
        line=self.stream.readline(128*2**20+1)
        if not line or len(line)>128*2**20 or not line.endswith(b'\n'):raise ValueError('complete bounded newline record required')
        self.hash.update(line);value=json.loads(line)
        if type(value) is not dict:raise ValueError('raw record inventory')
        return value,dict(offset=self.stream.tell(),prefix_sha256=self.hash.hexdigest(),last_line_sha256=hashlib.sha256(line).hexdigest())
    def finish(self,expected):
        for tail in iter(lambda:self.stream.read(BUFFER),b''):
            self.hash.update(tail)
            if tail.strip():raise ValueError('hidden raw tail')
        if fingerprint(self.before)!=fingerprint(self.path.stat()) or fingerprint(self.before)!=fingerprint(os.fstat(self.stream.fileno())):
            raise ValueError('raw identity changed during assessment')
        if self.stream.tell()!=expected['size'] or self.hash.hexdigest()!=expected['sha256']:raise ValueError('whole raw size/SHA')
        return self.hash.hexdigest()

"""Explicit transport of immutable admission metadata to a NEW audit root."""
import json
from pathlib import Path
from peri_elastic_normal import record,verify

def copy_admission(native_root,audit_root):
    source=Path(native_root)/'header_reference_admission.json'
    data=source.read_bytes();manifest=json.loads(data)
    verify([manifest['raw']])
    for name in ['initial_candidate.bin','initial.phy','contract.bin','describe_initial.attempt.json']:
        assert (Path(native_root)/name).is_file(),('required native input missing',name)
    target=Path(audit_root)/'header_reference_admission.json'
    with target.open('xb') as f:f.write(data)
    assert target.read_bytes()==data
    return [record(source),record(target)]

"""Bounded CLI assembly controls; no S4, fluid or experimental campaign."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--capsule', type=Path, required=True)
    parser.add_argument('--gel', type=Path, required=True)
    args = parser.parse_args()
    source = Path(__file__).resolve().parent.parent
    results = []
    with tempfile.TemporaryDirectory(prefix='peri-model-config-') as temporary:
        root = Path(temporary)

        def check(name, executable, text=None, expected=None, refuse=False, existing=None):
            directory = root / name
            directory.mkdir()
            config = existing
            if text is not None:
                config = directory / 'case.cfg'
                config.write_text(text, encoding='utf-8')
            command = [str(executable), '--pd-checkpoint-capabilities']
            if config is not None:
                command += ['--case-config', str(config)]
            completed = subprocess.run(command, cwd=directory, capture_output=True, text=True, timeout=10)
            assert not any((directory / p).exists() for p in ['output', 'input', 'reload', 'pd_checkpoints']), name
            if refuse:
                assert completed.returncode != 0, name
                assert completed.stderr.strip(), name
            else:
                assert completed.returncode == 0, (name, completed.stderr)
                value = json.loads(completed.stdout)
                for key, expected_value in (expected or {}).items():
                    assert value[key] == expected_value, (name, key, value)
            results.append(dict(name=name, passed=True, refused=refuse))

        for label, executable in [('capsule', args.capsule), ('gel', args.gel)]:
            for config in sorted((source / 'peri/configs/materials').glob('*.cfg')):
                tokens = dict(line.split('=', 1) for line in config.read_text().splitlines()
                              if line and not line.startswith('#'))
                check(label + '_' + config.stem, executable, existing=config,
                      expected={k: tokens[k] for k in ['material_model', 'hardening_model', 'fracture_model']})
            old = source / 'tests/3d_examples' / ('extra_3d_PD_' + ('capsule_compression' if label == 'capsule' else 'gel_compression')) / 'configs'
            for config in sorted(old.glob('*.cfg')):
                check(label + '_legacy_' + config.stem, executable, existing=config)
            check(label + '_default', executable)
            for index, text in enumerate([
                'material_model=unknown\n',
                'material_model=elastic\nfracture_model=principal_plastic_strain\n',
                'material_model=j2\nhardening_model=perfect\nhardening_modulus=10\n',
                'material_model=j2\nplastic_viscosity=1\n',
                'material_model=j2\nfracture_model=unknown\n',
                'material_model=drucker_prager\nhardening_model=perfect\n',
            ]):
                check(label + '_reject_' + str(index), executable, text=text, refuse=True)
        check('gel_reject_fractured_viscoelastic', args.gel,
              text='material_model=viscoelastic\nfracture_model=principal_stress\n', refuse=True)
        check('gel_reject_legacy_conflict', args.gel,
              text='material_model=j2\ngel_viscoelastic=true\n', refuse=True)
        check('capsule_reject_unimplemented_gel_history', args.capsule,
              text='material_model=viscoelastic\n', refuse=True)
        check('capsule_reject_legacy_conflict', args.capsule,
              text='material_model=j2\nshell_viscoplastic=true\n', refuse=True)
        for index, text in enumerate(['youngs_modulus=1e6\nyoungs_shell=2e6\n',
                                      'youngs_shell=2e6\nyoungs_modulus=1e6\n']):
            check('capsule_alias_conflict_' + str(index), args.capsule, text=text, refuse=True)

        # One new dry, tiny, intact test: at most a handful of PD steps, no SPH
        # evolution, no experimental fit, no old state/reviewer/producer reuse.
        directory = root / 'tiny_dry_integration'
        directory.mkdir()
        config = directory / 'case.cfg'
        text = (source / 'peri/configs/materials/j2_perfect.cfg').read_text()
        text += '\n'.join(['', 'capsule_radius_out=6e-5', 'shell_thickness=3e-5', 'shell_layers=3',
            'plate_radius=1e-4', 'plate_thickness=2e-5', 'initial_gap=0', 'compression_distance=2e-9',
            'compression_speed=0.1', 'gravity_g=0', 'with_liquid=false', 'output_frames=1']) + '\n'
        config.write_text(text, encoding='utf-8')
        completed = subprocess.run([str(args.capsule), '--case-config', str(config), '--pd-threads=1',
            '--pd-no-final-checkpoint', '--state_recording=false'], cwd=directory, capture_output=True, text=True, timeout=20)
        assert completed.returncode == 0, (completed.stdout, completed.stderr)
        assert 'completed' in (directory / 'output/run_status.txt').read_text()
        with (directory / 'output/diagnostics.csv').open() as stream:
            records = list(csv.DictReader(stream))
        assert records and 0 < int(records[-1]['pd_steps']) <= 4, records
        assert all(float(r['damage_max']) == 0 and int(r['numerical_failures']) == 0 for r in records)
        check('effective_config_roundtrip', args.capsule,
              existing=directory / 'output/effective_parameters.cfg',
              expected={'material_model': 'j2', 'hardening_model': 'perfect', 'fracture_model': 'none'})
        results.append(dict(name='tiny_dry_integration', passed=True, pd_steps=int(records[-1]['pd_steps'])))
    print(json.dumps(dict(controls=len(results), all_passed=True, cases=results), indent=2))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Check the gel executable's non-solver CLI and a forced pre-allocation stop."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile


def inspect(executable):
    executable = Path(executable).resolve(strict=True)
    checks = []
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)

        def invoke(*arguments, **kwargs):
            return subprocess.run([str(executable), *arguments], cwd=root, capture_output=True,
                                  text=True, timeout=10, **kwargs)

        description = invoke('--describe-model')
        if description.returncode != 0:
            raise ValueError(description.stderr)
        model = json.loads(description.stdout)
        if model['model'] != 'solid_sphere_no_separate_liquid' or model['calibrated']:
            raise ValueError('incorrect model identity or calibration claim')
        if model['max_particles'] != 1000 or model['default_threads'] > 8:
            raise ValueError('bounded execution limits changed')
        checks.append('model_description_without_particle_allocation')

        rejected = invoke('--describe-model', '--pd-threads', '9')
        if rejected.returncode == 0 or 'at most eight' not in rejected.stderr:
            raise ValueError('thread budget was not rejected')
        checks.append('thread_budget_rejected')

        invalid = root / 'invalid.cfg'
        invalid.write_text('poisson_ratio=0.5\n')
        rejected = invoke('--describe-model', '--case-config', str(invalid))
        if rejected.returncode == 0 or 'invalid particle' not in rejected.stderr:
            raise ValueError('singular material was not rejected before allocation')
        checks.append('invalid_material_rejected_before_allocation')

        # A guaranteed absent host path forces failure before SPHSystem exists.
        # This is a negative preflight control, not a small solver simulation.
        environment = dict(os.environ, PERI_HOST_STORAGE=str(root / 'absent-host-volume'))
        stopped = invoke(env=environment)
        if stopped.returncode != 125 or 'host_storage_unknown' not in stopped.stderr:
            raise ValueError('unverifiable host storage did not fail closed')
        if sorted(path.name for path in root.iterdir()) != ['invalid.cfg']:
            raise ValueError('preflight unexpectedly created solver output')
        checks.append('unknown_storage_stops_before_output_or_particles')
    return {'status': 'preflight_checks_passed', 'checks': checks, 'solver_simulation_performed': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(inspect(args.executable), indent=2))

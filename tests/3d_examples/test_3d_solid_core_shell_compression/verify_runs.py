"""Run and inspect the uncalibrated prototype; no experimental acceptance implied."""
import argparse
import bisect
import concurrent.futures
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inspect(folder, expected_source):
    summary = json.loads((folder / 'summary.json').read_text())
    params = dict(line.split('=', 1) for line in (folder / 'parameters.txt').read_text().splitlines())
    rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader((folder / 'history.csv').open())]
    assert rows and all(math.isfinite(v) for row in rows for v in row.values()), 'non-finite history'
    assert params['source_sha256'] == summary['source_sha256'] == expected_source, 'source mismatch'
    assert math.isclose(rows[-1]['time_s'], float(params['load_time']) + float(params['hold_time']), rel_tol=1e-10), 'incomplete duration'
    assert summary['status'] == 'COMPLETED_DEMONSTRATION'
    assert max(row['upper_N'] for row in rows) > 1e-10 and max(row['lower_N'] for row in rows) > 1e-10, 'no reciprocal plate loading'
    assert min(row['min_J'] for row in rows) > 0.2, 'invalid deformation'
    tail = [row for row in rows if row['time_s'] >= .8 * rows[-1]['time_s']]
    upper = statistics.mean(row['upper_N'] for row in tail)
    lower = statistics.mean(row['lower_N'] for row in tail)
    radius = float(params['diameter']) / 2
    if params['case'] == 'coated':
        radius -= float(params['thickness'])
    mass_ref = 4 * math.pi / 3 * radius**3 * float(params['core_density'])
    assert math.isclose(float(params['core_mass']), mass_ref, rel_tol=1e-12)
    if params['core_material'] == 'j2':
        assert max(row['max_plastic_alpha'] for row in rows) > 0, 'plastic update not exercised'
    return {
        'numerical_run_completed': True,
        'steps': int(rows[-1]['step']),
        'core_particles': int(params['core_particles']),
        'shell_particles': int(params['shell_particles']),
        'min_J': min(row['min_J'] for row in rows),
        'max_plastic_alpha': max(row['max_plastic_alpha'] for row in rows),
        'max_tie_gap_m': max(row['max_tie_gap_m'] for row in rows),
        'max_interface_force_residual': max(row['interface_force_residual'] for row in rows),
        'max_interface_power_residual_W': max(abs(row['interface_power_residual_W']) for row in rows),
        'last_20pct_time_mean_upper_N': upper,
        'last_20pct_time_mean_lower_N': lower,
        'mean_plate_imbalance': abs(upper-lower) / ((upper+lower)/2),
        'experimental_validation': False,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    exe = args.exe.resolve(strict=True)
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()):
        raise RuntimeError('verification output must be new or empty')
    out.mkdir(parents=True, exist_ok=True)
    executable_sha = sha(exe)
    source_sha = sha(Path(__file__).with_name('solid_capsule.cpp'))
    jobs = [
        ('bare_elastic', ['--case=bare']),
        ('coated_elastic', ['--case=coated']),
        ('bare_j2', ['--case=bare', '--material=j2']),
        ('coated_j2', ['--case=coated', '--material=j2']),
        ('coated_halfdt', ['--case=coated', '--dt-scale=0.5']),
    ]
    receipt = {'source_sha256': source_sha, 'executable_sha256': executable_sha, 'runs': {},
               'qualification': 'runtime_and_interface_checks_only', 'experimental_calibration': False,
               'hertz_spatial_contact_penalty_and_quasistatic_qualification': False,
               'complete_system_energy_balance': False}

    def run(job):
        name, options = job
        command = [str(exe), *options, '--output=' + str(out / name)]
        start = time.monotonic()
        with (out / (name + '.log')).open('w') as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, cwd=out, timeout=600)
        data = {'command': command, 'exit_code': result.returncode, 'wall_seconds': time.monotonic()-start}
        assert result.returncode == 0, f'{name} solver failed; see its log'
        data.update(inspect(out / name, source_sha))
        return name, data

    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for future in concurrent.futures.as_completed([pool.submit(run, job) for job in jobs]):
            name, data = future.result()
            receipt['runs'][name] = data
            (out / 'verification.json').write_text(json.dumps(receipt, indent=2))
            print(name, 'COMPLETED', flush=True)
    assert sha(exe) == executable_sha, 'executable changed during verification'
    assert sha(Path(__file__).with_name('solid_capsule.cpp')) == source_sha, 'source changed during verification'
    def history(name):
        return [{k: float(v) for k,v in r.items()} for r in csv.DictReader((out/name/'history.csv').open())]
    base, fine = history('coated_elastic'), history('coated_halfdt')
    fine_time = [r['time_s'] for r in fine]
    error = norm = 0
    for row in base:
        t = row['time_s']
        if t < 3e-5:
            continue
        j = min(max(bisect.bisect_right(fine_time, t), 1), len(fine)-1)
        a, b = fine[j-1], fine[j]
        w = (t-a['time_s']) / (b['time_s']-a['time_s'])
        f = .5 * (row['upper_N']+row['lower_N'])
        ff = .5 * ((1-w)*(a['upper_N']+a['lower_N']) + w*(b['upper_N']+b['lower_N']))
        error += (f-ff)**2
        norm += f*f
    receipt['coated_halfdt_relative_L2_mean_plate_force'] = math.sqrt(error/norm)
    (out/'verification.json').write_text(json.dumps(receipt, indent=2))
    print('Verification complete; this is not material calibration or G3 physical qualification.', flush=True)


if __name__ == '__main__':
    main()

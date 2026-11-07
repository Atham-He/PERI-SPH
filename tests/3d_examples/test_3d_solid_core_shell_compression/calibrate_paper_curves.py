"""Whole-curve, two-parameter calibration followed by a frozen-core shell run.

Standard-library-only driver. Every objective evaluation solves the complete
displacement history. For H=0 and fixed nu, geometry, mass, dimensionless
time-step/loading/damping controls and penalty factors, E scales force while
sigma_y/E controls shape. A final solve with the actual fitted parameters checks
this reduction. This is an effective calibration at the specified discretization,
not a claim of unique intrinsic material constants or spatial convergence.
"""
import argparse
import bisect
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def interpolate(points, x):
    xs = [p[0] for p in points]
    i = max(1, min(bisect.bisect_right(xs, x), len(xs) - 1))
    if not xs[0] <= x <= xs[-1]:
        raise ValueError('requested displacement outside digitized experimental curve')
    a, b = points[i - 1], points[i]
    return a[1] + (b[1] - a[1]) * (x - a[0]) / (b[0] - a[0])


def metrics(predicted, reference):
    peak = max(reference)
    errors = [a - b for a, b in zip(predicted, reference)]
    return dict(rmse_N=math.sqrt(sum(e * e for e in errors) / len(errors)),
                nrmse=math.sqrt(sum(e * e for e in errors) / len(errors)) / peak,
                max_abs_error_N=max(abs(e) for e in errors),
                max_abs_error_over_peak=max(abs(e) for e in errors) / peak,
                mean_signed_error_N=sum(errors) / len(errors))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ['exe', 'base-config', 'paper-bare', 'paper-coated', 'source', 'out']:
        parser.add_argument('--' + option, type=Path, required=True)
    parser.add_argument('--source-commit', required=True)
    parser.add_argument('--resume-control-change', action='store_true',
                        help='reuse completed solver histories only if solver/data/config hashes are unchanged')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'configs').mkdir(exist_ok=True)
    base = args.base_config.read_text()
    expected_hashes = dict(exe_sha256=digest(args.exe), source_sha256=digest(args.source),
                           base_config_sha256=digest(args.base_config),
                           paper_bare_sha256=digest(args.paper_bare),
                           paper_coated_sha256=digest(args.paper_coated),
                           driver_sha256=digest(__file__), source_commit=args.source_commit)
    receipt = args.out / 'launch_receipt.json'
    if receipt.exists():
        previous = json.loads(receipt.read_text())
        if previous != expected_hashes:
            assert args.resume_control_change, 'resume provenance mismatch'
            unchanged = set(expected_hashes) - {'driver_sha256', 'source_commit'}
            assert all(previous.get(k) == expected_hashes[k] for k in unchanged), 'solver/input changed on resume'
            archive = args.out / ('launch_receipt_control_' + digest(receipt)[:12] + '.json')
            if not archive.exists():
                archive.write_bytes(receipt.read_bytes())
            save(receipt, expected_hashes)
    else:
        save(receipt, expected_hashes)
    def paper(path):
        return sorted((float(r['strain_percent']), float(r['force_N']))
                      for r in csv.DictReader(path.open()))
    bare_paper, coated_paper = paper(args.paper_bare), paper(args.paper_coated)
    x = [.5 * i for i in range(1, 21)]
    reference = [interpolate(bare_paper, p) for p in x]
    coated_reference = [interpolate(coated_paper, p) for p in x]

    def simulate(name, E, sy, coated=False):
        config = args.out / 'configs' / (name + '.cfg')
        text = base + f'\n# Frozen explicit overrides for {name}\ncase={"coated" if coated else "bare"}\n'
        text += f'core-E={E:.17g}\nyield-stress={sy:.17g}\n'
        if config.exists():
            assert config.read_text() == text, 'configuration changed on resume'
        else:
            config.write_text(text)
        run = args.out / name
        started = time.monotonic()
        if not (run / 'summary.json').exists():
            print(f'START {name} E={E:.8g} sy={sy:.8g}', flush=True)
            with (args.out / (name + '.log')).open('w') as log:
                subprocess.run([str(args.exe), '--config=' + str(config), '--output=' + str(run)],
                               stdout=log, stderr=subprocess.STDOUT, check=True)
        summary = json.loads((run / 'summary.json').read_text())
        assert summary['status'] == 'COMPLETED_CONFIGURED_QUASISTATIC_CASE'
        assert summary['source_sha256'] == expected_hashes['source_sha256']
        assert summary['equilibrated_load_steps'] == 20
        params = dict(line.split('=', 1) for line in (run / 'parameters.txt').read_text().splitlines())
        assert float(params['core_E']) == E and float(params['yield_stress']) == sy
        assert float(params['hardening']) == 0 and float(params['core_nu']) == .3
        rows = list(csv.DictReader((run / 'equilibria.csv').open()))
        assert len(rows) == 20
        for i, row in enumerate(rows):
            assert all(math.isfinite(float(v)) for v in row.values())
            assert math.isclose(float(row['strain']) * 100, x[i], abs_tol=1e-9)
            for field, limit in [('force_residual', .02), ('plate_imbalance', .01), ('kinetic_ratio', 1e-5)]:
                assert float(row[field]) <= limit
        force = [(float(r['upper_N']) + float(r['lower_N'])) / 2 for r in rows]
        print(f'DONE {name} F10={force[-1]:.6g} N; {time.monotonic()-started:.1f} s', flush=True)
        return force, summary

    evaluations = []
    def evaluate(log_ratio):
        for result in evaluations:
            if abs(result['log10_yield_over_E'] - log_ratio) < 1e-9:
                return result
        ratio = 10 ** log_ratio
        name = f'fit_ratio_{len(evaluations):02d}'
        failure_file = args.out / (name + '.failure.json')
        if failure_file.exists():
            result = json.loads(failure_file.read_text())
            evaluations.append(result)
            return result
        try:
            values, summary = simulate(name, 1e9, 1e9 * ratio)
        except subprocess.CalledProcessError as error:
            result = dict(run=name, log10_yield_over_E=log_ratio, yield_over_E=ratio,
                          status='FAILED_NUMERICAL_EQUILIBRIUM', solver_exit_code=error.returncode)
            save(failure_file, result)
            evaluations.append(result)
            save(args.out / 'optimization_trace.json', evaluations)
            print(f'REJECT {name}: solver did not complete; no fit data used', flush=True)
            return result
        # Bound the two material parameters without changing numerical controls.
        Emin, Emax = max(2e7, 1e5 / ratio), min(2e10, 3e8 / ratio)
        scale = sum(a*b for a, b in zip(values[:16], reference[:16])) / sum(a*a for a in values[:16])
        scale = max(Emin / 1e9, min(Emax / 1e9, scale))
        predicted = [v * scale for v in values]
        result = dict(run=name, status='COMPLETED', log10_yield_over_E=log_ratio, yield_over_E=ratio,
                      E_Pa=1e9*scale, yield_stress_Pa=1e9*ratio*scale,
                      fit=metrics(predicted[:16], reference[:16]),
                      holdout=metrics(predicted[16:], reference[16:]),
                      scaled_force_N=predicted,
                      max_plastic_alpha=float(list(csv.DictReader(
                          (args.out / name / 'equilibria.csv').open()))[-1]['max_plastic_alpha']))
        evaluations.append(result)
        save(args.out / 'optimization_trace.json', evaluations)
        print(f'FIT ratio={ratio:.6g} E={result["E_Pa"]:.6g} sy={result["yield_stress_Pa"]:.6g} '
              f'fit_nrmse={result["fit"]["nrmse"]:.5f}', flush=True)
        return result

    # A bounded one-dimensional material-shape search, not a suite of prototypes.
    def objective(result):
        return result['fit']['nrmse'] if result['status'] == 'COMPLETED' else math.inf
    def adequate_fit(result):
        # Stop based only on the fitting window. Holdout NEVER selects parameters.
        return objective(result) <= .05 and result['fit']['max_abs_error_over_peak'] <= .10
    grid = [-3., -2., -1., -.5]
    adequate = None
    for value in [-2., -1., -3., -.5]:
        current = evaluate(value)
        if adequate_fit(current):
            adequate = current
            print('ADEQUATE FIT: stop material search; holdout has not selected parameters', flush=True)
            break
    if adequate is None:
        best_grid = min(range(len(grid)), key=lambda i: objective(evaluate(grid[i])))
        low, high = grid[max(0, best_grid - 1)], grid[min(len(grid) - 1, best_grid + 1)]
        golden = (math.sqrt(5) - 1) / 2
        a, b = high - golden * (high - low), low + golden * (high - low)
        ea, eb = evaluate(a), evaluate(b)
        for _ in range(5):
            if adequate_fit(ea) or adequate_fit(eb):
                adequate = ea if adequate_fit(ea) else eb
                break
            if objective(ea) < objective(eb):
                high, b, eb = b, a, ea
                a = high - golden * (high - low)
                ea = evaluate(a)
            else:
                low, a, ea = a, b, eb
                b = low + golden * (high - low)
                eb = evaluate(b)
    best = min(evaluations, key=objective)
    if not math.isfinite(objective(best)):
        raise RuntimeError('no material candidate completed equilibrium; calibration not accepted')
    E, sy = float(f'{best["E_Pa"]:.8g}'), float(f'{best["yield_stress_Pa"]:.8g}')
    actual, bare_summary = simulate('bare_calibrated', E, sy)
    fit = metrics(actual[:16], reference[:16])
    holdout = metrics(actual[16:], reference[16:])
    passed = fit['nrmse'] <= .10 and fit['max_abs_error_over_peak'] <= .20 and holdout['nrmse'] <= .15
    near_best = [r for r in evaluations if objective(r) <= objective(best) + .01]
    frozen = dict(status='ACCEPTED_EFFECTIVE_CALIBRATION' if passed else 'PROVISIONAL_BEST_FIT_NOT_ACCEPTED',
                  constitutive='native finite-strain J2 perfect plasticity', core_E_Pa=E,
                  core_yield_stress_Pa=sy, core_hardening_Pa=0, core_nu=.3, core_density_kg_m3=1355,
                  fit=fit, holdout=holdout, best_ratio_run=best['run'],
                  force_scaling_final_solve_max_difference_N=max(abs(a-b) for a, b in zip(actual, best['scaled_force_N'])),
                  nearly_equal_fit_E_range_Pa=[min(r['E_Pa'] for r in near_best), max(r['E_Pa'] for r in near_best)],
                  nearly_equal_fit_yield_range_Pa=[min(r['yield_stress_Pa'] for r in near_best), max(r['yield_stress_Pa'] for r in near_best)],
                  parameters_are_discretization_dependent=True, unique_intrinsic_material_constants_claimed=False,
                  identifiability_profile_completed=False,
                  search_stopped_for_adequate_fit=adequate is not None,
                  **expected_hashes)
    freeze_file = args.out / 'frozen_core.json'
    save(freeze_file, frozen)
    freeze_hash = digest(freeze_file)
    print(f'GATE1 {frozen["status"]}; frozen_core_sha256={freeze_hash}', flush=True)

    coated, coated_summary = simulate('coated_frozen_core', E, sy, coated=True)
    assert digest(freeze_file) == freeze_hash, 'core parameters changed during coated comparison'
    comparison = metrics(coated, coated_reference)
    result = dict(gate1_bare_calibration=frozen['status'],
                  gate2_coated_calculation='COMPLETED',
                  gate2_experimental_agreement='PASS' if comparison['nrmse'] <= .15 else 'FAIL',
                  coated_comparison=comparison, core_parameters_identical=True,
                  frozen_core_sha256=freeze_hash, shell_E_Pa=6e8,
                  shell_parameters_fitted=False, bare_summary=bare_summary, coated_summary=coated_summary)
    save(args.out / 'calibration_comparison_summary.json', result)
    with (args.out / 'comparison.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['strain_percent', 'compression_um', 'paper_bare_N', 'sim_bare_N',
                         'paper_coated_N', 'sim_coated_N', 'bare_is_fit'])
        writer.writerows((p, p*10, a, b, c, d, int(i < 16))
                         for i, (p, a, b, c, d) in enumerate(zip(x, reference, actual, coated_reference, coated)))
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()

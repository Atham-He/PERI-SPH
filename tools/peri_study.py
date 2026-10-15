#!/usr/bin/env python3
"""Sealed, sequential capsule time/space study; no partial run is convergence evidence."""
import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone

from peri_run import digest, parameters, verify, require_completion_audit


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if path.name == "study.json":
        fields = ("status", "active_case", "failed_case", "pending_case", "reason", "return_code", "source_git_head")
        summary = {key: value[key] for key in fields if key in value}
        summary["updated_at_utc"] = datetime.now(timezone.utc).isoformat()
        (path.parent / "status.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def plan(args):
    if not 0 < args.tolerance < 1:
        raise ValueError("study tolerance must be in (0,1)")
    root = args.root.resolve()
    root.mkdir(parents=True, exist_ok=False)
    config = parameters(args.config)
    if not math.isclose(config["capsule_radius_out"], 5e-4) or not math.isclose(config["compression_distance"], 5.8e-4):
        raise ValueError("full study requires the original 1 mm diameter and full 0.58 mm stroke")
    layers = config.get("shell_layers", 3)
    config["pd_horizon"] = config.get("pd_horizon", 3.015 * config["shell_thickness"] / layers)
    if layers != 3:
        raise ValueError("this three/four/five-layer plan requires a three-layer baseline")
    config["shell_layers"] = layers
    repo = Path(__file__).resolve().parents[1]
    frozen = root / "frozen_source"
    frozen.mkdir()
    for name in ("src", "cmake", "tests/3d_examples/extra_3d_PD_capsule_compression"):
        shutil.copytree(repo / name, frozen / name)
    for name in ("CMakeLists.txt", "PERI_SPH_FRAMEWORK.txt", "PERI_SPH_VALIDATION.txt", "PERI_SPH_VERSION.txt", "UPSTREAM_SYNC.md"):
        shutil.copy2(repo / name, frozen / name)
    (frozen / "tools").mkdir()
    for script in (repo / "tools").glob("peri_*.py"):
        shutil.copy2(script, frozen / "tools" / script.name)
    (frozen / "bin").mkdir()
    executable = frozen / "bin" / args.executable.name
    shutil.copy2(args.executable, executable)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo, capture_output=True, text=True)
    patch = subprocess.run(["git", "diff", "--binary", "HEAD"], cwd=repo, capture_output=True)
    (root / "source_changes.patch").write_bytes(patch.stdout)
    cases = []
    for name, group, step_factor, level in (
        ("baseline", "baseline", 1, layers),
        ("dt_half", "time", 0.5, layers),
        ("dt_quarter", "time", 0.25, layers),
        ("space_4", "space", 0.25, 4),
        ("space_5", "space", 0.25, 5),
    ):
        values = dict(config)
        values.update(pd_cfl=config.get("pd_cfl", 0.1) * step_factor,
                      acoustic_cfl=config.get("acoustic_cfl", 0.3) * step_factor,
                      shell_layers=level)
        target = root / (name + ".cfg")
        target.write_text("".join(f"{key}={value:.17g}\n" for key, value in sorted(values.items())), encoding="utf-8")
        cases.append({"name": name, "group": group, "config": target.name,
                      "config_sha256": digest(target), "shell_layers": level,
                      "spacing_m": config["shell_thickness"] / level,
                      "pd_horizon_m": config["pd_horizon"]})
    source_hashes = {str(p.relative_to(frozen)): digest(p) for p in sorted(frozen.rglob("*")) if p.is_file()}
    data = {"created_at_utc": datetime.now(timezone.utc).isoformat(), "status": "planned",
            "source_git_head": head.stdout.strip() if head.returncode == 0 else None,
            "source_changes_sha256": digest(root / "source_changes.patch"),
            "executable": str(executable.relative_to(root)), "source_sha256": source_hashes,
            "cases": cases, "force_relative_tolerance": args.tolerance,
            "energy_monitoring_tolerance": 0.05,
            "space_policy": "fixed physical PD horizon, fixed geometry/material/loading; shell/fluid/plate particle spacings refine together",
            "time_policy": "halve both solid and acoustic CFL; keep exchange substep cap fixed",
            "memory_policy": "conservative screening using measured baseline RSS times refinement ratio^6; profile refined cases on a larger node for an exact estimate"}
    save(root / "study.json", data)
    print(json.dumps({"root": str(root), "cases": cases, "status": "planned"}, indent=2))


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]


def relative_curve_error(a, b, field="total_Fy_N"):
    """Piecewise linear interpolation, integrated on the reference stroke grid."""
    if a[-1]["stroke_m"] != b[-1]["stroke_m"] and not math.isclose(a[-1]["stroke_m"], b[-1]["stroke_m"], rel_tol=1e-12):
        raise ValueError("curves do not cover the same full stroke")
    j = 0
    errors, reference, grid = [], [], []
    for row in b:
        x = row["stroke_m"]
        while j + 1 < len(a) - 1 and a[j+1]["stroke_m"] < x:
            j += 1
        u, v = a[j], a[j+1]
        fraction = (x-u["stroke_m"]) / (v["stroke_m"]-u["stroke_m"])
        value = u[field] + fraction * (v[field]-u[field])
        errors.append((value-row[field])**2); reference.append(row[field]**2); grid.append(x)
    def integrate(values):
        return sum(0.5*(u+v)*(y-x) for u,v,x,y in zip(values,values[1:],grid,grid[1:]))
    return math.sqrt(integrate(errors) / max(integrate(reference), 1e-300))


def compare_group(root, data, names, group):
    configs, curves, budgets, verification, audits = {}, {}, {}, {}, {}
    for name in names:
        run = root / "runs" / name
        verification[name] = verify(run)
        audits[name] = require_completion_audit(run)
        manifest = json.loads((run / "run_manifest.json").read_text(encoding="utf-8"))
        if manifest["executable_sha256"] != data["source_sha256"][data["executable"].removeprefix("frozen_source/")]:
            raise ValueError("different executable used in " + name)
        configs[name] = parameters(run / "output/effective_parameters.cfg")
        curves[name] = read_rows(run / "output/diagnostics.csv")
        budgets[name] = read_rows(run / "output/energy_budget.csv")
    allowed = {"pd_cfl", "acoustic_cfl"} if group == "time" else {"shell_layers"}
    base = configs[names[0]]
    for name in names[1:]:
        if configs[name].keys() != base.keys() or any(configs[name][key] != value for key,value in base.items() if key not in allowed):
            raise ValueError("uncontrolled parameter changed in " + name)
    reference = names[-1]
    result = {"group": group, "reference": reference, "comparisons": [], "mass_and_energy_checks": verification, "completion_audits": audits}
    def work(name):
        row = budgets[name][-1]
        return sum(row[key] for key in ("plate_contact_work_J", "plate_fluid_pressure_work_J", "plate_fluid_viscous_work_J"))
    for name in names[:-1]:
        peak = max(row["total_Fy_N"] for row in curves[name])
        refpeak = max(row["total_Fy_N"] for row in curves[reference])
        result["comparisons"].append({"case": name,
            "reaction_curve_relative_L2": relative_curve_error(curves[name], curves[reference]),
            "peak_reaction_relative_error": abs(peak-refpeak) / max(abs(refpeak),1e-300),
            "plate_work_relative_error": abs(work(name)-work(reference)) / max(abs(work(reference)),1e-300),
            "final_damage_mean_absolute_error": abs(curves[name][-1]["damage_mean"]-curves[reference][-1]["damage_mean"])})
    adjacent = result["comparisons"][-1]
    threshold = data["force_relative_tolerance"]
    result["adjacent_difference_within_tolerance"] = all(adjacent[key] <= threshold for key in
        ("reaction_curve_relative_L2", "peak_reaction_relative_error", "plate_work_relative_error", "final_damage_mean_absolute_error"))
    result["energy_monitoring_passed"] = all(v["mechanical_energy_residual_max_relative"] <= data["energy_monitoring_tolerance"] for v in verification.values())
    coarse=result["comparisons"][0]
    result["refinement_difference_trend"] = {
        key: coarse[key] + 1e-10 >= adjacent[key]
        for key in ("reaction_curve_relative_L2", "plate_work_relative_error")}
    result["accepted"] = (result["adjacent_difference_within_tolerance"] and
                          all(result["refinement_difference_trend"].values()) and
                          result["energy_monitoring_passed"])
    result["limits"] = "three-level differences with a provisional tolerance; not a formal asymptotic order or experimental validation"
    save(root / (group + "_comparison.json"), result)
    return result


def run(args):
    if args.memory_limit_gib <= 0 or args.baseline_memory_gib <= 0:
        raise ValueError("memory limits must be positive")
    root=args.root.resolve(); data=json.loads((root/"study.json").read_text(encoding="utf-8"))
    frozen=root/"frozen_source"
    for name, sha in data["source_sha256"].items():
        if digest(frozen/name) != sha:
            raise ValueError("sealed study source or executable changed: " + name)
    for case in data["cases"]:
        if digest(root/case["config"]) != case["config_sha256"]:
            raise ValueError("study config changed: "+case["name"])
        if case["group"] == "space":
            time_result=compare_group(root,data,["baseline","dt_half","dt_quarter"],"time")
            if not time_result["accepted"]:
                data.update(status="time_validation_failed", reason="time or energy gate failed; space runs were not launched")
                save(root/"study.json",data); return 2
        destination=root/"runs"/case["name"]
        if destination.exists():
            try:
                verify(destination)  # interrupted or failed simulations cannot silently be reused
                require_completion_audit(destination, snapshots=not args.no_states)
            except (OSError, ValueError, KeyError) as error:
                data.update(status="audit_failed", failed_case=case["name"], reason=str(error))
                save(root/"study.json", data); return 2
            continue
        ratio=case["shell_layers"]/data["cases"][0]["shell_layers"]
        estimate=args.baseline_memory_gib*1.05*ratio**6
        meminfo = Path("/proc/meminfo")
        available = math.inf
        if meminfo.exists():
            available = next(float(line.split()[1]) * 1024 / 2**30 for line in meminfo.read_text().splitlines() if line.startswith("MemAvailable:"))
        if estimate > min(args.memory_limit_gib, available):
            data.update(status="resource_limited", pending_case=case["name"], estimated_screening_memory_gib=estimate,
                        reason="conservative memory screen exceeded; no refined simulation was launched")
            save(root/"study.json",data); return 3
        data.update(status="running", active_case=case["name"]); save(root/"study.json",data)
        command=[sys.executable,str(frozen/"tools/peri_run.py"),"run","--executable",str(root/data["executable"]),
                 "--config",str(root/case["config"]),"--run-dir",str(destination)]
        if args.no_states: command.append("--no-states")
        rc=subprocess.run(command).returncode
        if rc:
            data.update(status="simulation_failed", failed_case=case["name"],return_code=rc)
            save(root/"study.json",data); return rc
        try:
            require_completion_audit(destination, snapshots=not args.no_states)
        except (OSError, ValueError, KeyError) as error:
            data.update(status="audit_failed", failed_case=case["name"], reason=str(error))
            save(root/"study.json", data); return 2
    compare_group(root,data,["baseline","dt_half","dt_quarter"],"time")
    spatial=compare_group(root,data,["dt_quarter","space_4","space_5"],"space")
    data.update(status="accepted" if spatial["accepted"] else "space_validation_failed")
    save(root/"study.json",data)
    return 0 if spatial["accepted"] else 2


def main():
    parser=argparse.ArgumentParser(description=__doc__); sub=parser.add_subparsers(dest="command",required=True)
    p=sub.add_parser("plan"); p.add_argument("--root",type=Path,required=True)
    p.add_argument("--config",type=Path,required=True); p.add_argument("--executable",type=Path,required=True)
    p.add_argument("--tolerance",type=float,default=0.05)
    r=sub.add_parser("run"); r.add_argument("--root",type=Path,required=True)
    r.add_argument("--memory-limit-gib",type=float,required=True)
    r.add_argument("--baseline-memory-gib",type=float,default=13.0)
    r.add_argument("--no-states",action="store_true")
    c=sub.add_parser("compare"); c.add_argument("--root",type=Path,required=True)
    c.add_argument("--group",choices=["time","space"],required=True)
    args=parser.parse_args()
    try:
        if args.command=="plan": plan(args); return 0
        if args.command=="run": return run(args)
        data=json.loads((args.root/"study.json").read_text(encoding="utf-8"))
        names=["baseline","dt_half","dt_quarter"] if args.group=="time" else ["dt_quarter","space_4","space_5"]
        print(json.dumps(compare_group(args.root,data,names,args.group),indent=2)); return 0
    except (ValueError,OSError,KeyError) as exc:
        print("PERI-SPH study rejected: "+str(exc),file=sys.stderr); return 1


if __name__=="__main__":
    raise SystemExit(main())

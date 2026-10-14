#!/usr/bin/env python3
"""Create isolated, reproducible PERI-SPH runs and verify their recorded state."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import uuid


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parameters(path):
    result = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            key, value = line.split("=", 1)
            result[key.strip()] = float(value.strip())
    return result


def status(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    result = {"status": lines[0]}
    result.update(line.split("=", 1) for line in lines[1:] if "=" in line)
    return result


def verify(run_dir, require_damage=False):
    output = run_dir / "output"
    record = status(output / "run_status.txt")
    if record["status"] != "completed":
        raise ValueError("simulation did not complete: " + record.get("reason", record["status"]))
    config = parameters(output / "effective_parameters.cfg")
    with (output / "diagnostics.csv").open(newline="", encoding="utf-8") as stream:
        rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]
    if len(rows) != int(config["output_frames"]) + 1:
        raise ValueError("incomplete diagnostic output")
    if not all(math.isfinite(value) for row in rows for value in row.values()):
        raise ValueError("diagnostics contain non-finite data")
    expected_time = config["compression_distance"] / config["compression_speed"]
    if not math.isclose(rows[-1]["time_s"], expected_time, rel_tol=1e-12, abs_tol=1e-16):
        raise ValueError("end time does not match prescribed stroke and speed")
    if not math.isclose(rows[-1]["stroke_m"], config["compression_distance"], rel_tol=1e-12, abs_tol=1e-16):
        raise ValueError("final stroke is incorrect")
    if any(b["time_s"] <= a["time_s"] for a, b in zip(rows, rows[1:])):
        raise ValueError("output time is not strictly increasing")
    if any(row["numerical_failures"] != 0 for row in rows):
        raise ValueError("numerically failed shell particles were recorded")
    if any(not (0 <= row["damage_mean"] <= row["damage_max"] <= 1) for row in rows):
        raise ValueError("invalid damage range")
    if any(abs(row["total_Fy_N"] - row["contact_Fy_N"] - row["fluid_Fy_N"]) > 1e-12 for row in rows):
        raise ValueError("reaction components do not reconcile")
    if require_damage and max(row["damage_max"] for row in rows) <= 0:
        raise ValueError("fracture regression did not produce damage")
    mass_error = 0.0
    if config["with_liquid"]:
        initial_mass = rows[0]["liquid_mass_kg"]
        if initial_mass <= 0:
            raise ValueError("missing liquid mass")
        mass_error = max(abs(row["liquid_mass_kg"] / initial_mass - 1) for row in rows)
        if mass_error > 1e-12:
            raise ValueError("liquid particle mass is not conserved")
        if max(row["fluid_density_error_max"] for row in rows) > config["max_density_error"]:
            raise ValueError("recorded fluid density exceeds configured limit")
        if max(row["fluid_speed_max_m_s"] for row in rows) > config["max_fluid_mach"] * config["c_f"]:
            raise ValueError("recorded Mach number exceeds configured limit")
    result = {
        "status": "verified",
        "diagnostic_rows": len(rows),
        "end_time_s": rows[-1]["time_s"],
        "stroke_m": rows[-1]["stroke_m"],
        "liquid_mass_relative_error": mass_error,
        "fluid_density_error_max": max(row["fluid_density_error_max"] for row in rows),
        "fluid_speed_max_m_s": max(row["fluid_speed_max_m_s"] for row in rows),
        "peak_total_reaction_N": max(row["total_Fy_N"] for row in rows),
        "acoustic_steps": int(rows[-1]["acoustic_steps"]),
        "pd_steps": int(rows[-1]["pd_steps"]),
        "scope": "output integrity and numerical limits; physical convergence and calibration require separate studies",
    }
    (run_dir / "verification.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, allow_nan=False))
    return result


def launch(args):
    executable = args.executable.resolve(strict=True)
    config = args.config.resolve(strict=True)
    run_dir = args.run_dir.resolve() if args.run_dir else args.run_root.resolve() / ("run-" + uuid.uuid4().hex)
    run_dir.mkdir(parents=True, exist_ok=False)
    shutil.copy2(config, run_dir / "case.cfg")
    repo = Path(__file__).resolve().parents[1]
    paths = [repo / "src", repo / "cmake", repo / "CMakeLists.txt", repo / "tools/peri_run.py"]
    paths += list((repo / "tests/3d_examples").glob("extra_3d_PD_*"))
    paths += [repo / name for name in ("src/shared/include/sphinxsys.h", "src/shared/particle_neighborhood/neighborhood.h", "src/shared/particle_neighborhood/neighborhood.cpp")]
    files = sorted({p for path in paths for p in (path.rglob("*") if path.is_dir() else [path]) if p.is_file()})
    target_case = repo / "tests/3d_examples" / executable.stem
    stale = [str(p.relative_to(repo)) for p in files if (p.is_relative_to(repo / "src/shared") or p.is_relative_to(repo / "src/for_3D_build") or p.is_relative_to(target_case)) and p.suffix in (".cpp", ".h", ".hpp", ".hxx") and p.stat().st_mtime > executable.stat().st_mtime + 0.001]
    if stale:
        raise ValueError("source is newer than the executable; rebuild before running: " + stale[0])
    hashes = {str(p.relative_to(repo)): digest(p) for p in files}
    head, patch = None, b""
    try:
        revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo, capture_output=True, text=True)
        if revision.returncode == 0:
            head = revision.stdout.strip()
            patch = subprocess.run(["git", "diff", "--binary", "HEAD"], cwd=repo, capture_output=True, check=True).stdout
    except FileNotFoundError:
        pass  # A distributed source archive need not include Git; source hashes remain recorded.
    (run_dir / "source_changes.patch").write_bytes(patch)
    with tarfile.open(run_dir / "source_snapshot.tar.gz", "w:gz") as archive:
        for p in files:
            archive.add(p, arcname=str(p.relative_to(repo)), recursive=False)
    metadata = {
        "started_at_utc": datetime.now(timezone.utc).isoformat(), "git_head": head,
        "executable": str(executable), "executable_sha256": digest(executable),
        "case_sha256": digest(config), "source_sha256": hashes,
    }
    (run_dir / "run_manifest.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    command = [str(executable), "--case-config", "case.cfg", "--state_recording=" + ("false" if args.no_states else "true")]
    with (run_dir / "run.log").open("w", encoding="utf-8") as log:
        process = subprocess.run(command, cwd=run_dir, stdout=log, stderr=subprocess.STDOUT)
    metadata.update(finished_at_utc=datetime.now(timezone.utc).isoformat(), return_code=process.returncode)
    (run_dir / "run_manifest.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    if args.command == "fault-test":
        record = status(run_dir / "output/run_status.txt")
        if process.returncode != 1 or record["status"] != "failed" or record.get("phase") != "accepted step checks" or "Velocity" not in record.get("reason", ""):
            raise ValueError("injected-NaN run failed for an unexpected reason or was accepted")
        if not list((run_dir / "output").glob("LiquidCore_ite_100000*.vtp")):
            raise ValueError("failure snapshot was not saved")
        print("Injected NaN was diagnosed, saved, and rejected with return code 1.")
        return 0
    if process.returncode != 0:
        print((run_dir / "run.log").read_text(encoding="utf-8")[-2000:], file=sys.stderr)
        return process.returncode
    verify(run_dir)
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("run", "fault-test"):
        run = commands.add_parser(name)
        run.add_argument("--executable", type=Path, required=True)
        run.add_argument("--config", type=Path, required=True)
        destination = run.add_mutually_exclusive_group(required=True)
        destination.add_argument("--run-dir", type=Path)
        destination.add_argument("--run-root", type=Path)
        run.add_argument("--no-states", action="store_true")
    check = commands.add_parser("verify")
    check.add_argument("--run-dir", type=Path, required=True)
    check.add_argument("--require-damage", action="store_true")
    args = parser.parse_args()
    try:
        if args.command == "verify":
            verify(args.run_dir.resolve(), args.require_damage)
            return 0
        return launch(args)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print("PERI-SPH run rejected:", error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

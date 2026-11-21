#!/usr/bin/env python3
"""Create isolated, reproducible PERI-SPH runs and verify their recorded state."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import uuid


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def source_snapshot_files(repo):
    """Collect authored inputs without walking generated solver output trees.

    Include the standalone Peri CMake entry and experiment reference profiles.
    Keep authored geometry/reload inputs; omit generated outputs and binaries.
    """
    repo = Path(repo).resolve()
    roots = [repo / name for name in ('src', 'cmake', 'peri', 'CMakeLists.txt')]
    roots += list((repo / 'tools').glob('peri_*.py'))
    roots += list((repo / 'tests/3d_examples').glob('extra_3d_PD_*'))
    excluded = {'CMakeFiles', '_deps', '__pycache__', '.git', '.venv'}
    generated_suffixes = {'.pyc', '.pyo', '.o', '.a', '.so', '.dll', '.exe', '.log'}
    files = set()
    for root in roots:
        if root.is_file():
            files.add(root)
        elif root.is_dir():
            def fail_walk(error):
                raise error
            for directory, subdirs, names in os.walk(root, onerror=fail_walk):
                subdirs[:] = [name for name in subdirs if name not in excluded]
                if root.parent == repo / 'tests/3d_examples' and Path(directory) == root:
                    subdirs[:] = [name for name in subdirs if name not in {'output', 'checkpoints', 'build'}]
                for name in names:
                    path = Path(directory) / name
                    if path.suffix.lower() not in generated_suffixes and path.is_file():
                        files.add(path)
    return sorted(files)


def parameters(path):
    result = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            key, value = line.split("=", 1)
            token = value.strip()
            if token in ("true", "false"):
                value = float(token == "true")
            else:
                try:
                    value = float(token)
                except ValueError:
                    value = token
            result[key.strip()] = value
    return result


def status(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    result = {"status": lines[0]}
    result.update(line.split("=", 1) for line in lines[1:] if "=" in line)
    return result


def require_completion_audit(run_dir, snapshots=False):
    """Warnings in a completed run cannot silently become convergence evidence."""
    from peri_audit import inspect
    audit = inspect(run_dir, snapshots=snapshots)
    (run_dir / "completion_audit.json").write_text(json.dumps(audit, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if not audit.get("simulation_completed") or audit["status"] != "sampled_checks_passed":
        reason = "; ".join(audit.get("warnings", [])) or audit.get("reason", audit["status"])
        raise ValueError("completion audit rejected run: " + reason)
    return audit


def verify(run_dir, require_damage=False, energy_tolerance=None):
    if energy_tolerance is not None and (not math.isfinite(energy_tolerance) or energy_tolerance <= 0):
        raise ValueError("energy tolerance must be finite and positive")
    output = run_dir / "output"
    record = status(output / "run_status.txt")
    if record["status"] != "completed":
        raise ValueError("simulation did not complete: " + record.get("reason", record["status"]))
    config = parameters(output / "effective_parameters.cfg")
    with (output / "diagnostics.csv").open(newline="", encoding="utf-8") as stream:
        rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]
    if len(rows) != int(config.get("diagnostic_frames", config["output_frames"])) + 1:
        raise ValueError("incomplete diagnostic output")
    if not all(math.isfinite(value) for row in rows for value in row.values()):
        raise ValueError("diagnostics contain non-finite data")
    expected_time = config["compression_distance"] / config["compression_speed"]
    if not math.isclose(rows[-1]["time_s"], expected_time, rel_tol=1e-12, abs_tol=1e-16):
        raise ValueError("end time does not match prescribed stroke and speed")
    if not math.isclose(rows[-1]["stroke_m"], config["compression_distance"], rel_tol=1e-12, abs_tol=1e-16):
        raise ValueError("final stroke is incorrect")
    if any(not math.isclose(row["stroke_m"], row["time_s"] * config["compression_speed"], rel_tol=1e-11, abs_tol=1e-16) for row in rows):
        raise ValueError("measured plate displacement does not match prescribed motion")
    if any(b["stroke_m"] <= a["stroke_m"] for a, b in zip(rows, rows[1:])):
        raise ValueError("measured plate stroke is not strictly increasing")
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
    # These budgets are independent of the sampled reaction CSV. The kinetic
    # identities use every actual solver kick; the mechanical residual remains
    # an assessment, because finite-strain hypoelasticity/fracture have no
    # complete thermodynamic potential in this model.
    def read_budget(name):
        with (output / name).open(newline="", encoding="utf-8") as stream:
            values = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]
        if len(values) != len(rows) or any(a["time_s"] != b["time_s"] for a, b in zip(values, rows)):
            raise ValueError(name + ": sample times do not match diagnostics")
        if not all(math.isfinite(x) for row in values for x in row.values()):
            raise ValueError(name + ": non-finite budget")
        return values
    mass_rows = read_budget("mass_budget.csv")
    energy_rows = read_budget("energy_budget.csv")
    inventory_error = max(abs(row[key]) for row in mass_rows for key in
        ("shell_mass_relative_error", "liquid_mass_relative_error", "max_particle_mass_relative_error"))
    if inventory_error > 1e-12:
        raise ValueError("particle inventory mass changed")
    if any(row["shell_count"] != mass_rows[0]["shell_count"] or row["liquid_count"] != mass_rows[0]["liquid_count"] for row in mass_rows):
        raise ValueError("particle inventory count changed")
    if any(abs(row["total_mass_kg"] - row["shell_mass_kg"] - row["liquid_mass_kg"]) > 1e-12 * row["total_mass_kg"] for row in mass_rows):
        raise ValueError("mass components do not reconcile")
    volume_mass_error = max(abs(row["density_volume_mass_residual_kg"]) / max(row["liquid_mass_kg"], 1e-300) for row in mass_rows)
    if volume_mass_error > 1e-12:
        raise ValueError("density-volume mass is inconsistent with particle mass")
    if any(row["liquid_mass_outside_domain_kg"] > 0 for row in mass_rows):
        raise ValueError("liquid particles left the computational domain")
    identity_error = max(abs(row[key]) for row in energy_rows for key in
        ("shell_kinetic_identity_residual_J", "fluid_kinetic_identity_residual_J"))
    kinetic_scale = max(max(abs(row[key]) for key in
        ("kinetic_J", "hourglass_signed_transfer_J", "pd_damping_signed_loss_J", "acoustic_inner_signed_loss_J", "acoustic_wall_signed_loss_J")) for row in energy_rows)
    if identity_error > max(1e-24, 1e-8 * kinetic_scale):
        raise ValueError("actual force-work does not reconcile with kinetic-energy change")
    max_closure = max(abs(row["mechanical_residual_relative"]) for row in energy_rows)
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
        "shell_particles": int(mass_rows[0]["shell_count"]),
        "liquid_particles": int(mass_rows[0]["liquid_count"]),
        "inventory_mass_relative_error": inventory_error,
        "shell_initial_mass_geometry_relative_error": mass_rows[0]["shell_mass_kg"] / (config["rho0_shell"] * 4 * math.pi / 3 * (config["capsule_radius_out"]**3 - (config["capsule_radius_out"]-config["shell_thickness"])**3)) - 1,
        "liquid_initial_mass_geometry_relative_error": mass_rows[0]["liquid_mass_kg"] / (config["rho0_f"] * 4 * math.pi / 3 * (config["capsule_radius_out"]-config["shell_thickness"])**3) - 1 if config["with_liquid"] else None,
        "density_volume_mass_relative_error": volume_mass_error,
        "kinetic_work_identity_max_absolute_J": identity_error,
        "mechanical_energy_residual_max_relative": max_closure,
        "mechanical_energy_residual_final_J": energy_rows[-1]["mechanical_residual_J"],
        "mechanical_energy_assessment": "within provisional 5% monitoring level" if max_closure <= 0.05 else "unclosed: exceeds provisional 5% monitoring level",
        "energy_model_limits": "shell elastic storage is a stress-compliance proxy; plastic work includes hardening and rate dependence; hourglass is signed transfer; fracture energy and finite-strain hypoelastic closure remain unresolved",
        "scope": "output integrity, mass inventory, actual kinetic-work identities and numerical limits; mechanical closure, convergence and calibration are reported separately",
    }
    (run_dir / "verification.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, allow_nan=False))
    if energy_tolerance is not None and max_closure > energy_tolerance:
        raise ValueError("mechanical energy residual exceeds requested tolerance")
    return result


def launch(args):
    executable = args.executable.resolve(strict=True)
    config = args.config.resolve(strict=True)
    run_dir = args.run_dir.resolve() if args.run_dir else args.run_root.resolve() / ("run-" + uuid.uuid4().hex)
    repo = Path(__file__).resolve().parents[1]
    files = source_snapshot_files(repo)
    target_case = repo / "tests/3d_examples" / executable.stem
    stale = [str(p.relative_to(repo)) for p in files if (p.is_relative_to(repo / "src/shared") or p.is_relative_to(repo / "src/for_3D_build") or p.is_relative_to(target_case)) and p.suffix in (".cpp", ".h", ".hpp", ".hxx") and p.stat().st_mtime > executable.stat().st_mtime + 0.001]
    if stale:
        raise ValueError("source is newer than the executable; rebuild before running: " + stale[0])
    # Do not create a seemingly dispatched run when source preflight already fails.
    run_dir.mkdir(parents=True, exist_ok=False)
    shutil.copy2(config, run_dir / "case.cfg")
    # Keep the exact executable even if the development build is replaced.
    frozen_executable = run_dir / executable.name
    shutil.copy2(executable, frozen_executable)
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
        "executable": str(executable), "frozen_executable": str(frozen_executable), "executable_sha256": digest(frozen_executable),
        "case_sha256": digest(config), "source_sha256": hashes,
        "state_recording": not args.no_states,
        "source_snapshot_scope": "authored source, Peri build/profile inputs and case assets; generated output/checkpoint/build trees excluded",
        "checkpoint_capabilities": checkpoint_capabilities(frozen_executable),
    }
    metadata["checkpoint_format"] = metadata["checkpoint_capabilities"].get("format", 0)
    return execute(args, run_dir, metadata, frozen_executable)


def checkpoint_capabilities(executable):
    # Old binaries reject this option. Probe inside an isolated directory so
    # their normal I/O initialization cannot touch any user output.
    with tempfile.TemporaryDirectory(prefix="peri-capabilities-") as directory:
        try:
            result = subprocess.run([str(executable), "--pd-checkpoint-capabilities"], cwd=directory,
                                    capture_output=True, text=True, timeout=15)
            value = json.loads(result.stdout)
            return value if result.returncode == 0 and value.get("format") == 1 else {}
        except (OSError, ValueError, subprocess.TimeoutExpired):
            return {}


def solver_alive(metadata):
    if sys.platform.startswith("linux") and metadata.get("linux_process_start_ticks"):
        try:
            fields = (Path("/proc")/str(metadata["solver_pid"])/"stat").read_text().rsplit(")",1)[1].split()
            return (fields[0] != "Z" and fields[19] == metadata["linux_process_start_ticks"] and
                    Path("/proc/sys/kernel/random/boot_id").read_text().strip() == metadata["linux_boot_id"])
        except FileNotFoundError:
            return False
    return False


def committed_checkpoints(directory):
    if (directory/"manifest.json").is_file():
        candidates=[directory]
    else:
        candidates=sorted(directory.glob("checkpoint-*"), reverse=True)
    result=[]
    for path in candidates:
        if path.is_symlink():continue
        try:
            data=json.loads((path/"manifest.json").read_text())
            if data["format"] == 1 and data["status"] == "complete" and (path/"state.bin").stat().st_size == data["state_bytes"]:
                result.append({"path": str(path), **data})
        except (OSError, ValueError, KeyError):
            continue
    return result


def resume(args):
    run_dir=args.run_dir.resolve(strict=True)
    metadata=json.loads((run_dir/"run_manifest.json").read_text())
    if solver_alive(metadata):
        raise ValueError("solver still exists; request a disk checkpoint pause before resuming")
    executable=run_dir/Path(metadata["frozen_executable"]).name
    if digest(executable)!=metadata["executable_sha256"] or digest(run_dir/"case.cfg")!=metadata["case_sha256"]:
        raise ValueError("sealed executable or configuration changed")
    if metadata.get("checkpoint_format")!=1:
        raise ValueError("this run has no supported full PD checkpoint capability")
    requested=args.checkpoint.resolve() if args.checkpoint else run_dir/"checkpoints"
    if not committed_checkpoints(requested):
        raise ValueError("no committed disk checkpoint available")
    args.no_states=not metadata["state_recording"]
    return execute(args,run_dir,metadata,executable,requested)


def execute(args, run_dir, metadata, frozen_executable, checkpoint=None):
    previous_controls=metadata.get("run_controls",{})
    if args.checkpoint_seconds is None:args.checkpoint_seconds=previous_controls.get("checkpoint_seconds",600)
    if args.checkpoint_steps is None:args.checkpoint_steps=previous_controls.get("checkpoint_steps",0)
    if not math.isfinite(args.checkpoint_seconds) or args.checkpoint_seconds < 0 or args.checkpoint_steps < 0 or (args.threads is not None and args.threads < 1):
        raise ValueError("invalid checkpoint interval or thread count")
    capable=metadata.get("checkpoint_format")==1
    if capable and args.threads is None:
        args.threads=metadata.get("run_controls",{}).get("threads") or metadata.get("checkpoint_capabilities",{}).get("default_threads")
    metadata["run_controls"]={"threads":args.threads,"checkpoint_seconds":args.checkpoint_seconds,"checkpoint_steps":args.checkpoint_steps}
    if not capable and (checkpoint is not None or args.checkpoint_steps or args.stop_after_checkpoint):
        raise ValueError("executable does not support full PD checkpoints")
    command=[str(frozen_executable),"--case-config","case.cfg","--state_recording="+("false" if args.no_states else "true")]
    if capable:
        command += ["--pd-checkpoint-seconds="+str(args.checkpoint_seconds),"--pd-checkpoint-steps="+str(args.checkpoint_steps)]
        if args.stop_after_checkpoint: command.append("--pd-stop-after-checkpoint")
        if args.threads is not None: command.append("--pd-threads="+str(args.threads))
        if checkpoint is not None: command.append("--pd-resume="+str(checkpoint))
    if checkpoint is not None and (run_dir/"output/resume_info.json").exists():
        old=run_dir/"output/resume_info.json"
        old.rename(old.with_name("resume_info-attempt-"+str(len(metadata.get("attempts",[])))+".json"))
    attempt={"started_at_utc":datetime.now(timezone.utc).isoformat(),"command":command,"resume_requested":str(checkpoint) if checkpoint else None}
    metadata.setdefault("attempts",[]).append(attempt)
    metadata.pop("live_audit_failure",None);metadata.pop("finished_at_utc",None);metadata.pop("return_code",None)
    metadata["status"]="restoring" if checkpoint else "running"
    def save_metadata():
        path=run_dir/"run_manifest.json";temporary=path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(metadata,indent=2)+"\n");temporary.replace(path)
    with (run_dir/"run.log").open("a" if checkpoint else "w",encoding="utf-8") as log:
        if checkpoint: log.write("\nRESUME ATTEMPT "+attempt["started_at_utc"]+"\n");log.flush()
        process=subprocess.Popen(command,cwd=run_dir,stdout=log,stderr=subprocess.STDOUT)
        metadata["solver_pid"]=process.pid
        if sys.platform.startswith("linux"):
            metadata["linux_process_start_ticks"]=(Path("/proc")/str(process.pid)/"stat").read_text().rsplit(")",1)[1].split()[19]
            metadata["linux_boot_id"]=Path("/proc/sys/kernel/random/boot_id").read_text().strip()
        attempt["solver_pid"]=process.pid;save_metadata()
        while True:
            try: process.wait(timeout=30);break
            except subprocess.TimeoutExpired:
                if checkpoint is not None and not (run_dir/"output/resume_info.json").exists():
                    live={"status":"restoring_checkpoint","solver_alive":True}
                elif (run_dir/"output/effective_parameters.cfg").exists():
                    from peri_audit import inspect
                    try: live=inspect(run_dir)
                    except (OSError,ValueError,KeyError) as error: live={"status":"failed","reason":str(error)}
                else: continue
                live["checked_at_utc"]=datetime.now(timezone.utc).isoformat()
                (run_dir/"live_audit.json").write_text(json.dumps(live,indent=2,allow_nan=False)+"\n")
                if live["status"]=="failed":
                    metadata["live_audit_failure"]=live["reason"];process.kill();process.wait();break
    attempt.update(finished_at_utc=datetime.now(timezone.utc).isoformat(),return_code=process.returncode)
    metadata.update(finished_at_utc=attempt["finished_at_utc"],return_code=process.returncode)
    if process.returncode==75:
        record=status(run_dir/"output/run_status.txt")
        checkpoints=committed_checkpoints(run_dir/"checkpoints")
        if record["status"]!="checkpointed" or not checkpoints:
            raise ValueError("checkpoint stop did not publish a complete state")
        metadata.update(status="paused",last_checkpoint=checkpoints[0]);save_metadata()
        from peri_audit import inspect
        live=inspect(run_dir);(run_dir/"live_audit.json").write_text(json.dumps(live,indent=2,allow_nan=False)+"\n")
        print(json.dumps({"status":"paused","run_dir":str(run_dir),"checkpoint":checkpoints[0]},indent=2))
        return 75
    metadata["status"]="completed" if process.returncode==0 else "failed";save_metadata()
    if args.command=="fault-test":
        record=status(run_dir/"output/run_status.txt")
        if process.returncode!=1 or record["status"]!="failed" or record.get("phase")!="accepted step checks" or "Velocity" not in record.get("reason",""):
            raise ValueError("injected-NaN run failed for an unexpected reason or was accepted")
        if not list((run_dir/"output").glob("LiquidCore_ite_100000*.vtp")):raise ValueError("failure snapshot was not saved")
        print("Injected NaN was diagnosed, saved, and rejected with return code 1.");return 0
    if process.returncode!=0:
        print((run_dir/"run.log").read_text()[-2000:],file=sys.stderr);return process.returncode
    verify(run_dir,energy_tolerance=args.energy_tolerance)
    final_audit=require_completion_audit(run_dir,snapshots=not args.no_states)
    (run_dir/"live_audit.json").write_text(json.dumps(final_audit,indent=2,allow_nan=False)+"\n")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("run", "fault-test"):
        run = commands.add_parser(name)
        run.add_argument("--energy-tolerance", type=float)
        run.add_argument("--executable", type=Path, required=True)
        run.add_argument("--config", type=Path, required=True)
        destination = run.add_mutually_exclusive_group(required=True)
        destination.add_argument("--run-dir", type=Path)
        destination.add_argument("--run-root", type=Path)
        run.add_argument("--no-states", action="store_true")
    resume_parser=commands.add_parser("resume")
    resume_parser.add_argument("--run-dir",type=Path,required=True)
    resume_parser.add_argument("--checkpoint",type=Path)
    resume_parser.add_argument("--energy-tolerance",type=float)
    for operational in (commands.choices["run"],commands.choices["fault-test"],resume_parser):
        operational.add_argument("--checkpoint-seconds",type=float)
        operational.add_argument("--checkpoint-steps",type=int)
        operational.add_argument("--stop-after-checkpoint",action="store_true")
        operational.add_argument("--threads",type=int)
    check = commands.add_parser("verify")
    check.add_argument("--run-dir", type=Path, required=True)
    check.add_argument("--require-damage", action="store_true")
    check.add_argument("--energy-tolerance", type=float)
    args = parser.parse_args()
    try:
        if args.command == "verify":
            verify(args.run_dir.resolve(), args.require_damage, args.energy_tolerance)
            return 0
        if args.command=="resume": return resume(args)
        return launch(args)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print("PERI-SPH run rejected:", error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

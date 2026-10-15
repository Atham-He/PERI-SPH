#!/usr/bin/env python3
"""Read running capsule budgets and stream VTPs without loading large XML trees."""
import argparse
from array import array
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import time
from xml.parsers import expat

from peri_run import parameters


def rows(path):
    """Ignore only a final uncommitted line while a writer is active."""
    if not path.exists():
        return []
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    if lines and not lines[-1].endswith("\n"):
        lines.pop()
    values = []
    for row in csv.DictReader(lines):
        if None in row or any(value is None for value in row.values()):
            raise ValueError(path.name + ": malformed committed row")
        value = {key: float(value) for key, value in row.items()}
        if not all(math.isfinite(v) for v in value.values()):
            raise ValueError(path.name + ": non-finite committed row")
        values.append(value)
    return values


def scan_vtp(path, require_full_precision=True):
    stack, fields, stored = [], {}, {}
    current = None
    point_count = None
    position_min, position_max = [math.inf]*3, [-math.inf]*3
    speed2, max_speed = 0., 0.

    def token(value):
        nonlocal speed2, max_speed
        x = float(value)
        if not math.isfinite(x):
            raise ValueError(path.name + ": non-finite " + current["name"])
        i = current["count"]
        current["count"] += 1
        current["min"] = min(current["min"], x)
        current["max"] = max(current["max"], x)
        if current["name"] in stored:
            stored[current["name"]].append(x)
        if current["name"] == "Position":
            axis = i % 3
            position_min[axis] = min(position_min[axis], x)
            position_max[axis] = max(position_max[axis], x)
        if current["name"] == "Velocity":
            speed2 += x*x
            if i % 3 == 2:
                max_speed = max(max_speed, math.sqrt(speed2)); speed2 = 0

    def text(data):
        if current is None:
            return
        value = current["tail"] + data
        parts = value.split()
        if value and not value[-1].isspace():
            current["tail"] = parts.pop() if parts else value
        else:
            current["tail"] = ""
        for part in parts:
            token(part)

    def start(name, attrs):
        nonlocal current, point_count
        parent = stack[-1] if stack else ""
        stack.append(name)
        if name == "Piece":
            point_count = int(attrs["NumberOfPoints"])
        if name != "DataArray":
            return
        if attrs.get("format", "ascii") != "ascii":
            raise ValueError("stream audit supports ASCII VTP only: " + path.name)
        label = attrs.get("Name", "unnamed")
        components = int(attrs.get("NumberOfComponents", "1"))
        expected = (point_count * components if parent in ("PointData", "Points") else
                    int(attrs.get("NumberOfTuples", "0"))*components if parent == "FieldData" else None)
        current = {"name": label, "type": attrs.get("type"), "count": 0, "expected": expected,
                   "min": math.inf, "max": -math.inf, "tail": "", "parent": parent}
        if label in ("Mass", "Density", "VolumetricMeasure", "TimeValue"):
            stored[label] = array("d")

    def end(name):
        nonlocal current
        if name == "DataArray":
            if current["tail"]:
                token(current["tail"])
            if current["expected"] is not None and current["count"] != current["expected"]:
                raise ValueError(path.name + ": wrong tuple count for " + current["name"])
            if current["count"]:
                fields[current["name"]] = {key: current[key] for key in ("type", "count", "min", "max", "parent")}
            current = None
        stack.pop()

    parser = expat.ParserCreate()
    parser.StartElementHandler = start; parser.EndElementHandler = end; parser.CharacterDataHandler = text
    with path.open("rb") as stream:
        while block := stream.read(65536):
            parser.Parse(block, False)
        parser.Parse(b"", True)
    if "Position" not in fields or fields["Position"]["parent"] != "Points":
        raise ValueError(path.name + ": missing particle coordinates")
    if "Density" in fields and fields["Density"]["min"] <= 0:
        raise ValueError(path.name + ": non-positive density")
    if "Mass" not in stored or not stored["Mass"] or min(stored["Mass"]) <= 0:
        raise ValueError(path.name + ": missing or non-positive particle Mass (possible output rounding)")
    if require_full_precision and any(v["type"] != "Float64" for v in fields.values()
                                     if v["type"].startswith("Float") and v["parent"] in ("Points", "PointData")):
        raise ValueError(path.name + ": floating fields are not Float64")
    if "NumericalFailureReason" in fields and fields["NumericalFailureReason"]["max"] != 0:
        raise ValueError(path.name + ": failed shell particles")
    result = {"file": str(path), "points": point_count, "fields_checked": len(fields),
              "mass_sum_kg": math.fsum(stored["Mass"]), "mass_min_kg": min(stored["Mass"]),
              "mass_max_kg": max(stored["Mass"]), "position_min_m": position_min,
              "position_max_m": position_max, "speed_max_m_s": max_speed,
              "time_s": stored["TimeValue"][0] if stored.get("TimeValue") else None}
    if "Density" in fields:
        result["density_min_kg_m3"] = fields["Density"]["min"]
        result["density_max_kg_m3"] = fields["Density"]["max"]
    if path.name.startswith("LiquidCore") and all(key in stored for key in ("Density", "VolumetricMeasure")):
        density_mass = math.fsum(r*v for r,v in zip(stored["Density"], stored["VolumetricMeasure"]))
        result["density_volume_mass_relative_error"] = abs(density_mass/result["mass_sum_kg"]-1)
        if result["density_volume_mass_relative_error"] > 1e-12:
            raise ValueError(path.name + ": density-volume mass inconsistent")
    return result


def acoustic_work_assessment(energy):
    """Grad-div components share one divergence; their signs are not passivity tests.

    Legacy wall work excludes prescribed rigid-wall input. Only a separately
    recorded boundary work permits testing the combined acoustic dissipation.
    Keep the raw signed work and every negative component visible for review.
    """
    observations, losses, flags = [], [], []
    boundary_recorded = all("acoustic_boundary_work_J" in row for row in energy)
    for row in energy:
        scale = max(abs(row["kinetic_J"]), abs(row["gravity_potential_change_J"]),
                    abs(sum(row[key] for key in ("plate_contact_work_J", "plate_fluid_pressure_work_J", "plate_fluid_viscous_work_J")) + row.get("acoustic_boundary_work_J",0)), 1e-30)
        tolerance = max(1e-24, 1e-3*scale)
        for key in ("acoustic_inner_signed_loss_J", "acoustic_wall_signed_loss_J"):
            if row[key] < -tolerance:
                observations.append({"time_s":row["time_s"],"term":key,"value_J":row[key]})
        net = row["acoustic_inner_signed_loss_J"] + row["acoustic_wall_signed_loss_J"]
        if boundary_recorded:
            loss = net + row["acoustic_boundary_work_J"]
            losses.append(loss)
            if loss < -tolerance:
                flags.append({"time_s":row["time_s"],"term":"acoustic_combined_dissipation_J","value_J":loss})
    return {"boundary_work_recorded":boundary_recorded,
            "dissipation_sign_checked":boundary_recorded,
            "signed_net_transfer_last_J":energy[-1]["acoustic_inner_signed_loss_J"]+energy[-1]["acoustic_wall_signed_loss_J"],
            "combined_dissipation_last_J":losses[-1] if losses else None,
            "negative_component_observations":observations,
            "interpretation":"legacy acoustic columns are signed work transfers; separate component signs do not establish energy injection; prescribed boundary work is required for a combined dissipation sign test"}, flags


def inspect(run_dir, snapshots=False, stale_seconds=180):
    output = run_dir / "output"
    config = parameters(output / "effective_parameters.cfg")
    all_rows = [rows(output / name) for name in ("diagnostics.csv", "mass_budget.csv", "energy_budget.csv")]
    times = set.intersection(*(set(row["time_s"] for row in values) for values in all_rows))
    if not times:
        return {"status": "initializing", "scope": "no complete diagnostic sample yet"}
    diagnostic, mass, energy = [[row for row in values if row["time_s"] in times] for values in all_rows]
    end_time = config["compression_distance"] / config["compression_speed"]
    progress = rows(output / "progress.csv")
    accepted_time = max(diagnostic[-1]["time_s"], progress[-1]["time_s"] if progress else 0)
    warnings = []
    for row in mass:
        for key in ("shell_count", "liquid_count"):
            if row[key] != mass[0][key]:
                raise ValueError("particle inventory count changed")
        for key in ("shell_mass_relative_error", "liquid_mass_relative_error", "max_particle_mass_relative_error"):
            if abs(row[key]) > 1e-12:
                raise ValueError("particle inventory mass changed")
        if row["liquid_mass_outside_domain_kg"] != 0:
            raise ValueError("liquid mass outside computational bounds")
        if abs(row["density_volume_mass_residual_kg"]) > max(1e-300, row["liquid_mass_kg"]*1e-12):
            raise ValueError("density-volume mass inconsistent")
    for row in diagnostic:
        if row["numerical_failures"] != 0 or not 0 <= row["damage_mean"] <= row["damage_max"] <= 1:
            raise ValueError("invalid shell state")
        if row["fluid_density_error_max"] > config["max_density_error"] or row["fluid_speed_max_m_s"] > config["max_fluid_mach"]*config["c_f"]:
            raise ValueError("fluid numerical limit exceeded")
    scale = max(max(abs(row[key]) for key in ("kinetic_J", "hourglass_signed_transfer_J", "pd_damping_signed_loss_J", "acoustic_inner_signed_loss_J", "acoustic_wall_signed_loss_J")) for row in energy)
    identity = max(abs(row[key]) for row in energy for key in ("shell_kinetic_identity_residual_J", "fluid_kinetic_identity_residual_J"))
    if identity > max(1e-24,1e-8*scale):
        raise ValueError("kinetic work identity failed")
    closure = max(abs(row["mechanical_residual_relative"]) for row in energy)
    if closure > 0.05:
        warnings.append("mechanical residual exceeds provisional 5% monitoring level")
    acoustic_assessment, negative_losses = acoustic_work_assessment(energy)
    for row in energy:
        row_scale = max(abs(row["kinetic_J"]),abs(row["gravity_potential_change_J"]),
                        abs(row["plate_contact_work_J"]+row["plate_fluid_pressure_work_J"]+row["plate_fluid_viscous_work_J"]+row.get("acoustic_boundary_work_J",0)),1e-30)
        for key in ("pd_damping_signed_loss_J","viscous_signed_loss_J","friction_signed_loss_J"):
            if row[key] < -max(1e-24,1e-3*row_scale):
                negative_losses.append({"time_s":row["time_s"],"term":key,"value_J":row[key]})
    if negative_losses:
        warnings.append("signed loss is negative beyond 0.1% of current mechanical scale; investigate energy injection")
    completed = (output/"run_status.txt").exists() and (output/"run_status.txt").read_text().splitlines()[0] == "completed"
    age = max(0,time.time()-(output/"progress.csv").stat().st_mtime) if progress else None
    if age is not None and age > stale_seconds and not completed:
        warnings.append("progress timestamp is stale; process liveness needs checking")
    alive = None
    manifest_file = run_dir/"run_manifest.json"
    if manifest_file.exists():
        manifest = json.loads(manifest_file.read_text())
        if "linux_process_start_ticks" in manifest and Path("/proc/sys/kernel/random/boot_id").exists():
            try:
                stat = Path("/proc")/str(manifest["solver_pid"])/"stat"
                alive = (stat.read_text().rsplit(")",1)[1].split()[19] == manifest["linux_process_start_ticks"] and
                         Path("/proc/sys/kernel/random/boot_id").read_text().strip() == manifest["linux_boot_id"])
            except FileNotFoundError:
                alive = False
            if not alive and not completed:
                warnings.append("recorded solver process is no longer alive")
    eta = None
    if len(progress) > 1:
        previous = progress[max(0,len(progress)-21)]
        wall_interval=progress[-1]["wall_s"]-previous["wall_s"]
        physical_interval=progress[-1]["time_s"]-previous["time_s"]
        if wall_interval > 0 and physical_interval > 0:
            eta=(end_time-accepted_time)*wall_interval/physical_interval
    result = {"status":"needs_attention" if warnings else "sampled_checks_passed",
              "checked_at_utc":datetime.now(timezone.utc).isoformat(),"simulation_completed":completed,
              "accepted_time_s":accepted_time,"completion_fraction":accepted_time/end_time,
              "budget_time_s":diagnostic[-1]["time_s"],"diagnostic_samples_checked":len(diagnostic),
              "solver_alive":alive,"estimated_remaining_wall_s_at_recent_rate":eta,
              "progress_age_s":age,"shell_mass_kg":mass[-1]["shell_mass_kg"],"liquid_mass_kg":mass[-1]["liquid_mass_kg"],
              "inventory_mass_relative_error":max(abs(row[key]) for row in mass for key in ("shell_mass_relative_error","liquid_mass_relative_error","max_particle_mass_relative_error")),
              "density_error_max":max(row["fluid_density_error_max"] for row in diagnostic),
              "mechanical_residual_max_relative":closure,"mechanical_residual_last_J":energy[-1]["mechanical_residual_J"],
              "kinetic_identity_max_absolute_J":identity,"negative_loss_flags":negative_losses,"warnings":warnings,
              "acoustic_work_assessment":acoustic_assessment,
              "scope":"sampled running diagnostics only; does not certify full-stroke completion, convergence, shell impermeability or calibrated fracture-energy closure"}
    if snapshots:
        body_specs=[("CapsuleShell","shell_mass_kg","shell_count")]
        if config["with_liquid"]: body_specs.append(("LiquidCore","liquid_mass_kg","liquid_count"))
        frame_sets = [{p.name.split("_ite_")[1].split(".")[0] for p in output.glob(body+"_ite_*.vtp")} for body,_,_ in body_specs]
        common = set.intersection(*frame_sets)
        if not common:
            raise ValueError("no common shell/liquid snapshot frame")
        if completed:
            expected={f"{i:010d}" for i in range(int(config["output_frames"])+1)}
            if any(frames != expected for frames in frame_sets):
                raise ValueError("snapshot cadence differs from output_frames")
        frame = sorted(common,key=int)[-1]
        summaries=[]
        for body,key,count_key in body_specs:
            path=output/(body+"_ite_"+frame+".vtp")
            with path.open("rb") as stream:
                stream.seek(max(0,path.stat().st_size-256)); tail=stream.read()
            if b"</VTKFile>" not in tail:
                result["snapshot_status"]="latest frame is still being written"; return result
            summary=scan_vtp(path)
            if summary["points"] != mass[0][count_key] or not math.isclose(summary["mass_sum_kg"],mass[0][key],rel_tol=1e-10):
                raise ValueError(body+": snapshot particle inventory differs from CSV")
            summaries.append(summary)
        result["snapshot_status"]="checked";result["snapshots"]=summaries
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir",type=Path,required=True)
    parser.add_argument("--snapshots",action="store_true")
    parser.add_argument("--output",type=Path,help="write a separate review without overwriting the sealed monitor's live audit")
    parser.add_argument("--strict",action="store_true",help="fail on monitoring warnings, for controlled regression cases")
    parser.add_argument("--stale-seconds",type=float,default=180)
    args=parser.parse_args()
    try:
        result=inspect(args.run_dir.resolve(),args.snapshots,args.stale_seconds)
    except (OSError,ValueError,KeyError,expat.ExpatError) as exc:
        result={"status":"failed","checked_at_utc":datetime.now(timezone.utc).isoformat(),"reason":str(exc)}
    destination=args.output if args.output is not None else args.run_dir/"live_audit.json"
    destination.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,allow_nan=False))
    return 1 if result["status"]=="failed" or (args.strict and result["status"]=="needs_attention") else 0


if __name__=="__main__":
    raise SystemExit(main())

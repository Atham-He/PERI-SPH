#!/usr/bin/env python3
"""Compare uninterrupted runs against fresh-process disk restores, including SIGKILL."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sys
import subprocess
import time
import uuid

FILES=("diagnostics.csv","mass_budget.csv","energy_budget.csv","UpperPlate_TotalReaction.dat")

def latest(run):
    return run/"checkpoints"/(run/"checkpoints/LATEST").read_text().strip()

def last_rows(run,name):
    with (run/"output"/name).open() as stream:return list(csv.DictReader(stream))

def command(executable,config,steps=0,stop=False,resume=False):
    args=[str(executable),"--case-config",str(config),"--state_recording=true","--pd-threads=1", "--pd-checkpoint-seconds=0", "--pd-checkpoint-steps="+str(steps)]
    if stop:args.append("--pd-stop-after-checkpoint")
    if resume:args.append("--pd-resume=checkpoints")
    return args

def run_command(args,directory,expected=0):
    directory.mkdir(parents=True,exist_ok=True)
    with (directory/"test.log").open("a") as log:
        log.write("COMMAND "+repr(args)+"\n");log.flush()
        result=subprocess.run(args,cwd=directory,stdout=log,stderr=subprocess.STDOUT,timeout=180)
    if result.returncode!=expected:raise ValueError(f"{directory}: expected {expected}, got {result.returncode}: "+(directory/"test.log").read_text()[-3000:])

def compare(reference,resumed):
    flags={name:(reference/"output"/name).read_bytes()==(resumed/"output"/name).read_bytes() for name in FILES}
    a=json.loads((latest(reference)/"manifest.json").read_text());b=json.loads((latest(resumed)/"manifest.json").read_text())
    flags["complete_physical_state_digest"]=a["physical_state_fnv1a64"]==b["physical_state_fnv1a64"]
    flags["step_counts"]=(a["acoustic_steps"],a["pd_steps"])==(b["acoustic_steps"],b["pd_steps"])
    if not all(flags.values()):raise ValueError("disk restore differs from uninterrupted run: "+json.dumps({"run":str(resumed),"checks":flags,"reference":a,"resumed":b}))
    return flags


def control_roundtrip(executable,config,root):
    """Exercise the real study/runner/signal control chain with a small fixture."""
    study=root/"control-study";frozen=study/"frozen_source";(frozen/"tools").mkdir(parents=True);(frozen/"bin").mkdir()
    repo=Path(__file__).resolve().parents[1]
    for script in (repo/"tools").glob("peri_*.py"):shutil.copy2(script,frozen/"tools"/script.name)
    shutil.copy2(executable,frozen/"bin"/executable.name);shutil.copy2(config,study/"baseline.cfg")
    sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    cases=[{"name":"baseline","group":"baseline","config":"baseline.cfg","config_sha256":sha(study/"baseline.cfg"),"shell_layers":3},
           {"name":"memory_guard","group":"time","config":"baseline.cfg","config_sha256":sha(study/"baseline.cfg"),"shell_layers":100}]
    data={"status":"planned","source_git_head":"control-regression-fixture","executable":"frozen_source/bin/"+executable.name,
          "source_sha256":{str(p.relative_to(frozen)):sha(p) for p in frozen.rglob("*") if p.is_file()},"cases":cases}
    (study/"study.json").write_text(json.dumps(data))
    with (study/"test.log").open("w") as log:
        process=subprocess.Popen([sys.executable,str(frozen/"tools/peri_study.py"),"run","--root",str(study),"--memory-limit-gib","1","--baseline-memory-gib","0.1"],stdout=log,stderr=subprocess.STDOUT)
        try:
            deadline=time.monotonic()+60;run=study/"runs/baseline"
            while time.monotonic()<deadline:
                if process.poll() is not None:raise ValueError("control fixture ended before pause")
                if (run/"run.log").exists() and "frame=1 " in (run/"run.log").read_text():break
                time.sleep(0.02)
            else:raise ValueError("control fixture never reached a nonzero time")
            control=[sys.executable,str(frozen/"tools/peri_control.py")]
            paused=subprocess.run(control+["pause","--study",str(study),"--timeout","60"],capture_output=True,text=True,timeout=65)
            if paused.returncode:raise ValueError("disk pause control failed: "+paused.stderr)
            pause_info=json.loads(paused.stdout);assert process.wait(timeout=10)==75
            saved_time=pause_info["checkpoint"]["time_s"];assert saved_time>0
            resumed=subprocess.run(control+["resume","--study",str(study),"--background"],capture_output=True,text=True,timeout=15)
            if resumed.returncode:raise ValueError("disk resume control failed: "+resumed.stderr)
            deadline=time.monotonic()+60
            while time.monotonic()<deadline:
                state=json.loads((study/"study.json").read_text())
                if state["status"]=="resource_limited":break
                if state["status"] in ("simulation_failed","audit_failed","resume_rejected"):raise ValueError("control resume rejected: "+json.dumps(state))
                time.sleep(0.05)
            else:raise ValueError("resumed study did not finish its small baseline")
            manifest=json.loads((run/"run_manifest.json").read_text());loaded=json.loads((run/"output/resume_info.json").read_text())
            if manifest["status"]!="completed" or len(manifest["attempts"])!=2 or loaded["time_s"]!=saved_time:
                raise ValueError("control did not resume the same nonzero disk state")
            return {"pause_status":pause_info["status"],"saved_time_s":saved_time,"restored_time_s":loaded["time_s"],
                    "new_process":manifest["attempts"][0]["solver_pid"]!=manifest["attempts"][1]["solver_pid"],"baseline_status":manifest["status"],"post_baseline_resource_gate":state["status"]}
        finally:
            if process.poll() is None:process.kill();process.wait()

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--executable",type=Path,required=True)
    parser.add_argument("--configs",type=Path,required=True);parser.add_argument("--root",type=Path,required=True)
    args=parser.parse_args();base=args.root.resolve();base.mkdir(parents=True,exist_ok=True);root=base/("run-"+uuid.uuid4().hex);root.mkdir()
    executable=root/args.executable.name;shutil.copy2(args.executable,executable)
    report={"executable_sha256":hashlib.sha256(executable.read_bytes()).hexdigest(),"cases":[],"scope":"small 3D verification; all state digest and scientific CSV byte comparisons; no full-scale convergence claim"}
    for case in ("smoke_liquid","smoke_fracture","smoke_checkpoint_plastic"):
        config=root/(case+".cfg")
        config.write_text((args.configs/(case+".cfg")).read_text().replace("capsule_radius_out=8e-5","capsule_radius_out=5e-5").replace("plate_radius=1.2e-4","plate_radius=8e-5").replace("compression_distance=2e-6","compression_distance=5e-7"))
        reference=root/(case+"-continuous")
        run_command(command(executable,config),reference)
        full=json.loads((latest(reference)/"manifest.json").read_text());split=root/(case+"-resumed")
        run_command(command(executable,config,max(1,full["acoustic_steps"]//2),True),split,75)
        saved=json.loads((latest(split)/"manifest.json").read_text())
        saved_damage=max(float(row["damage_max"]) for row in last_rows(split,"diagnostics.csv"))
        saved_plastic=max(float(row["plastic_work_J"]) for row in last_rows(split,"energy_budget.csv"))
        if case=="smoke_fracture" and saved_damage<=0:raise ValueError("saved fracture state had no broken bonds")
        if case=="smoke_checkpoint_plastic" and saved_plastic<=0:raise ValueError("saved plastic state had no plastic history")
        run_command(command(executable,config,stop=True,resume=True),split,75)
        roundtrip=json.loads((latest(split)/"manifest.json").read_text())
        if any(saved[key]!=roundtrip[key] for key in ("physical_state_fnv1a64","time_s","acoustic_steps","pd_steps","completed_frame")):
            raise ValueError("restored memory differs before the next integration step")
        for name in FILES:
            with (split/"output"/name).open("a") as stream:stream.write("999,uncommitted-tail\n")
        sentinel=split/"output/CapsuleShell_ite_1000001234.vtp";sentinel.write_text("orphaned future output")
        run_command(command(executable,config,resume=True),split)
        flags=compare(reference,split)
        if sentinel.exists() or not list((split/"output").glob("orphaned-after-resume-*/CapsuleShell_ite_1000001234.vtp")):
            raise ValueError("uncommitted future snapshots were not preserved out of the active output set")
        damage=max(row["damage_max"] for row in [{k:float(v) for k,v in x.items()} for x in last_rows(reference,"diagnostics.csv")])
        plastic=max(float(row["plastic_work_J"]) for row in last_rows(reference,"energy_budget.csv"))
        if case=="smoke_fracture" and damage<=0:raise ValueError("brittle checkpoint test did not break bonds")
        if case=="smoke_checkpoint_plastic" and plastic<=0:raise ValueError("plastic checkpoint test did not produce plastic history")
        report["cases"].append({"case":case,"saved_time_s":saved["time_s"],"checks":flags,"restore_before_advancing_digest_identical":True,"damage_max":damage,"plastic_work_J":plastic,"saved_damage_max":saved_damage,"saved_plastic_work_J":saved_plastic})
        if case=="smoke_liquid":
            crash=root/"liquid-crash-and-corrupt-latest";crash.mkdir()
            with (crash/"test.log").open("w") as log:
                process=subprocess.Popen(command(executable,config,max(1,full["acoustic_steps"]//3)),cwd=crash,stdout=log,stderr=subprocess.STDOUT)
                deadline=time.monotonic()+180
                while time.monotonic()<deadline:
                    if process.poll() is not None:raise ValueError("crash scenario finished before a nonzero checkpoint was captured")
                    try:
                        if json.loads((latest(crash)/"manifest.json").read_text())["time_s"]>0:break
                    except (OSError,ValueError):pass
                    time.sleep(0.005)
                else:process.kill();raise ValueError("nonzero checkpoint was not published")
                process.kill();process.wait();assert process.returncode==-9
            corrupt=latest(crash)/"state.bin"
            with corrupt.open("r+b") as stream:stream.truncate(corrupt.stat().st_size-17)
            run_command(command(executable,config,resume=True),crash)
            report["crash_corrupt_latest_fallback"]={"checks":compare(reference,crash),"killed_return_code":-9}
    report["operational_study_pause_resume"]=control_roundtrip(executable,root/"smoke_liquid.cfg",root)
    report["run_root"]=str(root)
    (root/"comparison.json").write_text(json.dumps(report,indent=2)+"\n")
    (base/"latest_comparison.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2));return 0

if __name__=="__main__":raise SystemExit(main())

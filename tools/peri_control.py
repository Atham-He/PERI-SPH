#!/usr/bin/env python3
"""Disk checkpoint control. pause exits the integrator cleanly and releases RAM."""
import argparse
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
from peri_run import committed_checkpoints,digest,solver_alive


def paths(root):
    data=json.loads((root/"study.json").read_text())
    name=data.get("active_case","baseline")
    return data,root/"runs"/name


def live_identity(run):
    metadata=json.loads((run/"run_manifest.json").read_text())
    if metadata.get("checkpoint_format")!=1:
        raise ValueError("legacy executable has no full PD checkpoint support; no signal was sent")
    if not solver_alive(metadata):raise ValueError("recorded solver is not alive")
    executable=run/Path(metadata["frozen_executable"]).name
    if digest(executable)!=metadata["executable_sha256"]:raise ValueError("sealed binary changed")
    actual=(Path("/proc")/str(metadata["solver_pid"])/"exe").resolve()
    if actual!=executable.resolve():raise ValueError("solver executable identity changed")
    return metadata


def request(root,pause,timeout):
    _,run=paths(root);metadata=live_identity(run)
    latest=run/"checkpoints/LATEST"
    before=latest.read_text() if latest.exists() else None
    os.kill(metadata["solver_pid"],signal.SIGUSR2 if pause else signal.SIGUSR1)
    deadline=time.monotonic()+timeout
    while time.monotonic()<deadline:
        changed=latest.exists() and latest.read_text()!=before
        candidates=committed_checkpoints(run/"checkpoints") if changed or pause else []
        if candidates:
            if not pause:
                return {"status":"checkpoint_committed","checkpoint":candidates[0],"solver_continues":True}
            study=json.loads((root/"study.json").read_text())
            state=(run/"output/run_status.txt").read_text().splitlines() if (run/"output/run_status.txt").exists() else []
            if state and state[0]=="checkpointed" and not solver_alive(metadata) and study["status"]=="paused":
                return {"status":"paused_on_disk","checkpoint":candidates[0],"safe_to_exit_processes":True,
                        "resume_command":[sys.executable,str(Path(__file__).resolve()),"resume","--study",str(root),"--background"]}
        if not solver_alive(metadata) and not changed:
            state=(run/"output/run_status.txt").read_text().splitlines() if (run/"output/run_status.txt").exists() else []
            if not (pause and state and state[0]=="checkpointed"):
                raise ValueError("solver exited before a new checkpoint was committed; previous generations remain available")
        time.sleep(0.25)
    raise ValueError("checkpoint request still pending; process was not killed; inspect run.log before shutting down")


def restart_study(root,background):
    data,run=paths(root)
    manifest=json.loads((run/"run_manifest.json").read_text())
    if solver_alive(manifest):raise ValueError("solver still exists; request a disk pause first")
    if manifest.get("checkpoint_format")!=1 or not committed_checkpoints(run/"checkpoints"):
        raise ValueError("no committed full PD checkpoint; refusing a silent restart from zero")
    controls=data.get("run_controls",{})
    command=[sys.executable,str(root/"frozen_source/tools/peri_study.py"),"run","--root",str(root),
             "--memory-limit-gib",str(controls.get("memory_limit_gib",15)),
             "--baseline-memory-gib",str(controls.get("baseline_memory_gib",13))]
    if controls.get("no_states"):command.append("--no-states")
    if not background:return subprocess.run(command).returncode
    with (root/"launcher.log").open("a") as log:
        process=subprocess.Popen(["nice","-n","10"]+command,cwd=root,stdin=subprocess.DEVNULL,
                                 stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    time.sleep(0.1)
    if process.poll() is not None:raise ValueError("resume controller exited; inspect launcher.log")
    stat=(Path("/proc")/str(process.pid)/"stat").read_text().rsplit(")",1)[1].split()
    record={"kind":"sealed_study_driver","pid":process.pid,"entry_script":command[1],
            "source_git_head":data["source_git_head"],"boot_id":Path("/proc/sys/kernel/random/boot_id").read_text().strip(),
            "identity":{"pid":process.pid,"start_ticks":stat[19],"cmdline":(Path("/proc")/str(process.pid)/"cmdline").read_bytes().replace(b"\0",b" ").decode().strip()},
            "resumed_at_utc":datetime.now(timezone.utc).isoformat()}
    (root/"launcher.pid").write_text(str(process.pid)+"\n")
    (root/"current_controller.json").write_text(json.dumps(record,indent=2)+"\n")
    return {"status":"disk_resume_launched","controller_pid":process.pid,"study":str(root)}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("action",choices=["status","checkpoint","pause","resume"])
    parser.add_argument("--study",type=Path,required=True);parser.add_argument("--timeout",type=float,default=900)
    parser.add_argument("--background",action="store_true");args=parser.parse_args();root=args.study.resolve(strict=True)
    try:
        if args.action=="status":
            data,run=paths(root);manifest=json.loads((run/"run_manifest.json").read_text())
            result={"study_status":data["status"],"active_case":data.get("active_case"),"solver_alive":solver_alive(manifest),
                    "checkpoint_format":manifest.get("checkpoint_format",0),"committed_checkpoints":committed_checkpoints(run/"checkpoints")}
        elif args.action=="resume":result=restart_study(root,args.background)
        else:result=request(root,args.action=="pause",args.timeout)
        if isinstance(result,int):return result
        print(json.dumps(result,indent=2));return 0
    except (OSError,ValueError,KeyError) as error:
        print("PERI-SPH control rejected: "+str(error),file=sys.stderr);return 1

if __name__=="__main__":sys.exit(main())

#!/usr/bin/env python3
"""Isolated, resumable cantilever runs. No third-party Python modules required.

Each run is a separate process/directory. A failed or warned run is retained and
never silently accepted. Existing results are reused only with matching inputs,
executable SHA256 and a successful completion marker.
"""
import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def cases(suite):
    base = dict(cells=4, length_ratio=12, scale=1, cycles=2, cfl=0.2,
                coefficient=0.5, horizon_ratio=3.015, amplitude_ratio=0.001,
                stabilization=2, load_mode=0, samples_per_cycle=160, threads=2,
                numerical_damping=True)
    result = {}
    def add(name, **changes):
        result[name] = dict(base, **changes)
    if suite == "pilot":
        add("pilot_new_n4", cycles=0.25)
        return result
    if suite == "performance":
        for threads in (1, 4, 8):
            add("pilot_n6_threads" + str(threads), cells=6, cycles=0.125, threads=threads)
        return result
    add("new_n4")
    add("legacy_n4", stabilization=1)
    add("off_n4", stabilization=0)
    if suite == "comparison":
        return result
    add("new_n3", cells=3)
    add("new_n6", cells=6)
    add("new_n8", cells=8)
    for n in (3, 6, 8):
        add("legacy_n" + str(n), cells=n, stabilization=1)
    add("new_dt010", cfl=0.1)
    add("new_dt005", cfl=0.05)
    add("new_small_amplitude", amplitude_ratio=0.0005)
    for method, ident in (("new", 2), ("legacy", 1)):
        add(method + "_scale0001", stabilization=ident, scale=0.001)
        add(method + "_scale10", stabilization=ident, scale=10)
    # Fixed physical horizon: delta/H=3.015/4; changes quadrature, not delta.
    add("new_fixed_delta_n3", cells=3, horizon_ratio=3.015 * 3 / 4)
    add("new_fixed_delta_n6", cells=6, horizon_ratio=3.015 * 6 / 4)
    for n in (4, 6, 8):
        add("static_new_n" + str(n), cells=n, load_mode=1, cycles=4)
    add("static_legacy_n4", stabilization=1, load_mode=1, cycles=4)
    # The static result is admissible only if the analyzer's equilibrium gates pass.
    add("static_legacy_n8", cells=8, stabilization=1, load_mode=1, cycles=2)
    return result


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--suite", choices=("pilot", "performance", "comparison", "full"), default="full")
    parser.add_argument("--only", nargs="*")
    parser.add_argument("--numerical-damping", action="store_true", help="enable the production pair damper for every selected case")
    parser.add_argument("--no-numerical-damping", action="store_true", help="explicit undamped diagnostic; may fail or grow spuriously")
    parser.add_argument("--threads", type=int, help="TBB threads per case; effective value is recorded")
    parser.add_argument("--jobs", type=int, default=1, help="independent simulation processes (default 1)")
    parser.add_argument("--timeout", type=float, default=7200)
    args = parser.parse_args()
    if args.jobs < 1 or (args.threads is not None and args.threads < 1):
        parser.error("threads and jobs must be positive")
    executable = args.executable.resolve()
    root = Path(__file__).resolve().parent.parent
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    selected = cases(args.suite)
    if args.numerical_damping:
        selected = {name: dict(parameters, numerical_damping=True) for name, parameters in selected.items()}
    if args.no_numerical_damping:
        if args.numerical_damping:
            parser.error("conflicting damping options")
        selected = {name: dict(parameters, numerical_damping=False) for name, parameters in selected.items()}
    if args.only:
        unknown = set(args.only) - set(selected)
        if unknown:
            parser.error("unknown cases: " + ", ".join(sorted(unknown)))
        selected = {name: selected[name] for name in args.only}
    if args.threads is not None:
        selected = {name: dict(parameters, threads=args.threads) for name, parameters in selected.items()}
    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True, check=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--short"], cwd=root, capture_output=True, text=True, check=True).stdout
    sources = list((root / "src/shared/peridynamics").glob("*.[ch]*"))
    sources += list((root / "tests/3d_examples/extra_3d_PD_cantilever").glob("*"))
    sources += [Path(__file__).resolve(), root / "tools/analyze_cantilever.py",
                root / "tools/test_cantilever_analysis.py"]
    provenance = dict(executable=str(executable), executable_sha256=digest(executable),
                      git_head=revision, git_status=dirty,
                      sources={str(p.relative_to(root)): digest(p) for p in sources if p.is_file()})
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    for source in sources:
        if source.is_file():
            frozen = output / "source_snapshot" / source.relative_to(root)
            frozen.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, frozen)
    patch = subprocess.run(["git", "diff", "--binary"], cwd=root, capture_output=True, check=True).stdout
    (output / "tracked_changes.patch").write_bytes(patch)
    manifest = {p.parent.name: json.loads(p.read_text()) for p in output.glob("*/result.json")}
    def run_case(name, parameters):
        directory = output / name
        expected = dict(parameters=parameters, executable_sha256=provenance["executable_sha256"])
        if directory.exists():
            previous = json.loads((directory / "input.json").read_text())
            if previous != expected:
                raise RuntimeError("refusing to mix changed executable or parameters in " + str(directory))
            if (directory / "result.json").is_file():
                cached = json.loads((directory / "result.json").read_text())
                print(name + ": retained " + cached["status"], flush=True)
                return name, cached
            raise RuntimeError("incomplete prior run; inspect and choose a fresh output directory: " + str(directory))
        directory.mkdir()
        (directory / "input.json").write_text(json.dumps(expected, indent=2) + "\n")
        config = directory / "case.cfg"
        config.write_text("".join(f"{k}={str(v).lower() if isinstance(v, bool) else v}\n" for k, v in parameters.items()))
        command = [str(executable), "--case-config=" + str(config), "--state_recording=false"]
        print(name + ": starting", flush=True)
        start = time.monotonic()
        with (directory / "run.log").open("w") as log:
            try:
                run = subprocess.run(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT, timeout=args.timeout)
                code = run.returncode
            except subprocess.TimeoutExpired:
                code = -999
        text = (directory / "run.log").read_text(errors="replace")
        warnings = [line for line in text.splitlines() if any(token in line.lower() for token in ("warning", "failed", "numerical failure"))]
        completed = code == 0 and (directory / "completion.json").is_file() and not warnings
        result = dict(status="completed" if completed else "failed", returncode=code,
                      elapsed_seconds=time.monotonic() - start, warnings=warnings, command=command)
        (directory / "result.json").write_text(json.dumps(result, indent=2) + "\n")
        print(name + ": " + result["status"] + f" ({result['elapsed_seconds']:.1f}s)", flush=True)
        return name, result
    # Only the coordinator writes the aggregate manifest; each worker owns one directory.
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending = [pool.submit(run_case, name, parameters) for name, parameters in selected.items()]
        for future in as_completed(pending):
            name, result = future.result()
            manifest[name] = result
            (output / ("manifest_" + args.suite + ".json")).write_text(json.dumps(manifest, indent=2) + "\n")
    return 0 if all(manifest[name]["status"] == "completed" for name in selected) else 1


if __name__ == "__main__":
    raise SystemExit(main())

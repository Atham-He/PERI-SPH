#!/usr/bin/env python3
"""Analyze peri_cantilever_study outputs; requires numpy, scipy, matplotlib.

Uses fitted fundamental frequency, full-history error, static equilibrium checks,
and pairwise refinements. A simulation completing is not a validation pass.
"""
import argparse
import csv
import json
import hashlib
import platform
import shutil
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar, least_squares


def table(path):
    return np.genfromtxt(path, delimiter=",", names=True, ndmin=1)


def frequency(t, y):
    # t is normalized by the analytical period: fitted f is omega_PD/omega_EB.
    if (len(t) < 4 or np.any(np.diff(t) <= 0) or not np.all(np.isfinite(y)) or
            np.linalg.norm(y - np.mean(y)) <= 1e-14 * max(np.linalg.norm(y), 1e-100)):
        return {"fit_valid": False}
    def fit(f, decay=0, residual_vector=False):
        envelope = np.exp(-decay * t)
        a = np.column_stack((envelope*np.sin(2 * np.pi * f * t), envelope*np.cos(2 * np.pi * f * t), np.ones_like(t)))
        coefficients = np.linalg.lstsq(a, y, rcond=None)[0]
        if residual_vector:
            return (y - a @ coefficients) / np.sqrt(np.mean(y * y))
        return np.mean((y - a @ coefficients) ** 2), coefficients
    # Wide enough to diagnose the legacy operator's scale-dependent stiffness.
    upper = min(20.0, 0.4 / np.min(np.diff(t)))
    grid = np.arange(0.25, upper + 0.0025, 0.005)
    i = int(np.argmin([fit(f)[0] for f in grid]))
    if i in (0, len(grid) - 1):
        return {"fit_valid": False}
    opt = minimize_scalar(lambda f: fit(f)[0], bounds=(grid[i - 1], grid[i + 1]), method="bounded", options={"xatol": 1e-12})
    damped = least_squares(lambda z: fit(z[0], z[1], True), [opt.x, 0], bounds=([0.25, -2], [upper, 2]),
                           xtol=1e-12, ftol=1e-12, gtol=1e-12)
    fitted_frequency, decay = damped.x
    residual, coefficients = fit(fitted_frequency, decay)
    relative_residual = float(np.sqrt(residual / np.mean(y * y)))
    # A fitted number is not evidence of a resolved, single fundamental mode.
    valid = bool(damped.success and relative_residual <= 0.05 and
                 0.2501 < fitted_frequency < upper - 0.0001 and abs(decay) < 1.999)
    return dict(fit_valid=valid, frequency_ratio=float(fitted_frequency),
                decay_per_EB_period=float(decay),
                undamped_frequency_ratio=float(np.hypot(fitted_frequency, decay/(2*np.pi))),
                fit_amplitude=float(np.hypot(*coefficients[:2])),
                fit_residual_relative=relative_residual)


def analyze(directory):
    result = json.loads((directory / "result.json").read_text())
    parameters = json.loads((directory / "input.json").read_text())["parameters"]
    out = dict(name=directory.name, status=result["status"], parameters=parameters, wall_seconds=result["elapsed_seconds"])
    if result["status"] != "completed":
        out["failure"] = (directory / "failure.txt").read_text().strip() if (directory / "failure.txt").exists() else result["warnings"]
        if (directory / "history.csv").is_file():
            partial = table(directory / "history.csv")
            if len(partial):
                out["last_valid_t_over_T"] = float(partial["t_over_T"][-1])
                if parameters["load_mode"] == 0 and partial["kinetic"][0] > 0:
                    out["max_kinetic_over_initial"] = float(max(partial["kinetic"]) / partial["kinetic"][0])
        return out
    h = table(directory / "history.csv")
    p = table(directory / "profiles.csv")
    meta = json.loads((directory / "metadata.json").read_text())
    completion = json.loads((directory / "completion.json").read_text())
    expected_time = meta["period_eb"] * parameters["cycles"]
    if (not np.all(np.isfinite(np.column_stack([h[n] for n in h.dtype.names]))) or
        not np.all(np.isfinite(np.column_stack([p[n] for n in p.dtype.names]))) or
        not np.isclose(h["time"][-1], expected_time, atol=0, rtol=1e-10) or
        not np.isclose(completion["time"], expected_time, atol=0, rtol=1e-10) or
        len(h) != int(np.ceil(parameters["cycles"] * parameters["samples_per_cycle"])) + 1 or
        np.any(h["failure_count"] != 0) or np.any(h["damaged_count"] != 0) or
        np.any(h["live_count"] != meta["particles"])):
        raise ValueError("invalid completed-run data: " + str(directory))
    out["particles"] = meta["particles"]
    out["directed_bonds"] = meta["directed_bonds"]
    out["steps"] = int(h["steps"][-1])
    out["max_strain_gradient"] = float(max(h["max_F_minus_I"]))
    out["max_clamp_drift_over_H"] = float(max(h["clamp_drift_over_H"]))
    # Ignore t=0 and low-force crossings when quoting relative force cancellation.
    interior = (h["t_over_T"] > 0.1) & (np.abs(h["modal_over_H"]) > 0.1 * parameters["amplitude_ratio"])
    out["hg_net_relative_max"] = float(max(h["hg_net_relative"][interior], default=0))
    out["hg_torque_relative_max"] = float(max(h["hg_torque_relative"][interior], default=0))
    out["internal_net_relative_max"] = float(max(h["internal_net_relative"][interior], default=0))
    if parameters["load_mode"] == 0:
        out["max_kinetic_over_initial"] = float(max(h["kinetic"]) / h["kinetic"][0])
        out["vibration_bounded"] = out["max_kinetic_over_initial"] < 2
        if not out["vibration_bounded"]:
            out["fit_valid"] = False
            out["quality_status"] = "rejected_energy_growth"
            return out
        out["quality_status"] = "bounded"
        if h["t_over_T"][-1] < 1.5:
            out["fit_valid"] = False
            return out
        out.update(frequency(h["t_over_T"], h["modal_over_H"]))
        if not out["fit_valid"]:
            out["quality_status"] = "bounded_but_single_mode_fit_rejected"
        out["tip_history_L2_relative"] = float(np.linalg.norm(h["tip_over_H"] - h["tip_eb_over_H"]) / np.linalg.norm(h["tip_eb_over_H"]))
        t_peak = h["t_over_T"][np.argmax(np.abs(h["modal_over_H"]))]
        profile = p[np.isclose(p["t_over_T"], t_peak, atol=1e-10, rtol=0)]
        b = 1.875104068711961
        sigma = (np.cosh(b) + np.cos(b)) / (np.sinh(b) + np.sin(b))
        x = profile["x_over_L"]
        phi = np.cosh(b*x)-np.cos(b*x)-sigma*(np.sinh(b*x)-np.sin(b*x))
        y = profile["mean_uy_over_H"]
        fitted = phi * (phi @ y) / (phi @ phi)
        out["mode_shape_L2_relative"] = float(np.linalg.norm(y - fitted) / np.linalg.norm(y))
    else:
        last = p[np.isclose(p["t_over_T"], h["t_over_T"][-1], atol=1e-10, rtol=0)]
        eb, tim, y = last["eb_over_H"], last["timoshenko_static_over_H"], last["mean_uy_over_H"]
        out["tip_error_eb"] = float(y[-1] / eb[-1] - 1)
        out["tip_error_timoshenko"] = float(y[-1] / tim[-1] - 1)
        out["static_shape_L2_relative"] = float(np.linalg.norm(y - eb) / np.linalg.norm(eb))
        tail = h[h["t_over_T"] >= h["t_over_T"][-1] - 0.25]
        out["tail_tip_drift_relative"] = float(np.ptp(tail["tip_over_H"]) / abs(y[-1]))
        out["final_acc_rms_over_gravity"] = float(h["free_acc_rms_normalized"][-1])
        out["final_vel_rms_normalized"] = float(h["free_vel_rms_normalized"][-1])
        out["reaction_error"] = float(h["reaction_y"][-1] / (-meta["load_per_length"] * meta["length"]) - 1)
        out["static_equilibrated"] = (out["tail_tip_drift_relative"] < 1e-3 and
                                     out["final_acc_rms_over_gravity"] < 1e-2 and
                                     out["final_vel_rms_normalized"] < 1e-3)
        out["static_reaction_balanced"] = abs(out["reaction_error"]) <= 0.01
        out["quality_status"] = ("not_equilibrated" if not out["static_equilibrated"] else
                                 "equilibrated" if out["static_reaction_balanced"] else
                                 "stationary_but_reaction_imbalanced")
    return out


def compare_history(root, a, b, normalized_amplitude=False):
    x, y = table(root/a/"history.csv"), table(root/b/"history.csv")
    pa = json.loads((root/a/"input.json").read_text())["parameters"]
    pb = json.loads((root/b/"input.json").read_text())["parameters"]
    if abs(x["t_over_T"][-1] - y["t_over_T"][-1]) > 1e-8:
        raise ValueError("incompatible final times for history comparison")
    ya = x["modal_over_H"] / (pa["amplitude_ratio"] if normalized_amplitude else 1)
    yb = np.interp(x["t_over_T"], y["t_over_T"], y["modal_over_H"]) / (pb["amplitude_ratio"] if normalized_amplitude else 1)
    return dict(a=a, b=b, modal_history_L2_relative=float(np.linalg.norm(ya-yb)/np.linalg.norm(ya)),
                modal_history_max_relative=float(np.max(np.abs(ya-yb))/np.max(np.abs(ya))))


def verification_checks(root, results, comparisons):
    """Explicit engineering targets, not publication-wide universal tolerances."""
    rows = {r["name"]: r for r in results}
    checks = {}
    for name in ("static_new_n4", "static_new_n8", "static_legacy_n4", "static_legacy_n8"):
        if "reaction_error" in rows.get(name, {}):
            r = rows[name]
            checks[name + "_reaction_balance_1_percent"] = dict(
                passed=r["static_equilibrated"] and r["static_reaction_balanced"],
                value=abs(r["reaction_error"]), limit=0.01, stationary=r["static_equilibrated"])
    if "new_n8" in rows and rows["new_n8"].get("fit_valid"):
        error = abs(rows["new_n8"]["frequency_ratio"] - 1)
        checks["finest_frequency_within_5_percent"] = dict(passed=error <= 0.05, value=error, limit=0.05)
    if "static_new_n8" in rows and "tip_error_eb" in rows["static_new_n8"]:
        r = rows["static_new_n8"]
        checks["finest_static_within_5_percent"] = dict(passed=r["static_equilibrated"] and abs(r["tip_error_eb"]) <= 0.05,
            value=abs(r["tip_error_eb"]), limit=0.05, equilibrated=r["static_equilibrated"])
    by_pair = {(c["a"], c["b"]): c for c in comparisons}
    key = ("new_dt010", "new_dt005")
    if key in by_pair:
        value = by_pair[key]["modal_history_L2_relative"]
        checks["time_history_change_below_0_1_percent"] = dict(passed=value <= 0.001, value=value, limit=0.001)
    keys = [("new_n4", "new_scale0001"), ("new_n4", "new_scale10")]
    if all(k in by_pair for k in keys):
        value = max(by_pair[k]["modal_history_max_relative"] for k in keys)
        topology = all(rows[k[0]]["particles"] == rows[k[1]]["particles"] and
                       rows[k[0]]["directed_bonds"] == rows[k[1]]["directed_bonds"] for k in keys)
        checks["geometric_similarity"] = dict(passed=value <= 1e-6 and topology, value=value, limit=1e-6, identical_topology=topology)
    names = ("new_n3", "new_n4", "new_n6", "new_n8")
    if all(rows.get(n, {}).get("fit_valid") for n in names):
        errors = [abs(rows[n]["frequency_ratio"]-1) for n in names]
        checks["frequency_error_decreases_with_refinement"] = dict(passed=all(a>b for a,b in zip(errors,errors[1:])), errors=errors)
        change = abs(rows["new_n8"]["frequency_ratio"] / rows["new_n6"]["frequency_ratio"] - 1)
        checks["finest_frequency_grid_independent_1_percent"] = dict(passed=change <= 0.01, value=change, limit=0.01)
    key = ("new_n4", "new_small_amplitude")
    if key in by_pair:
        value = by_pair[key]["modal_history_L2_relative"]
        checks["amplitude_linearity"] = dict(passed=value <= 0.001, value=value, limit=0.001)
    return checks


def text_report(results, comparisons, checks):
    dynamic_damping = {r["parameters"].get("numerical_damping", False) for r in results
                       if r["parameters"]["load_mode"] == 0}
    damping_text = ("含现有成对数值阻尼" if dynamic_damping == {True} else
                    "关闭成对数值阻尼" if dynamic_damping == {False} else "含不同阻尼配置，须逐算例核对")
    lines = ["PERI-SPH 悬臂梁验证：计算摘要", "", "这是验证结果，不把完成运行等同于物理验证通过。",
             "基准：L/H=12，方截面，E=10 MPa，rho=1000 kg/m³，nu=0.3，a/H=0.001。",
             f"动态结果为{damping_text}的单模态拟合频率；EB 是细长梁参考，不是三维非局部精确解。",
             "静态比较在最后粒子截面 x=L-dx/2，采用均布体载荷；先检查静力平衡，再解释挠度误差。",
             "细化分固定 delta/dx 与固定 delta/H 两组；没有拟合有效梁长或重新调参。", "", "算例结果："]
    for r in results:
        line = f"{r['name']}: {r['status']}"
        if r.get("fit_valid"):
            line += f" | omega/omega_EB={r['frequency_ratio']:.8f}, 频率偏差={(r['frequency_ratio']-1)*100:+.4f}%, 振型L2={r['mode_shape_L2_relative']*100:.4f}%"
        if "tip_error_eb" in r:
            line += f" | 静挠度偏差={r['tip_error_eb']*100:+.4f}%, 残余运动达标={r['static_equilibrated']}, 反力偏差={r['reaction_error']*100:+.4f}%, 末段漂移={r['tail_tip_drift_relative']*100:.6f}%, RMS加速度/g={r['final_acc_rms_over_gravity']:.5g}"
        if "quality_status" in r:
            line += " | " + r["quality_status"]
        if "failure" in r:
            line += " | " + str(r["failure"])
        lines.append(line)
    lines.extend(["", "验证门槛（本项工程目标，并非通用发表标准）："])
    for name, check in checks.items():
        lines.append(("PASS " if check["passed"] else "NOT MET ") + name + ": " + json.dumps(check, ensure_ascii=False))
    lines.extend(["", "归一化历史比较："])
    for c in comparisons:
        lines.append(f"{c['a']} → {c['b']}: L2={c['modal_history_L2_relative']:.9g}, 最大差/基准峰值={c['modal_history_max_relative']:.9g}")
    lines.extend(["", "判读边界：", "1. 解析误差下降说明具有细化趋势，不等于已经网格无关。",
                  "2. 原算子仅保留其力公式，其他求解流程使用当前代码；这不是整个冰射流历史程序的重跑。",
                  "3. 原公式保留 delta^-3，新公式为三维 delta^-5；系数均为0.5，无逐尺度标定。",
                  "4. 稳定化合力与力矩检查包含夹持粒子、使用力评价时刻坐标；不等同于整梁无外力。",
                  "5. 无数值阻尼诊断另行保留；主矩阵的成功不能推出新算子单独保证稳定。",
                  "6. 小应变弹性能未含稳定化能，不能据此宣称完整能量闭合。", "",
                  "解析来源：", "TU Graz: https://lampz.tugraz.at/~hadley/memm/mechanics/cbeam.php",
                  "Clemson: https://www.clemson.edu/cecas/departments/me/academics/graduate/qualifying_exams/media/solid_mechanics/solid-3.pdf"])
    return "\n".join(lines)+"\n"


def figures(root, results, dest):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial"],
                         "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 13,
                         "axes.spines.right": False, "axes.spines.top": False,
                         "axes.linewidth": 2.5, "legend.frameon": False})
    colors = {"new": "#305C83", "legacy": "#B66C45", "off": "#999999"}
    def save(fig, name):
        for ext in ("png", "svg", "pdf"):
            fig.savefig(dest/(name+"."+ext), dpi=220, bbox_inches="tight", facecolor="white")
        plt.close(fig)
    good = {r["name"]: r for r in results if r["status"] == "completed"}
    if "new_n4" not in good:
        return
    fig = plt.figure(figsize=(14, 8.6), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=(1.3, 1))
    ax = fig.add_subplot(grid[0, :])
    for name, label, family in (("new_n4", "New pair stabilization", "new"), ("legacy_n4", "Original stabilization", "legacy"), ("off_n4", "No stabilization", "off")):
        if name in good:
            h = table(root/name/"history.csv")
            ax.plot(h["t_over_T"], h["tip_over_H"]*1000, label=label, color=colors[family], lw=1.8)
    h = table(root/"new_n4"/"history.csv")
    ax.plot(h["t_over_T"], h["tip_eb_over_H"]*1000, "k--", lw=1.6, label="Euler–Bernoulli")
    ax.set(xlabel=r"Time $t/T_{EB}$", ylabel=r"Tip-section $u_y/H$ ($10^{-3}$)")
    damping_label = "with pair damping" if good["new_n4"]["parameters"].get("numerical_damping") else "without pair damping"
    ax.set_title(f"a  Small-amplitude vibration {damping_label} (H/Δx = 4)", pad=48)
    ax.legend(ncols=4, fontsize=10, loc="upper center", bbox_to_anchor=(0.5, 1.13))
    ax = fig.add_subplot(grid[1, 0])
    for method in ("new", "legacy"):
        scales, frequency_values = [], []
        for name in (method+"_scale0001", method+"_n4", method+"_scale10"):
            if name in good and good[name].get("fit_valid"):
                scales.append(good[name]["parameters"]["scale"]); frequency_values.append(good[name]["frequency_ratio"])
        ax.plot(scales, frequency_values, "o-", color=colors[method], label="New" if method == "new" else "Original")
    ax.axhline(1, color="k", ls="--", lw=1)
    ax.set(xscale="log", xlabel="Geometric scale s", ylabel=r"$\omega_{PD}/\omega_{EB}$", title="b  Geometric similarity")
    if any(r.get("frequency_ratio", 0) > 3 for r in good.values()):
        ax.set_yscale("log")
    for r in results:
        if r["name"] == "legacy_scale10" and not r.get("fit_valid"):
            ax.text(10, 0.20, "Original: failed" if r["status"] == "failed" else "Original: invalid fit",
                    transform=ax.get_xaxis_transform(), ha="right", va="bottom", fontsize=10, color=colors["legacy"])
    ax.legend(fontsize=11)
    ax = fig.add_subplot(grid[1, 1])
    for names, label, style in ((("new_n3", "new_n4", "new_n6", "new_n8"), "Fixed δ/Δx = 3.015", "o-"),
                              (("new_fixed_delta_n3", "new_n4", "new_fixed_delta_n6"), "Fixed δ/H = 0.75375", "s--")):
        rows = [good[name] for name in names if name in good and good[name].get("fit_valid")]
        ax.plot([1/r["parameters"]["cells"] for r in rows], [r["frequency_ratio"] for r in rows], style, label=label, color=colors["new"])
    legacy_rows = [good[n] for n in ("legacy_n3", "legacy_n4", "legacy_n6", "legacy_n8") if n in good and good[n].get("fit_valid")]
    if len(legacy_rows) > 1:
        ax.plot([1/r["parameters"]["cells"] for r in legacy_rows], [r["frequency_ratio"] for r in legacy_rows], "^-", color=colors["legacy"], label="Original: fixed δ/Δx")
    ax.axhline(1, color="k", ls=":", lw=1)
    ax.set(xlabel="Particle spacing Δx/H", ylabel=r"$\omega_{PD}/\omega_{EB}$", title="c  Spatial refinement")
    if any(r["name"] == "legacy_n3" and r["status"] == "failed" for r in results):
        ax.text(1/3, 0.94, "Original n=3: failed", transform=ax.get_xaxis_transform(),
                ha="left", va="top", fontsize=10, color=colors["legacy"])
    ax.invert_xaxis(); ax.legend(fontsize=10)
    save(fig, "cantilever_dynamic")
    if {"new_n4", "new_dt010", "new_dt005", "new_scale0001", "new_scale10"} <= set(good):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), layout="constrained")
        fine = table(root/"new_dt005"/"history.csv")
        for name, label, color in (("new_n4", "CFL 0.20 − 0.05", "#8CA7BE"), ("new_dt010", "CFL 0.10 − 0.05", colors["new"])):
            h = table(root/name/"history.csv")
            difference = h["modal_over_H"] - np.interp(h["t_over_T"], fine["t_over_T"], fine["modal_over_H"])
            axes[0].plot(h["t_over_T"], difference / good[name]["parameters"]["amplitude_ratio"], color=color, label=label)
        axes[0].set(xlabel=r"Time $t/T_{EB}$", ylabel="Modal difference / prescribed amplitude", title="a  Temporal refinement")
        axes[0].ticklabel_format(axis="y", style="sci", scilimits=(0,0)); axes[0].legend(fontsize=10)
        for name, label, style, color in (("new_n4", "s=1", "-", "#305C83"), ("new_scale0001", "s=0.001", "--", "#587D9F"), ("new_scale10", "s=10", ":", "#9AB1C7")):
            h = table(root/name/"history.csv")
            axes[1].plot(h["t_over_T"], h["modal_over_H"] * 1000, style, label=label, lw=2, color=color)
        axes[1].set(xlabel=r"Time $t/T_{EB}$", ylabel=r"Modal displacement / H ($10^{-3}$)", title="b  New operator: scale collapse")
        axes[1].legend(fontsize=10)
        save(fig, "cantilever_time_and_scale")
    statics = [good[n] for n in ("static_new_n4", "static_new_n6", "static_new_n8", "static_legacy_n4", "static_legacy_n8") if n in good]
    if statics:
        fig, axes = plt.subplots(1, 3, figsize=(18, 4.8), layout="constrained")
        for i, row in enumerate(statics):
            name = row["name"]; p = table(root/name/"profiles.csv"); h = table(root/name/"history.csv")
            p = p[np.isclose(p["t_over_T"], max(p["t_over_T"]), atol=1e-10, rtol=0)]
            label = ("New" if "new" in name else "Original") + ", n=" + str(row["parameters"]["cells"])
            color = {4:"#9AB1C7",6:"#587D9F",8:"#305C83"}[row["parameters"]["cells"]] if "new" in name else colors["legacy"]
            style = "--" if name == "static_legacy_n8" else "-"
            axes[0].plot(p["x_over_L"], p["mean_uy_over_H"]*1000, label=label, lw=1.8, color=color, ls=style)
            axes[1].plot(h["t_over_T"], h["tip_over_H"]*1000, label=label, lw=1.8, color=color, ls=style)
            meta = json.loads((root/name/"metadata.json").read_text())
            axes[2].plot(h["t_over_T"], -h["reaction_y"]/(meta["load_per_length"]*meta["length"]),
                         label=label, lw=1.8, color=color, ls=style)
        axes[0].plot(p["x_over_L"], p["eb_over_H"]*1000, "k--", label="Euler–Bernoulli")
        axes[0].plot(p["x_over_L"], p["timoshenko_static_over_H"]*1000, color="#777777", ls=":", label="With shear correction")
        axes[0].set(xlabel="x/L", ylabel=r"$u_y/H$ ($10^{-3}$)", title="a  Final static profile")
        axes[1].set(xlabel=r"Relaxation time $t/T_{EB}$", ylabel=r"Tip-section $u_y/H$ ($10^{-3}$)", title="b  Relaxation history")
        axes[2].axhline(1, color="k", ls=":", lw=1.5)
        axes[2].set(xlabel=r"Relaxation time $t/T_{EB}$", ylabel=r"Reaction ratio $-R_y/(qL)$", title="c  Global force balance")
        axes[0].legend(fontsize=9); axes[1].legend(fontsize=10)
        save(fig, "cantilever_static")


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument("root", type=Path)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    root = args.root; dest = args.output or root/"analysis"
    dest.mkdir(parents=True, exist_ok=True)
    results = [analyze(p.parent) for p in sorted(root.glob("*/result.json"))]
    good = {r["name"] for r in results if r["status"] == "completed"}
    comparisons = []
    pairs = [("new_n4", "new_dt010"), ("new_dt010", "new_dt005"),
             ("new_n3", "new_n4"), ("new_n4", "new_n6"), ("new_n6", "new_n8"),
             ("new_fixed_delta_n3", "new_n4"), ("new_n4", "new_fixed_delta_n6")]
    pairs += [("legacy_n3", "legacy_n4"), ("legacy_n4", "legacy_n6"), ("legacy_n6", "legacy_n8")]
    for method in ("new", "legacy"):
        pairs.extend([(method+"_n4", method+"_scale0001"), (method+"_n4", method+"_scale10")])
    for a, b in pairs:
        if a in good and b in good:
            comparisons.append(compare_history(root, a, b))
    if {"new_n4", "new_small_amplitude"} <= good:
        comparisons.append(compare_history(root, "new_n4", "new_small_amplitude", True))
    summary = dict(cases=results, comparisons=comparisons,
                   cautions=["EB is a slender-beam reference, not an exact 3-D nonlocal solution.",
                             "Fixed delta/dx refinement changes both quadrature and horizon.",
                             "Fixed delta/H refinement tests a nonlocal model at fixed horizon.",
                             "Legacy coefficient retains the original dimensional scaling; no per-scale retuning.",
                             "Completion and small strain do not establish analytical accuracy or convergence."])
    checks = verification_checks(root, results, comparisons)
    summary["verification_checks"] = checks
    (dest/"summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")
    (dest/"verification_report.txt").write_text(text_report(results, comparisons, checks), encoding="utf-8")
    import scipy
    import matplotlib
    script = Path(__file__).resolve()
    shutil.copy2(script, dest/"analyze_cantilever.py")
    (dest/"analysis_provenance.json").write_text(json.dumps(dict(python=platform.python_version(), numpy=np.__version__,
        scipy=scipy.__version__, matplotlib=matplotlib.__version__, script_sha256=hashlib.sha256(script.read_bytes()).hexdigest(),
        data_root=str(root.resolve())), indent=2)+"\n")
    keys = sorted(set().union(*(set(r)-{"parameters", "failure"} for r in results)))
    with (dest/"metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys, extrasaction="ignore")
        writer.writeheader(); writer.writerows(results)
    figures(root, results, dest)
    for row in results:
        print(row["name"], row["status"], {k: round(row[k], 7) for k in ("frequency_ratio", "tip_error_eb", "tail_tip_drift_relative", "final_acc_rms_over_gravity") if k in row})
    print("Comparisons:", json.dumps(comparisons, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""3 塔チェーン予測（v4）。"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation

import paths_v4 as paths
from column_chain_v4 import simulate_chain
from competitive_adsorption_v4 import Params
from evaluate_model_accuracy_v4 import metrics_at_obs, plot_adsorption_profile_gif, plot_overlay_sparse
from feed_composition_v4 import feed_for_batch
from params_columns_v4 import (
    Cin_from_process_row,
    compute_auto_t_end,
    flow_df_from_process,
    get_params_by_name,
)
from parameter_fitting_v4 import apply_fitted_vector, load_fitted_vector


def _parse_list(s: str):
    return [x.strip() for x in s.split(",") if x.strip()]


def _load_column_params(name: str) -> Params:
    p = get_params_by_name(name)
    if name == "R003":
        fit = load_fitted_vector(paths.RESULTS_DIR / "R003" / "latest" / "fitted_vector.json")
        return apply_fitted_vector(p, fit)
    if name == "R004A":
        fit = load_fitted_vector(paths.RESULTS_DIR / "R004A" / "latest" / "fitted_vector.json")
        return apply_fitted_vector(p, fit)
    if name == "R004B":
        # R004A 流用
        fit = load_fitted_vector(paths.RESULTS_DIR / "R004A" / "latest" / "fitted_vector.json")
        p = apply_fitted_vector(p, fit)
        out = paths.RESULTS_DIR / "R004B" / "latest"
        out.mkdir(parents=True, exist_ok=True)
        payload = {
            "based_on": "R004A",
            "q_total": p.q_total,
            "eps_b": p.eps_b,
            "k_hyd": p.k_hyd,
            "keys": fit.get("keys"),
            "values": fit.get("values"),
        }
        (out / "used_parameters.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
        (out / "readme.txt").write_text("R004B uses R004A fitted parameters (prediction only).\n", encoding="utf-8")
        return p
    raise ValueError(name)


def _obs_for_column(name: str, batch: int, suffix: str, proc: pd.DataFrame) -> pd.DataFrame:
    if name == "R003":
        return proc[proc["batch"] == batch].copy()
    if name == "R004A":
        df = pd.read_csv(paths.r004a_csv(suffix))
        return df[df["batch"] == batch].copy()
    if name == "R004B":
        df = pd.read_csv(paths.r004b_csv(suffix))
        return df[df["batch"] == batch].copy()
    return pd.DataFrame()


def _plot_chain_overlay(
    results_by_col: Dict[str, dict],
    species_list: Sequence[str],
    out_path: Path,
    feed_concentrations: Optional[Dict[str, float]] = None,
) -> None:
    cols = list(results_by_col.keys())
    n_r, n_c = len(cols), len(species_list)
    fig, axes = plt.subplots(n_r, n_c, figsize=(3.6 * n_c, 2.8 * n_r), sharex=False)
    if n_r == 1 and n_c == 1:
        axes = np.array([[axes]])
    elif n_r == 1:
        axes = np.array([axes])
    elif n_c == 1:
        axes = np.array([[ax] for ax in axes])

    for i, col in enumerate(cols):
        eff = results_by_col[col]["effluent"]
        obs = results_by_col[col]["observed"]
        for j, sp in enumerate(species_list):
            ax = axes[i, j]
            ax.plot(eff["t_min"], eff[f"C_{sp}_out"], color="#1f4e79", lw=1.8, label="model")
            col_o = f"C_{sp}_exp"
            if col_o in obs.columns:
                m = obs[col_o].notna()
                ax.scatter(obs.loc[m, "t_min"], obs.loc[m, col_o], s=18, color="#c45c26", label="obs", zorder=3)
            if feed_concentrations is not None:
                cin = float(feed_concentrations.get(sp, 0.0) or 0.0)
                if cin > 0:
                    ax.axhline(cin, color="gray", ls="--", lw=1.0, label=f"{sp} feed (Cin={cin:g})")
            ax.set_title(f"{col} / {sp}", fontsize=9)
            ax.grid(True, alpha=0.25)
            if i == n_r - 1:
                ax.set_xlabel("t [min]")
            if j == 0:
                ax.set_ylabel("C")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def _plot_chain_adsorption_profile_gif(
    results_by_col: Dict[str, dict],
    params_by_col: Dict[str, Params],
    out_path: Path,
    *,
    frame_interval: float = 10.0,
    gif_duration_ms: int = 800,
) -> None:
    cols = list(results_by_col.keys())
    # 共通時刻
    all_t = []
    for col in cols:
        all_t.append(results_by_col[col]["grid_df"]["t_min"].unique())
    t_union = np.sort(np.unique(np.concatenate(all_t)))
    sel = [t_union[0]]
    for ti in t_union[1:]:
        if ti - sel[-1] >= frame_interval - 1e-9:
            sel.append(ti)
    if sel[-1] != t_union[-1]:
        sel.append(t_union[-1])

    colors = {"VE1": "#2a6fbb", "VE2": "#3aa76d", "FA": "#d17a22"}
    n = len(cols)
    fig, axes = plt.subplots(n, 1, figsize=(7.0, 2.6 * n), sharex=False)
    if n == 1:
        axes = [axes]

    frames = []
    for ti in sel:
        for ax, col in zip(axes, cols):
            ax.clear()
            ax2 = ax.twinx()
            gdf = results_by_col[col]["grid_df"]
            # 最近傍時刻
            ts = gdf["t_min"].unique()
            t_near = ts[np.argmin(np.abs(ts - ti))]
            g = gdf[np.isclose(gdf["t_min"], t_near)]
            p = params_by_col[col]
            z = g["z_cm"].to_numpy(dtype=float)
            qstack = np.zeros_like(z)
            for sp in p.adsorbing_species:
                q = g[f"q_{sp}"].to_numpy(dtype=float)
                ax.fill_between(z, qstack, qstack + q, color=colors.get(sp, "gray"), alpha=0.55, label=f"q_{sp}")
                qstack = qstack + q
            for sp in p.adsorbing_species:
                ax2.plot(z, g[f"C_{sp}"].to_numpy(dtype=float), color=colors.get(sp, "gray"), lw=1.3)
            ax.set_xlim(0, p.L)
            ax.set_ylabel(f"{col}\nq")
            ax.set_title(f"{col}  t={ti:.1f} min (snap {t_near:.1f})", fontsize=9)
        axes[-1].set_xlabel("z [cm]")
        fig.tight_layout()
        fig.canvas.draw()
        frames.append(np.asarray(fig.canvas.buffer_rgba()).copy())
        # twin axes cleanup
        for ax in axes:
            for child in ax.figure.axes:
                if child not in axes and child.axes.transAxes == ax.transAxes:
                    pass
        # remove extra twin axes created each frame
        while len(fig.axes) > n:
            fig.delaxes(fig.axes[-1])
    plt.close(fig)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig2, axim = plt.subplots(figsize=(7.0, 2.6 * n))
    im = axim.imshow(frames[0])
    axim.axis("off")

    def _upd(i):
        im.set_data(frames[i])
        return (im,)

    ani = animation.FuncAnimation(fig2, _upd, frames=len(frames), blit=True)
    ani.save(out_path, writer=animation.PillowWriter(fps=max(1, int(1000 / gif_duration_ms))))
    plt.close(fig2)


def main():
    ap = argparse.ArgumentParser(description="3-column chain prediction (v4)")
    ap.add_argument("--t_end", type=float, default=None)
    ap.add_argument("--batches", type=str, default="2,3,4,6")
    ap.add_argument("--columns", type=str, default="R003,R004A,R004B")
    ap.add_argument("--data-suffix", type=str, default="")
    ap.add_argument("--label", type=str, default="chain")
    ap.add_argument(
        "--plot_new_species",
        action="store_true",
        help="chain_overlay に FAEE/ET/OH パネルを追加",
    )
    args = ap.parse_args()

    batch_ids = [int(x) for x in _parse_list(args.batches)]
    col_names = _parse_list(args.columns)
    proc = pd.read_csv(paths.PROCESS_CSV)

    params_list = [_load_column_params(c) for c in col_names]
    out_root = paths.OUTPUTS_DIR / args.label
    summary_rows = []

    for bid in batch_ids:
        dfp = proc[proc["batch"] == bid].copy()
        if dfp.empty:
            print(f"[skip] batch {bid}: no process data")
            continue
        # t_end
        obs_arrays = [dfp["t_min"]]
        for c in col_names:
            o = _obs_for_column(c, bid, args.data_suffix, proc)
            if not o.empty:
                obs_arrays.append(o["t_min"])
        t_end = float(args.t_end) if args.t_end is not None else compute_auto_t_end(*obs_arrays)

        plist = []
        for p in params_list:
            pp = p.copy()
            pp.t_end = t_end
            plist.append(pp)
        cin = Cin_from_process_row(dfp, plist[0].species_names)
        cin.update(feed_for_batch(bid))
        plist[0].C_in = {**plist[0].C_in, **cin}
        flow = flow_df_from_process(dfp)

        print(f"batch {bid}: t_end={t_end:.1f}, columns={col_names}")
        chain = simulate_chain(plist, v_T_series=flow, C_in_first=cin, n_cols=len(plist))

        batch_dir = out_root / f"batch_{bid}"
        results_by_col = {}
        params_by_col = {}
        for name, (grid_df, effluent, meta, _), pcol in zip(col_names, chain, plist):
            cdir = batch_dir / name
            cdir.mkdir(parents=True, exist_ok=True)
            effluent.to_csv(cdir / "model_output.csv", index=False)
            obs = _obs_for_column(name, bid, args.data_suffix, proc)
            obs.to_csv(cdir / "observed.csv", index=False)
            sps = [sp for sp in pcol.adsorbing_species if f"C_{sp}_exp" in obs.columns]
            if not sps:
                sps = list(pcol.adsorbing_species)
            for sp in sps:
                plot_overlay_sparse(
                    effluent,
                    obs,
                    sp,
                    cdir / f"overlay_{sp}.png",
                    feed_concentrations=cin,
                    title=f"batch{bid} {name} {sp}",
                )
            try:
                plot_adsorption_profile_gif(grid_df, pcol, cdir / "adsorption_profile.gif")
            except Exception as e:
                print(f"[warn] gif {name}: {e}")
            mets = metrics_at_obs(effluent, obs, sps, edge_tol=1.0)
            for sp, m in mets.items():
                summary_rows.append(
                    {
                        "batch": bid,
                        "column": name,
                        "species": sp,
                        "rmse": m["rmse"],
                        "mape": m["mape"],
                        "n_obs": m["n_obs"],
                        "t_end": t_end,
                    }
                )
            results_by_col[name] = {"effluent": effluent, "observed": obs, "grid_df": grid_df}
            params_by_col[name] = pcol

        plot_sps = ["VE1", "VE2", "FA"]
        if args.plot_new_species:
            plot_sps = ["VE1", "VE2", "FA", "FAEE", "ET", "OH"]
        _plot_chain_overlay(results_by_col, plot_sps, batch_dir / "chain_overlay.png", cin)
        try:
            _plot_chain_adsorption_profile_gif(
                results_by_col, params_by_col, batch_dir / "chain_adsorption_profile.gif"
            )
        except Exception as e:
            print(f"[warn] chain gif: {e}")

        # per-batch summary
        pd.DataFrame([r for r in summary_rows if r["batch"] == bid]).to_csv(
            batch_dir / "chain_summary.csv", index=False
        )

    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(out_root / "chain_summary.csv", index=False)
        print(f"Wrote {out_root / 'chain_summary.csv'}")
    print("Chain prediction done.")


if __name__ == "__main__":
    main()

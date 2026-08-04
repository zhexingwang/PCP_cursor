"""モデル精度評価・可視化（v3）。"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation

from competitive_adsorption_v4 import Params, simulate
from parameter_fitting_v4 import _obs_species, _rmse_at_obs

FeedProvider = Callable[[int, pd.DataFrame], Optional[Dict[str, float]]]
InletProvider = Callable[[int, pd.DataFrame], Optional[Any]]
FlowProvider = Callable[[int, pd.DataFrame], Optional[pd.DataFrame]]
TEndProvider = Callable[[int, pd.DataFrame], Optional[float]]


def mape(y_true: np.ndarray, y_pred: np.ndarray, eps: float = 1e-8) -> float:
    denom = np.maximum(np.abs(y_true), eps)
    return float(np.mean(np.abs((y_true - y_pred) / denom)) * 100.0)


def metrics_at_obs(
    effluent: pd.DataFrame,
    df_obs: pd.DataFrame,
    species: Sequence[str],
    edge_tol: float = 1.0,
) -> Dict[str, Dict[str, float]]:
    t_model = effluent["t_min"].to_numpy(dtype=float)
    out: Dict[str, Dict[str, float]] = {}
    for sp in species:
        col_m = f"C_{sp}_out"
        col_o = f"C_{sp}_exp"
        if col_m not in effluent.columns or col_o not in df_obs.columns:
            continue
        t_obs = df_obs["t_min"].to_numpy(dtype=float)
        y_obs = df_obs[col_o].to_numpy(dtype=float)
        mask = np.isfinite(t_obs) & np.isfinite(y_obs)
        t_obs, y_obs = t_obs[mask], y_obs[mask]
        if t_obs.size == 0 or t_model.size == 0:
            continue
        y_hat = np.interp(
            t_obs,
            t_model,
            effluent[col_m].to_numpy(dtype=float),
            left=np.nan,
            right=np.nan,
        )
        t_max = float(t_model[-1])
        near = (t_obs > t_max) & (t_obs <= t_max + edge_tol)
        y_hat = np.asarray(y_hat, dtype=float)
        y_hat[near] = float(effluent[col_m].iloc[-1])
        valid = np.isfinite(y_hat)
        if not np.any(valid):
            continue
        err = y_hat[valid] - y_obs[valid]
        out[sp] = {
            "rmse": float(np.sqrt(np.mean(err**2))),
            "mape": mape(y_obs[valid], y_hat[valid]),
            "n_obs": int(np.sum(valid)),
        }
    return out


def plot_overlay_sparse(
    effluent: pd.DataFrame,
    df_obs: pd.DataFrame,
    species: str,
    out_path: Path,
    *,
    feed_concentrations: Optional[Dict[str, float]] = None,
    title: Optional[str] = None,
) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    t_m = effluent["t_min"].to_numpy(dtype=float)
    y_m = effluent[f"C_{species}_out"].to_numpy(dtype=float)
    ax.plot(t_m, y_m, color="#1f4e79", lw=2.0, label="model")
    col = f"C_{species}_exp"
    if col in df_obs.columns:
        t_o = df_obs["t_min"].to_numpy(dtype=float)
        y_o = df_obs[col].to_numpy(dtype=float)
        m = np.isfinite(t_o) & np.isfinite(y_o)
        ax.scatter(t_o[m], y_o[m], color="#c45c26", s=28, zorder=3, label="observed")
    if feed_concentrations is not None:
        cin = float(feed_concentrations.get(species, 0.0) or 0.0)
        if cin > 0:
            ax.axhline(
                cin,
                color="gray",
                ls="--",
                lw=1.1,
                label=f"{species} feed (Cin={cin:g})",
            )
    ax.set_xlabel("t [min]")
    ax.set_ylabel(f"C_{species} [mmol/cm³]")
    ax.set_title(title or species)
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


def plot_adsorption_profile_gif(
    grid_df: pd.DataFrame,
    params: Params,
    out_path: Path,
    *,
    frame_interval: float = 10.0,
    gif_duration_ms: int = 800,
) -> None:
    ads = list(params.adsorbing_species)
    times = np.sort(grid_df["t_min"].unique())
    if times.size == 0:
        return
    # 間引き
    sel = [times[0]]
    for ti in times[1:]:
        if ti - sel[-1] >= frame_interval - 1e-9:
            sel.append(ti)
    if sel[-1] != times[-1]:
        sel.append(times[-1])

    colors = {"VE1": "#2a6fbb", "VE2": "#3aa76d", "FA": "#d17a22"}
    fig, ax = plt.subplots(figsize=(7.0, 3.8))
    ax2 = ax.twinx()

    def draw(ti: float):
        ax.clear()
        ax2.clear()
        g = grid_df[np.isclose(grid_df["t_min"], ti)]
        z = g["z_cm"].to_numpy(dtype=float)
        qstack = np.zeros_like(z)
        for sp in ads:
            q = g[f"q_{sp}"].to_numpy(dtype=float)
            ax.fill_between(z, qstack, qstack + q, color=colors.get(sp, "gray"), alpha=0.55, label=f"q_{sp}")
            qstack = qstack + q
        for sp in ads:
            c = g[f"C_{sp}"].to_numpy(dtype=float)
            ax2.plot(z, c, color=colors.get(sp, "gray"), lw=1.5, label=f"C_{sp}")
        ax.set_xlim(0, params.L)
        ax.set_xlabel("z [cm]")
        ax.set_ylabel("q [mmol/cm³-resin]")
        ax2.set_ylabel("C [mmol/cm³]")
        ax.set_title(f"t = {ti:.1f} min")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=7, ncol=2)

    frames = []
    for ti in sel:
        draw(ti)
        fig.canvas.draw()
        frames.append(np.asarray(fig.canvas.buffer_rgba()).copy())
    plt.close(fig)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    # pillow writer
    fig2, axim = plt.subplots(figsize=(7.0, 3.8))
    im = axim.imshow(frames[0])
    axim.axis("off")

    def _upd(i):
        im.set_data(frames[i])
        return (im,)

    ani = animation.FuncAnimation(fig2, _upd, frames=len(frames), blit=True)
    ani.save(out_path, writer=animation.PillowWriter(fps=max(1, int(1000 / gif_duration_ms))))
    plt.close(fig2)


def _evaluate_one_set(
    params: Params,
    dfb: pd.DataFrame,
    out_dir: Path,
    *,
    flow: Optional[pd.DataFrame] = None,
    cin_series: Optional[Any] = None,
    feed_concentrations: Optional[Dict[str, float]] = None,
    make_gif: bool = True,
) -> Dict[str, Any]:
    p = params.copy()
    if cin_series is not None:
        p.C_in_series = cin_series
    grid_df, effluent, meta, _ = simulate(p, v_T_series=flow)
    out_dir.mkdir(parents=True, exist_ok=True)
    effluent.to_csv(out_dir / "model_output.csv", index=False)
    dfb.to_csv(out_dir / "observed.csv", index=False)

    sps = _obs_species(dfb, p.adsorbing_species) or list(p.adsorbing_species)
    mets = metrics_at_obs(effluent, dfb, sps, edge_tol=1.0)
    feed = feed_concentrations if feed_concentrations is not None else dict(p.C_in)
    for sp in sps:
        plot_overlay_sparse(
            effluent,
            dfb,
            sp,
            out_dir / f"overlay_{sp}.png",
            feed_concentrations=feed,
            title=f"{out_dir.name} {sp}",
        )
    if make_gif:
        try:
            plot_adsorption_profile_gif(grid_df, p, out_dir / "adsorption_profile.gif")
        except Exception as e:
            print(f"[warn] gif failed: {e}")
    return {
        "metrics": mets,
        "effluent": effluent,
        "grid_df": grid_df,
        "meta": meta,
        "feed": feed,
    }


def evaluate_model_batches(
    base_params: Params,
    batches: Dict[int, pd.DataFrame],
    out_root: Path,
    *,
    inlet_provider: Optional[InletProvider] = None,
    flow_provider: Optional[FlowProvider] = None,
    feed_provider: Optional[FeedProvider] = None,
    t_end_provider: Optional[TEndProvider] = None,
    make_gif: bool = True,
) -> Dict[int, Dict[str, Any]]:
    results = {}
    for bid, dfb in batches.items():
        p = base_params.copy()
        if t_end_provider is not None:
            te = t_end_provider(bid, dfb)
            if te is not None:
                p.t_end = float(te)
        cin = {sp: float(dfb.iloc[0].get(f"C0_{sp}", p.C_in.get(sp, 0.0)) or 0.0) for sp in p.species_names}
        if feed_provider is not None:
            fed = feed_provider(bid, dfb)
            if fed:
                cin.update({k: float(v) for k, v in fed.items()})
        p.C_in = {**p.C_in, **cin}

        flow = flow_provider(bid, dfb) if flow_provider else None
        if flow is None and "v_T" in dfb.columns:
            from params_columns_v4 import append_flow_stop

            flow = append_flow_stop(dfb[["t_min", "v_T"]].dropna().copy())

        cin_series = None
        if inlet_provider is not None:
            inl = inlet_provider(bid, dfb)
            if isinstance(inl, dict) and "C_in_series" in inl:
                cin_series = inl["C_in_series"]
            elif inl is not None:
                cin_series = inl

        out_dir = out_root / f"batch_{bid}"
        results[bid] = _evaluate_one_set(
            p,
            dfb,
            out_dir,
            flow=flow,
            cin_series=cin_series,
            feed_concentrations=dict(p.C_in),
            make_gif=make_gif,
        )
    return results

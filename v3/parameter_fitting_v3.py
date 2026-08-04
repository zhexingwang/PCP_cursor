"""バッチ対応パラメータフィット（Differential Evolution）v3。"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution

from competitive_adsorption_v3 import Params, simulate

DEFAULT_BOUNDS: Dict[str, Tuple[float, float]] = {
    "q_total": (0.2, 2.5),  # mmol/cm3-solid（論文スケール）
    "eps_b": (0.25, 0.55),
    "eps_p": (0.30, 0.80),
    "k_hyd": (0.0, 5.0),
    "k_ads.VE1": (1e-3, 5.0),
    "k_ads.VE2": (1e-3, 5.0),
    "k_ads.FA": (1e-3, 5.0),
    "Keq_ads.VE1": (1.0, 500.0),
    "Keq_ads.VE2": (1.0, 500.0),
    "Keq_ads.FA": (1.0, 1000.0),
    "k_ex.FA->VE1": (1e-4, 5.0),
    "k_ex.FA->VE2": (1e-4, 5.0),
    "Keq_ex.FA->VE1": (0.1, 50.0),
    "Keq_ex.FA->VE2": (0.1, 50.0),
    "H.VE1": (0.5, 10.0),
    "H.VE2": (0.5, 10.0),
    "H.FA": (0.5, 10.0),
    "H.OH": (1.0, 20.0),
}

InletProvider = Callable[[int, pd.DataFrame], Optional[Dict[str, Any]]]
FlowProvider = Callable[[int, pd.DataFrame], Optional[pd.DataFrame]]
FeedProvider = Callable[[int, pd.DataFrame], Optional[Dict[str, float]]]
TEndProvider = Callable[[int, pd.DataFrame], Optional[float]]


def _obs_species(df: pd.DataFrame, candidates: Sequence[str]) -> List[str]:
    out = []
    for sp in candidates:
        if f"C_{sp}_exp" in df.columns and df[f"C_{sp}_exp"].notna().any():
            out.append(sp)
    return out


def _rmse_at_obs(
    effluent: pd.DataFrame,
    df_obs: pd.DataFrame,
    species: Sequence[str],
    edge_tol: float = 0.0,
) -> Dict[str, float]:
    t_model = effluent["t_min"].to_numpy(dtype=float)
    rmses = {}
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
        right = y_model_right(effluent[col_m].to_numpy(dtype=float), t_model, edge_tol)
        y_hat = np.interp(t_obs, t_model, effluent[col_m].to_numpy(dtype=float), left=np.nan, right=right)
        # edge_tol: モデル最終時刻 + tol 以内なら最終値を使う
        if edge_tol > 0:
            t_max = float(t_model[-1])
            near = (t_obs > t_max) & (t_obs <= t_max + edge_tol)
            y_hat = y_hat.copy()
            y_hat[near] = float(effluent[col_m].iloc[-1])
        valid = np.isfinite(y_hat)
        if not np.any(valid):
            continue
        err = y_hat[valid] - y_obs[valid]
        rmses[sp] = float(np.sqrt(np.mean(err**2)))
    return rmses


def y_model_right(y: np.ndarray, t: np.ndarray, edge_tol: float) -> float:
    # フィット目的関数は端許容なし（right=nan）
    return np.nan if edge_tol <= 0 else float(y[-1])


def _objective_batches(
    vec: np.ndarray,
    base: Params,
    batches: Dict[int, pd.DataFrame],
    target_species: Optional[Sequence[str]],
    inlet_provider: Optional[InletProvider],
    flow_provider: Optional[FlowProvider],
    feed_provider: Optional[FeedProvider],
    t_end_provider: Optional[TEndProvider],
) -> float:
    p = base.copy()
    p.update_from_vector(vec)
    rmses_all: List[float] = []
    for bid, dfb in batches.items():
        pb = p.copy()
        if t_end_provider is not None:
            te = t_end_provider(bid, dfb)
            if te is not None:
                pb.t_end = float(te)

        cin = {sp: float(dfb.iloc[0].get(f"C0_{sp}", pb.C_in.get(sp, 0.0)) or 0.0) for sp in pb.species_names}
        if feed_provider is not None:
            fed = feed_provider(bid, dfb)
            if fed:
                cin.update({k: float(v) for k, v in fed.items()})
        pb.C_in = {**pb.C_in, **cin}

        flow = None
        if flow_provider is not None:
            flow = flow_provider(bid, dfb)
        if flow is None and "v_T" in dfb.columns:
            from params_columns_v3 import append_flow_stop

            flow = append_flow_stop(dfb[["t_min", "v_T"]].dropna().copy())

        # inlet_provider: 時系列入口（R004A 等）。指定時は C_in_series を使う
        if inlet_provider is not None:
            inl = inlet_provider(bid, dfb)
            if inl is not None:
                if isinstance(inl, dict) and "C_in_series" in inl:
                    pb.C_in_series = inl["C_in_series"]
                elif isinstance(inl, (pd.DataFrame, dict)):
                    pb.C_in_series = inl

        _, effluent, _, _ = simulate(pb, v_T_series=flow)
        sps = list(target_species) if target_species is not None else _obs_species(dfb, pb.adsorbing_species)
        rmses = _rmse_at_obs(effluent, dfb, sps, edge_tol=0.0)
        if rmses:
            rmses_all.append(float(np.mean(list(rmses.values()))))
        else:
            rmses_all.append(1e3)
    return float(np.mean(rmses_all)) if rmses_all else 1e3


def fit_parameters_batches(
    base_params: Params,
    batches: Dict[int, pd.DataFrame],
    *,
    target_species: Optional[Sequence[str]] = None,
    inlet_provider: Optional[InletProvider] = None,
    flow_provider: Optional[FlowProvider] = None,
    feed_provider: Optional[FeedProvider] = None,
    t_end_provider: Optional[TEndProvider] = None,
    bounds: Optional[Dict[str, Tuple[float, float]]] = None,
    maxiter: int = 25,
    popsize: int = 12,
    seed: int = 42,
    workers: int = 1,
    polish: bool = False,
) -> Dict[str, Any]:
    keys = base_params.fit_keys()
    bmap = bounds or DEFAULT_BOUNDS
    bseq = [bmap.get(k, (0.0, 1.0)) for k in keys]

    history: List[Dict[str, Any]] = []

    def cb(xk, convergence):
        history.append(
            {
                "n": len(history) + 1,
                "x": [float(v) for v in xk],
                "convergence": float(convergence),
            }
        )
        return False

    t0 = time.time()
    result = differential_evolution(
        _objective_batches,
        bounds=bseq,
        args=(
            base_params,
            batches,
            target_species,
            inlet_provider,
            flow_provider,
            feed_provider,
            t_end_provider,
        ),
        maxiter=maxiter,
        popsize=popsize,
        seed=seed,
        workers=workers,
        polish=polish,
        updating="immediate",
        callback=cb,
        tol=1e-3,
    )
    elapsed = time.time() - t0
    best = base_params.copy()
    best.update_from_vector(result.x)
    out = {
        "keys": keys,
        "values": [float(v) for v in result.x],
        "vector": {k: float(v) for k, v in zip(keys, result.x)},
        "fun": float(result.fun),
        "nit": int(result.nit),
        "success": bool(result.success),
        "message": str(result.message),
        "elapsed_sec": float(elapsed),
        "history": history,
        "params": best,
    }
    for k, v in zip(keys, result.x):
        out[k] = float(v)
    return out


def save_fit_result(
    out_dir: Path,
    fit: Dict[str, Any],
    *,
    extra: Optional[Dict[str, Any]] = None,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "keys": fit["keys"],
        "values": fit["values"],
        **fit["vector"],
        "fun": fit["fun"],
        "nit": fit["nit"],
        "success": fit["success"],
        "message": fit["message"],
        "elapsed_sec": fit["elapsed_sec"],
    }
    if extra:
        payload.update(extra)
    path = out_dir / "fitted_vector.json"
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    (out_dir / "de_progress.json").write_text(
        json.dumps(fit.get("history", []), indent=2), encoding="utf-8"
    )
    return path


def load_fitted_vector(path: Path) -> Dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def apply_fitted_vector(params: Params, fitted: Dict[str, Any]) -> Params:
    p = params.copy()
    keys = fitted.get("keys")
    values = fitted.get("values")
    if keys and values and len(keys) == len(values):
        # fit_targets をキーに合わせて有効化
        targets = {k: False for k in p.fit_targets}
        for k in keys:
            targets[k] = True
            if k not in p.fit_targets:
                targets[k] = True
        p.fit_targets = {**p.fit_targets, **targets}
        p.update_from_vector(values)
        return p
    # 個別キー
    for k in ("q_total", "eps_b", "eps_p", "k_hyd"):
        if k in fitted:
            p._set_by_key(k, float(fitted[k]))
    return p

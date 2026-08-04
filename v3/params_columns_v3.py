"""塔別 Params 生成ユーティリティ（v3）。"""
from __future__ import annotations

from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd

from competitive_adsorption_v3 import Params

T_END_MARGIN_FACTOR = 1.0


def _set_fit_targets_q_eps(p: Params) -> None:
    targets = {k: False for k in p.fit_targets}
    targets["q_total"] = True
    targets["eps_b"] = True
    p.fit_targets = targets


def get_R003_params() -> Params:
    p = Params()
    p.d = 25.0
    p.L = 83.2
    _set_fit_targets_q_eps(p)
    return p


def get_R004A_params() -> Params:
    p = Params()
    p.d = 25.0
    p.L = 42.0
    _set_fit_targets_q_eps(p)
    return p


def get_R004B_params() -> Params:
    p = Params()
    p.d = 25.0
    p.L = 42.0
    _set_fit_targets_q_eps(p)
    return p


def get_params_by_name(name: str) -> Params:
    name = name.upper()
    if name == "R003":
        return get_R003_params()
    if name == "R004A":
        return get_R004A_params()
    if name == "R004B":
        return get_R004B_params()
    raise ValueError(f"unknown column: {name}")


def apply_dt_react_cap_all(params_list: Sequence[Params], dt_react_cap: float) -> None:
    for p in params_list:
        p.dt_react_cap = float(dt_react_cap)


def compute_auto_t_end(
    *t_min_arrays: Iterable[float],
    margin_factor: float = T_END_MARGIN_FACTOR,
) -> float:
    vals = []
    for arr in t_min_arrays:
        a = np.asarray(list(arr), dtype=float)
        if a.size:
            vals.append(float(np.nanmax(a)))
    if not vals:
        return 350.0
    return float(max(vals) * margin_factor)


def append_flow_stop(flow_df: pd.DataFrame, eps: float = 1e-9) -> pd.DataFrame:
    """工程終了直後に v_T=0 を付加し、延長時の通液継続を防ぐ。"""
    df = flow_df.copy().sort_values("t_min").reset_index(drop=True)
    if df.empty:
        return df
    t_last = float(df["t_min"].iloc[-1])
    row = {c: df[c].iloc[-1] for c in df.columns}
    row["t_min"] = t_last + eps
    row["v_T"] = 0.0
    return pd.concat([df, pd.DataFrame([row])], ignore_index=True)


def flow_df_from_process(df_batch: pd.DataFrame) -> pd.DataFrame:
    cols = ["t_min", "v_T"]
    out = df_batch[cols].dropna().copy()
    return append_flow_stop(out)


def Cin_from_process_row(df_batch: pd.DataFrame, species: Optional[Sequence[str]] = None) -> dict:
    """バッチ先頭行の C0_* から入口濃度 dict を作る。"""
    species = list(species or ["VE1", "VE2", "FA", "OH"])
    row = df_batch.iloc[0]
    cin = {}
    for sp in species:
        col = f"C0_{sp}"
        cin[sp] = float(row[col]) if col in df_batch.columns and pd.notna(row[col]) else 0.0
    return cin

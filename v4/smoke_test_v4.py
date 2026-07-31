#!/usr/bin/env python3
"""v4 動作検証（付録 F）。"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from competitive_adsorption_v4 import Params, simulate
from feed_composition_v4 import feed_for_batch, reload_feed_csv
from params_columns_v4 import append_flow_stop, get_R003_params


def _flow(t_end: float, v_T: float = 40.0) -> pd.DataFrame:
    t = np.arange(0.0, t_end + 1e-9, 5.0)
    return append_flow_stop(pd.DataFrame({"t_min": t, "v_T": np.full_like(t, v_T)}))


def test_v3_compat():
    """k_hyd=0・原料 FAEE/ET/水=0 で既存成分が有限・非負。"""
    p = get_R003_params()
    p.k_hyd = 0.0
    p.t_end = 80.0
    p.N = 12
    p.dt = 0.2
    p.sample_dt = 2.0
    p.C_in = {"VE1": 0.02635, "VE2": 0.012, "FA": 0.008, "OH": 0.0, "FAEE": 0.0, "ET": 0.0}
    _, e, _, _ = simulate(p, v_T_series=_flow(p.t_end))
    for sp in ("VE1", "VE2", "FA", "OH"):
        y = e[f"C_{sp}_out"].to_numpy(dtype=float)
        assert np.all(np.isfinite(y)), sp
        assert np.all(y >= -1e-12), sp
    print("[F.1] v3-compat (k_hyd=0, no feed water/FAEE) OK")


def test_hydrolysis_direction():
    p0 = get_R003_params()
    p0.k_hyd = 0.0
    p0.t_end = 200.0
    p0.N = 16
    p0.dt = 0.2
    p0.sample_dt = 2.0
    cin = {"VE1": 0.026, "VE2": 0.012, "FA": 0.008, "OH": 0.243, "FAEE": 0.774, "ET": 0.774}
    p0.C_in = cin
    p1 = p0.copy()
    p1.k_hyd = 1.0
    flow = _flow(200.0)
    _, e0, _, _ = simulate(p0, v_T_series=flow)
    _, e1, _, _ = simulate(p1, v_T_series=flow)
    # 出口積算（単純和）
    s0 = {sp: float(e0[f"C_{sp}_out"].sum()) for sp in ("FAEE", "FA", "ET")}
    s1 = {sp: float(e1[f"C_{sp}_out"].sum()) for sp in ("FAEE", "FA", "ET")}
    print(f"  FAEE sum {s0['FAEE']:.3f} -> {s1['FAEE']:.3f}")
    print(f"  FA   sum {s0['FA']:.3f} -> {s1['FA']:.3f}")
    print(f"  ET   sum {s0['ET']:.3f} -> {s1['ET']:.3f}")
    assert s1["FAEE"] < s0["FAEE"]
    assert s1["FA"] > s0["FA"]
    assert s1["ET"] > s0["ET"]
    print("[F.2] hydrolysis direction OK")

    # モル保存: FAEE+ET は反応で不変（同じ θ）
    sum0 = e0["C_FAEE_out"].to_numpy() + e0["C_ET_out"].to_numpy()
    sum1 = e1["C_FAEE_out"].to_numpy() + e1["C_ET_out"].to_numpy()
    # 時刻グリッドが同一前提
    diff = np.max(np.abs(sum0 - sum1))
    print(f"  max |(FAEE+ET)_0 - (FAEE+ET)_1| = {diff:.3e}")
    assert diff < 1e-6
    print("[F.3] mole invariant OK")

    for sp in ("VE1", "VE2", "FA", "OH", "FAEE", "ET"):
        y = e1[f"C_{sp}_out"].to_numpy()
        assert np.all(np.isfinite(y)) and np.all(y >= -1e-12)
    print("[F.4] numeric sanity OK")


def test_feed_composition():
    reload_feed_csv()
    f6 = feed_for_batch(6)
    assert abs(f6["FAEE"] - 0.665) < 1e-12
    assert abs(f6["ET"] - 0.665) < 1e-12
    assert abs(f6["OH"] - 0.243) < 1e-12
    f12 = feed_for_batch(12)
    assert abs(f12["FAEE"] - 0.774) < 1e-12
    import warnings

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        f2 = feed_for_batch(2)
        assert f2["FAEE"] == 0.0 and f2["ET"] == 0.0 and abs(f2["OH"] - 0.243) < 1e-12
        assert any("batch 2" in str(x.message) for x in w)
    print("[F.5] feed_composition_v4 OK")


def test_cli_help():
    import subprocess

    for script in (
        "run_fit_R003_v4.py",
        "run_fit_R004A_v4.py",
        "run_chain_predict_v4.py",
        "run_chain_visualize_v4.py",
    ):
        r = subprocess.run([sys.executable, str(ROOT / script), "--help"], capture_output=True)
        assert r.returncode == 0, script
    print("[F.6] CLI --help OK")


def main():
    test_v3_compat()
    test_hydrolysis_direction()
    test_feed_composition()
    test_cli_help()
    print("\nAll smoke tests passed.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""合成プロセスデータ生成（動作確認用）。"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from competitive_adsorption_v4 import simulate
from params_columns_v4 import append_flow_stop, get_R003_params, get_R004A_params, get_R004B_params
from column_chain_v4 import simulate_chain

BATCHES = {
    2: {"C0_VE1": 0.0240, "C0_VE2": 0.0110, "C0_FA": 0.0075, "v_T": 38.0, "t_end": 420.0},
    3: {"C0_VE1": 0.02635, "C0_VE2": 0.0120, "C0_FA": 0.0080, "v_T": 40.0, "t_end": 426.0},
    4: {"C0_VE1": 0.0270, "C0_VE2": 0.0130, "C0_FA": 0.0085, "v_T": 42.0, "t_end": 450.0},
    6: {"C0_VE1": 0.0255, "C0_VE2": 0.0115, "C0_FA": 0.0078, "v_T": 41.0, "t_end": 460.0},
    12: {"C0_VE1": 0.0260, "C0_VE2": 0.0125, "C0_FA": 0.0082, "v_T": 39.0, "t_end": 426.0},
}

TRUE_Q = {"R003": 0.55, "R004A": 0.50, "R004B": 0.50}  # mmol/cm3-solid
TRUE_EPS = {"R003": 0.40, "R004A": 0.41, "R004B": 0.41}


def _flow(t_end: float, v_T: float) -> pd.DataFrame:
    t = np.arange(0.0, t_end + 1e-9, 5.0)
    df = pd.DataFrame({"t_min": t, "v_T": np.full_like(t, v_T)})
    return append_flow_stop(df)


def main():
    rng = np.random.default_rng(0)
    proc_rows = []
    r4a_rows = []
    r4b_rows = []

    for bid, cfg in BATCHES.items():
        p3 = get_R003_params()
        p3.q_total = TRUE_Q["R003"]
        p3.eps_b = TRUE_EPS["R003"]
        p3.t_end = cfg["t_end"]
        p3.N = 16
        p3.dt = 0.2
        p3.sample_dt = 2.0
        cin = {
            "VE1": cfg["C0_VE1"],
            "VE2": cfg["C0_VE2"],
            "FA": cfg["C0_FA"],
            "OH": 0.0,
        }
        p3.C_in = {**p3.C_in, **cin}
        flow = _flow(cfg["t_end"], cfg["v_T"])

        p4a = get_R004A_params()
        p4a.q_total = TRUE_Q["R004A"]
        p4a.eps_b = TRUE_EPS["R004A"]
        p4a.t_end = cfg["t_end"]
        p4a.N = 16
        p4a.dt = 0.2
        p4a.sample_dt = 2.0

        p4b = get_R004B_params()
        p4b.q_total = TRUE_Q["R004B"]
        p4b.eps_b = TRUE_EPS["R004B"]
        p4b.t_end = cfg["t_end"]
        p4b.N = 16
        p4b.dt = 0.2
        p4b.sample_dt = 2.0

        chain = simulate_chain([p3, p4a, p4b], v_T_series=flow, C_in_first=cin)
        eff3, eff4a, eff4b = chain[0][1], chain[1][1], chain[2][1]

        # R003 dense observations every 10 min
        for t in np.arange(0, cfg["t_end"] + 1e-9, 10.0):
            row = {
                "batch": bid,
                "t_min": t,
                "v_T": cfg["v_T"],
                "C0_VE1": cfg["C0_VE1"],
                "C0_VE2": cfg["C0_VE2"],
                "C0_FA": cfg["C0_FA"],
                "C0_OH": 0.0,
            }
            for sp in ("VE1", "VE2", "FA"):
                y = float(np.interp(t, eff3["t_min"], eff3[f"C_{sp}_out"]))
                noise = 1.0 + float(rng.normal(0, 0.02))
                row[f"C_{sp}_exp"] = max(0.0, y * noise)
            proc_rows.append(row)

        # R004A sparse 4 points
        for t in np.linspace(cfg["t_end"] * 0.35, cfg["t_end"] * 0.95, 4):
            row = {"batch": bid, "t_min": float(t)}
            for sp in ("VE1", "VE2", "FA"):
                y = float(np.interp(t, eff4a["t_min"], eff4a[f"C_{sp}_out"]))
                # 仕様: R003 出口の 30-50% 程度 → ここではモデル出口に近い値
                row[f"C_{sp}_exp"] = max(0.0, y * (1.0 + float(rng.normal(0, 0.03))))
            r4a_rows.append(row)

        # R004B 1 point near end
        t = cfg["t_end"] * 0.98
        row = {"batch": bid, "t_min": float(t)}
        for sp in ("VE1", "VE2", "FA"):
            y = float(np.interp(t, eff4b["t_min"], eff4b[f"C_{sp}_out"]))
            row[f"C_{sp}_exp"] = max(0.0, y)
        r4b_rows.append(row)

    out = Path(__file__).resolve().parent
    pd.DataFrame(proc_rows).to_csv(out / "process_runs_batches.csv", index=False)
    pd.DataFrame(r4a_rows).to_csv(out / "r004a_observed_batches_test.csv", index=False)
    pd.DataFrame(r4b_rows).to_csv(out / "r004b_observed_batches_test.csv", index=False)

    # 実データ用テンプレ（空）
    pd.DataFrame(columns=["batch", "t_min", "C_VE1_exp", "C_VE2_exp", "C_FA_exp"]).to_csv(
        out / "r004a_observed_batches.csv", index=False
    )
    pd.DataFrame(columns=["batch", "t_min", "C_VE1_exp", "C_VE2_exp", "C_FA_exp"]).to_csv(
        out / "r004b_observed_batches.csv", index=False
    )
    # 便宜上、実データ CSV にもテストデータを入れておく（空だと本番パスが落ちるため）
    # → 仕様は空テンプレ。実データはユーザーが記入。テストは --data-suffix _test。
    print("Wrote synthetic CSVs to", out)


if __name__ == "__main__":
    main()

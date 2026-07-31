#!/usr/bin/env python3
"""v3 最小スモークテスト。"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from column_chain_v3 import simulate_chain
from competitive_adsorption_v3 import Params, effluent_to_Cin_series, simulate
from params_columns_v3 import append_flow_stop, get_R003_params, get_R004A_params, get_R004B_params


def test_cin_series_and_chain():
    p = get_R003_params()
    p.t_end = 30.0
    p.N = 10
    p.dt = 0.2
    p.sample_dt = 1.0
    assert p.C_in_series is None
    assert abs(p.C_in_at(10.0)["VE1"] - p.C_in["VE1"]) < 1e-15

    flow = append_flow_stop(pd.DataFrame({"t_min": [0.0, 30.0], "v_T": [40.0, 40.0]}))
    p2 = get_R004A_params()
    p2.t_end = 30.0
    p2.N = 10
    p2.dt = 0.2
    p2.sample_dt = 1.0
    p3 = get_R004B_params()
    p3.t_end = 30.0
    p3.N = 10
    p3.dt = 0.2
    p3.sample_dt = 1.0

    chain = simulate_chain([p, p2, p3], v_T_series=flow, C_in_first=p.C_in)
    cin = effluent_to_Cin_series(chain[0][1], p.species_names)
    # リレー濃度破壊なし（出口値と Cin 列が一致）
    for sp in p.adsorbing_species:
        d = np.max(np.abs(chain[0][1][f"C_{sp}_out"].to_numpy() - cin[sp].to_numpy()))
        assert d == 0.0, (sp, d)
    maxes = [ch[1]["C_VE1_out"].max() for ch in chain]
    assert maxes[0] >= maxes[1] >= maxes[2]
    print("[A.2/A.3] C_in_series + chain relay OK", maxes)


def test_cli():
    for script in (
        "run_fit_R003_v3.py",
        "run_fit_R004A_v3.py",
        "run_chain_predict_v3.py",
        "run_chain_visualize_v3.py",
    ):
        r = subprocess.run([sys.executable, str(ROOT / script), "--help"], capture_output=True)
        assert r.returncode == 0, script
    print("[A.4] CLI --help OK")


def main():
    test_cin_series_and_chain()
    test_cli()
    print("\nAll v3 smoke tests passed.")


if __name__ == "__main__":
    main()

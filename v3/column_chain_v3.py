"""3 塔リレー実行コア（v3）。"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import pandas as pd

from competitive_adsorption_v3 import Params, effluent_to_Cin_series, simulate

SimResult = Tuple[pd.DataFrame, pd.DataFrame, Dict[str, Any], List[str]]


def _effluent_to_Cin_series(effluent_df: pd.DataFrame, species_names: Sequence[str]) -> pd.DataFrame:
    return effluent_to_Cin_series(effluent_df, species_names)


def simulate_chain(
    params_list: Sequence[Params],
    v_T_series: Optional[Union[pd.DataFrame, Tuple]] = None,
    C_in_first: Optional[Dict[str, float]] = None,
    n_cols: Optional[int] = None,
) -> List[SimResult]:
    """R003 → R004A → R004B を順に simulate し、前塔出口を後塔入口にリレーする。"""
    n = len(params_list) if n_cols is None else min(int(n_cols), len(params_list))
    results: List[SimResult] = []
    next_Cin_series: Optional[pd.DataFrame] = None

    for i in range(n):
        p = params_list[i].copy()
        if v_T_series is not None:
            p.v_T_series = v_T_series

        if i == 0:
            if C_in_first is not None:
                p.C_in = {**p.C_in, **C_in_first}
            p.C_in_series = None
        else:
            assert next_Cin_series is not None
            p.C_in_series = next_Cin_series

        grid_df, effluent_df, meta, fields = simulate(p)
        results.append((grid_df, effluent_df, meta, fields))
        next_Cin_series = _effluent_to_Cin_series(effluent_df, p.species_names)

    return results

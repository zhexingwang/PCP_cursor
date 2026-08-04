"""原料組成（FAEE / ET / 水=OH）管理（v4）。"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

import paths_v4 as paths

# CSV に行が無いバッチ向けの固定近似値
FEED_FIXED_DEFAULTS: Dict[str, float] = {
    "FAEE": 0.0,
    "ET": 0.0,
    "OH": 0.243,
}

_cache: Optional[pd.DataFrame] = None


def _load_csv(path: Optional[Path] = None) -> pd.DataFrame:
    global _cache
    path = path or paths.FEED_COMPOSITION_CSV
    if _cache is not None and path == paths.FEED_COMPOSITION_CSV:
        return _cache
    if not Path(path).exists():
        df = pd.DataFrame(columns=["batch", "C0_FAEE", "C0_ET", "C0_OH"])
    else:
        # # コメント行をスキップ
        df = pd.read_csv(path, comment="#")
    if path == paths.FEED_COMPOSITION_CSV:
        _cache = df
    return df


def reload_feed_csv() -> None:
    global _cache
    _cache = None
    _load_csv()


def feed_for_batch(batch_id: int, csv_path: Optional[Path] = None) -> Dict[str, float]:
    """バッチの原料組成を返す。CSV 優先、無ければ FEED_FIXED_DEFAULTS。"""
    df = _load_csv(csv_path)
    out = dict(FEED_FIXED_DEFAULTS)
    hit = df[df["batch"] == int(batch_id)] if "batch" in df.columns else df.iloc[0:0]
    if hit.empty:
        warnings.warn(
            f"feed_composition: batch {batch_id} not in CSV — "
            f"using defaults FAEE={out['FAEE']}, ET={out['ET']}, OH={out['OH']}",
            stacklevel=2,
        )
        return out
    row = hit.iloc[0]
    for key, col in (("FAEE", "C0_FAEE"), ("ET", "C0_ET"), ("OH", "C0_OH")):
        if col in hit.columns and pd.notna(row[col]):
            out[key] = float(row[col])
    return out


def make_feed_provider(csv_path: Optional[Path] = None):
    def _provider(batch_id: int, dfb=None):
        return feed_for_batch(batch_id, csv_path=csv_path)

    return _provider

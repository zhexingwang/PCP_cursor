"""v3 パス定義。"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results_v3"
OUTPUTS_DIR = ROOT / "outputs_chain_v3"

PROCESS_CSV = DATA_DIR / "process_runs_batches.csv"


def r004a_csv(suffix: str = "") -> Path:
    return DATA_DIR / f"r004a_observed_batches{suffix}.csv"


def r004b_csv(suffix: str = "") -> Path:
    return DATA_DIR / f"r004b_observed_batches{suffix}.csv"


COLUMN_NAMES = ("R003", "R004A", "R004B")

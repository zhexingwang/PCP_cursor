#!/usr/bin/env python3
"""保存済みパラメータで 3 塔を再可視化（フィットなし）v3。"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description="Re-visualize chain with saved parameters (v4)")
    ap.add_argument("--t_end", type=float, default=None)
    ap.add_argument("--batches", type=str, default="6")
    ap.add_argument("--columns", type=str, default="R003,R004A,R004B")
    ap.add_argument("--data-suffix", type=str, default="")
    ap.add_argument("--label", type=str, default="visualize")
    ap.add_argument("--plot_new_species", action="store_true")
    args = ap.parse_args()

    cmd = [
        sys.executable,
        str(Path(__file__).resolve().parent / "run_chain_predict_v4.py"),
        "--batches",
        args.batches,
        "--columns",
        args.columns,
        "--data-suffix",
        args.data_suffix,
        "--label",
        args.label,
    ]
    if args.t_end is not None:
        cmd.extend(["--t_end", str(args.t_end)])
    if args.plot_new_species:
        cmd.append("--plot_new_species")
    raise SystemExit(subprocess.call(cmd))


if __name__ == "__main__":
    main()

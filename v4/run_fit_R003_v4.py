#!/usr/bin/env python3
"""R003 単塔フィット（v4: q_total + k_hyd）。"""
from __future__ import annotations

import argparse
import shutil
from datetime import date
from pathlib import Path

import pandas as pd

import paths_v4 as paths
from feed_composition_v4 import make_feed_provider
from params_columns_v4 import (
    Cin_from_process_row,
    compute_auto_t_end,
    flow_df_from_process,
    get_R003_params,
)
from parameter_fitting_v4 import apply_fitted_vector, fit_parameters_batches, save_fit_result
from evaluate_model_accuracy_v4 import evaluate_model_batches


def _parse_batches(s: str):
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def _dated_dir(root: Path) -> Path:
    base = root / date.today().isoformat()
    if not base.exists():
        return base
    i = 2
    while True:
        cand = root / f"{date.today().isoformat()}_{i}"
        if not cand.exists():
            return cand
        i += 1


def main():
    ap = argparse.ArgumentParser(description="Fit R003 competitive adsorption parameters (v4)")
    ap.add_argument("--t_end", type=float, default=None)
    ap.add_argument("--train_batches", type=str, default="6")
    ap.add_argument("--test_batches", type=str, default=None)
    ap.add_argument("--species", type=str, default=None)
    ap.add_argument("--maxiter", type=int, default=15)
    ap.add_argument("--popsize", type=int, default=10)
    args = ap.parse_args()

    df = pd.read_csv(paths.PROCESS_CSV)
    train_ids = _parse_batches(args.train_batches)
    test_ids = _parse_batches(args.test_batches) if args.test_batches else list(train_ids)
    train = {b: g.copy() for b, g in df.groupby("batch") if b in train_ids}
    test = {b: g.copy() for b, g in df.groupby("batch") if b in test_ids}
    if not train:
        raise SystemExit(f"no train batches found: {train_ids}")

    p = get_R003_params()
    target_species = [s.strip() for s in args.species.split(",")] if args.species else None
    if target_species is not None:
        print(f"Fit species : {target_species}")
    else:
        print("Fit species : auto")

    t_end_mode = "fixed" if args.t_end is not None else "auto_per_batch"
    t_end_map = {}
    if args.t_end is not None:
        p.t_end = float(args.t_end)
        t_end_desc = f"fixed {args.t_end}"
        t_end_provider = None
    else:
        for bid, dfb in {**train, **test}.items():
            t_end_map[bid] = compute_auto_t_end(dfb["t_min"])
        t_end_desc = ", ".join(f"b{b}:{t_end_map[b]:.1f}" for b in sorted(t_end_map))
        t_end_provider = lambda bid, dfb: t_end_map.get(bid)

    def flow_provider(bid, dfb):
        return flow_df_from_process(dfb)

    feed_provider = make_feed_provider()

    print(f"Train batches: {sorted(train)}")
    print(f"t_end mode   : {t_end_mode} ({t_end_desc})")
    print(f"Fit targets  : {p.fit_keys()}")
    fit = fit_parameters_batches(
        p,
        train,
        target_species=target_species,
        flow_provider=flow_provider,
        feed_provider=feed_provider,
        t_end_provider=t_end_provider,
        maxiter=args.maxiter,
        popsize=args.popsize,
    )
    print(f"DE done: fun={fit['fun']:.6g}, nit={fit['nit']}, keys={fit['keys']}, values={fit['values']}")

    out_dir = _dated_dir(paths.RESULTS_DIR / "R003")
    extra = {
        "column": "R003",
        "train_batches": train_ids,
        "test_batches": test_ids,
        "fit_species": target_species if target_species else "auto",
        "fit_species_desc": ",".join(target_species) if target_species else "auto (all observed)",
        "t_end_mode": t_end_mode,
        "t_end_desc": t_end_desc,
        "t_end": args.t_end,
    }
    save_fit_result(out_dir, fit, extra=extra)
    readme = [
        f"Column: R003",
        f"Train batches: {train_ids}",
        f"Test batches: {test_ids}",
        f"Fit species: {extra['fit_species_desc']}",
        f"t_end: {t_end_desc}",
        f"Objective: {fit['fun']}",
        f"Elapsed[s]: {fit['elapsed_sec']:.1f}",
        f"Keys: {fit['keys']}",
        f"Values: {fit['values']}",
    ]
    (out_dir / "readme.txt").write_text("\n".join(readme) + "\n", encoding="utf-8")

    latest = paths.RESULTS_DIR / "R003" / "latest"
    if latest.exists():
        shutil.rmtree(latest)
    shutil.copytree(out_dir, latest)
    print(f"Saved: {out_dir}")

    # 評価
    p_fit = fit["params"]
    evaluate_model_batches(
        p_fit,
        test,
        paths.OUTPUTS_DIR / "R003" / "test",
        flow_provider=flow_provider,
        feed_provider=feed_provider,
        t_end_provider=t_end_provider,
        make_gif=True,
    )
    print("Evaluation done.")


if __name__ == "__main__":
    main()

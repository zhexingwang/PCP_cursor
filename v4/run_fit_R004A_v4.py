#!/usr/bin/env python3
"""R004A フィット（R003 出口を入口にリレー）v4。"""
from __future__ import annotations

import argparse
import json
import shutil
from datetime import date
from pathlib import Path

import pandas as pd

import paths_v4 as paths
from competitive_adsorption_v4 import simulate
from feed_composition_v4 import feed_for_batch
from params_columns_v4 import (
    Cin_from_process_row,
    compute_auto_t_end,
    flow_df_from_process,
    get_R003_params,
    get_R004A_params,
)
from parameter_fitting_v4 import apply_fitted_vector, fit_parameters_batches, load_fitted_vector, save_fit_result
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
    ap = argparse.ArgumentParser(description="Fit R004A using R003 effluent as inlet (v4)")
    ap.add_argument("--t_end", type=float, default=None)
    ap.add_argument("--train_batches", type=str, default="6")
    ap.add_argument("--test_batches", type=str, default=None)
    ap.add_argument("--species", type=str, default=None)
    ap.add_argument("--data-suffix", type=str, default="")
    ap.add_argument("--maxiter", type=int, default=15)
    ap.add_argument("--popsize", type=int, default=10)
    args = ap.parse_args()

    r003_fit_path = paths.RESULTS_DIR / "R003" / "latest" / "fitted_vector.json"
    if not r003_fit_path.exists():
        raise SystemExit(f"FileNotFoundError: {r003_fit_path} — run run_fit_R003_v4.py first")

    proc = pd.read_csv(paths.PROCESS_CSV)
    r4a = pd.read_csv(paths.r004a_csv(args.data_suffix))
    train_ids = _parse_batches(args.train_batches)
    test_ids = _parse_batches(args.test_batches) if args.test_batches else list(train_ids)

    fitted_r003 = load_fitted_vector(r003_fit_path)
    p_r003 = apply_fitted_vector(get_R003_params(), fitted_r003)
    p_r4a = get_R004A_params()

    target_species = [s.strip() for s in args.species.split(",")] if args.species else None
    print(f"Fit species : {target_species if target_species else 'auto'}")

    # per-batch t_end
    t_end_map = {}
    for bid in sorted(set(train_ids + test_ids)):
        dfp = proc[proc["batch"] == bid]
        dfo = r4a[r4a["batch"] == bid]
        if args.t_end is not None:
            t_end_map[bid] = float(args.t_end)
        else:
            t_end_map[bid] = compute_auto_t_end(dfp["t_min"], dfo["t_min"])
    t_end_mode = "fixed" if args.t_end is not None else "auto_per_batch"
    t_end_desc = ", ".join(f"b{b}:{t_end_map[b]:.1f}" for b in sorted(t_end_map))

    # 事前に R003 出口を生成
    inlet_cache = {}
    for bid in sorted(set(train_ids + test_ids)):
        dfp = proc[proc["batch"] == bid].copy()
        if dfp.empty:
            continue
        pr = p_r003.copy()
        pr.t_end = t_end_map[bid]
        cin0 = Cin_from_process_row(dfp, pr.species_names)
        cin0.update(feed_for_batch(bid))
        pr.C_in = {**pr.C_in, **cin0}
        flow = flow_df_from_process(dfp)
        _, effluent, _, _ = simulate(pr, v_T_series=flow)
        from competitive_adsorption_v4 import effluent_to_Cin_series

        inlet_cache[bid] = effluent_to_Cin_series(effluent, pr.species_names)

    def inlet_provider(bid, dfb):
        return {"C_in_series": inlet_cache[bid]}

    def flow_provider_common(bid, dfb):
        dfp = proc[proc["batch"] == bid]
        return flow_df_from_process(dfp)

    def t_end_provider(bid, dfb):
        return t_end_map.get(bid)

    train = {b: r4a[r4a["batch"] == b].copy() for b in train_ids if (r4a["batch"] == b).any()}
    test = {b: r4a[r4a["batch"] == b].copy() for b in test_ids if (r4a["batch"] == b).any()}
    if not train:
        raise SystemExit("no R004A train observations")

    print(f"Train batches: {sorted(train)}")
    print(f"t_end mode   : {t_end_mode} ({t_end_desc})")
    fit = fit_parameters_batches(
        p_r4a,
        train,
        target_species=target_species,
        inlet_provider=inlet_provider,
        flow_provider=flow_provider_common,
        t_end_provider=t_end_provider,
        maxiter=args.maxiter,
        popsize=args.popsize,
    )
    print(f"DE done: fun={fit['fun']:.6g}, nit={fit['nit']}, values={fit['values']}")

    out_dir = _dated_dir(paths.RESULTS_DIR / "R004A")
    extra = {
        "column": "R004A",
        "train_batches": train_ids,
        "test_batches": test_ids,
        "based_on_R003": {
            k: fitted_r003.get(k) for k in ("q_total", "eps_b", "k_hyd", "keys", "values")
        },
        "fit_species": target_species if target_species else "auto",
        "fit_species_desc": ",".join(target_species) if target_species else "auto",
        "t_end_mode": t_end_mode,
        "t_end_desc": t_end_desc,
        "data_suffix": args.data_suffix,
    }
    save_fit_result(out_dir, fit, extra=extra)
    (out_dir / "readme.txt").write_text(
        "\n".join(
            [
                "Column: R004A",
                f"Train batches: {train_ids}",
                f"Fit species: {extra['fit_species_desc']}",
                f"t_end: {t_end_desc}",
                f"Objective: {fit['fun']}",
                f"Elapsed[s]: {fit['elapsed_sec']:.1f}",
                f"Values: {fit['values']}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    latest = paths.RESULTS_DIR / "R004A" / "latest"
    if latest.exists():
        shutil.rmtree(latest)
    shutil.copytree(out_dir, latest)
    print(f"Saved: {out_dir}")

    evaluate_model_batches(
        fit["params"],
        test,
        paths.OUTPUTS_DIR / "R004A" / "test",
        inlet_provider=inlet_provider,
        flow_provider=flow_provider_common,
        t_end_provider=t_end_provider,
        make_gif=True,
    )
    print("Evaluation done.")


if __name__ == "__main__":
    main()

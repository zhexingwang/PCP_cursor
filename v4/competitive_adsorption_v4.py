"""競争吸着カラムモデル v4（FAEE 加水分解 + C_in_series）。

論文: Competitive Adsorption Model for Process Design ...
      J. Chem. Eng. Japan, Vol.53 No.9, 2020 の定式化に基づく。
v3 を拡張し、非吸着種 FAEE/ET と液相加水分解 FAEE+OH→FA+ET を追加。
"""
from __future__ import annotations

from dataclasses import dataclass, field, fields
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import copy
import math

import numpy as np
import pandas as pd

SpeciesDict = Dict[str, float]
SeriesSpec = Union[
    pd.DataFrame,
    Dict[str, Tuple[np.ndarray, np.ndarray]],
]


def _default_species() -> List[str]:
    return ["VE1", "VE2", "FA", "OH", "FAEE", "ET"]


def _default_adsorbing() -> List[str]:
    return ["VE1", "VE2", "FA"]


def _default_nonadsorbing() -> List[str]:
    return ["FAEE", "ET"]


def _default_H() -> SpeciesDict:
    return {"VE1": 2.0, "VE2": 2.0, "FA": 2.03, "OH": 5.22, "FAEE": 1.0, "ET": 1.0}


def _default_k_ads() -> SpeciesDict:
    return {"VE1": 0.35, "VE2": 0.25, "FA": 0.80}


def _default_Keq_ads() -> SpeciesDict:
    return {"VE1": 80.0, "VE2": 45.0, "FA": 250.0}


def _default_k_ex() -> Dict[str, float]:
    # FA が VE を置換する交換反応
    return {"FA->VE1": 0.15, "FA->VE2": 0.12}


def _default_Keq_ex() -> Dict[str, float]:
    return {"FA->VE1": 4.0, "FA->VE2": 3.5}


def _default_C_in() -> SpeciesDict:
    return {
        "VE1": 0.02635,
        "VE2": 0.0120,
        "FA": 0.0080,
        "OH": 0.243,
        "FAEE": 0.0,
        "ET": 0.0,
    }


def _default_fit_targets() -> Dict[str, bool]:
    return {
        "q_total": True,
        "eps_b": False,
        "eps_p": False,
        "k_hyd": True,
        "k_ads.VE1": False,
        "k_ads.VE2": False,
        "k_ads.FA": False,
        "Keq_ads.VE1": False,
        "Keq_ads.VE2": False,
        "Keq_ads.FA": False,
        "k_ex.FA->VE1": False,
        "k_ex.FA->VE2": False,
        "Keq_ex.FA->VE1": False,
        "Keq_ex.FA->VE2": False,
        "H.VE1": False,
        "H.VE2": False,
        "H.FA": False,
        "H.OH": False,
    }


@dataclass
class Params:
    # カラム寸法
    d: float = 25.0  # cm
    L: float = 83.2  # cm
    N: int = 24

    # 充填・サイト
    eps_b: float = 0.40
    eps_p: float = 0.55
    q_total: float = 0.25  # mmol/cm3-resin

    # 時間積分
    t_end: float = 350.0  # min
    dt: float = 0.1  # min
    dt_react_cap: float = 0.1  # min（反応サブステップ上限）
    sample_dt: float = 1.0  # min（出力サンプリング間隔）

    # 分散（0 なら移流のみ・風上）
    D_ax: float = 0.5  # cm2/min

    # 化学種
    species_names: List[str] = field(default_factory=_default_species)
    adsorbing_species: List[str] = field(default_factory=_default_adsorbing)
    carrier_species: str = "OH"
    nonadsorbing_species: List[str] = field(default_factory=_default_nonadsorbing)

    H: SpeciesDict = field(default_factory=_default_H)
    k_ads: SpeciesDict = field(default_factory=_default_k_ads)
    Keq_ads: SpeciesDict = field(default_factory=_default_Keq_ads)
    k_ex: Dict[str, float] = field(default_factory=_default_k_ex)
    Keq_ex: Dict[str, float] = field(default_factory=_default_Keq_ex)

    # 入口・流量
    C_in: SpeciesDict = field(default_factory=_default_C_in)
    C_in_series: Optional[SeriesSpec] = None
    v_T: float = 40.0  # L/h
    v_T_series: Optional[Union[pd.DataFrame, Tuple[np.ndarray, np.ndarray]]] = None

    # 初期条件
    q0: Optional[SpeciesDict] = None  # None → 全吸着種 0

    # フィット対象フラグ（フォールバック。通常は params_columns で上書き）
    fit_targets: Dict[str, bool] = field(default_factory=_default_fit_targets)

    # 加水分解速度定数 [cm³/mmol/min]
    k_hyd: float = 0.0

    def copy(self) -> "Params":
        return copy.deepcopy(self)

    def area_cm2(self) -> float:
        return math.pi * (self.d * 0.5) ** 2

    def superficial_u(self, v_T_Lh: float) -> float:
        """体積流量 [L/h] → 空塔速度 [cm/min]。"""
        Q_cm3_min = v_T_Lh * 1000.0 / 60.0
        return Q_cm3_min / max(self.area_cm2(), 1e-12)

    def theta(self, sp: str) -> float:
        """液相ホールドアップ θ_i = ε_b + (1-ε_b) ε_p / H_i。"""
        H = float(self.H.get(sp, 1.0))
        H = max(H, 1e-12)
        return self.eps_b + (1.0 - self.eps_b) * self.eps_p / H

    def alpha(self) -> float:
        return 1.0 - self.eps_b

    def eps_liq(self) -> float:
        return self.eps_b + (1.0 - self.eps_b) * self.eps_p

    # ----- 時系列評価 -----
    def vT_at(self, t: float) -> float:
        ser = self.v_T_series
        if ser is None:
            return float(self.v_T)
        if isinstance(ser, pd.DataFrame):
            tt = ser["t_min"].to_numpy(dtype=float)
            vv = ser["v_T"].to_numpy(dtype=float)
        else:
            tt, vv = np.asarray(ser[0], dtype=float), np.asarray(ser[1], dtype=float)
        if len(tt) == 0:
            return float(self.v_T)
        return float(np.interp(t, tt, vv, left=vv[0], right=vv[-1]))

    def C_in_at(self, t: float) -> SpeciesDict:
        """入口濃度を時刻 t で評価。C_in_series があれば補間、なければ C_in 固定。"""
        out = {sp: float(self.C_in.get(sp, 0.0)) for sp in self.species_names}
        ser = self.C_in_series
        if ser is None:
            return out
        if isinstance(ser, pd.DataFrame):
            tt = ser["t_min"].to_numpy(dtype=float)
            for sp in self.species_names:
                if sp in ser.columns:
                    vv = ser[sp].to_numpy(dtype=float)
                    out[sp] = float(np.interp(t, tt, vv, left=vv[0], right=vv[-1]))
            return out
        # Dict[species, (t, v)]
        for sp, pair in ser.items():
            if sp not in out:
                continue
            tt = np.asarray(pair[0], dtype=float)
            vv = np.asarray(pair[1], dtype=float)
            if len(tt) == 0:
                continue
            out[sp] = float(np.interp(t, tt, vv, left=vv[0], right=vv[-1]))
        return out

    # ----- フィットベクトル -----
    def fit_keys(self) -> List[str]:
        return [k for k, on in self.fit_targets.items() if on]

    def _get_by_key(self, key: str) -> float:
        if key == "q_total":
            return float(self.q_total)
        if key == "eps_b":
            return float(self.eps_b)
        if key == "eps_p":
            return float(self.eps_p)
        if key == "k_hyd":
            return float(self.k_hyd)
        if "." in key:
            group, name = key.split(".", 1)
            d = getattr(self, group)
            return float(d[name])
        raise KeyError(key)

    def _set_by_key(self, key: str, value: float) -> None:
        if key == "q_total":
            self.q_total = float(value)
            return
        if key == "eps_b":
            self.eps_b = float(value)
            return
        if key == "eps_p":
            self.eps_p = float(value)
            return
        if key == "k_hyd":
            self.k_hyd = float(value)
            return
        if "." in key:
            group, name = key.split(".", 1)
            d = getattr(self, group)
            d[name] = float(value)
            return
        raise KeyError(key)

    def to_vector(self) -> np.ndarray:
        return np.array([self._get_by_key(k) for k in self.fit_keys()], dtype=float)

    def update_from_vector(self, vec: Sequence[float]) -> None:
        keys = self.fit_keys()
        if len(vec) != len(keys):
            raise ValueError(f"vector length {len(vec)} != n_keys {len(keys)}")
        for k, v in zip(keys, vec):
            self._set_by_key(k, float(v))


def reaction_rhs_ads(
    C: Dict[str, np.ndarray],
    q: Dict[str, np.ndarray],
    p: Params,
) -> Dict[str, np.ndarray]:
    """固相吸着量の反応速度 dq/dt（セル配列）。"""
    ads = list(p.adsorbing_species)
    n = next(iter(C.values())).shape[0]
    dq = {sp: np.zeros(n, dtype=float) for sp in ads}

    q_sum = np.zeros(n, dtype=float)
    for sp in ads:
        q_sum += q[sp]
    q_free = np.clip(p.q_total - q_sum, 0.0, None)

    Cres = {sp: p.H.get(sp, 1.0) * C[sp] for sp in p.species_names}
    C_OH = Cres[p.carrier_species]

    for sp in ads:
        kf = float(p.k_ads.get(sp, 0.0))
        Keq = max(float(p.Keq_ads.get(sp, 1.0)), 1e-12)
        kr = kf / Keq
        dq[sp] += kf * Cres[sp] * q_free - kr * q[sp] * C_OH

    # 交換反応 FA ↔ VE
    for key, kf in p.k_ex.items():
        if "->" not in key:
            continue
        a, b = key.split("->", 1)  # a displaces b (FA -> VE)
        if a not in ads or b not in ads:
            continue
        Keq = max(float(p.Keq_ex.get(key, 1.0)), 1e-12)
        kr = float(kf) / Keq
        rate = float(kf) * Cres[a] * q[b] - kr * Cres[b] * q[a]
        dq[a] += rate
        dq[b] -= rate

    return dq


def _transport_step(
    C: Dict[str, np.ndarray],
    u: float,
    dz: float,
    dt: float,
    p: Params,
    Cin: SpeciesDict,
    source: Dict[str, np.ndarray],
) -> None:
    """陽的移流（風上）+ 分散 + ソース。C を in-place 更新。"""
    N = p.N
    D = float(p.D_ax)
    for sp in p.species_names:
        c = C[sp]
        th = max(p.theta(sp), 1e-12)
        # 風上移流
        flux_in = np.empty(N, dtype=float)
        flux_in[0] = u * Cin.get(sp, 0.0)
        flux_in[1:] = u * c[:-1]
        flux_out = u * c
        adv = -(flux_out - flux_in) / dz

        disp = np.zeros(N, dtype=float)
        if D > 0.0 and N >= 3:
            disp[1:-1] = D * (c[2:] - 2.0 * c[1:-1] + c[:-2]) / (dz * dz)
            disp[0] = D * (c[1] - c[0]) / (dz * dz)
            disp[-1] = D * (c[-2] - c[-1]) / (dz * dz)

        c += dt * (adv + disp + source[sp]) / th
        np.maximum(c, 0.0, out=c)


def simulate(
    params: Params,
    v_T_series: Optional[Union[pd.DataFrame, Tuple[np.ndarray, np.ndarray]]] = None,
    C_in: Optional[SpeciesDict] = None,
    progress: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, Any], List[str]]:
    """カラム競争吸着シミュレーション。

    Returns
    -------
    grid_df : 時空間プロファイル（t_min, z_cm, C_*, q_*）
    effluent_df : 出口濃度時系列（t_min, C_*_out）
    meta : メタ情報
    fields : フィールド名リスト
    """
    p = params.copy()
    if v_T_series is not None:
        p.v_T_series = v_T_series
    if C_in is not None:
        p.C_in = {**p.C_in, **C_in}

    N = int(p.N)
    dz = p.L / N
    dt = float(p.dt)
    t_end = float(p.t_end)
    ads = list(p.adsorbing_species)
    species = list(p.species_names)
    alpha = p.alpha()

    C = {sp: np.zeros(N, dtype=float) for sp in species}
    q = {sp: np.zeros(N, dtype=float) for sp in ads}
    if p.q0:
        for sp, val in p.q0.items():
            if sp in q:
                q[sp][:] = float(val)

    t_hist: List[float] = []
    C_out_hist: Dict[str, List[float]] = {sp: [] for sp in species}
    # 間引き保存用
    snap_t: List[float] = []
    snap_C: Dict[str, List[np.ndarray]] = {sp: [] for sp in species}
    snap_q: Dict[str, List[np.ndarray]] = {sp: [] for sp in ads}
    next_sample = 0.0
    sample_dt = max(float(p.sample_dt), dt)

    t = 0.0
    # 仕様どおり while t < t_end + 1e-12; 浮動小数で最終点が欠けることあり
    while t < t_end + 1e-12:
        Cin = p.C_in_at(t)
        u = p.superficial_u(p.vT_at(t))

        # 反応サブステップ
        dt_left = dt
        dq_accum = {sp: np.zeros(N, dtype=float) for sp in ads}
        while dt_left > 1e-15:
            dtr = min(dt_left, float(p.dt_react_cap))
            dqdt = reaction_rhs_ads(C, q, p)
            # 加水分解（v3 では k_hyd=0。v4 で有効化）
            r_hyd = float(p.k_hyd) * C.get("FAEE", np.zeros(N)) * C.get("OH", np.zeros(N))
            S_hyd = p.eps_liq() * r_hyd

            for sp in ads:
                q[sp] += dtr * dqdt[sp]
                np.maximum(q[sp], 0.0, out=q[sp])
                dq_accum[sp] += dtr * dqdt[sp]

            # 液相ソース（反応由来）をこのサブステップで即時反映しない。
            # 移流と同時に扱うため、反応による液相変化は平均レートで後段へ。
            # ただし加水分解と吸着による液相消費/生成は反応サブステップでも更新する。
            for sp in ads:
                th = max(p.theta(sp), 1e-12)
                C[sp] -= (alpha / th) * dtr * dqdt[sp]
                np.maximum(C[sp], 0.0, out=C[sp])
            th_oh = max(p.theta(p.carrier_species), 1e-12)
            # 吸着で OH を液相へ放出
            C[p.carrier_species] += (alpha / th_oh) * dtr * sum(dqdt[sp] for sp in ads)
            np.maximum(C[p.carrier_species], 0.0, out=C[p.carrier_species])

            if "FAEE" in C and float(p.k_hyd) != 0.0:
                th_f = max(p.theta("FAEE"), 1e-12)
                th_et = max(p.theta("ET"), 1e-12) if "ET" in C else 1.0
                th_fa = max(p.theta("FA"), 1e-12)
                C["FAEE"] -= dtr * S_hyd / th_f
                C[p.carrier_species] -= dtr * S_hyd / th_oh
                C["FA"] += dtr * S_hyd / th_fa
                if "ET" in C:
                    C["ET"] += dtr * S_hyd / th_et
                for sp in ("FAEE", "ET", "FA", p.carrier_species):
                    if sp in C:
                        np.maximum(C[sp], 0.0, out=C[sp])

            dt_left -= dtr

        # 移流・分散（反応ソースは既に液相へ適用済みなので source=0）
        zero_src = {sp: np.zeros(N, dtype=float) for sp in species}
        _transport_step(C, u, dz, dt, p, Cin, zero_src)

        t += dt

        # 出口記録（毎ステップ）
        if t + 1e-12 >= next_sample or t >= t_end - 1e-12:
            t_hist.append(t)
            for sp in species:
                C_out_hist[sp].append(float(C[sp][-1]))
            snap_t.append(t)
            for sp in species:
                snap_C[sp].append(C[sp].copy())
            for sp in ads:
                snap_q[sp].append(q[sp].copy())
            next_sample += sample_dt

    z = (np.arange(N) + 0.5) * dz
    rows = []
    for i, ti in enumerate(snap_t):
        for j in range(N):
            row = {"t_min": ti, "z_cm": float(z[j])}
            for sp in species:
                row[f"C_{sp}"] = float(snap_C[sp][i][j])
            for sp in ads:
                row[f"q_{sp}"] = float(snap_q[sp][i][j])
            rows.append(row)
    grid_df = pd.DataFrame(rows)

    effluent = {"t_min": t_hist}
    for sp in species:
        effluent[f"C_{sp}_out"] = C_out_hist[sp]
    effluent_df = pd.DataFrame(effluent)

    meta = {
        "L": p.L,
        "d": p.d,
        "N": p.N,
        "t_end": p.t_end,
        "dt": p.dt,
        "q_total": p.q_total,
        "eps_b": p.eps_b,
        "k_hyd": p.k_hyd,
        "species_names": list(species),
        "adsorbing_species": list(ads),
    }
    field_names = ["t_min", "z_cm"] + [f"C_{sp}" for sp in species] + [f"q_{sp}" for sp in ads]
    return grid_df, effluent_df, meta, field_names


def effluent_to_Cin_series(effluent_df: pd.DataFrame, species_names: Sequence[str]) -> pd.DataFrame:
    """出口 DataFrame を後塔の C_in_series（t_min + 成分列）に変換。"""
    out = pd.DataFrame({"t_min": effluent_df["t_min"].to_numpy(dtype=float)})
    for sp in species_names:
        col = f"C_{sp}_out"
        if col in effluent_df.columns:
            out[sp] = effluent_df[col].to_numpy(dtype=float)
        else:
            out[sp] = 0.0
    return out

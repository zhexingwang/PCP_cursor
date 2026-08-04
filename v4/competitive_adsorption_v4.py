"""競争吸着カラムモデル v4（FAEE 加水分解 + C_in_series）。

論文: Hiromori et al., J. Chem. Eng. Japan, 53(9), 477–484 (2020)
      Eqs.(1)–(12), Tables 2–3 に準拠。

v3 からの拡張:
  - 非吸着種 FAEE / ET
  - 液相加水分解 FAEE + OH → FA + ET（k_hyd）
  - 原料水（入口 OH）対応
"""
from __future__ import annotations

from dataclasses import dataclass, field
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

# 論文 Table 3: k [dm3/mol/s] = [cm3/mmol/s] → ×60 で [cm3/mmol/min]
_K1_PER_MIN = 0.973e-2 * 60.0   # 0.5838
_K2_PER_MIN = 11.3e-2 * 60.0    # 6.78
_K3_PER_MIN = 4.25e-2 * 60.0    # 2.55


def _default_species() -> List[str]:
    return ["VE1", "VE2", "FA", "OH", "FAEE", "ET"]


def _default_adsorbing() -> List[str]:
    return ["VE1", "VE2", "FA"]


def _default_nonadsorbing() -> List[str]:
    return ["FAEE", "ET"]


def _default_H() -> SpeciesDict:
    # Table 2 + 非吸着種 H=1
    return {"VE1": 1.96, "VE2": 1.96, "FA": 2.03, "OH": 5.22, "FAEE": 1.0, "ET": 1.0}


def _default_k_ads() -> SpeciesDict:
    # j=1 (VE), j=2 (FA)
    return {"VE1": _K1_PER_MIN, "VE2": _K1_PER_MIN, "FA": _K2_PER_MIN}


def _default_Keq_ads() -> SpeciesDict:
    return {"VE1": 47.4, "VE2": 47.4, "FA": 412.0}


def _default_k_ex() -> Dict[str, float]:
    # j=3: FaH が VE を置換
    return {"FA->VE1": _K3_PER_MIN, "FA->VE2": _K3_PER_MIN}


def _default_Keq_ex() -> Dict[str, float]:
    return {"FA->VE1": 324.0, "FA->VE2": 324.0}


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

    # 充填・サイト（q_total は論文の solid 基準 [mmol/cm3-solid]）
    eps_b: float = 0.40
    eps_p: float = 0.226  # Table 推定値
    q_total: float = 1.10  # Fig.4 飽和付近

    # 時間積分
    t_end: float = 350.0  # min
    dt: float = 0.1  # min
    dt_react_cap: float = 0.05  # min（反応サブステップ上限）
    sample_dt: float = 1.0  # min

    # 粒子径（分散 Em = 2 u Rp）
    R_p: float = 0.023  # cm (= 0.23 mm, Table 1)
    # D_ax > 0 のとき Em を上書き。None/負なら Em=2uRp
    D_ax: Optional[float] = None

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

    # 入口・流量（v_T 入力は L/h）
    C_in: SpeciesDict = field(default_factory=_default_C_in)
    C_in_series: Optional[SeriesSpec] = None
    v_T: float = 40.0  # L/h
    v_T_series: Optional[Union[pd.DataFrame, Tuple[np.ndarray, np.ndarray]]] = None

    q0: Optional[SpeciesDict] = None
    fit_targets: Dict[str, bool] = field(default_factory=_default_fit_targets)

    # 加水分解速度定数 [cm3/mmol/min]
    k_hyd: float = 0.0

    def copy(self) -> "Params":
        return copy.deepcopy(self)

    def area_cm2(self) -> float:
        return math.pi * (self.d * 0.5) ** 2

    def vT_cm3_min(self, v_T_Lh: float) -> float:
        """[L/h] → [cm3/min]。"""
        return float(v_T_Lh) * 1000.0 / 60.0

    def interstitial_u(self, v_T_Lh: float) -> float:
        """空隙内線速度 u [cm/min] = (v_T/A) / ε_b。"""
        u_s = self.vT_cm3_min(v_T_Lh) / max(self.area_cm2(), 1e-12)
        return u_s / max(self.eps_b, 1e-12)

    def theta(self, sp: str) -> float:
        """Eq.(11) 左辺係数: ε_b + (1-ε_b) ε_p H_i。"""
        H = float(self.H.get(sp, 1.0))
        return self.eps_b + (1.0 - self.eps_b) * self.eps_p * H

    def alpha_solid(self) -> float:
        """Eq.(11) 固相係数: (1-ε_b)(1-ε_p)。"""
        return (1.0 - self.eps_b) * (1.0 - self.eps_p)

    def eps_liq(self) -> float:
        return self.eps_b + (1.0 - self.eps_b) * self.eps_p

    def Em(self, u: float) -> float:
        """Eq.(12): Em = 2 u Rp。D_ax 指定時はそれを使用。"""
        if self.D_ax is not None and self.D_ax >= 0:
            return float(self.D_ax)
        return 2.0 * float(u) * float(self.R_p)

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
        for sp, pair in ser.items():
            if sp not in out:
                continue
            tt = np.asarray(pair[0], dtype=float)
            vv = np.asarray(pair[1], dtype=float)
            if len(tt) == 0:
                continue
            out[sp] = float(np.interp(t, tt, vv, left=vv[0], right=vv[-1]))
        return out

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
            return float(getattr(self, group)[name])
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
            getattr(self, group)[name] = float(value)
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
    """論文 Eqs.(8)(9): 固相吸着速度 dq/dt。

    q_OH = q_total - Σ q_ads として従属的に用いる。
    """
    ads = list(p.adsorbing_species)
    n = next(iter(C.values())).shape[0]
    dq = {sp: np.zeros(n, dtype=float) for sp in ads}

    q_sum = np.zeros(n, dtype=float)
    for sp in ads:
        q_sum += q[sp]
    q_OH = np.clip(p.q_total - q_sum, 0.0, None)

    Cres = {sp: float(p.H.get(sp, 1.0)) * C[sp] for sp in p.species_names}
    C_OH = Cres[p.carrier_species]

    # 吸着 j=1,2: iH + S+(OH-) ⇌ S+(i-) + H2O
    for sp in ads:
        kf = float(p.k_ads.get(sp, 0.0))
        Keq = max(float(p.Keq_ads.get(sp, 1.0)), 1e-12)
        kr = kf / Keq
        dq[sp] += kf * q_OH * Cres[sp] - kr * q[sp] * C_OH

    # 交換 j=3: FaH + S+(VE-) ⇌ S+(Fa-) + VEH
    for key, kf in p.k_ex.items():
        if "->" not in key:
            continue
        a, b = key.split("->", 1)  # a=FA displaces b=VE
        if a not in ads or b not in ads:
            continue
        Keq = max(float(p.Keq_ex.get(key, 1.0)), 1e-12)
        kr = float(kf) / Keq
        rate = float(kf) * q[b] * Cres[a] - kr * q[a] * Cres[b]
        dq[a] += rate
        dq[b] -= rate

    return dq


def simulate(
    params: Params,
    v_T_series: Optional[Union[pd.DataFrame, Tuple[np.ndarray, np.ndarray]]] = None,
    C_in: Optional[SpeciesDict] = None,
    progress: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, Any], List[str]]:
    """カラム競争吸着シミュレーション（論文 Eq.(11)）。

    Returns
    -------
    grid_df, effluent_df, meta, fields
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
    alpha_s = p.alpha_solid()
    carrier = p.carrier_species

    C = {sp: np.zeros(N, dtype=float) for sp in species}
    q = {sp: np.zeros(N, dtype=float) for sp in ads}
    if p.q0:
        for sp, val in p.q0.items():
            if sp in q:
                q[sp][:] = float(val)

    t_hist: List[float] = []
    C_out_hist: Dict[str, List[float]] = {sp: [] for sp in species}
    snap_t: List[float] = []
    snap_C: Dict[str, List[np.ndarray]] = {sp: [] for sp in species}
    snap_q: Dict[str, List[np.ndarray]] = {sp: [] for sp in ads}
    next_sample = 0.0
    sample_dt = max(float(p.sample_dt), dt)

    t = 0.0
    while t < t_end + 1e-12:
        Cin = p.C_in_at(t)
        u = p.interstitial_u(p.vT_at(t))
        Em = p.Em(u)
        # 対流係数: Eq.(11) の -u ε_b ∂C/∂x
        u_eps = u * p.eps_b
        Em_eps = Em * p.eps_b

        dt_left = dt
        while dt_left > 1e-15:
            dtr = min(dt_left, float(p.dt_react_cap))
            dqdt = reaction_rhs_ads(C, q, p)

            # 固相更新
            for sp in ads:
                q[sp] += dtr * dqdt[sp]
                np.maximum(q[sp], 0.0, out=q[sp])
            # サイト総量を超えないようクリップ
            q_sum = np.zeros(N, dtype=float)
            for sp in ads:
                q_sum += q[sp]
            over = q_sum > p.q_total
            if np.any(over):
                scale = np.ones(N, dtype=float)
                scale[over] = p.q_total / np.maximum(q_sum[over], 1e-15)
                for sp in ads:
                    q[sp] *= scale

            # 固相速度（クリップ後の実効値で液相ソースを組む）
            dq_liq = {sp: dqdt[sp].copy() for sp in ads}
            dq_OH = -sum(dq_liq[sp] for sp in ads)

            # 加水分解（v4）: FAEE + OH → FA + ET
            r_hyd = float(p.k_hyd) * C.get("FAEE", np.zeros(N)) * C.get(carrier, np.zeros(N))
            S_hyd = p.eps_liq() * r_hyd

            # 移流（風上）+ 分散 + 固相ソース（Eq.11）
            for sp in species:
                c = C[sp]
                th = max(p.theta(sp), 1e-12)
                flux_in = np.empty(N, dtype=float)
                flux_in[0] = u_eps * Cin.get(sp, 0.0)
                flux_in[1:] = u_eps * c[:-1]
                flux_out = u_eps * c
                adv = -(flux_out - flux_in) / dz

                disp = np.zeros(N, dtype=float)
                if Em_eps > 0.0 and N >= 3:
                    disp[1:-1] = Em_eps * (c[2:] - 2.0 * c[1:-1] + c[:-2]) / (dz * dz)
                    disp[0] = Em_eps * (c[1] - c[0]) / (dz * dz)
                    disp[-1] = Em_eps * (c[-2] - c[-1]) / (dz * dz)

                if sp in ads:
                    solid_src = -alpha_s * dq_liq[sp]
                elif sp == carrier:
                    solid_src = -alpha_s * dq_OH
                else:
                    solid_src = 0.0

                # 加水分解ソース（仕様: θ ∂C/∂t = … − S_i）
                hyd = 0.0
                if float(p.k_hyd) != 0.0:
                    if sp == "FAEE" or sp == carrier:
                        hyd = -S_hyd
                    elif sp == "FA" or sp == "ET":
                        hyd = +S_hyd

                c += dtr * (adv + disp + solid_src + hyd) / th
                np.maximum(c, 0.0, out=c)

            dt_left -= dtr

        t += dt
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
        "eps_p": p.eps_p,
        "k_hyd": p.k_hyd,
        "species_names": list(species),
        "adsorbing_species": list(ads),
        "paper": "Hiromori et al., JCEJ 53(9) 2020",
    }
    field_names = ["t_min", "z_cm"] + [f"C_{sp}" for sp in species] + [f"q_{sp}" for sp in ads]
    return grid_df, effluent_df, meta, field_names


def effluent_to_Cin_series(effluent_df: pd.DataFrame, species_names: Sequence[str]) -> pd.DataFrame:
    out = pd.DataFrame({"t_min": effluent_df["t_min"].to_numpy(dtype=float)})
    for sp in species_names:
        col = f"C_{sp}_out"
        out[sp] = effluent_df[col].to_numpy(dtype=float) if col in effluent_df.columns else 0.0
    return out

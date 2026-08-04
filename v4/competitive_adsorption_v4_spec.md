# 競争吸着パイプライン v4 仕様書（FAEE 加水分解モデルの追加）

本仕様書は v3 の全機能・付録を取り込んだ自己完結版である。各節の見出しに出自タグを付けて v3 と v4 を区別する。

| タグ | 意味 |
|:---|:---|
| **【v4新規】** | v4 で新設した節 |
| **【v4変更】** | v3 の内容をベースに v4 で修正が入った節（変更点は節内に明記） |
| **【v3継承】** | v3 から内容そのまま（ファイル名等の `_v3` → `_v4` 置換のみ） |

---

## §0 v4 の位置づけ 【v4新規】

v2 までの仕様（§1〜§19 相当）は v2 フォルダ内の `competitive_adsorption_pipeline_spec.md` に記述されており、**引き続き有効**である。v3 の全仕様は本書に取り込み済みのため、v3 仕様書を別途参照する必要はない（v3 フォルダは凍結し、比較用ベースラインとして保持する）。

| 項目 | v3 | v4 |
|:-----|:---|:---|
| 化学種 | VE1, VE2, FA, OH（キャリア） | **+ FAEE, ET（非吸着種・液相のみ）** |
| 液相反応 | なし | **FAEE 加水分解（液相均一・不可逆）: FAEE + H₂O → FA + ET** |
| 原料中の水 | 考慮なし（入口 OH = 0） | **入口 OH = 0.243 mol/L（原料含有水）** |
| 原料組成の入力 | — | **固定デフォルト + `data/feed_composition_v4.csv` による外部入力（CSV 優先）** |
| フィット対象 | q_total のみ（eps_b は 0.40 固定） | **q_total + k_hyd**（k_hyd は R003 / R004A で独立にフィット。eps_b は 0.40 固定を継続） |
| 後方互換 | — | **k_hyd=0・原料 FAEE/ET/水=0 のとき v3 と出口濃度が完全一致**（付録 F.1 で検証済み） |
| フォルダ分離 | v2 を保持し v3 を別フォルダで開発 | **v3 を完全に保持し、v4 は別フォルダで開発** |

v3 はリファレンス実装として保持し、v4 で FAEE 加水分解による VE ピーク精度改善を狙う。

---

## §1 背景と課題（v4 開発の動機） 【v4新規】

### 1.1 観測された問題

ESX004 / ESX006（batch 4, 6）で R003・R004A をフィット（フィット対象は q_total のみ。R004A は VE2 のみ、R003 は全成分）し、ESX012（batch 12）を予測した結果（v3、通液時間をデータ上の流量存在範囲より約 2 時間延長して予測）:

- **R003 の VE1・VE2 出口濃度のピーク付近が乖離**する（モデルの山が低く・なだらか）。
- **R004A の VE1 もピークが乖離**する。
- R004A の VE1・VE2 の山はモデルではフィード濃度付近で頭打ちの平坦な形になるが、**現場の経験では濃度の山はもっとシャープ**であり、平坦な形は違和感がある。

該当図: `v3/outputs_chain_v3/chain_2026_7_ESX004_006_pre012_*/batch_12/chain_overlay.png`

### 1.2 課題の本質

VE のピークが平ら・低いのは、**FA による VE の置換（ロールアップ／オーバーシュート）の駆動力が不足している**ためと考えられる。現行モデルで VE を脱着させる経路は交換反応 (FA→VE1), (FA→VE2) のみであり、その駆動力は液相 FA 濃度である。フィード FA だけでは置換が足りない。

### 1.3 傍証

R003 の FA 出口観測は 180〜360 min の時間帯で既に約 0.004 mmol/cm³ 出ているのに対し、v3 モデルは 400 min 過ぎまでほぼ 0 である。この「FA の早期リーク」は、**塔内で FA が途中生成している**という仮説（§2）と整合的である。

---

## §2 仮説と追加モデル 【v4新規】

### 2.1 仮説

吸着原料液中の FAEE（脂肪酸エチルエステル）が、水と反応して遊離脂肪酸（FFA = 成分 FA）とエタノール（ET）を生成する（加水分解。これまで未知でありモデルに入っていなかった）:

```
FAEE + H₂O → FA + ET
```

- 水は「吸着反応で生成する水（= 既存キャリア OH、§4.4 参照）」と「原料中の含有水」の両方。
- 生成した FA が交換反応を介して吸着済み VE を追い出し、VE ピークをシャープに・高くする。

**メカニズムの因果の鎖:**

1. VE・FA が吸着 → その場で水（OH）が液相に放出される
2. 放出された水（+ 原料含有水）がフィード中の FAEE と反応し FA を生成
3. 生成した FA が交換反応 k_ex(FA→VE) を介して吸着済み VE を追い出す
4. VE 出口濃度がフィード濃度を超えて鋭く立ち上がる（オーバーシュートの強化）
5. 副次効果として FA 出口が中間時間帯から底上げされる（§1.3 の早期リークと整合するか検証可能）

水の生成場所が吸着帯であるため、FA の生成が VE 吸着帯を後ろから追いかける形になり、まさにオーバーシュートを強める方向に働く。

### 2.2 反応速度式（Phase 1: 不可逆）

```
r_hyd = k_hyd · C_FAEE · C_OH    [mmol/cm³-liquid/min]
```

- `k_hyd` [cm³/mmol/min] をフィット対象に追加（1 パラメータのみ）。
- **反応の場はバルク液均一反応**とする。理由: FAEE は OH 形イオン交換樹脂上での触媒反応が考えにくい。反応の空間的局在性は水濃度分布（吸着帯で生成）が自然に与える。
- 平衡定数は不要（不可逆近似）。エステル加水分解の平衡は EtOH 過剰下で左に押されるため、フィットされる `k_hyd` は平衡抑制込みの**見かけの速度定数**である。原料中の ET・水はバッチ間でほぼ一定のため、この近似はピーク形状再現の目的には十分機能する見込み。

### 2.3 収支への組み込み（実装形）

反応はバルク液＋粒子内液の全液相で起きるとし、ベッド体積あたりのソースを

```
S_hyd = ε_liq · r_hyd ,    ε_liq = ε_b + (1 − ε_b) ε_p
```

として、**同一の `S_hyd`** を各成分の液相収支 `θᵢ ∂Cᵢ/∂t = ⋯ − Sᵢ` に共通で与える:

| 成分 | ソース `Sᵢ` への追加 | 意味 |
|:---|:---|:---|
| FAEE | `+S_hyd` | 消費 |
| OH（水） | `+S_hyd`（既存の吸着由来放出 `−α Σ dq/dt` に加算） | 消費 |
| FA | `−S_hyd`（既存の吸着消費 `+α dq_FA/dt` に加算） | 生成 |
| ET | `−S_hyd` | 生成 |

同一の `S_hyd` を使うため **FAEE 消費量 = FA 生成量 = ET 生成量 = OH 消費量（モル保存）** が成立する（付録 F.3 で検証済み）。実装は `competitive_adsorption_v4.simulate()` 内。

### 2.4 Phase 2（可逆化・将来拡張）

Phase 1 で FA の生成が過剰・VE ピークが行き過ぎ等の兆候が出た場合は、逆反応（再エステル化）を追加する:

```
r_hyd = k_hyd ( C_FAEE · C_OH − (1/K_hyd) · C_FA · C_ET )
```

ET は v4 で既に化学種として実装済みのため、逆反応項と `K_hyd`（フィット対象 +1）を追加するだけで切り替えられる。

### 2.5 留意点

- **原料水の追加自体が既存挙動を変える**: 脱着の逆反応は液相 OH 濃度（樹脂相換算 H=5.22 倍）に駆動されるため、入口水 0.243 mol/L は加水分解と独立に破過を早める方向に働く。したがって **q_total と k_hyd は必ず同時にフィット**する（v3 の q_total を流用しない）。
- **識別可能性**: k_hyd は既存の k_ex, Keq_ex（固定値）とトレードオフになり得る。v4 では既存パラメータは固定のまま k_hyd のみ追加する。
- **数値安定性**: 陽的時間積分のため k_hyd が大きいと dt 制約が厳しくなる。DE 境界の上限（5.0）は dt_react_cap=0.1 を前提としており、フィット値が上限に張り付く場合は dt_react_cap を縮めて境界を緩める。
- **水消費の副作用**: 加水分解が水を消費すると脱着逆反応が弱まり吸着量が増える方向の副作用がある。挙動が不自然な場合は「水は過剰でほぼ消費されない」近似（OH を消費しない）への切替を検討する。

---

## §3 原料組成データ 【v4新規】

### 3.1 確定値（2026-07-29 打合せ）

| ロット | batch | FAEE [mol/L] | ET [mol/L] | 水 [mol/L] |
|:---|:---:|---:|---:|---:|
| ESX004 | 4 | 0.999 | 0.999 | 0.243 |
| ESX006 | 6 | 0.665 | 0.665 | 0.243 |
| ESX012 | 12 | 0.774 | 0.774 | 0.243 |

- FAEE と ET は各ロット同値（近似値）。
- 水は全ロット共通の近似値 0.243 mol/L（分析値は収集中。揃い次第 CSV を更新する）。
- 単位は mol/L = mmol/cm³（モデル内部単位と同じ。換算不要）。

### 3.2 固定値 / 外部入力の切替仕様

`feed_composition_v4.py` が管理する。値の優先順位:

1. `data/feed_composition_v4.csv` に該当バッチの行があれば **CSV 値（外部入力）**
2. 無ければ `FEED_FIXED_DEFAULTS`（固定近似値）: FAEE=0, ET=0, OH=0.243

- CSV に行が無いバッチ（例 batch 2, 3）は FAEE=ET=0（加水分解が起きない安全側）となり、**警告を表示**する。
- 分析値が揃ったら CSV を更新するだけで反映される（コード変更不要）。
- 各実行スクリプトは `feed_for_batch(batch_id)` の戻り値で R003 入口濃度 `C_in` を上書きする（`feed_provider` 機構、§5.1 参照）。プロセス CSV（`process_runs_batches.csv`）の C0_* 列より優先される（特に C0_OH=0 は 0.243 に上書きされる）。

### 3.3 CSV 書式

```
batch,C0_FAEE,C0_ET,C0_OH
4,0.999,0.999,0.243
6,0.665,0.665,0.243
12,0.774,0.774,0.243
```

`#` 始まりの行はコメント。列が欠けている成分は固定デフォルトにフォールバックする。

---

## §4 対象プロセスとモデル思想 【v4変更】

### 4.1 装置構成（v3 継承）

```
原料 → [R003: d=25, L=83.2] → [R004A: d=25, L=42.0] → [R004B: d=25, L=42.0] → 製品
       ↑                       ↑                       ↑
       既存モデル化            新規モデル化            新規モデル化（予測のみ）
```

- 全塔とも **同一樹脂**（d=25 cm 共通）
- R003 と R004A/R004B では **樹脂量（カラム長）が異なる**
- 直列配管・**v_T は全塔で共通**
- 各塔とも吸着工程前に洗浄工程があり、**初期吸着量 q0 = 0**

### 4.2 化学種一覧（v4 で拡張）

| 種 | 分類 | 液相 C | 固相 q | H | 入口濃度の由来 |
|:---|:---|:---:|:---:|---:|:---|
| VE1 | 吸着種（Toc-α + T3-α） | ○ | ○ | 2.0 | プロセス CSV `C0_VE1` |
| VE2 | 吸着種（Toc-β + T3-β） | ○ | ○ | 2.0 | プロセス CSV `C0_VE2` |
| FA | 吸着種（遊離脂肪酸） | ○ | ○ | 2.03 | プロセス CSV `C0_FA` |
| OH | キャリア（= 論文の H₂O 役、§4.4） | ○ | q_free として従属 | 5.22 | **v4: 原料組成（水 0.243）** |
| **FAEE** | **非吸着種（v4 新規）** | ○ | — | 1.0（固定） | **原料組成 CSV / 固定値** |
| **ET** | **非吸着種（v4 新規）** | ○ | — | 1.0（固定） | **原料組成 CSV / 固定値** |

非吸着種は `Params.nonadsorbing_species` で定義され、吸着・交換反応・サイト収支（q_free）に関与しない。移流・分散・加水分解反応のみを受ける。

### 4.3 塔ごとの劣化想定・フィット方針（v4 変更）

樹脂の種類は同じだが、塔ごとに通液履歴・劣化具合が異なる可能性があるため、**q_total は塔ごとに独立にフィット**する。v4 では加えて **k_hyd（加水分解速度定数）も塔ごとに独立にフィット**する（運用上の決定。物理的には同一液・同一温度なら同じ値のはずなので、両塔の値の整合性を妥当性チェックに使う。§12.3）。他のパラメータ（H, k_ads, Keq_ads, k_ex, Keq_ex）は樹脂物性として共通の固定値、eps_b は 0.40 固定。

R004B はフィットせず予測のみ（v3 継承）: R004B の出口は VE 漏れがほぼ無いことが運用上の前提であり、観測値も「ロットごと 1 点・ほぼゼロ」となるためフィットには不向き。R004B は **R004A と同じパラメータ（q_total, k_hyd）を使用して予測のみ**を行う。

### 4.4 キャリア成分 OH の定義と役割（v3 継承 + v4 追記）

本モデルは陰イオン交換樹脂（Diaion PA306S、活性サイト = 水酸化物イオン形）を対象とした論文
（*Competitive Adsorption Model for Process Design ...*, J. Chem. Eng. Japan, Vol.53 No.9, 2020）
の定式化に基づく。成分 `OH` は他成分（VE1/VE2/FA）とは性質が異なる **キャリア（carrier）成分** であり、以下のように定義・実装されている。

**イオン交換反応（論文 Eqs.(2)–(4)）**
- (2) VEH + S⁺(OH⁻) ⇌ S⁺(VE⁻) + H₂O
- (3) FaH + S⁺(OH⁻) ⇌ S⁺(Fa⁻) + H₂O
- (4) FaH + S⁺(VE⁻) ⇌ S⁺(Fa⁻) + VEH

サイト S⁺ は常に OH⁻ / VE⁻ / Fa⁻ のいずれか 1 つを対イオンとして保持する。
交換の量論は **1:1**（VE/Fa が 1 吸着すると OH⁻ が 1 放出される）。

**固相側（サイト占有）**
- `S⁺(OH⁻)` = まだ標的成分が付いていない **空き（基準）サイト**。VE/Fa が交換で入っていく先。
- 論文 Eq.(10) `q_total = q_OH + q_VE + q_Fa` に対応。
- コードでは `q_OH` を独立変数として持たず、**`q_free = q_total − q_VE − q_Fa`（= q_OH）** として従属的に計算する（`reaction_rhs_*` 内の `q_free`）。
- OH は専用の固相負荷 `q[OH]` や吸着定数を持たず、フィット対象にもならない。

**液相側（C_OH）**
- コードの液相 `OH` は **樹脂液相の OH⁻ 濃度**（論文 C_OH）に対応する。
- 反応(2)(3)で VE/Fa が吸着すると **1:1 で OH⁻ が液相へ放出**される（`simulate` 内のキャリアソース項 `src_carrier = -alpha · Σ dq`）。これは論文の「生成物（H₂O）＝逆反応を駆動する項」に対応する。
- 放出で増えた液相 OH⁻（`Cres[OH] = H[OH]·C[OH]`、既定 `H[OH]=5.22`）が、**逆反応（脱着）項 `kr · q[i] · Cres[OH]` を駆動**する。

**v4 追記: OH = 水としての役割拡張**
- 加水分解仮説（§2）における「水」は、このキャリア OH をそのまま流用する（独立の水種は追加しない）。
- v4 では OH に 2 つの供給源と 1 つの消費源が加わる:
  - 供給1: 吸着に伴う塔内生成（v3 と同じ）
  - 供給2: **原料含有水（入口 C_in[OH] = 0.243 mol/L、v4 新規）**
  - 消費: **加水分解 FAEE + OH → FA + ET（v4 新規）**
- 初期値の扱い（v3 継承）: 液相初期濃度は全成分ゼロ（案A）で開始。OH は「入口からの流入 + 吸着に伴う塔内生成」で供給される。

**一言まとめ**
> OH は「固相では空きサイトの基準状態（`q_free` で従属）」「液相では脱着を駆動する生成物 = 水」であり、v4 ではさらに「原料から流入し、加水分解で消費される水」の役割を持つ。

---

## §5 ワークフロー 【v4変更】

```
Step 1: R003 フィット
  入口: 実験 C0（process_runs_batches.csv）+ 原料組成（feed_composition_v4）★v4
  出口: R003 観測（密）
  対象: q_total, k_hyd（fit_targets=True）★v4: eps_b は 0.40 固定
  出力: results_v4/R003/<date>/fitted_vector.json

Step 2: R004A フィット（R003 のパラメータ固定）
  R003 のフィット結果（q_total, k_hyd）をロード → R003 を simulate
  R003 出口時系列（FAEE / ET / OH 含む全成分）を R004A の C_in_series として渡す ★v4
  対象: R004A の q_total, k_hyd（R003 とは独立の値）★v4
  出口観測: R004A 実測（疎 3〜4 点）
  出力: results_v4/R004A/<date>/fitted_vector.json

Step 3: 3 塔チェーン予測
  R003・R004A のフィット済みパラメータ（q_total, k_hyd）をロード
  R004B は R004A と同じパラメータを使用（予測専用）
  3 塔を順に simulate して各塔の出口を出力
  出力: outputs_chain_v4/<label>/batch_X/{R003,R004A,R004B}/

Step 4: 妥当性検証（手動）
  - R003・R004A の q_total・k_hyd を比較（樹脂・液が同じはずなのに大きく違えば再検討）
  - VE ピークの高さ・位置・シャープさの改善を確認 ★v4
  - FA 早期リーク（180〜360 min）の再現を確認 ★v4
  - R004B 予測の VE 出口濃度が「ほぼゼロ」になっていることを確認
```

### 5.1 feed_provider 機構（v4 新規）

フィット・評価の内部では、バッチごとに `feed_provider(batch_id, dfb) -> Dict[species, C_in]` コールバックが呼ばれ、C0_* 列から作った `C_in` を原料組成で上書きする。

- `parameter_fitting_v4.fit_parameters_batches(..., feed_provider=...)`
- `evaluate_model_accuracy_v4.evaluate_model_batches(..., feed_provider=...)`
- プロセス CSV に C0_FAEE / C0_ET 列が無くても、feed_provider が値を返せばエラーにならない。
- `inlet_provider`（時系列入口）を使うバッチ（R004A フィット等）には適用されない（後段塔の入口は前段塔出口で決まるため）。

### 5.2 切り分けケース（v4 新規・推奨）

加水分解の効果と原料水の効果を分離するため、必要に応じて次の中間ケースを実行する:

| ケース | 設定 | 目的 |
|:---|:---|:---|
| A: v3 相当 | k_hyd=0, 原料 FAEE/ET/水=0 | ベースライン（v3 と完全一致） |
| B: 水のみ | k_hyd=0, 原料水=0.243 | 原料水単独の影響（破過の早まり）を確認 |
| C: v4 フル | k_hyd フィット値, 原料組成フル | 本番 |

ケース A/B は `feed_composition_v4.csv` の値を一時的に 0 にする（または CSV の該当行を削除して固定デフォルトを編集する）ことで実行できる。

---

## §6 入口濃度のリレー方式（B1 採用） 【v3継承】

R004A の入口は、**R003 モデルの出口時系列**（実測値ではなくモデル予測）を用いる。理由:

- R003 のフィット結果を「予測モデル」として確定させる目的のため、その出力を素直に下流に渡すのが自然
- 実運用時には実測値が存在しない時間帯まで予測する必要があるため、モデル出口を用いる方式が一貫している

R003 のフィット精度が高いため、R004A 入口に渡される濃度は実測と概ね一致するという前提で設計する。万一 R003 モデルに大きな誤差がある場合は、R004A のフィット品質に影響が出るため、Step 4 の整合性チェックで検出する。

**v4 注記**: 新種 FAEE / ET / OH も effluent の `C_*_out` 列として自動的にリレーされる（`column_chain_v4._effluent_to_Cin_series` は `species_names` の全成分を変換する）。

---

## §7 ファイル構成 【v4変更】

### 7.1 v4 フォルダの構造

```
v4/
├── competitive_adsorption_v4_spec.md           本ドキュメント
│
├── competitive_adsorption_v4.py                v3 をコピーして拡張（FAEE/ET 種・加水分解・k_hyd）
├── parameter_fitting_v4.py                     v3 をコピーして拡張（feed_provider・k_hyd 境界）
├── evaluate_model_accuracy_v4.py               v3 をコピーして拡張（feed_provider）
├── params_columns_v4.py                        v3 をコピーして拡張（fit_targets = q_total + k_hyd）
├── column_chain_v4.py                          v3 をコピー（機能変更なし）
├── feed_composition_v4.py                      ★新規: 原料組成（FAEE/ET/水）管理
├── paths_v4.py                                 v3 をコピーして拡張（FEED_COMPOSITION_CSV 追加）
│
├── run_fit_R003_v4.py                          R003 単塔フィット（q_total + k_hyd）
├── run_fit_R004A_v4.py                         R003 結果を用いた R004A フィット
├── run_chain_predict_v4.py                     3 塔チェーン予測（--plot_new_species 追加）
├── run_chain_visualize_v4.py                   保存済みパラメータで再可視化
├── smoke_test_v4.py                            ★新規: 動作検証スクリプト（付録 F）
│
├── data/
│   ├── process_runs_batches.csv                R003 入口・出口（v3 からコピー）
│   ├── feed_composition_v4.csv                 ★新規: ロット別 FAEE / ET / 水 濃度
│   ├── r004a_observed_batches.csv              R004A 出口観測（v3 からコピー）
│   ├── r004a_observed_batches_test.csv         R004A 人工テストデータ
│   ├── r004b_observed_batches.csv              R004B 出口観測（v3 からコピー）
│   └── r004b_observed_batches_test.csv         R004B 人工テストデータ
│
├── results_v4/
│   ├── R003/<date>/                            R003 フィット結果
│   ├── R004A/<date>/                           R004A フィット結果
│   └── R004B/latest/                           R004B 使用パラメータ記録（フィットなし）
│
└── outputs_chain_v4/
    ├── R003/, R004A/                           フィット時の train/test 評価出力
    └── <label>/batch_X/{R003,R004A,R004B}/     チェーン予測の overlay/gif/CSV
```

### 7.2 v3 ファイルとの対応

| v3 ファイル | v4 での扱い |
|:-----------|:-----------|
| `competitive_adsorption_v3.py` | コピー → `competitive_adsorption_v4.py`（FAEE/ET・加水分解・k_hyd 拡張） |
| `parameter_fitting_v3.py` | コピー → `parameter_fitting_v4.py`（feed_provider・k_hyd 境界） |
| `evaluate_model_accuracy_v3.py` | コピー → `evaluate_model_accuracy_v4.py`（feed_provider） |
| `params_columns_v3.py` | コピー → `params_columns_v4.py`（`_set_fit_targets_q_khyd`） |
| `column_chain_v3.py` | コピー → `column_chain_v4.py`（機能変更なし） |
| `paths_v3.py` | コピー → `paths_v4.py`（FEED_COMPOSITION_CSV 追加） |
| `run_fit_R003_v3.py` 他 run 系 | コピー → `run_*_v4.py`（feed 注入・k_hyd 対応） |
| （対応なし） | `feed_composition_v4.py`, `smoke_test_v4.py`, `data/feed_composition_v4.csv` は新規 |

v3 のファイルは**一切変更しない**。v4 はあくまで独立したコピーから派生する。

---

## §8 競争吸着モデル本体の拡張（C_in_series） 【v3継承】

### 8.1 v2 の `Params` 入口濃度

```python
# v2: 固定値のみ
C_in: Dict[str, float] = field(default_factory=lambda: {"VE1":..., "VE2":..., "FA":...})
```

### 8.2 v3 での拡張（v4 でも同じ）

`v_T_series` と同じ流儀で、時系列入口濃度を受けられる。

```python
C_in: Dict[str, float] = ...                    # 固定値
C_in_series: Optional[Dict[str, ...]] = None    # 時系列
```

`C_in_series` が指定されている成分は各時刻 t で補間して使用、未指定の成分は `C_in` 固定値を使用。

### 8.3 内部 API

```python
def C_in_at(t: float) -> Dict[str, float]:
    """入口濃度を時刻 t で評価する。C_in_series があれば補間、なければ C_in 固定値を返す。"""
```

`simulate()` 内で入口境界条件を組む際に毎ステップ評価される。

---

## §9 塔別パラメータ生成（params_columns_v4.py） 【v4変更】

```python
def get_R003_params() -> Params:   # d=25.0, L=83.2
def get_R004A_params() -> Params:  # d=25.0, L=42.0
def get_R004B_params() -> Params:  # d=25.0, L=42.0
```

樹脂物性（H, k_ads, Keq_ads, k_ex, Keq_ex）は全塔で同一の `Params` デフォルトを採用。寸法のみ塔ごとに上書きする。`q_total`, `k_hyd` はフィットで決まるため、生成時はデフォルト値を入れておき、フィット結果を後から上書きする。

> **フィット対象（fit_targets）の決定場所【重要】**
> 実際のフィット対象は、`get_*_params()` 内で呼ばれる `params_columns_v4._set_fit_targets_q_khyd()` が決定する。この関数が `p.fit_targets` を丸ごと作り直し、**`q_total` と `k_hyd` の 2 つだけを True**（他はすべて False。eps_b も False = 0.40 固定）に設定する。
>
> ```python
> def _set_fit_targets_q_khyd(p: Params) -> None:
>     targets = {k: False for k in p.fit_targets}
>     targets["q_total"] = True
>     targets["k_hyd"] = True
>     targets["eps_b"] = False
>     p.fit_targets = targets
> ```
>
> `competitive_adsorption_v4.py` の `Params.fit_targets` デフォルト値は「`Params()` を単体生成し、上書きも通さない場合」のみ効くフォールバックであり、通常フローでは**必ず上書きされる**。フィット対象を変えたい場合はこの関数（または呼び出し側で `p.fit_targets` を上書き）を編集すること。

**k_hyd の DE 境界**（`parameter_fitting_v4.DEFAULT_BOUNDS`）: `(0.0, 0.05)` [cm³/mmol/min]。VE ピーク高さ・時刻への感度が大きく、運用上は 0 近傍が妥当なため上限を小さく設定。フィット値が上限に張り付く場合のみ上限を緩める。

---

## §10 3 塔リレー実行コア（column_chain_v4.py） 【v3継承】

### 10.1 アルゴリズム

```
入力: params_list = [P_R003, P_R004A, P_R004B]
      v_T_series  = 全塔共通の流量時系列
      C_in_first  = R003 の入口濃度（実験 C0 + 原料組成）

Step 1: P_R003 を simulate → R003 出口濃度時系列 = effluent_R003
Step 2: P_R004A の C_in_series = effluent_R003 → simulate → effluent_R004A
Step 3: P_R004B の C_in_series = effluent_R004A → simulate → effluent_R004B

戻り値: 各塔の (grid_df, effluent_df, meta, fields) のリスト
```

各塔のシミュレーションは独立に行うため、塔ごとに `t_end`, `dt`, `N` を変えても良い。デフォルトは全塔で揃える。

### 10.2 入力モードの選択

```python
simulate_chain(params_list, ..., n_cols=3)  # 3 塔通し
simulate_chain(params_list[:2], ...)         # R003+R004A まで
simulate_chain(params_list[:1], ...)         # R003 のみ
```

中間で打ち切れるため、フィットの各段階で必要分だけ実行可能。

---

## §11 出力構造 【v4変更】

### 11.1 各塔の出力（v3 継承）

```
outputs_chain_v4/<label>/batch_X/{R003,R004A,R004B}/
  ├── overlay_VE1.png         実測 vs モデル（R003 密、R004A 疎、R004B 1 点）
  ├── overlay_VE2.png         + 原料濃度の水平点線（feed Cin）も併記 (§11.3)
  ├── overlay_FA.png
  ├── adsorption_profile.gif  カラム内プロファイル動画（吸着種のみ。FAEE/ET は非表示）
  ├── model_output.csv        モデル出口濃度時系列（v4: FAEE / ET / OH 列も含む）
  └── observed.csv            実測値（あれば）
```

R004B はフィットしないが、出力構造は他塔と揃える。

### 11.2 チェーン全体の集約（v3 継承 + v4 拡張）

```
outputs_chain_v4/<label>/batch_X/
  ├── chain_summary.csv               各塔の RMSE まとめ（全成分）
  ├── chain_overlay.png               3 塔出口濃度の一覧オーバーレイ図
  └── chain_adsorption_profile.gif    3 塔のカラム内プロファイルを縦並びで表示する動画
```

- `chain_overlay.png` は 3 塔 × 3 成分（9 サブプロット）の俯瞰図。**v4: `--plot_new_species` 指定時は FAEE / ET / OH のパネルを追加し 3 塔 × 6 成分（18 サブプロット）になる**（新種は観測が無いためモデル線と feed 点線のみ）。
- `chain_adsorption_profile.gif` は 3 塔の吸着プロファイルを縦並び・**全塔同一 snap_time** でまとめた統合動画（v3 継承）。

### 11.3 overlay 図への原料濃度（feed Cin）水平点線（v3 継承）

すべての overlay 図に、R003 入口濃度を灰色の破線（`color="gray", ls="--"`）として水平に描画する。`Cin > 0` のときだけ描画。凡例は `"<species> feed (Cin=<value>)"`。完全破過の到達レベルを視覚的に確認するのが狙い。

**実装場所**:
- `evaluate_model_accuracy_v4.plot_overlay_sparse(..., feed_concentrations=...)`
- `run_chain_predict_v4._plot_chain_overlay(..., feed_concentrations=...)`

### 11.4 chain_adsorption_profile.gif（3 塔縦並び動画、v3 継承）

| 仕様 | 値 |
|:---|:---|
| サブプロット | `n_col` 行 × 1 列（R003 が上、R004A 中、R004B が下） |
| 各サブプロットの軸 | x: カラム位置 [cm]（塔ごとに L が異なるため独立スケール） |
| 同一フレーム内の時刻 | **全塔で同じ `snap_time`**（時間軸統一） |
| `frame_interval` | 既定 10 min |
| `gif_duration_ms` | 1 フレームあたり 800 ms |
| 描画内容 | 左軸: 吸着量 q（積み上げ面）、右軸: 液相濃度 C（線）。**v4: 吸着種のみ**（FAEE/ET は q を持たないため非表示） |

**実装場所**: `run_chain_predict_v4._plot_chain_adsorption_profile_gif`

### 11.5 フィット結果（v4 変更）

```
results_v4/
  ├── R003/<date>/
  │   ├── fitted_vector.json       q_total, k_hyd の値（keys/values + 個別キー）★v4
  │   ├── de_progress.json         DE 収束履歴
  │   └── readme.txt               バッチ番号 / 計算時間 / 境界条件
  ├── R004A/<date>/
  │   └── （同上。based_on_R003 に R003 の k_hyd も記録）★v4
  └── R004B/
      └── latest/used_parameters.json   R004A から流用したパラメータ（k_hyd 含む）★v4
```

`results_v4/R003/latest/`, `results_v4/R004A/latest/` に最新版コピーを保持し、可視化スクリプトから安定参照する（v3 継承）。

---

## §12 想定する妥当性チェック 【v4変更】

### 12.1 樹脂物性の整合性（v3 継承）

R003 と R004A は同じ樹脂のはずなので、`q_total` の値は近いことが期待される。大きく違えば:

- **q_total 差大**: 塔ごとに劣化進行度が違う（先住塔 R003 の方が劣化大の可能性）
- モデルが捉えきれていない物理（粒径分布、流路偏流など）

### 12.2 R004B 予測の物理的妥当性（v3 継承）

R004B 出口の VE1, VE2 が「ほぼゼロ」になっているか確認。濃度が高く出る場合は R004A の樹脂劣化・流量過多・モデル構造の限界を検討。

### 12.3 k_hyd の塔間整合性（v4 新規）

k_hyd は R003 / R004A で独立にフィットするが、物理的には同一液・同一温度なら同じ値のはず。フィット後に両塔の値が桁で一致しているかを確認する。大きく違う場合:

- 反応がバルク均一でない（樹脂近傍で促進/抑制されている）可能性
- k_hyd が他の誤差（入口濃度誤差、R003 モデル誤差の伝播）を吸収している可能性

### 12.4 v4 固有の検証項目（v4 新規）

- **VE ピーク**: R003 の VE1/VE2、R004A の VE1 のピーク高さ・位置・シャープさが v3 より改善するか（§1.1 の課題）
- **FA 早期リーク**: R003 の FA 出口が観測の早期リーク（180〜360 min の約 0.004）を再現するようになるか（仮説の独立した裏付け、§1.3）
- **切り分け**: 「水のみ（k_hyd=0）」ケース（§5.2 B）で原料水単独の影響量を確認し、加水分解の寄与と分離する
- **バッチ間差**: FAEE フィード濃度の違い（0.999 / 0.665 / 0.774）でバッチ間の挙動差が説明できるか

### 12.5 物質収支（v3 継承 + v4 注記）

各塔の (入量) − (出量) − (吸着量) ≈ 0 の検算。v4 では FA の収支に加水分解生成分が加わることに注意（FAEE 消費量 = FA 生成量、付録 F.3 の不変量で確認可能）。

---

## §13 今後の拡張余地（実装スコープ外） 【v4変更】

| 項目 | 概要 |
|:-----|:-----|
| **Phase 2: 加水分解の可逆化** | 逆反応項 + K_hyd をフィット対象に追加（§2.4）。ET は実装済みのため小改修で可能 |
| **原料濃度の時系列入力** | 現在はバッチごと固定値。C_in_series 機構は既にあるため、原料組成の時間変化にも拡張可能 |
| **原料分析値への更新** | 水・ET の分析値が揃い次第 `feed_composition_v4.csv` を更新（コード変更不要） |
| 塔の本数増減 | column_chain_v4 は可変長リスト対応（v3 継承） |
| パラメータ共通化フィット | 「樹脂物性は共通」「k_hyd は全塔共通」という制約を入れた同時フィット（オプション） |
| バックフラッシュ・再生工程 | 現在は 1 工程内のシミュレーションのみ。再生サイクルは別途検討 |

---

## §14 観測データの切替（実データ ↔ 人工テストデータ） 【v3継承】

### 14.1 ファイル構成

| ファイル | 役割 |
|:--------|:------|
| `data/r004a_observed_batches.csv` | 実データ |
| `data/r004a_observed_batches_test.csv` | 人工テストデータ（コード動作確認用） |
| `data/r004b_observed_batches.csv` | 実データ |
| `data/r004b_observed_batches_test.csv` | 人工テストデータ |

### 14.2 切替方法（コマンドライン引数 `--data-suffix`）

```bash
# 人工テストデータで動作確認
python run_fit_R004A_v4.py --data-suffix _test
python run_chain_predict_v4.py --data-suffix _test

# 実データで本番実行（デフォルト）
python run_fit_R004A_v4.py
python run_chain_predict_v4.py
```

### 14.3 人工テストデータの設計方針

| 塔 | 設計根拠 |
|:---|:---------|
| R004A | R003 出口実観測の 30〜50% 程度に低減。各バッチ 3〜4 点 |
| R004B | VE 漏れなし方針を反映した微小値。バッチ後半 1 点 |

人工データの目的は**コード動作確認**であり、フィット結果の物理的正しさは保証しない。

---

## §15 実装ステップ（チェックリスト） 【v4新規】

- [x] v4 フォルダ作成 + v3 全ファイルコピー・`_v3`→`_v4` リネーム
- [x] `competitive_adsorption_v4.py`: FAEE / ET 種追加（nonadsorbing_species）
- [x] `competitive_adsorption_v4.py`: 加水分解ソース項（S_hyd、モル保存形）
- [x] `competitive_adsorption_v4.py`: k_hyd フィットキー対応（_get/_set_by_key, fit_targets）
- [x] `feed_composition_v4.py` + `data/feed_composition_v4.csv`（固定値/CSV 切替）
- [x] `parameter_fitting_v4.py`: feed_provider 機構 + DEFAULT_BOUNDS["k_hyd"]
- [x] `evaluate_model_accuracy_v4.py`: feed_provider 機構
- [x] `params_columns_v4.py`: `_set_fit_targets_q_khyd`（q_total + k_hyd）
- [x] `run_fit_R003_v4.py`: feed 注入・k_hyd 初期値・snapshot 拡張
- [x] `run_fit_R004A_v4.py`: R003 k_hyd 反映・feed 注入・全成分リレー
- [x] `run_chain_predict_v4.py`: k_hyd 反映・feed 注入・`--plot_new_species`
- [x] スモークテスト（`smoke_test_v4.py`、付録 F）
- [ ] **次フェーズ**: R003 実データフィット（q_total + k_hyd、ESX004/006）
- [ ] **次フェーズ**: R004A 実データフィット → ESX012 チェーン予測 → v3 結果と比較
- [ ] **次フェーズ**: 切り分けケース（水のみ）の実行と影響量の記録
- [ ] **次フェーズ**: VE ピーク・FA 早期リークの改善評価（§12.4）

---

## 付録 A: v3 時点の動作検証記録 【v3継承】

以下は v3 実装フェーズで実施したスモークテストの記録である（v3 時点のもの。v4 の検証は付録 F）。

### A.1 v2 後方互換性

| テスト項目 | 結果 |
|:----------|:-----|
| `Params(C_in_series=None)` で `C_in_at(t)` が固定値を返す | ✅ OK |
| 同じ入口条件を v2 と v3 で `simulate()` した結果の出口濃度差 | ✅ **0.00e+00**（完全一致） |
| `inlet_provider=None`, `flow_provider=None` での目的関数値 | ✅ 既存 v2 と同値 |

### A.2 C_in_series 拡張機能

| テスト項目 | 結果 |
|:----------|:-----|
| `Dict[species, (t_array, v_array)]` 形式の指定 | ✅ 線形補間動作 |
| `pd.DataFrame(columns=[t_min, VE1, VE2, FA])` 形式の指定 | ✅ 動作 OK |
| 一部成分のみ時系列、残りは固定値 `C_in` の併用 | ✅ 期待通り |
| 補間範囲外は端値で外挿 | ✅ 動作 OK |

### A.3 3 塔リレーの整合性

| テスト項目 | 結果 |
|:----------|:-----|
| 各塔の `simulate()` 完了 | ✅ OK |
| **R003 出口 vs R004A 入口の濃度差** | ✅ **0.00e+00**（リレー時に濃度破壊なし） |
| 物理的妥当性（下流塔ほど VE1 max が小さい） | ✅ R003 (4.6e-6) → R004A (2.0e-9) → R004B (8.4e-13) |

### A.4 スクリプト起動

v3 の 4 スクリプトすべて `--help` 起動 ✅

### A.5 人工テストデータでの通し動作（Step1→2→3 全通し）

| テスト項目 | 結果 |
|:----------|:-----|
| `run_fit_R003_v3.py --train_batches 6 --t_end 460` 完走 | ✅ 約 3 分 |
| `run_fit_R004A_v3.py --train_batches 6 --t_end 460 --data-suffix _test` 完走 | ✅ 約 3 分 |
| `run_chain_predict_v3.py --batches 6 --t_end 460 --data-suffix _test` 完走 | ✅ 約 25 秒 |
| 全 9 サブプロットで物理的妥当な破過カーブ | ✅ |
| `chain_summary.csv` に各塔の RMSE 記録 | ✅ |

### A.6 実データ batch3 単独フィット（v3 本番動作検証）

| テスト項目 | 結果 |
|:----------|:-----|
| `run_fit_R003_v3.py --train_batches 3 --t_end 460` 完走 | ✅ 4m 28s, DE 7 反復 |
| `run_fit_R004A_v3.py --train_batches 3 --t_end 460` 完走 | ✅ 8m 6s, DE 15 反復 |
| `run_chain_predict_v3.py --batches 3 --t_end 460` 完走 | ✅ 約 46 秒 |
| R003 RMSE | ✅ VE1=0.010, VE2=0.008, FA=0.032 |
| R004A RMSE | ✅ VE1=0.025, VE2=0.008, FA=0.003 |
| R004B RMSE（観測 1 点） | ✅ VE1=0.0003, VE2=1e-6, FA=0.002 |
| 物理的妥当性: R003→R004A→R004B でピーク時刻が約 100 min ずつ遅延 | ✅ |
| 物理的妥当性: VE1 オーバーシュートが Cin に収束（完全破過） | ✅ |

### A.7 補助図表機能

| テスト項目 | 結果 |
|:----------|:-----|
| `overlay_*.png` に feed Cin の水平点線が表示 | ✅ |
| `chain_overlay.png` の全サブプロットで feed 点線併記 | ✅ |
| `Cin <= 0` の場合は点線スキップ | ✅ |
| `chain_adsorption_profile.gif` が batch ごとに 1 個生成、全塔同一 snap_time | ✅ |

### A.8 v3 時点で未検証だった項目

| 項目 | 内容 | v4 での状況 |
|:-----|:-----|:-----|
| 多バッチ学習（実データ） | `--train_batches 4,6` 等での共通パラメータ学習 | v3 後期に ESX004/006 学習で実施済み（本書 §1.1 の予測はこの結果） |
| batch 単独間のばらつき評価 | 各 batch 単独フィット結果の比較 | 未実施 |
| fit_targets 拡張 | q_total/eps_b 以外を可変にした場合の収束性 | v4 で k_hyd を追加（本書の主題） |

---

## 付録 B: 使い方ガイド 【v4変更】

### B.1 全体フロー

```
[初回セットアップ]
  ├─ data/process_runs_batches.csv    : R003 入口・出口（v3 から流用済み）
  ├─ data/feed_composition_v4.csv     : 原料組成（FAEE/ET/水）を確認・更新 ★v4
  ├─ data/r004a_observed_batches.csv  : R004A 観測値
  └─ data/r004b_observed_batches.csv  : R004B 観測値

[Step 1] R003 単独フィット（q_total + k_hyd）
  python run_fit_R003_v4.py --train_batches 4,6
  → results_v4/R003/<date>/fitted_vector.json

[Step 2] R004A フィット（R003 の出口を入口として使う。q_total + k_hyd を独立フィット）
  python run_fit_R004A_v4.py --train_batches 4,6 --species VE2
  → results_v4/R004A/<date>/fitted_vector.json

[Step 3] 3 塔チェーン予測（R004B は R004A 流用）
  python run_chain_predict_v4.py --batches 12 --plot_new_species
  → outputs_chain_v4/chain/batch_12/{R003,R004A,R004B}/...
  → outputs_chain_v4/chain/chain_summary.csv

[再可視化] 保存済みパラメータで描画のみ（フィットしない）
  python run_chain_visualize_v4.py --batches 12
```

### B.2 各スクリプトの引数

#### B.2.1 `run_fit_R003_v4.py`

| 引数 | 説明 | デフォルト |
|:-----|:-----|:---------|
| `--t_end` | シミュレーション終了時間 [min] | 未指定時はバッチごと自動（付録 C） |
| `--train_batches` | 学習用バッチ番号（カンマ区切り） | `6` |
| `--test_batches` | 検証用バッチ番号 | train と同じ |
| `--species` | フィット対象の成分（付録 D） | 未指定=観測がある全成分 |

```bash
# 多バッチ学習（ESX004 + ESX006）
python run_fit_R003_v4.py --train_batches 4,6

# VE1, VE2 のみを目的関数に含める（FA を除外）
python run_fit_R003_v4.py --train_batches 6 --species VE1,VE2

# t_end を延長して通水停止後の挙動も含める
python run_fit_R003_v4.py --train_batches 6 --t_end 800
```

#### B.2.2 `run_fit_R004A_v4.py`

| 引数 | 説明 | デフォルト |
|:-----|:-----|:---------|
| `--t_end` | シミュレーション終了時間 [min] | 未指定時はバッチごと自動 |
| `--train_batches` | 学習用バッチ番号 | `6` |
| `--test_batches` | 検証用バッチ番号 | train と同じ |
| `--data-suffix` | 観測 CSV のサフィックス | `""`（実データ） |
| `--species` | フィット対象の成分（付録 D） | 未指定=観測がある全成分 |

```bash
# R004A を VE2 のみでフィット（現行運用）
python run_fit_R004A_v4.py --train_batches 4,6 --species VE2

# 人工テストデータで動作確認
python run_fit_R004A_v4.py --train_batches 6 --data-suffix _test
```

実行に先立ち `results_v4/R003/latest/fitted_vector.json` が必要（先に Step 1 を実施）。

#### B.2.3 `run_chain_predict_v4.py`

| 引数 | 説明 | デフォルト |
|:-----|:-----|:---------|
| `--t_end` | シミュレーション終了時間 [min] | 未指定時はバッチごと自動 |
| `--batches` | 予測対象バッチ番号 | `2,3,4,6` |
| `--columns` | 予測する塔 | `R003,R004A,R004B` |
| `--data-suffix` | 観測 CSV のサフィックス | `""` |
| `--label` | 出力フォルダのラベル | `chain` |
| `--plot_new_species` | ★v4: chain_overlay に FAEE / ET / OH パネルを追加 | off |

```bash
# ESX012 予測（通液時間 2 時間延長の例）
python run_chain_predict_v4.py --batches 12 --t_end 546 --label chain_pre012

# 新種（FAEE / ET / 水）の挙動も俯瞰図で確認
python run_chain_predict_v4.py --batches 12 --plot_new_species

# 別ラベルで結果を残す（過去結果を保護）
python run_chain_predict_v4.py --label chain_2026Q3_run1
```

### B.3 多バッチ学習について（v3 継承）

`--train_batches` をカンマ区切りで複数指定するだけ。`parameter_fitting_v4._objective_batches` が各バッチを個別に simulate し、バッチ×成分ごとの RMSE を単純平均して目的関数とする。つまり「全ロットで共通の (q_total, k_hyd)」を 1 回の DE 最適化で求める。

| 注意点 | 説明 |
|:-----|:-----|
| 計算時間 | バッチ数に比例 |
| 観測点の偏り | 各バッチの RMSE が等重みで平均（点数による重み付けなし） |
| 物理整合性 | 全バッチで同じパラメータを仮定。ロット間差が大きいと精度低下 |
| **v4: FAEE 濃度差** | バッチごとの FAEE フィード濃度差は feed_composition_v4 が自動反映するため、k_hyd は全バッチ共通で良い |

推奨手順: 1 バッチで確認 → 複数バッチに拡大 → 全バッチ本番（v3 継承）。

### B.4 出力ディレクトリ構造

```
results_v4/
├── R003/
│   ├── 2026-XX-XX/                       実行日付フォルダ（同日複数回は _2, _3）
│   │   ├── fitted_vector.json            q_total, k_hyd の値
│   │   ├── de_progress.json              DE 収束履歴
│   │   └── readme.txt                    実行条件（境界・バッチ・経過時間）
│   └── latest/fitted_vector.json         最新版コピー（自動更新）
├── R004A/（同上）
└── R004B/latest/used_parameters.json     R004A から流用（k_hyd 含む）

outputs_chain_v4/
├── R003/, R004A/                          フィット時の train/test 評価出力
├── chain/                                 run_chain_predict_v4.py の出力
│   ├── batch_X/
│   │   ├── R003/ R004A/ R004B/            overlay_*.png / gif / model_output.csv
│   │   ├── chain_overlay.png              3 塔俯瞰図（--plot_new_species で 6 成分）
│   │   └── chain_adsorption_profile.gif
│   └── chain_summary.csv                  全バッチ・全塔・全成分の RMSE/MAPE
└── <label>/                               --label で個別保存
```

### B.5 トラブルシューティング（v3 継承 + v4 追加）

| 症状 | 対処 |
|:-----|:-----|
| `FileNotFoundError: results_v4/R003/latest/fitted_vector.json` | Step 1 (`run_fit_R003_v4.py`) を先に実行する |
| R004A 観測点が 0 | `data/r004a_observed_batches.csv` を確認、または `--data-suffix _test` |
| フィット結果が境界に張り付く | `parameter_fitting_v4.DEFAULT_BOUNDS` の境界を緩める |
| 数値振動が出る | `params_columns_v4.apply_dt_react_cap_all()` で `dt_react_cap` を 0.05 や 0.03 に縮める |
| **k_hyd が上限 0.05 に張り付く** ★v4 | `DEFAULT_BOUNDS["k_hyd"]` の上限を緩める（必要なら dt_react_cap も縮小） |
| **`[feed] batch X: C0_FAEE が…無いため固定値 0.0 を使用します` 警告** ★v4 | そのバッチの原料組成が未登録。`data/feed_composition_v4.csv` に行を追加する（加水分解を効かせたい場合） |
| R004B 予測が VE 漏れあり | R004A のフィット結果が想定外の可能性。overlay 図と物理整合性を確認 |

---

## 付録 C: シミュレーション終了時間 `t_end` の自動決定（バッチ個別） 【v3継承】

### C.1 背景

バッチごとに運転時間が異なるため、固定 t_end では長いバッチの打ち切り・短いバッチの過剰計算が起きる。

### C.2 方針

- **`--t_end` 未指定時**: 各バッチの `t_min` 範囲に合わせて **バッチごとに個別に** `t_end` を自動決定する。
- **`--t_end` 指定時**: その値を **全バッチ一律** で使用（自動決定は無効）。

### C.3 自動決定式

```
base   = max( プロセスデータ t_min 最大, 観測データ t_min 最大 )
t_end  = base × T_END_MARGIN_FACTOR      (既定 1.0 = マージン無し)
```

- `T_END_MARGIN_FACTOR = 1.0`（`params_columns_v4.py`）＝マージン無し。工程終了（通液停止）の時刻で描画を切る。
- 工程終了後の内部緩和を可視化したい場合はこの定数を 1.05 などに増やす。増やした区間は C.3.1 の通り流量 0 として扱われる。

### C.3.1 工程終了後の流量の扱い（通液停止）

`vT_at` は流量時系列の範囲外を両端値で外挿する。v4 では「運転時間を延長したらどうなるか」の仮想計算を主用途とするため、**延長区間の既定は末尾流量のまま通液継続**（`prepare_flow_series(..., beyond="last")`）とする。平均値固定（`mean`）や旧仕様の通液停止（`stop`）は `--flow_beyond` で切替可能。実装は `params_columns_v4.prepare_flow_series`。fit・predict・visualize すべてで適用する。

- 観測点は `t_min.max` 以前にしか無いため、フィットの目的関数（RMSE）には影響しない。
- 有効な t_min が 1 つも無い場合は `None` を返し、`Params` デフォルト（350 min）へフォールバックする。

**v4 注記**: 通液停止後（流量 0 区間）も加水分解反応は進行する（反応のみ進行、の「反応」に含まれる）。`--t_end` 延長予測（例: ESX012 の 2 時間延長）ではこの区間の FA 生成・VE 置換が計算される。

### C.4 実装

| 場所 | 内容 |
|:-----|:-----|
| `params_columns_v4.compute_auto_t_end(*t_min_arrays, margin_factor=...)` | 複数 t_min 配列の最大 × マージンを返す共通ヘルパー |
| `params_columns_v4.prepare_flow_series(flow_df, beyond=...)` | 延長区間の流量: last（既定）/ mean / stop |
| `parameter_fitting_v4.fit_parameters_batches(..., t_end_provider=None)` | バッチごとに `t_end_provider(bid, dfb)` を評価し `p.t_end` に反映 |
| `evaluate_model_accuracy_v4.evaluate_model_batches(..., t_end_provider=None)` | 評価でも同様 |
| `run_fit_R003_v4.py` / `run_fit_R004A_v4.py` / `run_chain_predict_v4.py` | 未指定時にバッチ別 t_end を自動算出 |

### C.5 記録

- `fitted_vector.json` に `t_end_mode`（`"fixed"` / `"auto_per_batch"`）と `t_end_desc`（バッチ別 t_end の一覧）を記録。
- `readme.txt` の `t_end` 行にも同じ説明を出力する。

---

## 付録 D: フィット対象成分の選択（`--species`） 【v4変更（軽微）】

### D.1 背景（v3 での追加経緯）

観測濃度のスケールが成分間で大きく異なり（FA は VE1 の約 1/20 の希薄域）、目的関数（成分等重みの絶対 RMSE 平均）への FA の寄与が小さい。「希薄で当てはめの意味が薄い成分を目的関数から外したい」という運用要望から、フィット対象成分を実行時に選択できる機能を追加した。

### D.2 仕様

- 対象スクリプト: `run_fit_R003_v4.py` / `run_fit_R004A_v4.py`
- `--species VE1,VE2` のようにカンマ区切りで指定。指定成分だけを**目的関数（RMSE）**に用いる。
- 未指定時は観測列 `C_*_exp` がある全成分を自動採用（後方互換）。
- **評価図・metrics は指定に関わらず全観測成分を表示**する。
- 観測列が無い成分名は警告して無視。有効成分 0 個ならエラーで停止。

> **注意（重み付けではない）**: 「含める／含めない」の二択であり、成分ごとの連続的な重み付けはできない。

### D.3 使用方法

**指定できる成分は `VE1` / `VE2` / `FA` の 3 種類**（観測列 `C_*_exp` を持つ成分）。

**v4 注記**: FAEE / ET / OH は観測列を持たないため指定できない（指定しても警告して無視される）。これらの成分はフィットの目的関数に入らないが、モデル内では常に計算される。

| 成分名 | 内容 |
|:-------|:-----|
| `VE1` | ビタミンE（Toc-α + T3-α） |
| `VE2` | ビタミンE（Toc-β + T3-β） |
| `FA`  | 遊離脂肪酸（競合成分。v4 では加水分解による生成分も含む） |

```bash
python run_fit_R003_v4.py  --train_batches 6 --species VE1,VE2
python run_fit_R004A_v4.py --train_batches 6 --species VE2
```

### D.4 実装

| 場所 | 内容 |
|:-----|:-----|
| `parameter_fitting_v4.fit_parameters_batches(..., target_species=None)` | 目的関数に含める成分リスト |
| `run_fit_R003_v4.py` / `run_fit_R004A_v4.py` | `--species` を解析して `target_species` に渡す |

### D.5 記録

`fitted_vector.json` に `fit_species` / `fit_species_desc`、`readme.txt` に `Fit species` 行を出力。

---

## 付録 E: 工程終了時刻ちょうどの観測が評価から欠落する問題と対応（端の許容） 【v3継承】

### E.1 症状

R004B は観測が「工程終了時刻ちょうどの 1 点」しか無いため、その 1 点が評価から外れると R004B 全体が `chain_summary.csv` から欠落する。

### E.2 原因

1. `simulate` の時間積分で浮動小数点誤差が累積し、工程終了 t_end ちょうどのサンプルが記録されないことがある（モデル最終時刻が t_end − dt になる）。
2. 評価の補間が `np.interp(..., right=np.nan)` のため、モデル範囲外（= t_end ちょうど）の観測が NaN 扱いで除外されていた。

### E.3 対応（評価側で端の許容）

評価時の補間を共通ヘルパー **`evaluate_model_accuracy_v4.interp_obs_to_model()`** で行う。観測時刻がモデル最終時刻を**わずか（許容幅 `tol` 以内）**超える場合は端へスナップしてモデル端値と比較する。`tol` 未指定時は「モデル時間刻みの中央値 × 1.5」。許容幅を超える真の範囲外観測は従来どおり NaN（除外）。

### E.4 実装

| 場所 | 内容 |
|:-----|:-----|
| `evaluate_model_accuracy_v4.interp_obs_to_model(t_obs, t_mod, y_mod, tol=None, tol_factor=1.5)` | 端の許容つき補間ヘルパー |
| `evaluate_model_accuracy_v4._evaluate_one_set` / `plot_overlay_sparse` | metrics・overlay の補間に使用 |
| `run_chain_predict_v4._eval_overlay_one_column` | `chain_summary.csv` の補間に使用 |

### E.5 注意（フィット目的関数は対象外）

本対応は**評価側のみ**。フィットの目的関数の補間は従来どおり端の許容なしで、工程終了 t_end ちょうどの観測はフィットには使われない（挙動を変えないため据え置き）。

---

## 付録 F: v4 動作検証記録 【v4新規】

実装フェーズ（2026-07-29）で `smoke_test_v4.py` により実施。実行方法: `python smoke_test_v4.py`

### F.1 v3 完全互換

k_hyd=0・原料 FAEE/ET/水=0 の条件で、v3 と v4 の `simulate()` 出口濃度（VE1/VE2/FA/OH）を比較。

| テスト項目 | 結果 |
|:----------|:-----|
| 出口濃度の最大絶対差 | ✅ **0.00e+00**（完全一致） |

（v4 は成分数が増えるが、FAEE/ET がゼロのとき既存成分の時間刻み・更新式は v3 と同一のため、ビット単位で一致する）

### F.2 加水分解 ON の方向確認

原料 FAEE=ET=0.774, 水=0.243、k_hyd=1.0 vs 0.0 の比較（R003 相当・t_end=200 min）。

| テスト項目 | 結果 |
|:----------|:-----|
| FAEE 出口積算の減少（消費） | ✅ 1004.8 → 481.8 |
| FA 出口積算の増加（生成） | ✅ 0.007 → 286.3 |
| ET 出口積算の増加（生成） | ✅ 1004.8 → 1527.8 |

### F.3 モル保存の不変量

FAEE と ET は同じ θ を持ち消費量=生成量のため、**(C_FAEE + C_ET) は反応の有無に依らず同一**になるはず。

| テスト項目 | 結果 |
|:----------|:-----|
| (C_FAEE + C_ET) 出口時系列の k_hyd=0/1 間の最大差 | ✅ 2.2e-15（丸め誤差レベル） |

### F.4 数値健全性

| テスト項目 | 結果 |
|:----------|:-----|
| 全出口濃度に NaN / 負値なし | ✅ |

### F.5 feed_composition_v4

| テスト項目 | 結果 |
|:----------|:-----|
| batch 6 → FAEE=ET=0.665, OH=0.243（CSV 値） | ✅ |
| batch 12 → FAEE=ET=0.774（CSV 値） | ✅ |
| batch 2（CSV に無い）→ FAEE=ET=0（デフォルト）+ 警告表示、OH=0.243 | ✅ |

### F.6 実行スクリプト起動・チェーン通し

| テスト項目 | 結果 |
|:----------|:-----|
| `run_fit_R003_v4.py` / `run_fit_R004A_v4.py` / `run_chain_predict_v4.py` / `run_chain_visualize_v4.py` の `--help` 起動 | ✅ 全て exit 0 |
| `run_chain_predict_v4.py --batches 6 --t_end 120 --plot_new_species`（v3 フィット結果を一時流用）完走 | ✅ 約 15 秒 |
| chain_overlay.png が 3 塔 × 6 成分（VE1/VE2/FA/FAEE/ET/OH）で描画 | ✅ |
| FAEE / ET が非吸着のため速く破過（~50 min）、OH（H=5.22）は遅れて破過し、下流塔へ正しくリレー | ✅ |

※ この通し確認に使った一時的な results_v4（v3 パラメータのコピー）と outputs_chain_v4 は削除済み。**v4 の本番フィット（q_total + k_hyd）は未実施**であり、§15 の次フェーズとして実行する（原料水の追加により v3 の q_total は流用できない。§2.5）。

### F.7 まだ検証されていない項目

| 項目 | 内容 |
|:-----|:-----|
| **実データフィット** | R003 / R004A の q_total + k_hyd フィット（ESX004/006）と収束性 |
| **ESX012 予測の改善評価** | VE ピーク・FA 早期リークが v3 より改善するか（§12.4） |
| **切り分けケース** | 水のみ（k_hyd=0）の影響量の定量化（§5.2） |
| **k_hyd の塔間整合性** | R003 と R004A のフィット値比較（§12.3） |

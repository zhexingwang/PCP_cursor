# 競争吸着パイプライン v3 仕様書（三塔モデル化）

## §0 v3 の位置づけ

v2 までの仕様（§1〜§19）は v2 フォルダ内の `competitive_adsorption_pipeline_spec.md` に記述されており、**引き続き有効**である。本仕様書はその上で、**v3 における追加・変更点**を記述する。

| 項目 | v2 | v3 |
|:-----|:---|:---|
| モデル対象 | R003 単塔 | **R003 → R004A → R004B の 3 塔直列** |
| 入口濃度 | 固定値 `C_in: Dict[str, float]` | **時系列 `C_in_series` 対応**（前塔出口を後塔入口にリレー） |
| フィット対象 | 単塔の任意パラメータ | 各塔ごと独立フィット（R003, R004A） |
| 検証データ | R003 出口（密） | R003 出口（密）+ R004A 出口（疎 3〜4 点）+ R004B 出口（1 点） |
| フォルダ分離 | — | **v2 を完全に保持し、v3 は別フォルダで開発** |

v2 はリファレンス実装として保持し、v3 で 3 塔チェーン機能を完成させる。

---

## §1 対象プロセスとモデル思想

### 1.1 装置構成

```
原料 → [R003: d=25, L=83.2] → [R004A: d=25, L=42.0] → [R004B: d=25, L=42.0] → 製品
       ↑                       ↑                       ↑
       既存モデル化            新規モデル化            新規モデル化（予測のみ）
```

- 全塔とも **同一樹脂**（d=25 cm 共通）
- R003 と R004A/R004B では **樹脂量（カラム長）が異なる**
- 直列配管・**v_T は全塔で共通**
- 各塔とも吸着工程前に洗浄工程があり、**初期吸着量 q0 = 0**

### 1.2 塔ごとの劣化想定

樹脂の種類は同じだが、塔ごとに通液履歴・劣化具合が異なる可能性があるため、**`q_total` と `eps_b` は塔ごとに独立にフィット**する。他のパラメータ（`H`, `k_ads`, `Keq_ads`, `k_ex`, `Keq_ex`）は樹脂物性として共通の固定値を採用する。

### 1.3 R004B はフィットせず予測のみ

R004B の出口は VE 漏れがほぼ無いことが運用上の前提であり、観測値も「ロットごと 1 点・ほぼゼロ」となるためフィットには不向き。R004B は **R004A と同じパラメータを使用して予測のみ**を行い、結果が「ほぼゼロ」になることを妥当性確認の根拠とする。

### 1.4 キャリア成分 OH の定義と役割

本モデルは陰イオン交換樹脂（Diaion PA306S、活性サイト = 水酸化物イオン形）を対象とした論文
（*Competitive Adsorption Model for Process Design ...*, J. Chem. Eng. Japan, Vol.53 No.9, 2020）
の定式化に基づく。v2/v3 のコードはこの論文を参考に実装している。成分 `OH` は他成分（VE1/VE2/FA）とは
性質が異なる **キャリア（carrier）成分** であり、以下のように定義・実装されている。

**イオン交換反応（論文 Eqs.(2)–(4)）**
- (2) VEH + S⁺(OH⁻) ⇌ S⁺(VE⁻) + H₂O
- (3) FaH + S⁺(OH⁻) ⇌ S⁺(Fa⁻) + H₂O
- (4) FaH + S⁺(VE⁻) ⇌ S⁺(Fa⁻) + VEH

サイト S⁺ は常に OH⁻ / VE⁻ / Fa⁻ のいずれか 1 つを対イオンとして保持する。
交換の量論は **1:1**（VE/Fa が 1 吸着すると OH⁻ が 1 放出される）。

**固相側（サイト占有）**
- `S⁺(OH⁻)` = まだ標的成分が付いていない **空き（基準）サイト**。VE/Fa が交換で入っていく先。
- 論文 Eq.(10) `q_total = q_OH + q_VE + q_Fa` に対応。
- コードでは `q_OH` を独立変数として持たず、移項した形で **`q_free = q_total − q_VE − q_Fa`（= q_OH）** として従属的に計算する
  （`reaction_rhs_*` 内の `q_free`）。
- したがって OH はサイトを占有するが、**専用の固相負荷 `q[OH]` や吸着定数 `k_ads[OH]`/`Keq_ads[OH]` は持たず、フィット対象にもならない**。
  OH は「基準となる空きサイトそのもの」であり、その量は他成分の補数として決まる。

**液相側（C_OH）**
- コードの液相 `OH` は **樹脂液相の OH⁻ 濃度**（論文 C_OH）に対応する。
- 反応(2)(3)で VE/Fa が吸着すると **1:1 で OH⁻ が液相へ放出**される
  （`simulate` 内のキャリアソース項 `src_carrier = -alpha · Σ dq`）。これは論文の「生成物（H₂O）＝逆反応を駆動する項」に対応する。
- 放出で増えた液相 OH⁻（`Cres[OH] = H[OH]·C[OH]`、既定 `H[OH]=5.22`）が、
  **逆反応（脱着）項 `kr · q[i] · Cres[OH]` を駆動**する（液相 OH⁻ が高いほど脱着が進み、吸着平衡を抑制）。

**初期値の扱い**
- 液相初期濃度は全成分ゼロ（案A）で開始するため、OH も塔内初期濃度 0。
  OH は「入口からの流入（チェーンでは前段塔出口 OH）」＋「吸着に伴う塔内生成」で供給されるため、0 始まりで物理的な不足は生じない。
- OH は固相吸着量を独立に持たない（上記のとおり `q_free` で表現）ため、`q0[OH]` という概念は無い。

**一言まとめ**
> OH は「固相では空きサイトの基準状態（量は `q_free` で従属的に決まる）」であり、
> 「液相では交換で放出される OH⁻ 濃度（＝逆反応=脱着を駆動する生成物項）」として実装されている。

---

## §2 ワークフロー

```
Step 1: R003 フィット
  入口: 実験 C0（process_runs_batches.csv）
  出口: R003 観測（密）
  対象: q_total, eps_b（fit_targets=True）
  出力: results_v3/R003/<date>/fitted_vector.json

Step 2: R004A フィット（R003 のパラメータ固定）
  R003 のフィット結果をロード → R003 を simulate
  R003 出口時系列を R004A の C_in_series として渡す
  対象: R004A の q_total, eps_b
  出口観測: R004A 実測（疎 3〜4 点）
  出力: results_v3/R004A/<date>/fitted_vector.json

Step 3: 3 塔チェーン予測
  R003・R004A のフィット済みパラメータをロード
  R004B は R004A と同じパラメータを使用（予測専用）
  3 塔を順に simulate して各塔の出口を出力
  出力: outputs_chain_v3/<train|test>/batch_X/{R003,R004A,R004B}/

Step 4: 妥当性検証（手動）
  - R003・R004A のフィット結果を表で並べて q_total・eps_b を比較
  - 樹脂物性が同じはずなのに大きく違えば、データ・モデルの再検討
  - R004B 予測の VE 出口濃度が「ほぼゼロ」になっていることを確認
```

---

## §3 入口濃度のリレー方式（B1 採用）

R004A の入口は、**R003 モデルの出口時系列**（実測値ではなくモデル予測）を用いる。理由:

- R003 のフィット結果を「予測モデル」として確定させる目的のため、その出力を素直に下流に渡すのが自然
- 実運用時には実測値が存在しない時間帯まで予測する必要があるため、モデル出口を用いる方式が一貫している

R003 のフィット精度が高いため、R004A 入口に渡される濃度は実測と概ね一致するという前提で設計する。万一 R003 モデルに大きな誤差がある場合は、R004A のフィット品質に影響が出るため、Step 4 の整合性チェックで検出する。

---

## §4 ファイル構成

### 4.1 v3 フォルダの構造

```
v3/
├── competitive_adsorption_v3_spec.md           本ドキュメント
│
├── competitive_adsorption_v3.py                v2 をコピーして拡張（C_in_series 対応）
├── parameter_fitting_v3.py                     v2 をコピーして拡張
├── evaluate_model_accuracy_v3.py               v2 をコピーして拡張
│
├── params_columns_v3.py                        新規：塔別 Params の生成
├── column_chain_v3.py                          新規：3 塔リレー実行コア
│
├── run_fit_R003_v3.py                          R003 単塔フィット
├── run_fit_R004A_v3.py                         R003 結果を用いた R004A フィット
├── run_chain_predict_v3.py                     3 塔チェーン予測
├── run_chain_visualize_v3.py                   保存済みパラメータで再可視化
│
├── data/
│   ├── process_runs_batches.csv                R003 入口・出口（v2 からコピー）
│   ├── r004a_observed_batches.csv              R004A 出口観測（実データ用・空テンプレ）
│   ├── r004a_observed_batches_test.csv         R004A 人工テストデータ（コード動作確認用）
│   ├── r004b_observed_batches.csv              R004B 出口観測（実データ用・空テンプレ）
│   └── r004b_observed_batches_test.csv         R004B 人工テストデータ（コード動作確認用）
│
├── results_v3/
│   ├── R003/<date>/                            R003 フィット結果
│   ├── R004A/<date>/                           R004A フィット結果
│   └── R004B/<date>/                           R004B 使用パラメータ記録（フィットなし）
│
└── outputs_chain_v3/
    ├── train/batch_X/{R003,R004A,R004B}/       各塔の overlay/gif/CSV
    │   ├── overlay_VE1.png  overlay_VE2.png  overlay_FA.png
    │   ├── adsorption_profile.gif
    │   ├── model_output.csv
    │   └── observed.csv
    ├── train/batch_X/chain_summary.csv         3 塔の RMSE まとめ
    ├── train/batch_X/chain_overlay.png         3 塔出口の一覧オーバーレイ
    └── test/batch_X/                           （同上）
```

### 4.2 既存スクリプトとの関係

| v2 ファイル | v3 での扱い |
|:-----------|:-----------|
| `competitive_adsorption_v6_revised_v2.py` | v3 にコピー → `competitive_adsorption_v3.py`（C_in_series 拡張） |
| `parameter_fitting_v2_batches.py` | v3 にコピー → `parameter_fitting_v3.py` |
| `evaluate_model_accuracy_v2_batches.py` | v3 にコピー → `evaluate_model_accuracy_v3.py` |
| `run_fit_from_process_data_v2.py` | v3 では `run_fit_R003_v3.py` として再構成 |
| `run_visualize_only_v2.py` | v3 では `run_chain_visualize_v3.py` として 3 塔版に拡張 |

v2 のファイルは**一切変更しない**。v3 はあくまで独立したコピーから派生する。

---

## §5 競争吸着モデル本体の拡張（C_in_series）

### 5.1 既存 v2 の `Params` 入口濃度

```python
# v2: 固定値のみ
C_in: Dict[str, float] = field(default_factory=lambda: {"VE1":..., "VE2":..., "FA":...})
```

### 5.2 v3 での拡張

`v_T_series` と同じ流儀で、時系列入口濃度を受けられるようにする。

```python
# v3: 固定値 OR 時系列を選択可能
C_in: Dict[str, float] = ...                    # 既存（固定値）
C_in_series: Optional[Dict[str, ...]] = None    # 新規（時系列）
```

`C_in_series` が指定されていれば各時刻 t で補間して使用、未指定なら従来通り `C_in` を固定値として使用。**v2 の挙動は完全に保たれる**。

### 5.3 内部 API

```python
def C_in_at(t: float) -> Dict[str, float]:
    """入口濃度を時刻 t で評価する。C_in_series があれば補間、なければ C_in 固定値を返す。"""
```

`simulate()` 内で入口境界条件を組む際に `C_in[sp]` を直接参照していた箇所を `C_in_at(t)[sp]` に置き換える。

---

## §6 塔別パラメータ生成（params_columns_v3.py）

```python
def get_R003_params() -> Params:
    p = Params()
    p.d = 25.0
    p.L = 83.2
    return p

def get_R004A_params() -> Params:
    p = Params()
    p.d = 25.0
    p.L = 42.0
    return p

def get_R004B_params() -> Params:
    p = Params()
    p.d = 25.0
    p.L = 42.0
    return p
```

樹脂物性（H, k_ads, Keq_ads, k_ex, Keq_ex）は全塔で同一の `Params` デフォルトを採用。寸法のみ塔ごとに上書きする。

`q_total`, `eps_b` はフィットで決まるため、生成時はデフォルト値を入れておき、フィット結果を後から `update_from_vector()` で上書きする。

> **フィット対象（fit_targets）の決定場所【重要】**
> 実際のフィット対象は、上記 `get_*_params()` 内で呼ばれる `params_columns_v3._set_fit_targets_q_eps()` が決定する。この関数が `p.fit_targets` を丸ごと作り直し、**`q_total` と `eps_b` の 2 つだけを True**（他はすべて False）に設定する。
>
> ```python
> def _set_fit_targets_q_eps(p: Params) -> None:
>     targets = {k: False for k in p.fit_targets}
>     targets["q_total"] = True
>     targets["eps_b"] = True
>     p.fit_targets = targets
> ```
>
> `competitive_adsorption_v3.py` の `Params.fit_targets` デフォルト値は「`Params()` を単体生成し、かつ `_set_fit_targets_q_eps` も上書きも通さない場合」のみ効くフォールバックであり、R003 / R004A / チェーンの通常フローでは**使われない（必ず上書きされる）**。フィット対象を変えたい場合は `Params` のデフォルトではなく `_set_fit_targets_q_eps`（または呼び出し側で `p.fit_targets` を上書き）を編集すること。

---

## §7 3 塔リレー実行コア（column_chain_v3.py）

### 7.1 アルゴリズム

```
入力: params_list = [P_R003, P_R004A, P_R004B]
      v_T_series  = 全塔共通の流量時系列
      C_in_first  = R003 の入口濃度時系列（実験 C0）

Step 1: P_R003 を simulate(v_T_series, C_in=C_in_first)
        → R003 出口濃度時系列 = effluent_R003

Step 2: P_R004A の C_in_series = effluent_R003
        P_R004A を simulate()
        → R004A 出口濃度時系列 = effluent_R004A

Step 3: P_R004B の C_in_series = effluent_R004A
        P_R004B を simulate()
        → R004B 出口濃度時系列 = effluent_R004B

戻り値: 各塔の (grid_df, effluent_df, meta, fields) のリスト
```

各塔のシミュレーションは独立に行うため、塔ごとに `t_end`, `dt`, `N` を変えても良い。デフォルトは全塔で揃える。

### 7.2 入力モードの選択

```python
simulate_chain(params_list, ..., n_cols=3)  # 3 塔通し
simulate_chain(params_list[:2], ...)         # R003+R004A まで
simulate_chain(params_list[:1], ...)         # R003 のみ
```

中間で打ち切れるため、フィットの各段階で必要分だけ実行可能。

---

## §8 出力構造

### 8.1 各塔の出力

```
outputs_chain_v3/<label>/batch_X/{R003,R004A,R004B}/
  ├── overlay_VE1.png         実測 vs モデル（R003 密、R004A 疎、R004B 1 点）
  ├── overlay_VE2.png         + 原料濃度の水平点線（feed Cin）も併記 (§8.4)
  ├── overlay_FA.png
  ├── adsorption_profile.gif  カラム内プロファイル動画
  ├── model_output.csv        モデル出口濃度時系列
  └── observed.csv            実測値（あれば）
```

R004B はフィットしないが、出力構造は他塔と揃える。

### 8.2 チェーン全体の集約

```
outputs_chain_v3/<label>/batch_X/
  ├── chain_summary.csv               各塔の RMSE まとめ（全成分）
  ├── chain_overlay.png               3 塔出口濃度の一覧オーバーレイ図
  └── chain_adsorption_profile.gif    3 塔のカラム内プロファイルを縦並びで表示する動画
```

- `chain_overlay.png` は 3 塔 × 3 成分（9 サブプロット）の俯瞰図。R003 → R004A → R004B の濃度低下が一目で分かる。各サブプロットに **原料濃度の水平点線** (§8.4) も併記。
- `chain_adsorption_profile.gif` は 3 塔の `adsorption_profile.gif` を縦並びにまとめた統合動画。**全塔とも同じ snap_time** で 1 フレーム化するため、上流→下流の吸着帯伝搬を直接比較できる (§8.4)。

### 8.3 overlay 図への原料濃度（feed Cin）水平点線

すべての overlay 図（各塔の `overlay_*.png` および `chain_overlay.png`）に、**R003 入口濃度** (`Params.C_in[species]`) を灰色の破線として水平に描画する。

| 項目 | 仕様 |
|:---|:---|
| 線色・スタイル | `color="gray", ls="--", lw=1.0–1.2` |
| 表示条件 | `Cin > 0` のときだけ描画（R004A 評価時など、入口がほぼ 0 の場合は自動でスキップ） |
| 凡例ラベル | `"<species> feed (Cin=<value>)"`（例: `VE1 feed (Cin=0.02635)`） |

**狙い**: 完全破過の到達レベル（モデル予測曲線がこの点線に漸近すれば吸着サイト飽和）を視覚的に確認できる。

**実装場所**:
- `evaluate_model_accuracy_v3.plot_overlay_sparse(..., feed_concentrations=...)`：単塔 overlay 図
- `run_chain_predict_v3._plot_chain_overlay(..., feed_concentrations=...)`：チェーン俯瞰図
- 引数 `feed_concentrations` は `Dict[str, float]`（成分→濃度）。`None` 時は v2 と同じ挙動（点線なし）。

### 8.4 chain_adsorption_profile.gif（3 塔縦並び動画）

各バッチ単位で 1 個生成される統合動画。

| 仕様 | 値 |
|:---|:---|
| サブプロット | `n_col` 行 × 1 列（R003 が上、R004A 中、R004B が下） |
| 各サブプロットの軸 | x: カラム位置 [cm]（塔ごとに L が異なるため独立スケール） |
| 同一フレーム内の時刻 | **全塔で同じ `snap_time` を使う**（時間軸統一） |
| `frame_interval` | 既定 10 min（`_plot_chain_adsorption_profile_gif` の引数で変更可） |
| `gif_duration_ms` | 1 フレームあたり 800 ms |
| 描画内容 | 各塔とも左軸: 吸着量 q（積み上げ面）、右軸: 液相濃度 C（線） |
| 出力場所 | `outputs_chain_v3/<label>/batch_X/chain_adsorption_profile.gif` |

**実装場所**: `run_chain_predict_v3._plot_chain_adsorption_profile_gif`

**狙い**: 上流塔で破過した成分が下流塔のサイトに順次充填されていく様子を、同時刻の輪切りで直接比較できる。各塔単独 GIF (§8.1) は従来通り維持されており、こちらと併用可能。

### 8.5 フィット結果

```
results_v3/
  ├── R003/<date>/
  │   ├── fitted_vector.json       q_total, eps_b の値
  │   ├── de_progress.json         DE 収束履歴
  │   └── readme.txt               バッチ番号 / 計算時間 / 境界条件
  ├── R004A/<date>/
  │   └── （同上）
  └── R004B/<date>/
      ├── used_parameters.json     R004A から流用したパラメータ
      └── readme.txt               予測実行日時・参照元
```

`results_v3/R003/latest/`, `results_v3/R004A/latest/` に最新版コピーを保持し、可視化スクリプトから安定参照する。

---

## §9 想定する妥当性チェック

### 9.1 樹脂物性の整合性

R003 と R004A は同じ樹脂のはずなので、`q_total` と `eps_b` の値は近いことが期待される。大きく違えば次の可能性を検討:

- **q_total 差大**: 塔ごとに劣化進行度が違う（先住塔 R003 の方が劣化大の可能性）
- **eps_b 差大**: 充填密度の違い、または充填量の管理誤差
- **両方差大**: モデルが捉えきれていない物理（粒径分布、流路偏流など）

### 9.2 R004B 予測の物理的妥当性

R004B 出口の VE1, VE2 が「ほぼゼロ」になっているか確認。万一濃度が高く出る場合は:

- R004A の樹脂が想定より劣化している（fit_targets を増やして再フィット）
- 流量が多すぎる（実運用条件の見直し）
- モデル構造の限界

### 9.3 物質収支

各塔の (入量) - (出量) - (吸着量) ≈ 0 が成立しているか、`chain_summary.csv` で検算する仕組みを入れる。

---

## §10 今後の拡張余地（実装スコープ外）

| 項目 | 概要 |
|:-----|:-----|
| 塔の本数増減 | column_chain_v3 は 3 塔固定ではなく可変長リスト対応にしておく |
| パラメータ共通化フィット | 「樹脂物性は共通」という強い制約を入れた同時フィット（オプション） |
| 同時フィット切替 | R003+R004A を 4 変数（q_total×2, eps_b×2）で同時最適化 |
| バックフラッシュ・再生工程 | 現在は 1 工程内のシミュレーションのみ。再生サイクルは別途検討 |

---

## §11 観測データの切替（実データ ↔ 人工テストデータ）

### 11.1 ファイル構成

| ファイル | 役割 |
|:--------|:------|
| `data/r004a_observed_batches.csv` | 実データ用（空テンプレ。実値入手後に記入） |
| `data/r004a_observed_batches_test.csv` | 人工テストデータ（コード動作確認・スモークテスト用） |
| `data/r004b_observed_batches.csv` | 実データ用（空テンプレ） |
| `data/r004b_observed_batches_test.csv` | 人工テストデータ |

### 11.2 切替方法（コマンドライン引数 `--data-suffix`）

各実行スクリプトに `--data-suffix` 引数を用意し、ファイル名末尾に付与する文字列を指定する。

```bash
# 人工テストデータで動作確認（_test サフィックス付きのファイルを読む）
python run_fit_R004A_v3.py --data-suffix _test
python run_chain_predict_v3.py --data-suffix _test

# 実データで本番実行（サフィックスなし、デフォルト）
python run_fit_R004A_v3.py
python run_chain_predict_v3.py
```

スクリプト内では下記のように構築する：

```python
parser.add_argument("--data-suffix", default="",
                    help="観測データファイル名のサフィックス（例: '_test' で _test.csv を読む）")
args = parser.parse_args()
suffix = args.data_suffix
r4a_csv = f"data/r004a_observed_batches{suffix}.csv"
r4b_csv = f"data/r004b_observed_batches{suffix}.csv"
```

### 11.3 人工テストデータの設計方針

| 塔 | 設計根拠 |
|:---|:---------|
| R004A | R003 出口（v2 実観測値）の 30〜50% 程度に低減。観測時刻は R003 のタイムラインに合わせ、各バッチ 3〜4 点 |
| R004B | VE 漏れなし方針を反映し、VE1/VE2 は 1e-4 オーダー、FA は 1e-3 オーダーの微小値。観測点はバッチ後半 1 点 |

人工データの目的は**コード動作確認**であり、フィット結果が物理的に正しいことは保証しない。実データ入手後に上書きすること。

---

## §12 実装ステップ（チェックリスト）

- [x] v3 フォルダ作成 + v2 コアファイルコピー
- [x] R004A・R004B 観測データの空テンプレ CSV 作成
- [x] R004A・R004B 人工テストデータ作成（v2 実データを参考に）
- [x] `competitive_adsorption_v3.py` に `C_in_series` 拡張（v2 完全互換確認済み）
- [x] `params_columns_v3.py`（塔別 Params 生成）作成
- [x] `column_chain_v3.py`（3 塔リレー実行コア）作成（リレー濃度差 0 確認済み）
- [x] `parameter_fitting_v3.py` に inlet_provider / flow_provider フック追加
- [x] `evaluate_model_accuracy_v3.py` に inlet_provider / flow_provider フック追加
- [x] `run_fit_R003_v3.py` 作成
- [x] `run_fit_R004A_v3.py` 作成
- [x] `run_chain_predict_v3.py`（chain_summary.csv + chain_overlay.png 含む）作成
- [x] `run_chain_visualize_v3.py` 作成
- [ ] **次フェーズ**: R003 単塔フィット動作検証（v2 と同等の結果を得る）
- [ ] **次フェーズ**: 人工テストデータでの 3 塔チェーン動作検証
- [ ] **次フェーズ**: 実データ入手後、R004A フィット → 3 塔予測の本番運用

---

## 付録 A: 動作検証済み事項

実装フェーズで実施したスモークテストの結果を記録する。各テストは v2 から v3 への移行で
**既存挙動を壊していないこと**、および**新規機能（時系列入口濃度・3 塔リレー）が正しく動作すること**
を確認するためのもの。

### A.1 v2 後方互換性

`competitive_adsorption_v3.py` は v2 のコピーから派生したため、`C_in_series=None` の場合は
v2 と完全に同じ動作になることを検証済み。

| テスト項目 | 結果 |
|:----------|:-----|
| `Params(C_in_series=None)` で `C_in_at(t)` が固定値を返す | ✅ OK |
| 同じ入口条件を v2 と v3 で `simulate()` した結果の出口濃度差 | ✅ **0.00e+00**（完全一致） |
| `parameter_fitting_v3._objective_batches` に `inlet_provider=None`, `flow_provider=None` を渡したときの目的関数値 | ✅ 既存 v2 と同値 |

### A.2 C_in_series 拡張機能

| テスト項目 | 結果 |
|:----------|:-----|
| `Dict[species, (t_array, v_array)]` 形式の指定 | ✅ 線形補間動作 |
| `pd.DataFrame(columns=[t_min, VE1, VE2, FA])` 形式の指定 | ✅ 動作 OK |
| 一部成分のみ時系列、残りは固定値 `C_in` の併用 | ✅ 期待通り |
| 補間範囲外は端値で外挿（`np.interp(left=v[0], right=v[-1])`） | ✅ 動作 OK |
| `simulate()` 内部で時系列入口を反映できる | ✅ 動作 OK |

### A.3 3 塔リレーの整合性

3 塔チェーン (`R003 → R004A → R004B`) の最小スモークテスト（共通 Params, t_end=30 min, N=10）:

| テスト項目 | 結果 |
|:----------|:-----|
| 各塔の `simulate()` 完了 | ✅ OK |
| **R003 出口 vs R004A 入口の濃度差** | ✅ **0.00e+00**（リレー時に濃度破壊なし） |
| 物理的妥当性（下流塔ほど VE1 max が小さい） | ✅ R003 (4.6e-6) → R004A (2.0e-9) → R004B (8.4e-13) |

### A.4 スクリプト起動

| スクリプト | `--help` 起動 |
|:----------|:--------------|
| `run_fit_R003_v3.py` | ✅ |
| `run_fit_R004A_v3.py` | ✅ |
| `run_chain_predict_v3.py` | ✅ |
| `run_chain_visualize_v3.py` | ✅ |

### A.5 人工テストデータでの通し動作（Step1→2→3 全通し）

人工テストデータ (`*_test.csv`, `--data-suffix _test`) を使った全パイプライン検証。

| テスト項目 | 結果 |
|:----------|:-----|
| `run_fit_R003_v3.py --train_batches 6 --t_end 460` 完走 | ✅ 約 3 分 |
| `run_fit_R004A_v3.py --train_batches 6 --t_end 460 --data-suffix _test` 完走 | ✅ 約 3 分 |
| `run_chain_predict_v3.py --batches 6 --t_end 460 --data-suffix _test` 完走 | ✅ 約 25 秒 |
| 全 9 サブプロット (`chain_overlay.png`) で物理的妥当な破過カーブ | ✅ |
| 各塔の `adsorption_profile.gif` 生成 | ✅ |
| `chain_summary.csv` に各塔の RMSE が記録 | ✅ |

### A.6 実データ batch3 単独フィット（本番動作検証）

実観測データ（`プロセスデータ処理/data/*batch3*.xlsx` から自動抽出した CSV）で
Step1→2→3 を batch3 単独実行した結果。

| テスト項目 | 結果 |
|:----------|:-----|
| `run_fit_R003_v3.py --train_batches 3 --t_end 460` 完走 | ✅ 4m 28s, DE 7 反復 |
| `run_fit_R004A_v3.py --train_batches 3 --t_end 460` 完走（実 CSV） | ✅ 8m 6s, DE 15 反復 |
| `run_chain_predict_v3.py --batches 3 --t_end 460` 完走（実 CSV） | ✅ 約 46 秒（GIF 含む） |
| R003 RMSE | ✅ VE1=0.010, VE2=0.008, FA=0.032 |
| R004A RMSE | ✅ VE1=0.025, VE2=0.008, FA=0.003 |
| R004B RMSE（観測 1 点） | ✅ VE1=0.0003, VE2=1e-6, FA=0.002 |
| 物理的妥当性: R003→R004A→R004B でピーク時刻が約 100 min ずつ遅延 | ✅ |
| 物理的妥当性: VE1 オーバーシュートが Cin に収束（完全破過） | ✅ |

### A.7 補助図表機能

| テスト項目 | 結果 |
|:----------|:-----|
| `overlay_*.png` に feed Cin の水平点線（灰色破線）が表示 | ✅ |
| `chain_overlay.png` の全 9 サブプロットで feed 点線併記 | ✅ |
| `Cin <= 0` の場合は点線スキップ（R004A 評価時の入口≈0 を誤表示しない） | ✅ |
| `chain_adsorption_profile.gif` が batch ごとに 1 個生成 | ✅ |
| 同 GIF で R003/R004A/R004B が縦並び、全塔とも同じ snap_time | ✅ |

### A.8 まだ検証されていない項目

| 項目 | 内容 |
|:-----|:-----|
| **多バッチ学習（実データ）** | `--train_batches 2,3,4,6` で R003 / R004A をまとめて学習し、共通 (q_total, eps_b) を求める運用 |
| **batch 単独間のばらつき評価** | 各 batch ごとの単独フィット結果を比較し、樹脂劣化進行度の塔・ロット間差を定量化 |
| **fit_targets 拡張** | `q_total`/`eps_b` 以外（例: `H`, `Keq_ads`）を可変にした場合の収束性 |

---

## 付録 B: 使い方ガイド

### B.1 全体フロー

```
[初回セットアップ]
  ├─ data/process_runs_batches.csv : R003 入口・出口（v2 から流用済み）
  ├─ data/r004a_observed_batches.csv : R004A 観測値を記入（実データ）
  └─ data/r004b_observed_batches.csv : R004B 観測値を記入（実データ）

[Step 1] R003 単独フィット
  python run_fit_R003_v3.py --train_batches 2,4,6
  → results_v3/R003/<date>/fitted_vector.json

[Step 2] R004A フィット（R003 の出口を入口として使う）
  python run_fit_R004A_v3.py --train_batches 2,4,6
  → results_v3/R004A/<date>/fitted_vector.json

[Step 3] 3 塔チェーン予測（R004B は R004A 流用）
  python run_chain_predict_v3.py --batches 2,3,4,6
  → outputs_chain_v3/chain/batch_X/{R003,R004A,R004B}/...
  → outputs_chain_v3/chain/chain_summary.csv

[再可視化] 保存済みパラメータで描画のみ（フィットしない）
  python run_chain_visualize_v3.py --batches 6
```

### B.2 各スクリプトの引数

#### B.2.1 `run_fit_R003_v3.py`

| 引数 | 説明 | デフォルト |
|:-----|:-----|:---------|
| `--t_end` | シミュレーション終了時間 [min] | Params のデフォルト（350 min） |
| `--train_batches` | 学習用バッチ番号（カンマ区切り） | `6` |
| `--test_batches` | 検証用バッチ番号（カンマ区切り） | train と同じ |
| `--species` | フィット対象の成分（カンマ区切り、例 `VE1,VE2`）。詳細は付録 D | 未指定=観測がある全成分 |

```bash
# 単一バッチ学習
python run_fit_R003_v3.py --train_batches 6

# VE1, VE2 のみを目的関数に含める（FA を除外）
python run_fit_R003_v3.py --train_batches 6 --species VE1,VE2

# 多バッチ学習（共通パラメータを 1 回で同定）
python run_fit_R003_v3.py --train_batches 2,4,6

# 学習と検証を分ける
python run_fit_R003_v3.py --train_batches 2,4 --test_batches 6

# t_end を延長して通水停止後の挙動も含める
python run_fit_R003_v3.py --train_batches 6 --t_end 800
```

#### B.2.2 `run_fit_R004A_v3.py`

| 引数 | 説明 | デフォルト |
|:-----|:-----|:---------|
| `--t_end` | シミュレーション終了時間 [min] | Params のデフォルト |
| `--train_batches` | 学習用バッチ番号 | `6` |
| `--test_batches` | 検証用バッチ番号 | train と同じ |
| `--data-suffix` | 観測 CSV のサフィックス | `""`（実データ） |
| `--species` | フィット対象の成分（カンマ区切り、例 `VE1,VE2`）。詳細は付録 D | 未指定=観測がある全成分 |

```bash
# 多バッチで R004A を学習（実データ）
python run_fit_R004A_v3.py --train_batches 2,4,6

# VE1, VE2 のみを目的関数に含める（FA を除外）
python run_fit_R004A_v3.py --train_batches 6 --species VE1,VE2

# 人工テストデータで動作確認
python run_fit_R004A_v3.py --train_batches 6 --data-suffix _test
```

実行に先立ち `results_v3/R003/latest/fitted_vector.json` が必要（先に Step 1 を実施）。

#### B.2.3 `run_chain_predict_v3.py`

| 引数 | 説明 | デフォルト |
|:-----|:-----|:---------|
| `--t_end` | シミュレーション終了時間 [min] | Params のデフォルト |
| `--batches` | 予測対象バッチ番号 | `2,3,4,6` |
| `--columns` | 予測する塔 | `R003,R004A,R004B` |
| `--data-suffix` | 観測 CSV のサフィックス | `""` |
| `--label` | 出力フォルダのラベル | `chain` |

```bash
# 全 4 バッチ・3 塔フル予測
python run_chain_predict_v3.py

# R003+R004A までの 2 塔予測
python run_chain_predict_v3.py --columns R003,R004A

# 人工テストデータで通し動作確認
python run_chain_predict_v3.py --data-suffix _test --label chain_test

# 別ラベルで結果を残す（過去結果を保護）
python run_chain_predict_v3.py --label chain_2026Q2_run1
```

### B.3 多バッチ学習について

**v3 は単一バッチ学習も多バッチ学習も同じスクリプトで対応**。`--train_batches` をカンマ区切りで複数指定するだけ。

#### 仕組み

`parameter_fitting_v3._objective_batches` が以下を行う：

1. 各バッチについて個別に `simulate()` を実行
2. バッチごとに観測値とモデル予測値の **RMSE** を計算
3. 全バッチの RMSE を**単純平均**して目的関数値とする

つまり「全ロットで共通の (q_total, eps_b)」を 1 回の DE 最適化で求める。

#### 注意点

| 項目 | 説明 |
|:-----|:-----|
| **計算時間** | バッチ数に比例（例: 5 バッチで 5 倍） |
| **観測点の偏り** | 各バッチの RMSE が等しい重みで平均（観測点数による重み付けなし） |
| **物理整合性** | 全バッチで同じ樹脂物性を仮定するため、ロット間で挙動差が大きいとフィット精度が低下 |
| **目的関数の解釈** | 値は「平均 RMSE」であり、絶対値の良さは観測点ごとの誤差で判断 |

#### 推奨手順

1. **まず 1 バッチで確認**：`--train_batches 6` などで 1 バッチだけフィットし、結果が物理的に妥当か確認
2. **小さい複数バッチで試行**：`--train_batches 2,6` などで 2 バッチに広げる
3. **全バッチで本番**：`--train_batches 2,4,6` などで全ロット学習

### B.4 出力ディレクトリ構造

```
results_v3/
├── R003/
│   ├── 2026-06-08/                       実行日付フォルダ
│   │   ├── fitted_vector.json
│   │   ├── de_progress.json
│   │   └── readme.txt
│   ├── 2026-06-08_2/                     同日に複数回実行した場合は _2, _3
│   └── latest/
│       └── fitted_vector.json            最新版コピー（自動更新）
├── R004A/
│   └── (同上)
└── R004B/
    └── latest/
        └── used_parameters.json          R004A から流用したパラメータの記録

outputs_chain_v3/
├── chain/                                run_chain_predict_v3.py の出力
│   ├── batch_2/
│   │   ├── R003/
│   │   │   ├── overlay_VE1.png           feed Cin の水平点線つき (§8.3)
│   │   │   ├── overlay_VE2.png
│   │   │   ├── overlay_FA.png
│   │   │   ├── adsorption_profile.gif    各塔単独 GIF（従来通り）
│   │   │   ├── model_output.csv
│   │   │   └── observed.csv
│   │   ├── R004A/
│   │   ├── R004B/
│   │   ├── chain_overlay.png             3 塔×3 成分の俯瞰図 + feed 点線
│   │   └── chain_adsorption_profile.gif  3 塔縦並び・時間軸統一の統合 GIF (§8.4)
│   ├── batch_3/, batch_4/, batch_6/
│   └── chain_summary.csv                 全バッチ・全塔・全成分の RMSE/MAPE まとめ
├── visualize/                            run_chain_visualize_v3.py のデフォルト出力先
└── chain_2026Q2_run1/                    --label で個別保存
```

### B.5 トラブルシューティング

| 症状 | 対処 |
|:-----|:-----|
| `FileNotFoundError: results_v3/R003/latest/fitted_vector.json` | Step 1 (`run_fit_R003_v3.py`) を先に実行する |
| `R004A 観測点が 0` | `data/r004a_observed_batches.csv` に値が入っているか確認、または `--data-suffix _test` で人工データを試す |
| フィット結果が境界に張り付く | `parameter_fitting_v3.DEFAULT_BOUNDS` の該当パラメータの境界を緩める |
| 数値振動が出る | `params_columns_v3.apply_dt_react_cap_all()` で `dt_react_cap` を 0.05 や 0.03 に縮める |
| R004B 予測が VE 漏れあり | R004A のフィット結果が想定外（樹脂劣化大）の可能性。R004A の overlay 図と物理整合性を確認 |

---

## 付録 C: シミュレーション終了時間 `t_end` の自動決定（バッチ個別）

### C.1 背景

従来は `--t_end` 未指定時に `Params` のデフォルト（350 min）を全バッチ一律で使用していた。
バッチごとに運転時間が異なる（例: batch2≒460 min、batch6≒420 min）ため、350 固定では
**長いバッチは途中で打ち切られ**（後半の観測点が RMSE から脱落）、**短いバッチは運転終了後まで
過剰計算**する問題があった。

### C.2 方針

- **`--t_end` 未指定時**: 各バッチの `t_min` 範囲に合わせて **バッチごとに個別に** `t_end` を自動決定する。
- **`--t_end` 指定時**: 従来通り、その値を **全バッチ一律** で使用（自動決定は無効）。

複数バッチを学習・予測する場合も、全バッチ共通の 1 値（最大や最小）にはせず、
**各バッチが自分自身のデータ長を使う**（バッチ個別）。これは fit / predict とも
バッチごとに独立して `simulate()` しているため自然かつ精度上も適切なため。

### C.3 自動決定式

```
base   = max( プロセスデータ t_min 最大, 観測データ t_min 最大 )
t_end  = base × T_END_MARGIN_FACTOR      (既定 1.0 = マージン無し)
```

- 基準は「そのバッチ・その塔で有効データが存在する範囲」をカバーするよう、
  プロセス（`process_runs_batches.csv` の流量時系列）と観測（各塔の `C_*_exp`）の
  **t_min 最大の大きい方**を採用する。
- マージン `T_END_MARGIN_FACTOR = 1.0`（`params_columns_v3.py`）＝**マージン無し**。
  `t_end = t_min.max` ちょうどとし、工程終了（通液停止）の時刻で描画を切る（末尾流量の外挿区間を作らない）。
- 工程終了後の内部緩和（交換反応の進行）を可視化したい場合は、この定数を 1.05 などに増やす。
  増やした区間は下記 C.3.1 の通り**流量 0**として扱われる。

### C.3.1 工程終了後の流量の扱い（通液停止）

`t_min.max` は物理的に「工程終了（通液停止）」を意味する。`Params.vT_at` は流量時系列の
範囲外を末尾値のまま外挿するため、`t_min.max` を超える時刻（マージンを付けた場合や
`--t_end` で延長した場合）に「最後の流量のまま通液が続く」状態になってしまう
（実データの末尾 `v_T` は 0 ではなく ~30–43 L/h）。

これを避けるため、流量時系列の末尾（`t_min.max` の直後）に **`v_T = 0` の点を付加**し、
工程終了後は**流量 0**（移流・分散ゼロ、反応のみ進行）として扱う。実装は
`params_columns_v3.append_flow_stop(flow_df)`。各実行スクリプトが `simulate` へ渡す流量時系列に
適用する（fit・predict・visualize すべて）。

- 既定のマージン無し（1.0）では `t_end = t_min.max` なので、この流量 0 区間は実質的に描画されない
  （安全策として、`--t_end` 延長時やマージンを付けた場合にのみ効く）。
- 観測点は `t_min.max` 以前にしか無いため、**フィットの目的関数（RMSE）には影響しない**。
- 有効な t_min が 1 つも無い場合は `None` を返し、呼び出し側で `Params` デフォルト（350 min）へフォールバックする。

### C.4 実装

| 場所 | 内容 |
|:-----|:-----|
| `params_columns_v3.compute_auto_t_end(*t_min_arrays, margin_factor=...)` | 複数 t_min 配列の最大 × マージンを返す共通ヘルパー。`T_END_MARGIN_FACTOR` も同モジュールで定義。 |
| `params_columns_v3.append_flow_stop(flow_df)` | 流量時系列の末尾（工程終了直後）に `v_T=0` を付加し、工程終了後を通液停止として扱う。fit（R003 の `flow_provider`、R004A の `flow_provider_common` / R003 出口生成）・predict の流量作成に適用。 |
| `parameter_fitting_v3.fit_parameters_batches(..., t_end_provider=None)` | 目的関数 `_objective_batches` 内で、バッチごとに `t_end_provider(bid, dfb)` を評価し `p.t_end` に反映。`None` を返すとグローバル値のまま。 |
| `evaluate_model_accuracy_v3.evaluate_model_batches(..., t_end_provider=None)` | 評価（overlay/GIF/CSV）でも同様にバッチごとに `t_end` を反映。 |
| `run_fit_R003_v3.py` | 未指定時、`dfb["t_min"]` から自動決定する `t_end_provider` を fit と評価に渡す。 |
| `run_fit_R004A_v3.py` | 未指定時、プロセスと R004A 観測の t_min から `t_end_map`（バッチ→t_end）を作成。R003 出口生成・R004A フィット・評価すべてで同じ per-batch t_end を使用。 |
| `run_chain_predict_v3.py` | 未指定時、プロセスと予測対象塔の観測の t_min からバッチ専用 `t_end` を算出して 3 塔に適用。 |

### C.5 記録

- フィット結果 JSON (`fitted_vector.json`) に `t_end_mode`（`"fixed"` / `"auto_per_batch"`）と
  `t_end_desc`（人間可読な説明・バッチ別 t_end の一覧）を追加。従来の `t_end`（スカラー）は
  代表値（自動時は train バッチの最大）として維持。
- `readme.txt` の `t_end` 行にも同じ説明を出力する。

---

## 付録 D: フィット対象成分の選択（`--species`）

### D.1 背景（追加経緯）

R004A のフィット結果を検討する中で、次の点が確認された。

- 観測濃度のスケールが成分間で大きく異なる（batch6 R004A で VE1 ≈ 0.05–0.08、VE2 ≈ 0.01–0.03、
  **FA ≈ 0.0017–0.0039 と VE1 の約 1/20 の希薄域**）。
- 目的関数は「(バッチ×成分) ごとの絶対 RMSE の単純平均（成分等重み）」であり、FA は希薄ゆえ
  目的関数への寄与が小さい（R004A の例で約 4%）。すなわち **FA はフィット結果にほとんど影響しない**。
- 「希薄で当てはめの意味が薄い成分（FA など）を目的関数から外したい」という運用要望が生じた。

そこで、**どの成分を目的関数に含めるかを実行時に選択できる機能**を追加した。従来から
`parameter_fitting_v3` の関数は `target_species` 引数で成分を限定できたが、実行スクリプトからは
`None`（＝観測がある全成分を自動採用）固定で、選択手段が公開されていなかった。本機能はこれを
コマンドライン引数 `--species` として公開したもの。

### D.2 仕様

- 対象スクリプト: **`run_fit_R003_v3.py` / `run_fit_R004A_v3.py`**（フィットを行う 2 本）。
  `run_chain_predict_v3.py` はフィットしないため対象外。
- `--species VE1,VE2` のように**カンマ区切り**で指定。指定した成分だけを**目的関数（RMSE）**の計算に用いる。
- **未指定時は従来通り** `None` = 観測列 `C_*_exp` がある全成分を自動採用（後方互換）。
- **評価図・metrics は指定に関わらず全観測成分を表示**する。
  → フィット対象から外した成分（例 FA）が、フィット後にどう振る舞うかを overlay で確認できる。
- 妥当性チェック:
  - 観測列が無い成分名を指定した場合は警告して無視。
  - 有効な成分が 0 個になった場合は明確なエラーで停止。

> **注意（重み付けではない）**: 本機能は「含める／含めない」の二択であり、成分ごとの連続的な重み付けは
> できない（含めた成分は等重みの RMSE 平均）。VE を相対的に重視するといった重み付けが必要な場合は
> 目的関数側の拡張が別途必要。

### D.3 使用方法

**指定できる成分は `VE1` / `VE2` / `FA` の 3 種類**（観測列 `C_*_exp` を持つ成分）。
キャリア成分 `OH` は観測列を持たずフィット対象外のため指定できない（指定しても警告して無視）。

| 成分名 | 内容 |
|:-------|:-----|
| `VE1` | ビタミンE（Toc-α + T3-α） |
| `VE2` | ビタミンE（Toc-β + T3-β） |
| `FA`  | 遊離脂肪酸（競合成分） |

```bash
# VE1, VE2 のみをフィット対象にする（FA を目的関数から除外）
python run_fit_R003_v3.py  --train_batches 6 --species VE1,VE2
python run_fit_R004A_v3.py --train_batches 6 --species VE1,VE2

# 単一成分だけを対象にする（VE1 のみ）
python run_fit_R003_v3.py  --train_batches 6 --species VE1

# 3 種類すべてを明示指定（未指定時と同じ挙動）
python run_fit_R003_v3.py  --train_batches 6 --species VE1,VE2,FA

# 未指定（従来どおり、観測がある全成分を対象）
python run_fit_R003_v3.py  --train_batches 6
```

### D.4 実装

| 場所 | 内容 |
|:-----|:-----|
| `parameter_fitting_v3.fit_parameters_batches(..., target_species=None)` | 目的関数 `_objective_batches` に含める成分リスト。`None` なら観測列がある全成分を自動採用（従来からの機能）。 |
| `run_fit_R003_v3.py` / `run_fit_R004A_v3.py` | `--species` を解析して `target_species` に渡す。観測列（`C_*_exp`）と突き合わせて検証（未知成分は警告、空ならエラー）。評価（`evaluate_model_batches`）へは渡さず、図・metrics は全観測成分を表示。 |

### D.5 記録

- フィット結果 JSON (`fitted_vector.json`) に `fit_species`（成分リスト、未指定時は `"auto"`）と
  `fit_species_desc`（人間可読な説明）を追加。
- `readme.txt` に `Fit species` 行を出力。
- 実行時に標準出力へ `Fit species : ...` を表示。

---

## 付録 E: 工程終了時刻ちょうどの観測が評価から欠落する問題と対応（端の許容）

### E.1 症状

複数ロット学習・予測（例: ESX012 と 004 を学習し ESX012 を予測）の `chain_summary.csv` で、
**R004B の RMSE/MAPE が計算されず、結果に行が出力されない**事象が発生した。
R004B は観測が「工程終了時刻ちょうどの 1 点」（例: batch12 で `t_min=426`）しか無いため、
その 1 点が評価から外れると R004B 全体が欠落する。R004A でも同時刻の 1 点が落ち、
`n_obs` が 1 つ少なくなっていた（例: 8 点 → 7 点）。

### E.2 原因

1. `simulate` は `while t < t_end + 1e-12` で `t += dt` と時間積分する。
   dt の加算で**浮動小数点誤差が累積**し（例: `t=425.9000000000283`）、次ステップが
   ガード `t_end + 1e-12` をわずかに超えるため、**工程終了 t_end ちょうど（例 426）の
   サンプルが記録されない**。結果、モデル出力の最終時刻は `t_end - dt`（例 425.9）になる。
2. 評価は観測時刻へモデル値を `np.interp(..., right=np.nan)` で補間しており、
   モデル範囲外（= t_end ちょうど）の観測は **NaN 扱いで比較対象から除外**されていた。
3. R004B は該当の 1 点しか観測が無いため、除外により metrics が空になり、
   `chain_summary.csv` に行自体が出力されなかった。

### E.3 対応（案A: 評価側で端の許容）

評価時の補間を共通ヘルパー **`evaluate_model_accuracy_v3.interp_obs_to_model()`** に置き換えた。
観測時刻がモデル最終時刻を**わずか（許容幅 `tol` 以内）**超える場合は、刻み誤差による
最終サンプル欠落とみなして**端（最小/最大時刻）へスナップ**し、モデル端値と比較する。
`tol` 未指定時は「モデル時間刻みの中央値 × 1.5」（＝約 1 ステップ分の欠落を吸収）。
許容幅を超える真の範囲外観測は、従来どおり NaN（除外）とする。

これにより、工程終了 t_end ちょうどの観測（R004B の唯一点など）が
モデル最終時刻（t_end−dt）の値と正しく比較され、RMSE/MAPE が計算されるようになる。

### E.4 実装

| 場所 | 内容 |
|:-----|:-----|
| `evaluate_model_accuracy_v3.interp_obs_to_model(t_obs, t_mod, y_mod, tol=None, tol_factor=1.5)` | 端の許容つき補間ヘルパー（新規）。範囲外でも許容幅内なら端へスナップ。 |
| `evaluate_model_accuracy_v3._evaluate_one_set` / `plot_overlay_sparse` | metrics 計算・overlay 図の補間をこのヘルパーに置換（R003 / R004A フィット評価に反映）。 |
| `run_chain_predict_v3._eval_overlay_one_column` | チェーン予測の metrics（`chain_summary.csv`）の補間をこのヘルパーに置換。 |

### E.5 注意（フィット目的関数は対象外）

本対応は**評価（metrics・overlay）側のみ**に適用している。フィットの目的関数
（`parameter_fitting_v3._objective_batches`）の補間は従来どおり `right=np.nan`（端の許容なし）で、
**工程終了 t_end ちょうどの観測はフィットには使われない**。

- R004B は予測専用（フィットしない）ため、今回の症状は本対応で完全に解消する。
- R003 / R004A で「t_end ちょうどの観測点もフィットに含めたい」場合は、フィット目的関数側にも
  同じ端の許容を適用する拡張が別途必要（現状は未適用。挙動を変えないため据え置き）。

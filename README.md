# PCP_cursor — 競争吸着パイプライン

陰イオン交換樹脂カラムの競争吸着プロセスをモデリングするリポジトリです。

## 構成

| フォルダ | 内容 |
|:---|:---|
| `v3/` | 3 塔直列（R003→R004A→R004B）チェーンモデル |
| `v4/` | v3 + FAEE 加水分解（FAEE + H₂O → FA + ET） |

各フォルダは自己完結しており、仕様書・コード・データを含みます。

## セットアップ

```bash
pip install -r requirements.txt
```

## v3 の使い方（要約）

```bash
cd v3
python smoke_test_v3.py

# 人工テストデータでの通し
python run_fit_R003_v3.py --train_batches 6 --t_end 200 --maxiter 8
python run_fit_R004A_v3.py --train_batches 6 --t_end 200 --data-suffix _test --maxiter 8
python run_chain_predict_v3.py --batches 6 --t_end 200 --data-suffix _test --label chain_test
```

詳細は `v3/competitive_adsorption_v3_spec.md` を参照。

## v4 の使い方（要約）

```bash
cd v4
python smoke_test_v4.py

python run_fit_R003_v4.py --train_batches 6 --t_end 200 --maxiter 8
python run_fit_R004A_v4.py --train_batches 6 --t_end 200 --data-suffix _test --maxiter 8
python run_chain_predict_v4.py --batches 6 --t_end 120 --data-suffix _test --plot_new_species
```

詳細は `v4/competitive_adsorption_v4_spec.md` を参照。

## データ

- `data/process_runs_batches.csv` … R003 入口・出口
- `data/r004a_observed_batches.csv` … R004A 実測（テンプレ）
- `data/r004*_observed_batches_test.csv` … 動作確認用人工データ
- `v4/data/feed_composition_v4.csv` … ロット別 FAEE/ET/水

# wood-bridge — ダビンチアーチ橋 バルサ模型

5C 設計製図（2026年後期）の橋梁模型（支間 500 mm、バルサ材）の設計作業リポジトリ。
構造形式はダビンチアーチ（下路式・タイドアーチ）。現段階は **SVG 簡易図面 → AutoCAD 作図** の前段として、要件整理と簡易デザインを行っている。

## 構成

| パス | 内容 |
| --- | --- |
| [docs/requirements.md](docs/requirements.md) | 課題条件・構造コンセプト・ベース案の幾何・未決定事項 |
| [docs/structural-audit.md](docs/structural-audit.md) | 構造監査（FEM・座屈・形式比較・CP 評価） |
| [docs/design-exploration.md](docs/design-exploration.md) | ダビンチを維持した形状探索（CP 最大化）と推奨案 |
| [scripts/audit_fem.py](scripts/audit_fem.py) | 監査用 2D フレーム FEM（標準ライブラリのみ） |
| [scripts/explore_davinci.py](scripts/explore_davinci.py) | ダビンチ構成の総当たり探索 |
| [scripts/make_draft.py](scripts/make_draft.py) | 設計パラメータ（正本）から図面 SVG と確認用 STL を生成 |
| [work/bridge-draft.svg](work/bridge-draft.svg) | 簡易一般図（側面・平面・断面・部材表・照査結果） |
| [work/bridge-draft.stl](work/bridge-draft.stl) | 形状確認用 3D モデル（mm、Z 上向き） |
| [MEMORY.md](MEMORY.md) | セッションごとの作業記録 |

## 使い方

```sh
python3 scripts/make_draft.py   # Python 3.9 以上、標準ライブラリのみ
```

`work/` の SVG・STL を再生成し、照査結果・部材表・断面組合せ比較を標準出力に表示する。
寸法を変えるときはスクリプト冒頭の「設計パラメータ」を編集して再実行する。
SVG はブラウザで、STL は macOS のプレビュー等の 3D ビューアで確認できる。

## 参照

- 5C 設計製図 実施要領（R8年度、2026-10-14 版 PDF）— 課題条件の原本（リポジトリ未収録）

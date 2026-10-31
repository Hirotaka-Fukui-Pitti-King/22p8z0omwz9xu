# Nirei (2015) — 実装タスクと実行順序

## 目的と完成条件

Nirei (2015), *An Interaction-Based Foundation of Aggregate Investment Fluctuations*, Section 3 と Table 2 の Benchmark simulation を **Julia だけで**実装・実行する。完成条件は **110,000 社 × 1,100 期間 × 10 試行を実際に計算**し、初期 100 期間を除いた9モーメント、10試行の標準誤差、Table 2 との比較表、再実行に必要な Julia 環境・seed と時系列を保存すること。Issue を作ること自体は実装完了ではない。

本計画は [README.md](README.md)、[ALGORITHM.md](ALGORITHM.md)、[policy.json](policy.json) に従う。各 Issue には対応する `tasks/Txx.json` があり、個別作業・成果物・依存関係・完了判定を記録する。

## Issue とタスク

| 順 | Issue | 作業 | 主要成果物 | 先行タスク | タスク JSON |
|---|---|---|---|---|---|
| T01 | [#1](https://github.com/Hirotaka-Fukui-Pitti-King/22p8z0omwz9xu/issues/1) | Julia 環境、定数 A、定常状態 | `Project.toml`, `Manifest.toml`, `src/calibration.jl` | なし | [T01.json](tasks/T01.json) |
| T02 | [#2](https://github.com/Hirotaka-Fukui-Pitti-King/22p8z0omwz9xu/issues/2) | Appendix A の線形合理的期待均衡 | `src/expectations.jl` | T01 | [T02.json](tasks/T02.json) |
| T03 | [#3](https://github.com/Hirotaka-Fukui-Pitti-King/22p8z0omwz9xu/issues/3) | 企業分布、初期資本、生産性乱数 | `src/firms.jl` | T01 | [T03.json](tasks/T03.json) |
| T04 | [#4](https://github.com/Hirotaka-Fukui-Pitti-King/22p8z0omwz9xu/issues/4) | 有限企業の固定点全探索、ES1 | `src/equilibrium.jl` | T01, T03 | [T04.json](tasks/T04.json) |
| T05 | [#5](https://github.com/Hirotaka-Fukui-Pitti-King/22p8z0omwz9xu/issues/5) | 期待形成、ES1、投資ショック、マクロ・企業状態の逐次計算 | `src/simulation.jl`, `scripts/run_benchmark.jl` | T02, T03, T04 | [T05.json](tasks/T05.json) |
| T06 | [#6](https://github.com/Hirotaka-Fukui-Pitti-King/22p8z0omwz9xu/issues/6) | 実際の Benchmark 実行、モーメント、Table 2 比較 | `src/moments.jl`, `results/` | T05 | [T06.json](tasks/T06.json) |

## 依存関係と進行

```text
T01 Julia environment, A, steady state
 ├── T02 Linear rational expectations ─────┐
 └── T03 Firm initialization and RNG ──┐   │
                                      T04  │
                                      ES1  │
                                        └───┴─> T05 Integrated simulator
                                                     |
                                                     v
                                              T06 Full benchmark
```

**推奨実行順序**: T01 → T02 と T03（並行可能）→ T04 → T05 → T06。T04 は T03 後に実施できるため T02 と並行してもよい。Issue は作業を区切るためのものであり、途中段階で「検証だけ」の Issue や無関係なタスクを増やさない。

## 各タスクの実質的な引き渡し条件

**T01**: Julia 標準・公式パッケージ資料を確認し、依存環境を固定して、論文の定常状態と A を計算できる。T02 と T03 が同じ構造パラメータを参照できる。

**T02**: Appendix A の行列 B,D と安定固有解から、期首集計資本を入力すると近似的な期待投資と要素価格が計算できる。情報集合の時点を記録する。

**T03**: Appendix C に基づく `lambda_i`、`a_i,t`、一様ギャップ、初期 `k_i,0` が生成でき、後続期の企業状態と RNG を引き渡せる。企業ごとの `lambda_i` は固定。

**T04**: (13)、(14)、(16) から非線形の企業投資均衡を全閾値の昇順走査で計算し、**ES1** を選択する。単純な fictitious tâtonnement の ES2 を使用しない。

**T05**: T02–T04 を一つの Julia シミュレーションループに統合し、企業資本と `Y,C,X,L,K` の時系列を逐次生成する。厳密な有限企業資本と線形化資本の差は隠さない。

**T06**: 110,000 社 × 1,100 期 × 10 試行を走らせ、5つの対数標準偏差、産出との4相関、10試行間の標準誤差を算出し、論文と比較する。実行できなかった場合、Issue は完了扱いにしない。

## 統一的な実装原則

- **Julia-only**。Python、R、MATLAB、C++ 等で先に書いたコードを Julia へ翻訳しない。実装時には [Julia](https://julialang.org/)、[Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/)、[Pkg](https://pkgdocs.julialang.org/v1/) と各パッケージの一次資料を参照する。
- 期待均衡は連続体・一次近似で、企業投資は有限 N の非線形固定点で解く。資本をすべて線形化しない。
- Benchmark は Section 3.4 の ES1 を使用する。fictitious tâtonnement と ES2 は理論上の別の均衡選択である。
- 論文にある式、論文が指定しない初期化・tie-breaking・近似の接続を分けて記録する。結果に合わせて恣意的な smoothing、clip、外生集計ショックを加えない。
- Julia の関数、具象型配列、再利用する作業メモリを用い、固定点探索の目標計算量はソート込みで `O(N log N)` / 期間とする。
- テストや数値チェックは本体実装に組み込み、検証専用の大規模タスク群へ分岐しない。失敗時は当該 Issue 内で原因に直結する修正を完了する。
- Issue を閉じるときは、**具体的なコミット、実行コマンド、生成物、未解決の制限**を記載する。文書の追加だけでコードを実装済みと扱わない。

## 統計的再現対象

| モーメント | Benchmark |
|---|---:|
| SD(log Y) | 2.28% |
| SD(log C) | 0.86% |
| SD(log X) | 6.60% |
| SD(log L) | 3.06% |
| SD(log K) | 0.26% |
| Corr(Y,C) | 0.822 |
| Corr(Y,X) | 0.979 |
| Corr(Y,L) | 0.795 |
| Corr(Y,K) | −0.018 |

最終タスクの完了は、**理論上の計算可能性ではなく実際に出力した比較表**で判断する。論文値と異なるときは差を表示し、推測で「再現成功」としない。

## 現在のステータス

2026-11-01 に6件の GitHub Issue と各タスクの JSON を作成。**全タスクは実装前（open）**であり、Julia コード・実行結果はまだ存在しない。進捗の source of truth は個別の GitHub Issue と実際のコミット・成果物とする。

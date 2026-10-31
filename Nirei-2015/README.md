# Nirei (2015) — Aggregate Investment Fluctuations: Julia Replication

## 研究対象と目的

Makoto Nirei (2015), “An Interaction-Based Foundation of Aggregate Investment Fluctuations,” *Theoretical Economics* 10(3), 953–985, DOI: [10.3982/TE1611](https://doi.org/10.3982/TE1611) の Section 3 および Table 2 を、**Julia で再現する**ための研究・実装仕様を保存するフォルダです。

研究の核心は、企業固有の生産性ショック、不可分な企業別投資、戦略的補完性から、外生的な集計ショックを与えずに内生的な集計投資変動が生じる仕組みを再現することです。最終的な完成条件は **Benchmark calibration、企業数 110,000、1,100 期 × 10 試行を実行し、初期 100 期を除く Table 2 の標準偏差と相関を出力すること**です。

これは実装前の仕様文書です。**このコミットでは Julia ソース、実行結果、Table 2 の再現達成を主張しません。**

## 文書

- [ALGORITHM.md](ALGORITHM.md): 論文の均衡条件、定常状態、期待均衡の対数線形化、企業資本の初期化、固定点の全探索、ES1、内生投資ショック、マクロ変数の逐次更新、Monte Carlo モーメントの計算、論文で未指定の実装事項。
- [policy.json](policy.json): プログラミング言語、Julia-native 実装、公式資料の優先順位、モデルの解釈、数値計算の条件、完了条件、禁止事項、成果物仕様を機械可読形式で固定。

## 数値アルゴリズムの全体像

1. **事前計算**：構造パラメータと企業分布を設定し、定数 A、定常状態、および Appendix A の 2×2 線形合理的期待均衡を解く。安定固有値と saddle-path 係数を一度だけ求める。
2. **試行の初期化**：各企業の不可分性 λ_i を一度だけ生成・固定し、生産性 a_i,0 と資本ギャップ s_i,0 ~ U[0,1) から企業資本ベクトルを構成する。
3. **各期の期待形成**：期首の集計資本から近似的な期待資本と期待投資を計算する。これは連続体近似による expectation system であり、有限企業の実現値とは区別する。
4. **企業ショックと非連続投資**：a_i,t+1 を独立に生成し、式 (16) の投資閾値を計算する。有限企業の式 (13)、(14) の固定点を解く。
5. **均衡選択**：Section 3.4 に従って、実現投資と期待投資の対数乖離の絶対値が最小になる **ES1** を採用する。Section 2.4 の単純な fictitious tâtonnement は **ES2** に対応し、Benchmark replication の ES1 の代わりに使わない。
6. **集計ショックと遷移**：ε_t = log(X_t) − log(X_t^e) を計算し、企業別資本は非線形の正確な式で、マクロ変数は Appendix A の一次近似で更新する。
7. **統計量**：1,100 期のうち最初の 100 期を破棄。10 試行の各モーメントを計算し、その平均と Monte Carlo standard errors を Table 2 と比較する。

期待側の連続体近似と、実現側の有限 N の和を同じものとして扱わないことが重要です。とくに資本の近似的な線形遷移式と、企業資本の正確な集計値を区別します。

## Benchmark（論文 Table 1 / Appendix C）

| 量 | 値 |
|---|---:|
| α | 0.26 |
| δ | 0.02（四半期） |
| β | 0.99 |
| σ | 15 |
| η | 11 |
| ψ | 1 |
| ζ | 2 |
| λ_i | N(1.028, 0.004²) を平均 ±2 標準偏差で切断。各試行で企業ごとに一度だけ生成 |
| log a_i,t | i.i.d. 正規、標準偏差 0.0005。論文指定のモーメントで平均を正規化 |
| N | 110,000 |
| 各試行 | 1,100 期、最初の 100 期を破棄 |
| 試行数 | 10 |

## Table 2 Benchmark の再現ターゲット

| 指標 | 論文 |
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

異なる乱数列による標本誤差があるため、元の時系列のビット単位一致ではなく、論文の Monte Carlo モーメントとの整合を成功基準とします。乱数 seed、使用した Julia のバージョン、依存関係のバージョン、実装上の選択を保存します。

## Julia-only 実装の厳格な方針

**モデルのシミュレーション、均衡計算、数値積分、データ出力、集計統計量は Julia で実装します。Python、R、MATLAB、C++ 等で先に実装して Julia へ翻訳することは禁止します。** 数値手法の設計時点から Julia の型推論、列優先配列、インプレース更新、標準ライブラリ、パッケージ環境を前提とします。

実装前と依存関係の追加前に、[Julia](https://julialang.org/)、[Julia Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/)、[Pkg.jl](https://pkgdocs.julialang.org/v1/)、各パッケージの公式ドキュメント・リポジトリを参照します。

想定する依存対象は、標準ライブラリの LinearAlgebra、Random、Statistics、Printf、TOML、Distributions.jl（切断正規分布）、QuadGK.jl（一次元積分）です。再現可能な乱数ストリームを必要とする場合は StableRNGs.jl を検討します。不要なパッケージは追加しません。依存関係は Nirei-2015 内の Project.toml / Manifest.toml に固定します。

パフォーマンス上は、大規模ループを関数に閉じ込め、具象型配列と事前確保した作業領域を使い、反応関数の切替点をソートして一回の走査で固定点を列挙します。毎候補で企業全体を再集計する O(N²) 実装を採用しません。最適化は結果の正しさを犠牲にしない範囲で行います。

## 作業方針

**完成する 1 本の実行可能な replication を優先します。** 検証のための独立したミニプロジェクト、無関係な文献調査、新規モデル拡張、診断専用機能の増殖は行いません。必要なチェック（有限値、正の集計投資、固定点の区間整合、資本集計、シードの記録）は本体の数値計算に含めます。エラーが出たときだけ原因に直結する修正を行い、Table 2 の再現まで進めます。

## 出典と優先順位

一次資料は [Nirei (2015) 本文](https://doi.org/10.3982/TE1611)、とくに Section 2.3–2.4、Section 3.2–3.5、Appendices A–C です。論文で明記された式・選択ルールと、未指定箇所を埋める実装上の取り決めを [ALGORITHM.md](ALGORITHM.md) で明示的に分けます。前後のチャット回答よりも原論文の式を優先します。


## Issue #1 — Julia calibration のローカル実行

Issue #1 の Julia ソースは \`src/Nirei2015.jl\`、\`src/calibration.jl\`、\`scripts/run_calibration.jl\` にあります。Julia の数値計算は Windows PC で実行します。**Julia 未実行時点では Manifest.toml と calibration.toml をリポジトリに用意していません。**

クローン済みのリポジトリのルートディレクトリから Windows PowerShell で次を実行します。

\`\`\`powershell
git pull --ff-only origin main
Set-Location .\Nirei-2015
julia --version
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/run_calibration.jl
Get-Content .\results\calibration.toml
\`\`\`

最初の \`Pkg.instantiate()\` により、Julia/Pkg が \`Manifest.toml\` を生成します。パッケージ登録情報と依存バージョンが固定されます。次のスクリプトは、式 (16) の A、定常状態、埋め込み均衡残差を計算し、\`results/calibration.toml\` を作成します。初期段階では PC 上の実行結果とログを確認してから GitHub に Manifest と結果を追加し Issue #1 を閉じます。

このステップでの計算は企業別ショックや一般均衡シミュレーションを一切含みません。Issue #2 以降は、この Julia パッケージの \`BenchmarkParameters\` と \`CalibrationResult\` を直接再利用します。


## Issue #2 — Appendix A linear rational expectations（Julia）

追加ソース：\`src/expectations.jl\`、\`scripts/run_expectations.jl\`。既存の \`CalibrationResult\` を再利用し、Appendix A (36)–(37) の \`B\`、\`D\`、\`F = B \ D\`、二つの固有値と安定固有ベクトルを計算します。

**論文の情報時点**：\`compute_expected_equilibrium(solution, calibration, K_prior, mu_prior)\` は \(\widetilde K_{-1},\widetilde\mu_{-1}\) を *別々の predetermined inputs* として受け取り、\(E_{-1}[\widetilde K_0,\widetilde\mu_0,\widetilde X_0,\widetilde C_0,\widetilde L_0,\widetilde Y_0,\widetilde w_0,\widetilde R_0]\) を一次近似で返します。各期の実現 \(\widetilde\mu\) が自動的に安定多様体上に戻るとは仮定しません。\`compute_saddle_path_expectations\` は「初期状態が安定固有空間にある」と明示的に仮定する例示的インターフェースだけです。

**近似の意味**：\`*_next\` の level は条件付きの期待対数偏差の指数変換であり、厳密な水準の条件付き期待値ではありません。実現資本と期待資本を混同せず、後続の Issue #5 では投資ショック後の内生的な限界効用の推移を扱います。

Windows PowerShell で、クローン済み \`Nirei-2015\` ディレクトリから：

\`\`\`powershell
git status --short
git pull --ff-only origin main
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/run_expectations.jl
Get-Content .\results\expectations.toml
\`\`\`

\`Project.toml\` に \`LinearAlgebra\` 標準ライブラリを追加しました。既存の \`Manifest.toml\` を手動で作り直さず、Julia/Pkg に環境を同期させます。\`run_expectations.jl\` は定常点と期首資本を1%増加させた安定経路上の例で係数と方程式残差を \`results/expectations.toml\` に書き込みます。Julia のローカル実行と結果確認が完了するまで GitHub Issue #2 は open に保ちます。


## Issue #3 — 企業異質性・初期資本・固有生産性（Julia）

実装：\`src/firms.jl\` と \`scripts/run_firms.jl\`。Section 2.2、2.6、3.2–3.3、Appendix C の \(\lambda_i\)、\(a_{i,t}\)、\(s_{i,0}\) を既存の \`CalibrationResult\` で初期化します。外生集計ショックや企業間相関は追加しません。

- \`firm_rngs(seed, replication)\` は4つの明示的な \`Random.Xoshiro\` ストリームを生成します。ストリームごとのシードを記録しますが、異なるシードが数学的に独立であることや Julia のバージョンをまたぐビット一致は保証しません。
- \`initialize_firms(calibration, rngs)\` は、\(\lambda_i\) を切断正規分布から一度だけ生成し、\(\log a_{i,0}\) をモーメント正規化した正規分布から生成し、\(s_{i,0}\sim U[0,1)\) とします。初期資本は \(k_{i,0}=b_{i,0}\bar K\lambda_i^{s_{i,0}}\) です。
- \`aggregate_capital(state, calibration)\` は原論文 (13) に対応する有限 \(N\) の厳密なべき乗和を計算します。**有限 \(N\) の \(K_0\) を \(\bar K\) に強制一致させません。** 初期ギャップの一様性は \(\bar K\) を基準にした構成上の条件であり、実現 \(K_0\) を基準に計算したギャップが区間外に出ることがあります。その差と該当企業数をログに記録します。有限標本スケーリングは原論文 Appendix C では明示されておらず、この方法は実装上の選択です。
- \`draw_productivity!(state, rngs, calibration)\` は事前確保した \`next_productivity\` に次期の iid 生産性を生成し、現在の資本、現在の生産性、固定された \(\lambda_i\) を変更しません。\`advance_productivity!\` は資本更新後の時点移動で配列参照を交換します。Issue #4 は \`capital\` と \`next_productivity\` を直接利用できます。
- 本段階では企業の投資方策、有限企業の固定点・ES1、マクロ逐次シミュレーションは実装しません。

Windows PowerShell で、クローン済み \`Nirei-2015\` ディレクトリから次を実行します。

\`\`\`powershell
git status --short
git pull --ff-only origin main
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
julia --project=. scripts/run_firms.jl
Get-Content .\results\firms.toml
\`\`\`

今回は \`Project.toml\` に Julia 標準ライブラリ \`Random\`、\`Statistics\` を追加したため、従来の \`Manifest.toml\` との同期には \`Pkg.resolve()\` が必要です。\`Pkg.update()\` は不要です。実行成功後は \`Manifest.toml\` と新規生成した \`results/firms.toml\` のみを GitHub に追加してください。Issue #3 はローカル実行確認まで open とします。


## Issue #4 — Finite-N fixed-point enumeration and ES1

\`src/equilibrium.jl\` implements the exact finite-N investment reaction function in Nirei (2015) equations (10), (13), (14), (16), with the **Section 3.4 ES1** rule choosing the minimum absolute log deviation from expected investment among fixed points with positive aggregate net investment. It does not use ES2 or fictitious tâtonnement.

\`initialize_equilibrium_workspace(state)\` preallocates O(N) buffers. \`enumerate_fixed_points!(workspace,state,calibration)\` sorts 2N switch events and sweeps them once, incrementally updating the productivity-weighted capital power sum and net investment. Intervals are **(left,right]**, matching the paper's equality conditions. Identical switch values are grouped and processed only after the interval ending at that value is considered. No tolerance admits otherwise invalid fixed points. \`select_es1(workspace,expected_X,expected_K)\` excludes and counts nonpositive-X fixed points, for which the log criterion is undefined; exact ties prefer the fixed point closest to expected capital and then the smaller K. These last two policies are *numerical implementation choices*, not claims about the author's unpublished code.

\`solve_equilibrium!(workspace,state,calibration,expected_X,K_prior)\` enumerates all candidates, selects ES1, computes the new firm capital in a reusable buffer, recomputes the selected candidate's aggregate capital and investment exactly once, checks residuals and only then swaps the firm's capital and productivity to the next period. For the standalone one-period script, expectations are evaluated at a *specified saddle-path initial state*; the API nevertheless accepts arbitrary expected_X from the two-state Issue #2 expectation system.

In Windows PowerShell, within the existing \`Nirei-2015\` directory:

\`\`\`powershell
git status --short
git pull --ff-only origin main
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/run_equilibrium.jl
Get-Content .\results\equilibrium.toml
\`\`\`

The script generates \`results/equilibrium.toml\`, with equilibrium counts, ES1 choice, endogenous investment shock, fixed-point/accounting residuals, timing, Julia version and RNG seed. No manifest or numerical results are fabricated; Issue #4 remains open until the Windows run has been inspected.

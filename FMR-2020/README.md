# FMR-2020 — Household Heterogeneity and the Transmission of Foreign Shocks

## 研究対象と目的

de Ferra, Sergio; Mitman, Kurt; Romei, Federica (2020), **“Household Heterogeneity and the Transmission of Foreign Shocks,”** *Journal of International Economics* 124, 103303, [DOI: 10.1016/j.jinteco.2020.103303](https://doi.org/10.1016/j.jinteco.2020.103303)。

**目的は、異質家計・外貨建て債務・名目価格硬直性を持つ小国開放経済（HANKSOME）の非線形一般均衡をJuliaで独立実装し、論文の図表・定常値・ショック経路を数値的に再現することです。** 最初の完了条件はFigure 1–5の再現です。その後にFigure 6–14および厚生分析を対象にします。「図が似ている」だけを再現成功の条件にしません。

- 原論文（open access）：https://doi.org/10.1016/j.jinteco.2020.103303
- 著者公開Replication Code（MATLAB）：https://github.com/kurtmitman/HANKSOME
- 著者README：https://github.com/kurtmitman/HANKSOME/blob/main/Reamde.txt
- 著者コードのライセンス：MIT（原著者の著作権表示は出典に残す。コードの引用・再利用時は条件を遵守）。
- 当該研究の作業ディレクトリ：FMR-2020/。隣接するNirei-2015/とは独立したJulia環境にします。

**Stage 1（2026-11-02、著者EGM方式のWindows実行確認済み）**：本人のWindows 11 / Julia 1.12.6で17 Juliaファイルの構文検査、2,420/2,420件のユニットテスト、3定常均衡の求解、細グリッドEuler診断が成功しました。原著と同じ「内生EGMグリッドから直接細グリッド消費をPCHIP評価」する方式です。較正値は \(\beta=0.9832249176502228\)。各定常均衡の相対資産市場残差は \(9.1\times10^{-9}\) 以下、相対財市場残差は \(9.1\times10^{-11}\) 以下です。**細グリッドEuler最大残差は0.00347–0.00417で、事前の\(10^{-6}\)精度目標には未達**です。著者方式に切り替えた後のグリッド収束は未検証であり、Stage 1の全精度ゲートはまだ満たしません。Juliaの実行は研究者本人のPCに限定します。

## ファイル

| ファイル | 役割 |
|---|---|
| [README.md](README.md) | 研究目的、実装順、完了基準、状態 |
| [ALGORITHM.md](ALGORITHM.md) | 均衡条件、EGM、定常分布、非線形移行、突然の信用収縮、レバレッジ再評価、固定・変動相場制、IRF、監査項目 |
| [policy.json](policy.json) | 機械可読な実装・検証方針、段階別受入基準、採用候補ライブラリ、既知の照合課題 |

Stage 1のJuliaソースは `src/`、検証は `test/`、PC実行は `scripts/` にあります。Project.tomlは用意しました。Manifest.tomlは本人PCで生成されるまで未作成、results/は実行スクリプトで作成されます。著者のMATLABコードからの逐語翻訳ではなく、原論文の定常均衡条件に基づいて独立実装しています。

## 推奨する実装と検証の順序

1. **論文と著者実装の仕様監査**。式・正規化・時点・資産の通貨単位、変数の意味、パラメータを突合し、著者コードの不整合候補を記録する。原論文と公開コードが食い違う場合は、原論文準拠解と著者コード追跡解を混同しない。
2. **独立したJulia定常家計ソルバー**。Rouwenhorst、粗・細資産グリッド、借入制約を含むEGM、固定された価格の下での最適政策をテスト。
3. **定常分布と定常一般均衡**。資産遷移作用素、分布の質量保存、定常分布、純資産市場清算、β較正と異なるBでの定常均衡。
4. **非線形Perfect-Foresight Transition**。外生B経路、資本・労働の全期間経路を未知数とする市場清算、家計Backward EGM、分布Forward Simulation。Figure 1–3と比較。
5. **ショック時の外貨建て債務再評価**。40四半期後の状態を継承し、Gross Portfolioに対する予期しない交易条件変化を家計分布に反映。再評価の前後の総質量・Gross Assets・Liabilitiesを検証。
6. **相場制度別の非線形均衡**。変動相場制では国内価格安定化、固定相場制ではRotemberg価格設定と財市場清算を含める。Figure 4–5を比較。
7. **拡張実験**。22%信用供給収縮・5四半期半減期、国内財建て/外貨建て債務、低・中・高レバレッジ、Welfare/CEVを実装してFigure 6–14を比較。

各段階は前段のテストと収束記録が通過するまで「完了」にしません。問題を迂回するために予想経路・分布を固定してHANKをRANK/TANK化することは禁止します。

## 数値アルゴリズムの概要

- **定常状態**：EGMで消費・貯蓄政策関数 → 資産×個別所得のMarkov作用素 → 定常分布 → 集計資産と資本・対外純資産供給の市場清算。
- **移行経路（変動相場）**：外側で資本Kと労働Lの経路（T=200なら2T−1=399変数）を推測 → 生産・投資・q・収益率・交易条件を計算 → 家計を後ろ向きに解く → 分布を前向きに更新 → 資産・労働市場残差をゼロにする。
- **移行経路（固定相場）**：K、L、交易条件の経路（3T−1=599変数）を未知数とし、国内インフレとRotemberg NKPC、財市場残差を追加。
- **MIT Shock**：新しい外生経路は時点0で一度だけ予期せず観測され、その後の経路は完全予見。最初の流入拡大から40四半期後に別の予期せぬ信用収縮が発生する実験では、分布と資本の中間状態を必ず継承する。
- **外貨負債**：交易条件の下落で外貨建て負債の国内財建て価値が上昇する。再評価後の純資産を質量保存型補間で細かい資産グリッドへ再配分する。

詳細な順番、式、時点・配列の対応は [ALGORITHM.md](ALGORITHM.md) に記載します。

## Juliaを唯一の実装言語とする

実装は**Juliaによる独立した設計**とし、Python等でモデルを作成してJuliaへ翻訳することを禁止します。MATLAB公開コードは、アルゴリズム・仕様・結果を調査するための**参照資料**であり、行単位の機械的移植を開発方針にしません。Juliaの配列・型・多重ディスパッチ・関数バリア・疎行列・線形代数・メモリレイアウトに合った設計を採用します。

[Julia](https://julialang.org/) と [Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/) を一次資料とし、パッケージの利用前には**公式ドキュメント、GitHubリポジトリ、Juliaバージョン互換性とAPI**を確認します。

- 高負荷のカーネルは関数内に置き、型安定性を確保する。グローバルな可変状態を避け、preallocation、in-place更新、列優先走査を原則とする。
- 正規化、時点、資産グリッド、閾値を型付き構造体と明示的な関数引数で管理し、暗黙の共有状態をなくす。
- SparseArrays、LinearAlgebra、Statistics、Testなど標準ライブラリを優先する。NonlinearSolve.jl、LinearSolve.jl、Roots.jl、Interpolations.jlなどは**候補**であり、採用・バージョン・UUIDは環境検証後に確定する。
- 不等間隔資産グリッドでのEGM補間は単調性・境界挙動を数値テストで保証する。MATLAB pchipと意味が異なる補間を黙って置換しない。
- Pkgによる独立環境を用い、生成されたProject.toml/Manifest.tomlを記録する。依存関係のUUIDと互換性を推測して書き込まない。
- Juliaコード作成後は、利用可能な環境でMeta.parseall、パッケージ読み込み、最小実行、ユニット/スモークテスト、CI、結果照合を順に実施する。未実施の項目は必ず未検証と報告する。
- Windows PowerShellからjulia -eを使う場合、Pkg.add("Package")をそのまま埋め込まず、原則Pkg.add(string(:Package))の形式を使う。

## 実行環境と役割分担（必須）

**Juliaの実行と数値計算は、研究者本人のWindows PCで行います。** アシスタント側はGitHubへのJuliaネイティブ実装・修正、論文/コードの照合、静的レビュー、実行手順・PowerShellコマンドの提示、返されたログ・出力ファイルの分析を担当します。クラウド環境やアシスタント環境でJuliaコードを勝手に実行せず、GitHub ActionsによるJulia実行も本人の明示的な依頼がない限り設定・起動しません。

作業サイクルは「GitHubへ変更をコミット → 本人がPCで `git pull --ff-only` → 本人がJuliaで構文検査・依存解決・スモークテスト・本計算を実行 → コンソールログ・残差・出力を共有 → アシスタントが原因分析・修正をコミット」です。研究者PCでの実行結果が未共有なら**実行未検証・構文未検証・数値未検証**を維持します。報告されたエラーの修正後は関連コード全体を再点検し、本人PCで再検証します。

`Project.toml`、`Manifest.toml`、テストとスクリプトが揃った段階で初めて具体的な実行コマンドを確定します。許可されないPC操作やリモート実行を前提にしません。CIが存在しても実行結果を確認できない場合は成功と報告しません。

## Stage 1 ローカル検証ログと今回の修正

Windows 11 / Julia 1.12.6で14ファイルの`Meta.parseall`、スモークテスト、68/68件のユニットテストが成功したことを確認しました。一方、`run_steady_states.jl` と `grid_refinement.jl` は、粗いグリッドの消費政策を補間した細かいグリッドで、翌期資産が借入制約を約1e-5下回る問題のため失敗しました。

根本原因となった補間対象を、消費政策から**実行可能な貯蓄政策**へ変更し、細かいグリッドの消費を予算制約から再計算するコードを追加しました。違反を黙ってclampする処理ではありません。追加した回帰テストも含めて**この修正のPCでの再実行は未検証**です。詳細は[Stage 1監査](docs/STAGE1_AUDIT.md)を参照してください。

## 定常利子率ブラケットの修正（2026-11-02）

初回の借入制約違反修正はWindows 11 / Julia 1.12.6で157件のテストとグリッド感度分析まで成功しました。しかし、3定常均衡の利子率root finderが、根から遠い上端候補で家計貯蓄の最大資産グリッド逸脱（29688.016202460338）を検出しました。初回の二分法が均衡付近ではなく両端を先に評価していたためです。

利子率の探索を**基準定常利子率の近くから始めて、許容される候補間に残差の符号変化を見つけてから二分法で解く**方式へ変更しました。経済的に実行不能な上端候補だけを探索点から除外し、家計貯蓄を上限へ投影することはしません。今回の修正はまだPCでの実行・構文検証をしていないため、均衡収束は未検証です。詳しくは[Stage 1監査](docs/STAGE1_AUDIT.md)を参照してください。

## Stage 1 精度検証（3定常均衡の収束後）

すでに\(\beta\)較正と3定常均衡の資産・財市場が高精度で解けています。しかし、粗グリッドで計算した最大Euler残差は\(10^{-6}\)を上回るため、単純にGate 1を通過済みとすることはできません。

- `scripts/stage1_accuracy_audit.jl` は3つの定常均衡を再計算し、各細グリッド状態でEuler/KKT残差を独立評価し、確率質量で加重した平均・99パーセンタイル・最大誤差の位置・資産上限近傍の確率質量をCSVに保存します。最悪点が事実上観測されない状態か、正の質量を持つ状態かを判定します。
- `scripts/stage1_recalibrated_grid.jl` は97×385点と193×769点でそれぞれ割引因子と3定常均衡の利子率を**再較正**します。固定価格下の資産供給の比較（`grid_refinement.jl`）とは別実験です。385×1537点も実行する場合はPowerShellで`$env:FMR_LARGE_GRID = '1'`と指定します（大規模計算は任意）。
- 両スクリプトは今後の実行で正しく動作するか未検証です。既存の数値結果を変更するものではなく、検証用の追加計算です。

## 原著EGM方式の実行結果（2026-11-02）

Julia 1.12.6のWindows PC実行ログとGitHub保存済みの `results/stage1_euler_summary.csv`、`results/stage1_steady_states.csv` を確認しました。従来の貯蓄再補間方式から著者の内生EGMグリッド直接評価へ移行した後でも、`run_steady_states.jl` は3定常均衡に収束しています。

| 定常均衡 | B | r | 相対資産市場残差 | 細グリッド最大Euler残差 | 確率加重平均 |
|---|---:|---:|---:|---:|---:|
| Initial | 0 | 0.011279045284627722 | 1.2165000828531331e-9 | 0.0039822264227957075 | 1.4930139853175105e-5 |
| Calibrated | -4.34324989100106 | 0.010609825647710607 | 9.094922339427287e-9 | 0.004174528166204521 | 1.7034186774983578e-5 |
| Final | -11.132131056640317 | 0.009397361811040686 | 5.770314550993793e-9 | 0.0034732478454411048 | 1.970957640442883e-5 |

借入制約に近い低資産・低生産性状態には正の確率質量があり、最大Euler誤差を「分布で生じない状態だから問題ない」として無視しません。一方、経済の集計均衡は解けています。今回の結果で**定常GEソルバーは動作確認済み、Euler最大残差の精度基準および著者方式によるグリッド収束は未確認・未達**と判定します。

初回の補間方式比較（179件のテストと比較CSV）は切替前の暫定実装での実験です。現在の採用方式は著者方式のみです。

## 原著EGMの計算順序に統一（2026-11-02）

Replicationでは、原著の数値アルゴリズムを優先します。著者の`solve_R.m`および`solve_EGM_EL_open_tg.m`と同様、収束したEGMの内生資産ノット `a_today` と内生消費 `c_s` から、細グリッドの消費を直接PCHIP評価します。借入制約がbindingの領域は翌期資産を借入下限に置き、予算制約で消費を計算します。非制約領域の貯蓄も予算制約から計算します。粗グリッド貯蓄政策の再補間と方法選択オプションは廃止しました。

この変更は**アルゴリズムを著者に合わせたもの**であり、以前の均衡値をこの新方式で確認したことを意味しません。新コードの構文・実行・定常均衡結果はWindows PCで再検証します。Juliaコード自体はJuliaで独立設計しており、MATLABコードの機械的翻訳ではありません。

## Stage 1 — PC実行コマンド（Windows PowerShell）

以下は手元のGit cloneが `C:\\github\\22p8z0omwz9xu` にある場合の例です。ローカル変更があれば`git switch`前に退避・確認してください。

```powershell
Set-Location C:\github\22p8z0omwz9xu
git fetch origin
git switch research/fmr-2020-stage1-steady-state
git pull --ff-only
Set-Location FMR-2020
julia --version
julia --project=. scripts/check_syntax.jl
julia --project=. scripts/smoke_test.jl
julia --project=. test/runtests.jl
julia --project=. scripts/run_steady_states.jl
julia --project=. scripts/grid_refinement.jl
julia --project=. scripts/stage1_accuracy_audit.jl
julia --project=. scripts/stage1_recalibrated_grid.jl
```

最初に構文検査を実行し、失敗した場合はそこで止めてログを共有してください。パッケージ依存を含む実行環境を確定する場合は`julia --project=. -e 'using Pkg; Pkg.instantiate()'`を別途実行します。Stage 1のホットパスはBaseのみを利用し、外部パッケージのAPI・UUIDを推測して導入していません。実際のProject/Manifestは本人の実行環境で確定し、Manifestを生成した場合は記録してください。

`results/stage1_steady_states.csv`、`results/stage1_diagnostics.txt`が3定常均衡実行の出力です。PowerShellのログ、失敗時のスタックトレース、これらの出力を共有してください。本人PCの検証前にコードの成功・均衡解の収束を保証するものではありません。

[Stage 1 数式・較正監査](docs/STAGE1_AUDIT.md) を参照してください。


## Stage 2: Nonlinear transition implementation (source committed; PC test pending)

Stage 2 branch: research/fmr-2020-stage2-flexible-transition, based on Stage 1 branch.

The implementation contains the complete Julia source for the original nonlinear foreign-borrowing expansion with T=200 quarters and 399 K/L unknowns: transition_prices.jl, transition_household.jl, transition_distribution.jl, transition_equilibrium.jl, transition_diagnostics.jl. It reconstructs all prices and investment, performs backward EGM, moves household distributions forward for each aggregate guess, solves the full nonlinear market-clearing system, independently verifies residuals, and exports Figures 1–3 as SVG.

Official numerical reference: authors' Main.m, MIT_transition.m, Transition/solve_transition_bis.m, Transition/backsolve_egm.m, Figures.m. Refer to docs/STAGE2_AUDIT.md for exact indexing and departures in numerical optimizer. Stage 1 steady-state code was not modified.

**Critical source discrepancy:** The authors' transition labor condition omits the 1/varphi power present in the stationary model's union FOC (varphi=0.5). The Julia implementation follows the published author transition code to enable numerical comparison, but explicitly reports the implied stationary-limit mismatch. Do not silently correct the original author result or call a failed stationary-limit test a pass.

**Trade price:** With theta=1, the scalar trade-price equation used by the authors has an exact positive analytical root. The implementation solves this same economic equation directly; the general-equilibrium unknown vector remains the authors' K/L path.

**Figure 1 caveat:** The published Figure 1 includes both the credit expansion and a later unexpected credit contraction. Stage 2 only implements the expansion, and the resulting SVG is explicitly labeled partial; Figures 2 and 3 cover the expansion.

Windows PowerShell run order (execution must take place only on the owner's PC):

    Set-Location C:\github\22p8z0omwz9xu
    git fetch origin
    git switch research/fmr-2020-stage2-flexible-transition
    git pull --ff-only
    Set-Location FMR-2020
    julia --project=. scripts/check_syntax.jl
    julia --project=. test/runtests.jl
    julia --project=. scripts/validate_transition.jl
    julia --project=. scripts/run_transition_expansion.jl

Source, syntax, tests, 200-quarter convergence, and SVG output of these new Julia files are **not yet verified on the owner's PC**. A successful Stage 1 run must not be represented as Stage 2 validation. The full run writes stage2_transition_expansion.csv, stage2_equilibrium_diagnostics.csv, stage2_solver_log.txt, and figure_1.svg through figure_3.svg under results/ only after convergence; an unsuccessful root writes a failure log and raises an exception.

## 成功基準と公開物

次をすべて満たす場合のみ「Figure 1–5再現成功」と記録します。

1. 仕様表と原著/著者実装の差分・選択した解釈が明記されている。
2. 定常政策関数・分布・市場清算の独立テストが成功し、グリッド拡張に対する頑健性が示される。
3. 変動・固定両相場制でEuler/借入制約、分布質量、無裁定、資産・労働・財市場、価格設定と終端の検証を通る。
4. Figure 1–5の系列と原図・著者データ/参照解を、変数定義・基準経路・タイミングを揃えて比較する。閾値は結果を見て後から緩めず、事前に記録する。
5. コミットID、Juliaバージョン、Manifest、実行コマンド、実行環境、誤差表、図、失敗した検証も保存する。

Figure 6–14とCEVには**別の完了判定**を設けます。原著の全結果を再現したとは先に宣言しません。

## 一次資料と実装候補の公式ドキュメント

- [原論文](https://doi.org/10.1016/j.jinteco.2020.103303) — Section 2–7, Table 1, Figures 1–14.
- [著者の公開コード](https://github.com/kurtmitman/HANKSOME) — Main.m、MIT_transition.m、Steady_State_Ayagari/、Transition/、Figures.m。
- [Julia言語](https://julialang.org/), [Julia Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/), [Pkg environments](https://pkgdocs.julialang.org/v1/environments/).
- [NonlinearSolve.jl](https://docs.sciml.ai/NonlinearSolve/stable/), [LinearSolve.jl](https://docs.sciml.ai/LinearSolve/stable/), [Roots.jl](https://juliamath.github.io/Roots.jl/stable/), [Interpolations.jl](https://juliamath.github.io/Interpolations.jl/stable/).

本ディレクトリの文書は研究上の実装計画であり、原著者によるJulia実装の検証・承認を意味しません。


### 2026-11-02: Windows Stage 2 first-run failure and targeted correction

Owner Windows 11 / Julia 1.12.6 result: scripts/check_syntax.jl printed PARSE OK on 26 files, but its previous parser did not reject returned parse-error expressions. test/runtests.jl passed 2,420 Stage 1 tests, then Stage 2 test parsing failed at test/test_stage2.jl:53 because view(...,:,:,end) uses end outside indexing syntax.

The eight-period smoke completed backward EGM and forward distribution, including mass error 1.4432899320127035e-14, but its stationary candidate was NOT an asset-clearing GE (beta was imposed 0.9832). Its 0.009541840016043829 asset and 0.0012405672215367591 labor residuals are therefore not proof of a bad dynamic model.

The 200-period simulation finished three stationary GE solves but FAILED before the first root iteration: EGM saving at a far-right fine asset node was 29759.54040678395, above the upper distribution grid. Unlike the original author's backsolve_egm.m, the Julia backward household block erroneously imposed the distribution upper asset bound on all fine-grid policy nodes, including nodes with zero population mass.

The correction follows author EGM: every consumption stays strictly positive and savings satisfies the lower borrowing bound. No upper-grid bound is placed on policy evaluations at unoccupied high-wealth nodes. The forward distribution applies the upper-grid condition to all **positive-mass** states and throws an error for any true out-of-grid decision. Nothing is silently clipped or projected outside the distribution grid. A direct regression covers both occupied and unoccupied off-grid choices.

The parser now rejects error/incomplete AST nodes, the test replaces invalid end by explicit size(...,3), and the stationary-limit smoke calibrates a true stationary GE before constructing a constant-B path. The known author labor exponent inconsistency is preserved, not silently altered.

PATCH STATUS: committed but Julia syntax, full tests, eight-quarter true-stationary test, 399-variable nonlinear convergence and generated SVGs remain UNVERIFIED on the owner's PC.


### 2026-11-02: Windows 200-quarter nonlinear solve succeeded; terminal goods check isolated

Windows 11 / Julia 1.12.6 **owner-provided local execution log** for the previous commit:

- 26 Julia files passed improved Meta.parseall syntax checks.
- 2,420/2,420 stationary tests and 30/30 Stage 2 transition tests passed.
- True stationary GE constant-B eight-quarter limit: max scaled GE residual 1.2327219356400576e-7, max asset 6.362651475204854e-8, max labor 1.2327219356400576e-7, distribution mass deviation 7.327471962526033e-15.
- Full T=200, 399-variable nonlinear root **converged** in 7 iterations with 409 residual evaluations, max scaled 399-equation residual 3.642975265893966e-7.
- Independent asset market maximum level gap 9.364358355412605e-7; labor maximum 3.642975265893966e-7; capital accumulation gap 0; maximum household budget gap 1.3304388854170658e-16; maximum distribution-mass error 1.2212453270876722e-14.
- Existing independent goods-market maximum reported **0.2633831568791398**, but it pooled **all t=1:200** rather than only the interior dates for which individual household savings and consumption are part of the original author's transition closure.
- Maximum dynamic Euler/KKT violation 0.005496693743120638 (mass weighted mean 1.9816793256669015e-5). This does not meet the former strict 1e-6 maximum accuracy target.
- Terminal capital discrepancy abs(K[200]-K_final)=0.42899657754234255, terminal labor discrepancy 0.0071350035667030465. Both must be explicitly reported; they are finite-horizon approximations, not hidden by the root solver.
- CSV and SVG files were emitted before the original script discovered a failing goods diagnostic. Those files are **unaccepted provisional outputs**, not verified replications.

**Confirmed source-level endpoint mismatch:** the original author `backsolve_egm.m` aggregates household consumption and saving at dates 1:TT-1 and the original `solve_transition_bis.m` sets `c_t=[c_t';c_t(end)]`, reusing consumption of TT-1 for the last labor FOC. The T-th consumption policy and distribution are terminal evaluation objects, not an additional independent market-clearing equation. The Julia diagnostic mistakenly included T-th household consumption in the maximum independent goods-market residual, imposing a condition beyond the authors' 399 equations.

**Correction committed:** preserve the independent household-consumption goods check for **interior dates 1:199 only**; separately output `terminal_goods_gap_unconstrained`, `terminal_home_budget_gap_unconstrained`, and `terminal_wealth_gap_unconstrained`. Record `max_goods_market_period` and `max_interior_home_budget_gap` to determine if the 0.263383 gap was all at the endpoint. Do NOT loosen any tolerance, adjust the labor FOC, or change household EGM. Generate new CSV/SVG only after the independent interior GE acceptance gate passes. Regression test deliberately corrupts the T-th reported aggregate consumption and confirms the interior goods residual is unchanged, while the terminal diagnostic responds.

**Verification status:** this newest diagnostic patch is not yet syntax-checked or executed on the owner's Windows PC. It is still unknown whether the maximum interior goods residual meets 1e-4; the new code will report its period. The T=200 root convergence reported above is established for the previous commit, not the patched diagnostic.

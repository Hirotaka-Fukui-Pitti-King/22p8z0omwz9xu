# FMR-2020 — Household Heterogeneity and the Transmission of Foreign Shocks

## 研究対象と目的

de Ferra, Sergio; Mitman, Kurt; Romei, Federica (2020), **“Household Heterogeneity and the Transmission of Foreign Shocks,”** *Journal of International Economics* 124, 103303, [DOI: 10.1016/j.jinteco.2020.103303](https://doi.org/10.1016/j.jinteco.2020.103303)。

**目的は、異質家計・外貨建て債務・名目価格硬直性を持つ小国開放経済（HANKSOME）の非線形一般均衡をJuliaで独立実装し、論文の図表・定常値・ショック経路を数値的に再現することです。** 最初の完了条件はFigure 1–5の再現です。その後にFigure 6–14および厚生分析を対象にします。「図が似ている」だけを再現成功の条件にしません。

- 原論文（open access）：https://doi.org/10.1016/j.jinteco.2020.103303
- 著者公開Replication Code（MATLAB）：https://github.com/kurtmitman/HANKSOME
- 著者README：https://github.com/kurtmitman/HANKSOME/blob/main/Reamde.txt
- 著者コードのライセンス：MIT（原著者の著作権表示は出典に残す。コードの引用・再利用時は条件を遵守）。
- 当該研究の作業ディレクトリ：FMR-2020/。隣接するNirei-2015/とは独立したJulia環境にします。

**現在の状態（2026-11-01）**：アルゴリズム・再現方針・検証計画を文書化済み。Juliaソースコード、Project.toml、Manifest.toml、数値出力・図、CIは未作成です。Julia構文解析・実行・IRF一致は未検証です。著者MATLABコードについても、この作業で実行成功を確認したわけではありません。

## ファイル

| ファイル | 役割 |
|---|---|
| [README.md](README.md) | 研究目的、実装順、完了基準、状態 |
| [ALGORITHM.md](ALGORITHM.md) | 均衡条件、EGM、定常分布、非線形移行、突然の信用収縮、レバレッジ再評価、固定・変動相場制、IRF、監査項目 |
| [policy.json](policy.json) | 機械可読な実装・検証方針、段階別受入基準、採用候補ライブラリ、既知の照合課題 |

将来の実装時には、このディレクトリ内に独立したProject.toml、Manifest.toml、src/、test/、scripts/、results/を置きます。未作成のファイルや実行コマンドは存在するものとして扱いません。

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

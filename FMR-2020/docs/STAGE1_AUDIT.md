# Stage 1 — 原著との定常均衡照合と検証保留事項

## 対象と出典

- de Ferra, Mitman, and Romei (2020), JIE 124, 103303, Section 2–5, Table 1: https://doi.org/10.1016/j.jinteco.2020.103303
- 公開著者コードの Main.m: https://github.com/kurtmitman/HANKSOME/blob/main/Main.m
- 公開著者コードの calibrate_beta_open_tg.m: https://github.com/kurtmitman/HANKSOME/blob/main/Steady_State_Ayagari/calibrate_beta_open_tg.m
- calibrate_R_open_tg.m: https://github.com/kurtmitman/HANKSOME/blob/main/Steady_State_Ayagari/calibrate_R_open_tg.m
- solve_R.m / solve_EGM_EL_open_tg.m: https://github.com/kurtmitman/HANKSOME/tree/main/Steady_State_Ayagari
- ComputeDistHist_open.m: https://github.com/kurtmitman/HANKSOME/blob/main/Steady_State_Ayagari/ComputeDistHist_open.m
- Julia official language and performance: https://julialang.org/ ; https://docs.julialang.org/en/v1/manual/performance-tips/ ; https://pkgdocs.julialang.org/v1/environments/

この文書に記載する原著設定はソースを静的に読んだものです。実行・再現成功を意味しません。

## 経済学的な定常閉鎖

定常状態を家計の資産選択・生産・外部資産供給で閉じます。家庭の純資産供給は A^s=Σμg で、供給量は K+B。外部純資産 B<0 は国外からの純借入です。投資I=δK、q=1、国内財価格硬直性の調整費用=0（定常インフレ1）。企業の定常価格では mc=(ε_f−1)/ε_f、資本レンタルr+δ、賃金w=mc(1−α)Y/L、利益=Y/ε_f。組合の賃金FOCを用いてLを決定する。初期の財需要を決めるため、外国の輸出需要D x^(−θ*)と集計予算の国内財消費C_Hを一致させる。βまたはrの根探索の**それぞれの試行**で家計EGMと分布を解き直す。

CES需要からc_F/c_H=[(1−χ)/χ]x^θ、消費支出= [1+((1−χ)/χ)x^(θ−1)] c_H。今回のγ=θ=1では消費のEulerはβ(1+r)E c_H^(−γ)。これを一般のγ・θに拡張する場合は限界効用と組合の労働条件を再導出し、デフォルト以外の設定を無断で使わない。

## 元ソースで確認した較正情報

| 対象 | 元コードの設定・計算 | Stage 1 |
|---|---|---|
| 生産性 | nx=7, rho=.97, innovation std=.84 sqrt(1−rho²), level mean normalised to 1 | 直接対応。ここで0.84は**log生産性の定常標準偏差**であり、innovation stdとは異なる |
| 資本・GDP比 | capital per capita=(-NFA + net wealth); NFA=-6.5×.69、wealth=26.2、GDP=12.65/4 | K/Y ≈ 9.703（四半期GDP比）。対象の単位を保持 |
| 外貨純資産/GDP | (-6.5×.69)/(12.65/4) | 約−1.418（コードでの較正）。原論文Table 1の他の実験設定と混同しない |
| 定常資本 | L=1、K=(K/Y)^(1/(1−α)) | 同じ正規化 |
| 割引因子 | beta_guess=0.9832、fsolveで残差A−K−Bをゼロにする | 固定せず、根探索の参照値とする |
| 財配分 | ω=.4, χ=.6, θ=γ=1 | c_H消費と消費支出の単位に分ける |
| 需要シフターD | baseline Y−C_H−δK | 初期のx=1を基準に構築 |
| 3定常 | B_initial=0, B_calibrated=観測較正, B_final=拡大B経路の最終値 | それぞれrと家計分布を再計算 |
| グリッド | 原著のフルグリッド193 / 769、z∈[−.2,1]、z³ | フルグリッドを保存。Stage 1の定常家計計算は借入制約以上の部分を使用 |

## 独立実装と著者コードの違い（意図的）

- 著者の定常分布はsparse matrixの支配的固有ベクトルで計算する。Julia版は**mass-preserving, matrix-free iteration**を用い、政策が一定のとき資産ノットの上下インデックスと補間ウェイトを前計算する。数値的にμ=Q'μを満たすか独立診断で評価する。
- MATLABのpchipを単純なcubic splineに置き換えることはせず、JuliaネイティブのFritsch–Butland / Hermite形状保持補間を実装する。MATLABの出力との逐点一致は**未検証**。
- MATLABのグリッドには、外貨再評価のために借入制約を下回る点がある。JuliaのStage 1では全体のグリッドを残しつつ、定常の**翌期資産・定常分布は制約以上の活性集合でのみ計算**する。ショック時の負の期首資産はStage 2の家計問題として別に実装する。
- Julia版では根のブラケットを明示し、異常値のclippingや収束失敗の黙殺を避ける。βとrをroot finderで解き、解を再評価して検証する。
- 第1段階では外部Juliaパッケージに依存しない。Project.tomlのUUIDは本研究用に新規生成したプロジェクト固有IDで、外部パッケージUUIDを推測したものではない。

## 今後、Windows PCで確認する必要がある点

1. Meta.parseallが全Juliaファイルで成功すること。Julia moduleとTestの読み込みが成功すること。
2. 実際のEGM / invariant distributionが既定tolで収束し、βとrの初期ブラケット内に根が存在すること。
3. 1e-12のEGM update criterionは原著の設定値だが、Euler residualと同義ではない。補間誤差、KKT不等式を別に報告する。
4. bisectionのβ初期範囲[.980,.985]、r初期範囲[-.005,.055]（実際の上限は1/β−1−0.0001まで制限）は探索開始範囲であり、実行時に根が括れていないと報告された場合に、モデルと符号を監査してから拡大する。
5. 小グリッド・大グリッド間のasset supply, interest ratesの収束、stationary tail mass、最上位資産ノットへの質量滞留がないこと。
6. 著者Main.mと原論文Table 1で見られるB/GDPの相違が、論文の別実験・表記単位・コードの実験設定のいずれかを原論文と対応づけて説明できること。現時点では未確定。
7. β/r較正後、所得・配当の合計、財市場・資産市場・労働組合FOC・企業要素価格の整合が数値的に成立すること。

**報告ステータス**：JuliaコードはGitHubに追加済み、構文・読み込み・実行・数値結果は**本人PCでの実行ログ待ち（未検証）**。論文全体のReplicationではない。

## 2026-11-02: User-PC regression, fine-grid policy fix

Windows Julia 1.12.6 result supplied by the user:
- **Syntax**: Meta.parseall for 14 files PASS.
- **Smoke**: 1436 EGM updates (reported diff 9.997052075050306e-10), Euler 0.00110884601889777, KKT 0, budget 2.2021586079190927e-17; distribution mass 0.9999999999999907, 1565 iterations, 9.721244087490533e-10. Smoke PASS.
- **Tests**: 68/68 PASS.
- **Full stationary calibration**: FAIL at solve_household fine-grid policy with computed savings −9.862788178216508e-6 below the borrowing limit.
- **Grid refinement**: FAIL at same policy interpolation with savings −1.0261942101019716e-5 below the borrowing limit.
- **Diagnosis from source**: previously PCHIP-interpolated coarse consumption and then inferred savings from the budget. Even if coarse savings respects the bound, consumption interpolation can imply savings **below** that bound at interior fine-grid knots.
- **Correction committed**: interpolate the *coarse savings policy* with shape-preserving Hermite instead; derive fine-grid consumption exactly from budget, and fail on any remaining infeasible saving/consumption (no tolerance relaxation, no clipping). Report max fine budget error and minimum fine saving. Add direct regression on (61,241) and (193,769) grids at beta=0.980.
- **Status**: patch is NOT YET verified on the owner's PC. EGM solution, distribution, beta calibration and all three steady-state equilibrium outcomes remain unverified after the patch. The new savings interpolation differs from author code's consumption interpolation, so its numerical effect on asset demand and the original paper's calibration must be checked by grid refinement and comparison.


## 2026-11-02: second user-PC run, rate-bracketing failure

The borrowing-limit interpolation repair was verified on the owner's Windows
PC with Julia 1.12.6:
- Meta.parseall: 14/14 files successful.
- Smoke EGM/distribution: successful, 1436 and 1565 iterations respectively.
- Unit tests: 157/157 passed.
- Grid refinement: completed without violating borrowing constraints.
  Asset demands: 61/241=25.64759695956763; 97/385=25.401319414404544;
  193/769=25.246929207836487. These are fixed-beta/price evaluations,
  **not equilibria**, and convergence to a grid limit has not been established.
- Full three-steady-state run: beta calibration progressed, but the first
  interest-rate root solve failed with "fine-grid policy exceeds distribution
  grid upper bound"; reported chosen savings 29688.016202460338.
- Source-level cause: the current bisection evaluated the fixed high endpoint
  near 1/beta−1, where high-wealth households can prefer savings outside
  the finite asset grid. This is not evidence that an interior equilibrium root
  does not exist. Artificially truncating their policy would distort asset
  market residuals.

**Committed repair**: root bracketing anchored at ref.r, scanning nearby
candidate rates for a genuine sign change and only skipping the specifically
identified non-feasible *upper-grid* exploratory points. Other economic and
numerical failures propagate. Root computation uses ordinary bracketed
bisection only after a valid pair has been found. If none is feasible, raise
an error reporting the explored interval and number of excluded candidates.
No savings clipping and no looser tolerances. Add mock-root regression and
stage-by-stage progress diagnostics.

**Status**: latest patch syntax and execution UNVERIFIED pending Windows-PC
rerun; the accepted 157 tests refer to the previous commit. No claim about
the existence or values of the three stationary equilibria is made yet.


## 2026-11-02: user-PC full stationary equilibrium successfully computed

Source of results: pasted Windows PowerShell console log, Julia 1.12.6, based on
branch head \`4ba084bb8a02e53017f704b87c2a8ae9be9714f5\`. The numerical outputs
below are user-provided verified run logs, **not** estimates from an assistant runtime.

- Meta.parseall: 14/14 files parsed.
- Unit tests: 161/161 passed; prior smoke test and grid-refinement run also passed.
- Calibrated beta: 0.9832249340415002; beta scaled asset residual
  2.4759552789121424e-9; beta root iterations 24.
- Reference NFA/Y: -1.418181818181818; K/Y: 9.702766798418972.

| Steady state | B | r | Unscaled asset residual | Scaled asset residual | Scaled goods residual | Maximum coarse Euler residual |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Initial | 0 | 0.011279042995809365 | -5.774051814455561e-8 | 2.039225110212321e-9 | 5.699578568478156e-11 | 0.00026237630682234236 |
| Calibrated | -4.34324989100106 | 0.010609825647710607 | 6.848970457440373e-8 | 2.699427274785832e-9 | 8.257698680921228e-11 | 0.00026706681849458303 |
| Final | -11.132131056640317 | 0.009397365702031897 | -9.265521327961324e-8 | 4.3842010122971905e-9 | 3.4551611688112584e-11 | 0.00027831416802568147 |

All reported stationary-distribution masses were within approximately
1.4e-14 of one. The above establishes **stationary market clearing on the
author-code-style grid**, but not the stronger Euler/grid convergence
requirements in policy.json.

### Decision

G2 market-clearing numerical solve: **demonstrated** for three stationary
states. Full Stage 1 acceptance: **pending** because policy.json sets maximum
unconstrained Euler residual threshold 1e-6, whereas all three observed maxima
are near 2.7e-4. No threshold relaxation is authorized. The fixed-beta/price
asset supply sensitivity across (61,241), (97,385), (193,769) grids is
25.64759695956763, 25.401319414404544, 25.246929207836487. Recalibrated
GE grid convergence has not been established.

### Accuracy investigations committed for owner's PC

- scripts/stage1_accuracy_audit.jl: Fine-grid Euler gap and complementarity
  residual at each state; record stationary mass at worst error, mass-weighted
  mean and mass-weighted 99th percentile, and upper asset-grid tail mass.
- scripts/stage1_recalibrated_grid.jl: Re-run beta and r calibration for
  alternative coarse/fine grids, separately reporting three stationary
  equilibria for each grid.

No equations or tolerances were changed merely to make the acceptance gate
pass. Julia source in these newly added audit files remains **syntax and
execution unverified until owner's Windows-PC logs are provided**.

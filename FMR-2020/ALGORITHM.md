# FMR-2020 — 数値アルゴリズム仕様書（Julia独立実装用）

**対象**：de Ferra, Mitman & Romei (2020), *Household Heterogeneity and the Transmission of Foreign Shocks*, JIE 124, 103303。原論文 https://doi.org/10.1016/j.jinteco.2020.103303 。著者公開コード https://github.com/kurtmitman/HANKSOME 。

**文書の位置づけ**：原論文の均衡式、著者MATLABコードで確認できる実際の計算順序、今後採用するJuliaの独立実装・検証設計を区別して記述する。ここに書いた数値収束・Julia実行はまだ実施していない。変数の時点の切り方は新しいJuliaコードで一義的に定義し、MATLABの1始まり添字と「ショック後時点」を同一視しない。

## 0. 全体の依存関係と原著のコード

1. パラメータ・Rouwenhorst近似・資産グリッドを生成する。出典：Main.m、Exogenous_files/rouwenhorst.m。
2. 与えられた定常価格に対しEGMで家計消費・翌期資産政策を解く。出典：solve_EGM_EL_open_tg.m、solve_R.m。
3. 政策と所得遷移から定常分布・集計資産を計算する。出典：ComputeDistHist_open.m。
4. β較正と外国信用残高B別の定常一般均衡を解く。出典：calibrate_beta_open_tg.m、calibrate_R_open_tg.m。
5. 初期B=0から信用拡大の完全予見移行経路を解く。出典：MIT_transition.m、solve_transition_bis.m、backsolve_egm.m。
6. 約40四半期後のK、μ、前期のqと交易条件を継承し、信用収縮の新しい将来経路に対する完全予見移行を解く。
7. 外貨建てGross Liabilitiesの再評価を新しい期首家計分布へ反映する。レバレッジの配分も変える。
8. 変動相場制と固定相場制を比較する。固定相場制はsolve_transition_fixed.mでRotemberg NKPCと財市場を追加。
9. Figures.mでFigure 1–5。Figure 6–14/CEVは別途実装・確認が必要。

**核心**：非線形transitionの外側では集計数量経路を探索し、内側では（終端→現在）家計政策関数のBackward EGMと（現在→終端）家計分布のForward Propagationを実行して市場清算残差を返す。SSJの一次線形IRFや集計RANKモデルへの置換ではない。

## 1. モデル：定義と単位

家計状態は現在純資産aと個別労働生産性s。個別生産性は7状態Markov chain、遷移行列Π、定常平均E[s]=1。家計分布μ_t(a,s)は内生的な状態変数。生産Y_t=Z_t K_t^α L_t^(1−α)、投資I_t=K_(t+1)−(1−δ)K_t。

二財CESの価格指数の下で国内財消費をc_H、外国財消費をc_F、交易条件をx_t=P_H,t/P_F,tとする。χは国内財の支出ウェイト、ω=1−χ（コードのomega）。個別最適財配分は c_F/c_H=[ω/(1−ω)] x^θ。したがって国内財単位の消費支出は m(x)c_H、m(x)=1+[ω/(1−ω)]x^(θ−1) である。論文較正θ=1、γ=1（CES指数のθ→1の極限はCobb–Douglasとして評価し、0/0を数値評価しない）。

家計予算：m(x_t)c_H,t(a,s)+a'=(1+r^a_t)a+y_t(s)+d_t、a'≥a_min。資産の意味と収益率の時点を固定し、ショック時の再評価後の期首aをこの式に渡す。家計の労働時間は各家計が資産に応じて個別に自由選択するのではなく、組合による共通賃金設定から集計Lを求める（論文Section 2.1）。共通配当の分配ルールも原著・コードと照合する。

企業：国内財のマークアップ ε/(ε−1)、資本のレンタル r^K_t、実質賃金w_t。投資費用 Ψ_t=(ψ/2)K_t (I_t/K_t−δ)^2、q_t=1+ψ(I_t/K_t−δ)。資本の無裁定：1+r^a_(t+1)=[r^K_(t+1)+(1−δ)q_(t+1)]/q_t。

外国部門：外国の国内財需要C_H,t^*=D x_t^(−θ*)。外国純資産B_(t+1)は、外生的に与えるホーム財建て外国信用供給経路に対応する。資産市場：A^s_(t+1)=q_t K_(t+1)+B_(t+1)。国内財資源制約はY_t=C_H,t+C_H,t^*+I_t+Ψ_t+RotembergCosts_t。変動相場（国内財価格安定）ではRotembergCosts=0、固定相場では価格調整費用とNKPCが必要。

「外国債務の額」と「純外国資産B」の符号は厳密に区別する（外国からの借入増加はBが負方向）。粗・細グリッドのaは純資産であって、Gross AssetsとGross Liabilitiesは後段の追加配分である。

## 2. グリッド・個別所得過程

著者Main.mで確認した初期技術設定：
- n_s=7、個別所得AR(1)持続性ρ_s=0.97、Rouwenhorstを使用。正規化して定常平均s=1。
- 資産の粗グリッド193点、細グリッド769点。著者はzを−0.2から1まで等間隔に取り、a(z)=a_min+(a_max−a_min)z^3を使用。z=0（借入制約a_min）がグリッドに含まれる。
- a_min≈−10^(−7) K/Y × Y、a_maxはMain.mの定義から生成。非負に投影したグリッドに書き換えない。初期外貨負債再評価後にa<a_minとなり得るので、その状態を表すグリッドが必要。
- 定常家計EGMのtol=10^(−12)、最大5,000反復。移行期間T=200。これらは著者コードの**設定値**であり、独立実装の達成済み残差ではない。

Juliaでは資産グリッドを昇順に保持し、単調性・重複ノット・外挿・借入制約の位置を明示的に検証する。細グリッドの2点への質量配分に使う上下インデックスとウェイトを前計算できる箇所ではキャッシュする。

## 3. Stationary household EGM

ある定常価格(r,w,x,d,L)およびβの下で、粗グリッド上に消費政策関数c_H^(n)(a,s)を推測する。

(1) 各次期資産a'および現在所得状態sについて
    M(a',s)=β(1+r) Σ_(s') Π_(s,s') [c_H^(n)(a',s')]^(−γ)
を計算する。論文/著者のθ=1/γの制約下で、国内財消費のEulerのこの表現を用いる。

(2) 内生グリッドの消費：c_endo(a',s)=M(a',s)^(−1/γ)。

(3) 当期資産を予算式から逆算：
    a_endo(a',s)=[m(x)c_endo(a',s)+a'−y(s)−d]/(1+r)。

(4) 制約領域ではa'=a_minとし、
    c_con(a,s)=[(1+r)a+y(s)+d−a_min]/m(x)。
    非制約領域は(a_endo,c_endo)を元のaグリッドに単調補間する（著者はpchip）。補間・外挿・負の消費を検知し、黙って無理なclipをしない。消費を下限値に丸める著者コードの処理が数値結果に与える影響も検証する。

(5) ||c_new−c_old||∞が設定許容値以下になるまで、必要に応じてdamping付きで反復する。政策関数とa'=g(a,s)=(1+r)a+y(s)+d−m(x)c_H(a,s)を保存する。

(6) 得た粗グリッド政策を細グリッドに評価し、a'≥a_min、Euler残差とbinding領域のEuler不等式を検証する。最適化の誤差と補間の誤差を分ける。

## 4. Stationary distributionと集計

細グリッド上のg(a_i,s_j)が隣接a_ℓとa_(ℓ+1)の間なら、重みλ=[a_(ℓ+1)−g]/[a_(ℓ+1)−a_ℓ]を使い、翌期所得s_kの確率Π_(j,k)と掛ける：

Q[(i,j),(ℓ,k)] += λ Π_(j,k)、
Q[(i,j),(ℓ+1,k)] += (1−λ) Π_(j,k)。

Qは行確率遷移行列（確率行和1）の疎行列。定常分布はQ' μ̄=μ̄、sum(μ̄)=1。著者ComputeDistHist_open.mではsparseとeigs(Q',1)を利用する。JuliaではまずSparseArraysによる遷移作用素とpower iterationを独立テストし、その後に固有値法が必要かを判断する。

集計資産A^s=Σ_(a,s) μ̄(a,s)g(a,s)、国内消費C_H=Σ μ̄ c_H。政策関数微分からMPCを評価するときは微分の変数（資産aか一時的所得か）と単位を明記する。数値的質量保存を全反復で確認する。

## 5. 定常一般均衡とβ較正

Main.mはまず較正済みのK/Y、L等からK,Y、要素価格、rおよび利益を構築する。マークアップを考慮した定常要素価格：
    r^K=[(ε−1)/ε] α Y/K、r=r^K−δ、
    w=[(ε−1)/ε](1−α)Y/L。

関数 calibrate_beta_open_tg はβを仮定 → EGМ → 定常分布 → 資産A^s を計算し、残差 Fβ(β)=A^s−(K+B) を評価する。Main.mのbeta_guess=0.9832は初期推測値であり、最終解は根探索で得る。

次に関数calibrate_R_open_tgではβを固定し、外国純資産Bを与えて、rを未知数とする。r→K/L,w→集計労働・輸出・交易条件→EGM→定常分布→F_r(r;B)=A^s(r;B)−K(r;B)−B=0を解く。

保存する定常均衡は「B=0の初期定常」「信用拡大が続いた場合の最終定常」「外貨収縮後に向かう較正定常」の3つ。これらのB、政策、分布、q、xを混同しない。原論文Table 1とMain.mで較正方法/パラメータ値の差を別表で監査する。

## 6. 最初の信用拡大ショックの経路

著者Main.m：params.cur_acc=0.074、rho_shock=0.98、T=200。配列j=1,...,Tで、
    B_boom[j]=−0.074 Y_initial Σ_(k=0)^(j−1) (0.98)^k。
これを「外国からの資金流入に伴い純外国資産が減少する」経路として扱う。定常経常収支やGDPの期間単位（四半期・年率表示）を混同しない。

経済主体は初期時点に、これまで予期していなかったB_boomの全経路を知る。以降は完全予見。これにより「元の定常→外国借入拡大後定常」への非線形 transition を求める（論文Section 5–6、Figure 1–3）。

## 7. 変動相場制：外側399変数の移行均衡

T=200、初期K_0所与。未知ベクトルX=(K_1,...,K_(T−1), L_0,...,L_(T−1))、長さ2T−1=399。変動相場政策は国内財価格安定化を実現し、実質限界費用mc=(ε−1)/ε。

一回の関数評価 F_float!(R,X,inputs,workspace) の概念的順序：

1. XからK_t,L_tを展開、終端K_Tなどは定常状態を利用する（最終期の資産清算と終端政策の扱いを原著コードと一致させる）。
2. Y_t=Z_t K_t^α L_t^(1−α)。
3. I_t=K_(t+1)−(1−δ)K_t、Ψ_t=(ψ/2)K_t(I_t/K_t−δ)^2、q_t=1+ψ(I_t/K_t−δ)。終端q=1。
4. r^K_t=mc αY_t/K_t、w_t=mc(1−α)Y_t/L_t。
5. 1+r_t^a=[r^K_t+(1−δ)q_t]/q_(t−1)。前期qを外生初期条件から渡し、t=0の既保有資産の収益と予想t+1の収益を混同しない。
6. 各期の集計予算上の国内財消費C_H^budget(x)と、国内財資源制約C_H^resource(x)=Y−I−Ψ−D x^(−θ*)を一致させ、正のx_tを求める。著者は各期に一変数fsolveを呼ぶが、Juliaでは括り付きの正値根探索を優先して評価。
7. 家計終端政策c_T=c_finalを読み、Backward EGMを実行。
8. 初期分布μ_0を読み、Forward Distributionでμ_(t+1)を更新、A^s_(t+1)、C^s_H,tを集計。
9. R_K,t=(A^s_(t+1)−B_(t+1))/q_t−K_(t+1)、R_L,t=L^s_t−L_t を格納。
10. 全期間のRを返し、外側の非線形root solverでF_float(X)=0を求める。残差評価は常に内部の家計最適化と分布を含む。

著者のsolve_transition_bis.mでは価格・集計消費の内側で交易条件を決め、backsolve_egm.mがBackward/Forwardを行う。transition_start.matは「初期推測値」にすぎず、モデルの真値として固定しない。

## 8. 非定常Backward EGM：価格経路所与

終端政策c_T(a,s)を最終定常政策で与える。t=T−1,...,0で
    M_t(a',s)=β(1+r^a_(t+1)) Σ_(s')Π_(s,s') [c_(t+1)(a',s')]^(−γ)、
    c_endo,t=M_t^(−1/γ)、
    a_endo,t=[m(x_t)c_endo,t+a'−y_t(s)−d_t]/(1+r^a_t)。
借入制約領域ではa'=a_min。単調補間と境界判定からc_t(a,s)、g_t(a,s)を得て、粗・細グリッドに保存する。

注意：終端条件は「最終定常状態に到達した」とみなす有限Tの近似。T=200→300に変更して初期30四半期のIRFが安定しているかを検証する。

## 9. 非定常Forward Distribution：実現した分布の更新

前段が出力した政策関数を使い、μ_0（初回は初期定常分布）から前向きに
    μ_(t+1)=T_(g_t,Π)(μ_t)
を計算する。質量保存・非負性・期末純資産の市場価値を検証する。著者コードbacksolve_egm.mは細グリッドの翌期資産を隣接2ノットへ分割し、Πと掛けて更新する。Juliaではベクトル化に固執せず、列優先・事前確保・局所的な加算ループで実装する。並列化は同一配列への加算の競合/非決定性を避ける。

## 10. 外側の一般均衡ソルバー

変動相場の残差の数は2T−1。Newton・準Newtonなどの具体的な選択は性能と頑健性テストで決める。初期経路は初期/最終定常状態を補間して用い、実験を進めると前の均衡経路をwarm startにする。非正のK,L,x、不適切なq、借入制約違反では明示的に失敗を返し、黙って値をclipして均衡を偽装しない。

根探索の終了フラグと独立に、max normの（単位調整済み）市場残差、経済量の有限性・正値、家計Euler残差、資源制約、終端誤差を計算する。許容値はpolicy.jsonに事前記録した「提案値」と「原著の値」を区別する。

## 11. 第二のMIT Shock：40四半期後の予期せぬ外国信用収縮

著者MIT_transition.mのshock_start_t=41（MATLAB 1-based）。初期信用拡大の移行経路上にある家計分布とK、直前のqと交易条件を継承する。実装における「shock直前」と「shock時点」の配列対応を単体テストする。

著者コードはj=1,...,Tについて
    B_stop[j]=B_cal+(B_boom[40]−B_cal)(1−(j−1)/(T−1))^4
を作る（コードはshock_b_agg(shock_start_t-1)を参照）。これは信用残高のジャンプではなく、将来の資金流入経路の予期せぬ変更として扱われる。初期分布はdist_tran(:,:,41)、Kはk_agg(41)、qと交易条件は40番目の要素を渡す。定義を踏まえoff-by-oneを必ず検証する。

新たな終端定常均衡、B_stopと inherited state に対する非線形移行を再計算する。ここを初期定常分布に戻すと論文実験を再現できない。

## 12. Gross Portfolioと外貨建て負債再評価

家計の純資産aはGross Capital k_gross(a)と（負の）外国通貨建て債務b_F(a)の和：a=k_gross(a)+b_F(a)。公開コードは純資産政策だけを動的最適化し、Gross Portfolio配分を外生的な関数で割り当てる。これを多資産家計の最適ポートフォリオ解と誤解しない。

概念的な配分：k_gross(a)=max(λ a+κ,0)、b_F(a)=a−k_gross(a)。元コードのk_const_hatが0（初期）から1/16（実験）になるケースは、低純資産家計にも正のGross Assets・Gross Liabilitiesを割り当てる。

交易条件x_old→x_newが予期せず変わると、外貨負債は
    b_F,new(a)=b_F,old(a) × x_old/x_new
に再評価され、
    a_new(a)=k_gross(a)+b_F,new(a)。
外国通貨建て債務b_F<0、交易条件下落なら国内財建て負債絶対値が増える。

「再評価後の期首分布」は、各元グリッド点a_iに付随する確率質量μ_0(a_i,s)を新しい位置a_new(a_i)へ、隣接2グリッド点の線形ウェイトで保守的に移す。生産性sはこの瞬間にはMarkov遷移しない。負のa_newが表現可能か、グリッド外に質量が出ないかを確認する。元コードの条件チェックが十分か独立検証する。

**一般均衡フィードバック**：x_newは外側の推測K,Lから内生的に決まるため、各残差評価ごとに a_new→μ_new→C,A→市場残差 が変わる。再評価を一度だけ固定してはいけない。

## 13. 固定相場制：599変数・価格硬直性

固定名目為替レートe=constant、よって国内財インフレπ_H,t=x_t/x_(t−1)。未知X_fixed=(K_1,...,K_(T−1), L_0,...,L_(T−1), x_0,...,x_(T−1))、長さ3T−1=599。

企業のRotemberg条件（論文Section 2.2）
    0=(1−ε)+ε mc_t−ζ(π_t−1)π_t+[ζ/(1+r^a_(t+1))](Y_(t+1)/Y_t)(π_(t+1)−1)π_(t+1)。
したがって終端条件を置きmc_tを後ろ向きに計算し、w_t、r^K_t、r^a_t、企業利益を求める。価格調整資源コストは(ζ/2)(π_t−1)^2Y_tで、必ず国内財資源制約から差し引く。

この価格経路を家計のBackward EGMへ渡し、再評価後の初期分布からForward Simulationする。残差は
    R_K,t=(A^s_(t+1)−B_(t+1))/q_t−K_(t+1)、
    R_L,t=L^s_t−L_t、
    R_C,t=Y_t−C^s_H,t−D x_t^(−θ*)−I_t−Ψ_t−(ζ/2)(π_t−1)^2Y_t。
一度の外側評価で3系列の残差を計算し、F_fixed(X)=0を解く。固定相場の内生価格硬直性を柔軟価格のmcに置き換えることは禁止。

## 14. Fig 1–5の作図とIRF定義

Figure 1はBの拡大と予期せぬ反転（著者は年次GDPの%表示のため四半期GDP比を0.25倍している）。Figure 2–3はB=0の初期定常水準を基準とする初回移行。Figure 4–5は40四半期後の突然の反転で、著者Figures.mは時点近傍の**単一の移行経路水準**をベースに図示する項目が多い。時点ごとの「予期せぬショックなし経路」からの差を新規計算する場合、著者のプロットと混同しない。

ログIRF：100(log X_shock,t−log X_base)、投資の通常パーセント差：100(I_shock,t/I_base−1)、利子率は水準・percentage points。図の横軸はQuarterであり、第一の対外信用拡大図と第二のSudden Stop図の原点を区別する。総消費・国内財消費・外国財消費のラベルを厳密に区別する。

Figure 1–5について「元図表示の近似」「原著保存済みMATファイルとの数値一致」「経済モデルを独立に解いて同じ数値になる」を別の証拠水準として記録する。最重要なのは後者。

## 15. Figure 6–14と厚生（後続段階）

Section 7は定常経済から22%の外国信用収縮、半減期5四半期で戻るmean-reverting MIT shockを用い、国内財建て債務/外貨建て債務、変動/固定相場、低/中/高レバレッジを比較する。著者READMEで保証される公開コードの直接再現範囲はFigure 1–5であり、Figure 6–14やCEVの再現には追加実装が必要。

ショック経路ρ=2^(−1/5)は「半減期5四半期」から計算する候補であり、原著のショック定義との一致を明示的に確認してから採用する。すべての国別・世帯別厚生比較は家計価値関数のショック経路/反実仮想経路を別々に評価し、CRRA=1ならlog効用のCEV式の符号と基準消費を検証する。論文Section 6の0.81%/0.34% CEVは一次資料の報告値であり、未実行の推定結果として扱わない。

## 16. 原著コードと論文の照合課題（修正前の監査）

以下は公開MATLABソースを読んで発見した**不整合候補**。実行による誤差特定・修正済みという意味ではない。

1. 労働供給条件：Main.mのvarphi=0.5を用いる定常式は一般にL=[((1−ω)(ε_w−1)/ε_w)w/(χ_dis C_H)]^(1/varphi)だが、Transition/solve_transition_bis.mとsolve_transition_fixed.mのlab_primeには1/varphi累乗が見られない。労働不効用の指数と同じ式か要検証。
2. 固定相場制の資源制約：solve_transition_fixed.mのnet_cに価格調整費用の二重負号が見られる。投資調整費用で配列中の期別K_tではなくループ残存のtを参照する箇所が疑われる。数式から再導出し、固定・変動の資源制約を別々にテストする。
3. 家計労働所得：backsolve_egm.mの所得計算が一定マークアップの係数に依存し、固定相場でNKPCから求めた賃金・限界費用をそのまま渡していないように見える。賃金総額の家計収入との一致をチェック。
4. 時点・配列：shock_start_t=41とdist_tran(:,:,41)、k_agg(41)、q_agg(40)、tot(40)の継承がモデルの時点定義に合っているかを独立テストする。
5. 定常分布の評価：Eigenvectorの負値clipやinterpolation outside gridが結果を変えるか確認。失敗時にMATLABのkeyboardへ入る箇所はJulia版では例外と診断記録にする。
6. 較正の差：Table 1のB（GDP比−2）とMain.mのtarget.nfa_gdpの算出値が一見異なる。これは実験、四半期換算、変更された較正のどれに由来するかを確定するまで勝手に修正しない。
7. 専用ファイルを使った高速化：Main.mのbeta_guess=0.9832、transition_start.matのkl_pathは便利な初期推測値。保存済み参照解を「収束済みだから正しい」と断定しない。
8. 技術的な演算順・終端：200四半期目の投資・q・配当の終端近似を論文式と照合し、境界時点の意味を固定する。

正しい原著の数式と、論文図を直接生み出した公開コードの「再現」を区別するため、修正を行うときはbaseline-author / baseline-equationsの比較表を必ず保存する。

## 17. Juliaのネイティブ設計（翻訳禁止）

- JuliaコードはJuliaで**新規に設計**する。Python、MATLAB等のモデルの逐語的翻訳ではなく、経済学の数式と独立テストから構築。著者MATLABコードは仕様・参照解を確認するためだけに読む。
- 型・モジュール：パラメータ、グリッド、定常解、価格経路、家計政策、分布、残差用workspaceをそれぞれ明示的な型付き構造体で管理。型不安定なDict/Anyをホットループから排除。
- Function barriers：計算は実装関数内に置く。グローバル可変状態なし。配列長と時点を関数引数で指定し、notebook的な暗黙状態を禁止。
- 内部メソッド：rouwenhorst、stationary_egm!、distribution_step!、steady_state_residual!、transition_household_backward!、transition_distribution_forward!、flexible_residual!、fixed_residual!、portfolio_revalue!等の責務を独立させる。ここでの名称は設計候補であり、未作成API。
- Allocation：Matrix{Float64}は資産×所得（第1添字＝資産）で列優先走査。配列の再確保を避け、in-place!・preallocation・cacheを優先。フル密行列5383×5383を無条件に生成しない。定常MarkovはSparseArraysあるいはmatrix-free作用素を使う。
- 数値計算：LinearAlgebra、SparseArrays、Statistics、Testなどの標準ライブラリ優先。個別の正値交易条件のrootはRoots.jl、全期間非線形systemはNonlinearSolve.jl、線形系はLinearSolve.jl、補間はInterpolations.jl等を**候補**として公式資料・API・Juliaバージョンとの互換性を実環境で確認した後に決める。数学的に必要な単調・形状保持補間とパッケージが提供する補間の一致は未確認なので無条件に採用しない。
- 非線形外側：初期は明示Jacobian/数値微分を使った小規模テスト、性能上必要ならJacobian-vector/準Newton/continuationを比較する。ADは補間分岐と借入制約の非微分点を考慮して正当化し、実測前に高速と断定しない。
- Pkg：独立環境のProject.toml/Manifest.tomlをJuliaのPkgから生成しコミット。UUID・バージョンを推測で作らない。再現CIはWindowsとLinux、利用バージョンを明記する。
- コード品質：予約語を識別子に使用しない、Meta.parseallによる構文検証、依存導入・usingによるロード確認、小規模テスト、全スモーク、性能測定、CIの順。実行不可なら「構文未検証・実行未検証」と表示。
- Julia公式資料：https://julialang.org/ 、https://docs.julialang.org/en/v1/manual/performance-tips/ 、https://pkgdocs.julialang.org/v1/environments/ 。候補：https://docs.sciml.ai/NonlinearSolve/stable/ 、https://docs.sciml.ai/LinearSolve/stable/ 、https://juliamath.github.io/Roots.jl/stable/ 、https://juliamath.github.io/Interpolations.jl/stable/ 。

## 18. Verification gatesとDefinition of Done

**Gate 1（定常家計）**：Rouwenhorst確率行和=1、所得定常平均=1、EGMのEuler・制約・予算・単調性、細グリッドの補間境界。定常分布μ≥0、sum μ=1、||Q'μ−μ||∞が設定許容値未満。

**Gate 2（定常GE）**：β・r根探索後のA−K−B、企業価格、賃金条件、国内財市場が事前規定の許容値未満。初期・較正・最終定常の比較表。

**Gate 3（初回移行）**：T=200の各期資産・労働市場、投資・資本蓄積、qと無裁定、家計Backward/Forward、終端、T延長/資産グリッド倍化。Figure 1–3と系列比較。

**Gate 4（突然の信用収縮）**：継承時点、外国信用新経路、Gross Portfolio集計、再評価後家計資産・質量保存、再評価の符号、φ>0等の入力域を検証。変動相場のFigure 4–5を比較。

**Gate 5（固定相場）**：Rotemberg NKPC残差、所得フロー、資源制約（価格調整コストの符号を含む）、資産・労働市場、CEV以前の実物・価格IRFを比較。論文に相違があれば未解決として報告する。

**Gate 6（full paper）**：Figure 6–14、低/中/高leverage、債務通貨建て、22%/5四半期半減期実験、Welfare/CEVの家計別・集計別比較。これはGate 5と別の完了条件。

各Gateは、利用したコミットSHA、Julia実行環境（Project/Manifest）、実行ログ、数値残差、原著比較表、画像ファイル名、既知の不一致、未確認点を保存する。論文との差分を説明できないまま「再現済み」と書かない。

## 19. 後続の論文・アルゴリズム参考文献

- de Ferra, Mitman, Romei (2020), JIE 124, 103303, Sections 2–7, Figures 1–14（研究対象）。
- Carroll (2006), *The Method of Endogenous Gridpoints for Solving Dynamic Stochastic Optimization Problems*（EGM）。
- Rouwenhorst (1995), discrete approximation to persistent income processes（Markov離散化）。
- Boppart, Krusell, Mitman (2018), *Exploiting MIT Shocks in Heterogeneous-Agent Economies: The Impulse Response as a Numerical Derivative*（MIT shock）。
- Auclert, Bardóczy, Rognlie, Straub (2021), *Using the Sequence-Space Jacobian to Solve and Estimate Heterogeneous-Agent Models*（代替の線形化手法。今回の非線形IRFと混同しない）。

**現時点の到達点**：数式・公開コードの読解に基づく詳細仕様と監査計画を作成した。実装・構文解析・実行・結果再現は今後の作業である。

# Nirei (2015) — 数値計算アルゴリズム詳細仕様

> **2026-11-01: 検証は中断・未再現。** 本文は当初の設計資料です。実際の検証結果、未解決問題、統合履歴は [最終アーカイブ](docs/REPLICATION_ARCHIVE_2026-11-01.md) を優先してください。

## 0. スコープ、根拠、完了条件

この仕様は Nirei (2015), *Theoretical Economics* 10, 953–985 の **Section 3、Appendices A–C、Table 2 の Benchmark replication** を実装するためのものである。主たる均衡条件は本文 (7)–(19)、期待均衡の線形化は Appendix A (26)–(37)、企業の投資閾値と均衡選択は (10)、(13)、(14)、(16) と Section 3.4 である。

- **論文に指定あり**: モデル式、企業の離散選択、連続体近似による期待形成、ES1、キャリブレーション、N、T、burn-in、10試行。
- **再現実装上の選択**: 初期有限標本の構成、乱数 seed、数値積分、固定点の効率的な完全列挙、有限企業の厳密資本集計と近似的期待均衡の時間更新の接続、ゼロ投資の扱い、同点処理、数値許容誤差。これらは著者の未公開コードと同一だとは主張しない。

主目的は **11万社 × 1,100期 × 10試行を Julia で実行し、Table 2 の9モーメントおよび10試行間の standard error を得ること**。別個の「検証プロジェクト」を開始せず、本体を完成させる。未実行のものを「再現済み」と称しない。

## 1. 時間と記号の定義

実装では次の時点を採用する。

- 期首: 既存資本 k_{i,t}、集計 K_t、前期に決まった限界効用の状態が存在。
- 当期の期待形成: 有限企業の正確な状態を所与として連続体期待システムを使い、X_t^e と予想要素価格を計算する。論文が E_{t-1}[X_t] と書く情報集合に対応させる。
- 当期企業投資: 新たな a_{i,t+1} が観察され、次期資本 k_{i,t+1} と集計 K_{t+1} が同時に決まる。
- 当期マクロ実現: X_t が確定し、内生投資ショック ε_t を用いて C_t、Y_t、L_t、w_t、R_t を計算する。
- 期末: k_{i,t+1}、a_{i,t+1}、K_{t+1} を保存する。

\(\widetilde Z_t=\log(Z_t/\bar Z)\)。\(R_t\) はグロス金利、\(c_L=(1-1/\eta)(1-\alpha)\)、\(\rho=(1-1/\eta)\alpha/(1-c_L)\)。全企業は同じ \(\delta\) を持つ。企業別 \(\lambda_i\) は時間を通して固定。

**重要な添字上の注意**: 原論文式 (13) の印刷された \(a_{i,t}\) は式 (15)、Section 3.4 と次期資本の生産性加重集計の説明では \(a_{i,t+1}\) になっている。実装では次期実現生産性を使い、差異をログに記録する。文献の原文が完全に一致していると偽らない。

## 2. 構造・キャリブレーション

| parameter | benchmark |
|---|---:|
| alpha | 0.26 |
| delta | 0.02 |
| beta | 0.99 |
| sigma | 1.5 |
| eta | 11 |
| psi | 1 |
| zeta | 2 |
| E(lambda_i) | 1.028 |
| SD before truncation(lambda_i) | 0.004 |
| lambda_i truncation | [1.020, 1.036] |
| SD(log a_i,t) | 0.0005 |
| N | 110000 |
| periods per run | 1100 |
| burn-in | 100 |
| Monte Carlo runs | 10 |

\(\lambda_i\) は \(N(1.028,0.004^2)\) を平均 ±2 標準偏差で切断した分布から各試行の開始時に生成し、試行中固定。\(a_{i,t}\) は企業・時点を通じて i.i.d.、平均を \(E[a^{\rho/\{\alpha(1-\rho)\}}]=1\) に設定する。

\(q=\rho/[\alpha(1-\rho)]\) と置くと \(\log a\sim N(-q\sigma_a^2/2,\sigma_a^2)\) により \(E[a^q]=1\) が成立する。これは Appendix C の指定を lognormal distribution のパラメータに書き下したものである。

\(c_L=(1-1/\eta)(1-\alpha)\)、\(\rho=(1-1/\eta)\alpha/(1-c_L)\) を最初に計算する。有限 \(N\) の実現集計と連続体期待形成を区別する。

## 3. 定数 A の数値積分

論文式 (16):

\[
b_{i,t+1}=A\left[a_{i,t+1}^{\rho/\alpha}
\frac{\lambda_i^\rho-1}{\lambda_i-1}\right]^{1/(1-\rho)}
\]

\[
A=\left\{
E_{\lambda,a}\left[
\left(\frac{\lambda^\rho-1}{(\lambda-1)^\rho}\right)^{1/(1-\rho)}
\frac{a^{\rho/[\alpha(1-\rho)]}}{\rho\log\lambda}
\right]\right\}^{-1/\rho}.
\]

生産性と \(\lambda\) の独立性・前節の正規化を使い、\(A\) は切断正規分布の \(\lambda\) 上の **一変数積分**にする。QuadGK.jl の積分と Distributions.jl の切断正規 pdf を想定する。積分範囲は [1.020,1.036]。差分 \(\lambda^\rho-1\) は \(\exp(\rho\log\lambda)-1\) を使う数値安定な実装（Julia の \`expm1\`、\`log1p\`）を優先する。A は一度だけ計算し、各期では再計算しない。

### 注意

\(\lambda\) の下限 1.020 で \(\lambda(1-\delta)=0.9996<1\) となるため、論文の理論節で用いる厳密な \(\lambda(1-\delta)>1\) は分布の下端付近で満たされない。この境界上の理論条件と Benchmark calibration の間の軽微な齟齬を隠さない。論文の分布を勝手に切り直さない。

## 4. Steady state

家計効用:

\[
U(C,L)=\frac{C^{1-\sigma}(1-\psi L^\zeta)^{1-\sigma}}{1-\sigma}.
\]

計算に用いる定常関係:

\[
\bar R=\beta^{-1},\quad a_x=\frac{\bar X}{\bar Y}=\frac{\delta\rho(1-c_L)}{\bar R-1+\delta},\quad a_c=1-a_x,
\]

\[
\bar K/\bar Y=\frac{1-c_L}{A^{1-\rho}(\bar R-1+\delta)},\qquad
\bar L=\left[\frac{c_L}{\psi(c_L+\zeta a_c)}\right]^{1/\zeta}.
\]

\[
\bar Y=\bar K^\alpha\bar L^{1-\alpha},\quad
\bar C=a_c\bar Y,\quad \bar X=a_x\bar Y,\quad
\bar w=c_L\bar Y/\bar L,
\]

\[
\bar\mu=\bar C^{-\sigma}(1-\psi\bar L^\zeta)^{1-\sigma}.
\]

まず \(K/Y\) と \(\bar L\) が求まり、\(\bar Y=(K/Y)^{\alpha/(1-\alpha)}\bar L\)、\(\bar K=(K/Y)\bar Y\) により水準が求まる。論文 Table 1 の丸め値 \(c_L\simeq0.67\)、\(a_c\simeq0.84\) と比較するが、計算には丸め値を代入しない。

## 5. Appendix A の線形期待均衡

\(a_R=\bar R/(\bar R-1+\delta)\)、
\[
\eta_L=\left[\left(2-\frac1\sigma\right)\frac{c_L}{a_c}+\zeta-1\right]^{-1},\qquad
A_1=1-\alpha-\frac{c_L(\sigma-1)}{\sigma}.
\]

Appendix A 式 (36) と (37) から、論文の時点で
\[
B\begin{pmatrix}E_{-1}\widetilde K_0\\E_{-1}\widetilde\mu_0\end{pmatrix}
=D\begin{pmatrix}\widetilde K_{-1}\\\widetilde\mu_{-1}\end{pmatrix}.
\]

\[
B=\begin{pmatrix}
a_x/\delta & -a_c/\sigma-A_1/[\sigma(\alpha+\eta_L^{-1})]\\
0 & \eta_L/[\sigma(1+\alpha\eta_L)]+a_R/(1-\alpha)
\end{pmatrix},
\]

\[
D=\begin{pmatrix}
\alpha+a_x(1-\delta)/\delta+A_1\alpha/(\alpha+\eta_L^{-1}) &0\\
1/(1+\alpha\eta_L)&a_R/(1-\alpha)
\end{pmatrix}.
\]

\(F=B^{-1}D\)。Julia では、可能な限り明示的な逆行列演算を避け、\`F = B \ D\` として線形方程式を解く。固有分解の安定固有値 \(\xi_s\)（絶対値が1未満）と固有ベクトル \(v_s=(v_K,v_\mu)^\top\) から \(h=v_\mu/v_K\) を計算する。Benchmark の現行実装では \(\xi_s\approx0.96379\)、\(h\approx-0.65037\) である（2015年原論文に印刷された数値ではなく、実装に基づく計算値）。これらは論文が印刷した数値ではない。

**時点を厳密に扱う**: 期首資本 \(K_t\) と、その期首資本に対応する saddle-path 限界効用状態を入力とする場合、次期資本期待値は \(\xi_s\widetilde K_t\)。Appendix A で \(\widetilde K_{-1}\) が predetermined state であることと取り違えない。

\[
E\widetilde K_{t+1}=\xi_s\widetilde K_t,\qquad
E\widetilde X_t=\frac{\xi_s-(1-\delta)}{\delta}\widetilde K_t.
\]

ただし \(E\widetilde X_t=E[\log(X_t/\bar X)\mid\mathcal F_{t-1}]\) であり、厳密な \(\log E[X_t]\) と等しいわけではない。一次近似の水準期待値として \(X_t^e=\bar X\exp(E\widetilde X_t)\) を使う。

同様に、(33),(27),(28),(30),(31) から
\[
E\widetilde L_t=\frac{a_x E\widetilde X_t-\alpha\widetilde K_t-(a_c/\sigma)E\widetilde\mu_t}{A_1},
\]
\[
E\widetilde C_t=\frac{-E\widetilde\mu_t+(\sigma-1)(c_L/a_c)E\widetilde L_t}{\sigma},
\]
\[
E\widetilde Y_t=\alpha\widetilde K_t+(1-\alpha)E\widetilde L_t,
\]
\[
E\widetilde w_t=E\widetilde C_t+(\zeta-1+c_L/a_c)E\widetilde L_t.
\]

金利は原文の Appendix A (31) または (32) から一貫した時点添字で求めること。\(\widetilde R_t=\widetilde\mu_{t-1}-\widetilde\mu_t\) は **実現**した条件付きオイラー関係をそのまま置いたものではない。実現した限界効用ショックを含む場合は、実際の情報集合と条件付き期待を守り、\(R_t\) を過去の限界効用と当期の期待限界効用から計算する。目的は Table 2 に必要な5系列（Y,C,X,L,K）であり、金利の補助的な出力はその後でよい。

## 6. 有限企業経済の初期化

各 Monte Carlo 試行に対して:
1. 同定可能な seed と生成器を確保。
2. \(\lambda_i\) を切断正規分布から生成、時間を通して固定。
3. \(a_{i,0}\) を生産性分布から独立生成。
4. \(s_{i,0}\overset{iid}{\sim}U[0,1)\) を生成。
5. \(b_{i,0}=A[a_{i,0}^{\rho/\alpha}(\lambda_i^\rho-1)/(\lambda_i-1)]^{1/(1-\rho)}\) とし、暫定的に \(k_{i,0}=\bar K b_{i,0}\lambda_i^{s_{i,0}}\) を構成。
6. 有限標本の集計 \(K_0=[N^{-1}\sum_i(a_{i,0}^{1/\alpha}k_{i,0})^\rho]^{1/\rho}\) を記録。

**厳密性の注意**: この暫定初期化は、有限 \(N\) の \(K_0\) と \(\bar K\) の一致を厳密に保証しない。必要なら全 \(k_{i,0}\) を共通比率でスケーリングし、初期の aggregate capital を \(\bar K\) に正規化する。ただしその場合も \(s_i\) の解釈との整合を記録する。論文は初期 \(s_i\) の分布のみ指定し、有限標本のスケール合わせを詳細に指定していない。元の著者の初期化手順だと断言しない。

## 7. 逐次シミュレーション: t=0,...,1099

各期は **期待 → 生産性 → 投資均衡 → ES1 → ショック → マクロ変数 → 状態更新** の順序を変えない。

### 7.1 期待形成

期首の \(K_t\) から \(\widetilde K_t=\log(K_t/\bar K)\) を計算し、Section 5 の係数で期待投資を返す:
\[
X_t^e=\bar X\exp\left(\frac{\xi_s-1+\delta}{\delta}\widetilde K_t\right).
\]
期首の合理的期待資本基準は、連続体近似の厳密な式 (18) に沿って
\[
K_{t+1}^*=(1-\delta)K_t+\frac{A^{\rho-1}}{\rho}X_t^e.
\]
線形化された \(\bar K\exp(\xi_s\widetilde K_t)\) との違いは一次近似によるものであり、両者を無条件に同一視しない。

### 7.2 企業別生産性

毎期 \(a_{i,t+1}\) を Section 2 の分布から新しく生成。共通 aggregate productivity shock は生成しない。

### 7.3 投資閾値

\[
b_i=A\left[a_{i,t+1}^{\rho/\alpha}\frac{\lambda_i^\rho-1}{\lambda_i-1}\right]^{1/(1-\rho)},\quad d_i=(1-\delta)k_{i,t}.
\]
任意の候補 \(K>0\) に対する政策は
\[
\kappa_i(K)=
\begin{cases}
-1 & K\le d_i/(\lambda_i b_i),\\
0 & d_i/(\lambda_i b_i)<K\le d_i/b_i,\\
1 & K>d_i/b_i.
\end{cases}
\]
境界での厳密な等号は原論文の式 (2) の選択ルールに合わせること。

### 7.4 集計反応関数と固定点

\[
\Phi_t(K)=\left[\frac1N\sum_i(a_{i,t+1}^{1/\alpha}d_i\lambda_i^{\kappa_i(K)})^\rho\right]^{1/\rho}.
\]
均衡は \(K=\Phi_t(K)\)。

これを**全切替点の昇順ソート**によって正確に列挙する。各企業の2閾値:
\[
\ell_i=d_i/(\lambda_i b_i),\quad u_i=d_i/b_i.
\]
初めに全企業を -1 状態にして、集計べき乗和
\[
S=N^{-1}\sum_i(a_{i,t+1}^{1/\alpha}d_i/\lambda_i)^\rho
\]
を計算する。閾値を昇順に走査し、\(\ell_i\) の通過で
\[
\Delta S_i=N^{-1}(a_{i,t+1}^{1/\alpha}d_i)^\rho(1-\lambda_i^{-\rho})
\]
を加え、\(u_i\) の通過で
\[
\Delta S_i=N^{-1}(a_{i,t+1}^{1/\alpha}d_i)^\rho(\lambda_i^\rho-1)
\]
を加える。各区間で反応関数は \(S^{1/\rho}\) の一定値なので、それが区間内にある場合に固定点として記録する。複数同一閾値は一括処理し、端点の等号は政策の厳密な規則に従う。

これにより時間計算量は \(O(N\log N)\)、作業メモリは \(O(N)\)。企業を全固定点候補で毎回再集計する \(O(N^2)\) の実装は禁止する。

**これは ES1 を計算するために採用する Julia-native 数値アルゴリズムであり、論文が採用した特定の実装方法と断言しない。**

### 7.5 ES1

各固定点 \(K^{(j)}\) に対応する企業投資と集計投資:
\[
x_i^{(j)}=d_i(\lambda_i^{\kappa_i(K^{(j)})}-1),\quad X^{(j)}=N^{-1}\sum_i x_i^{(j)}.
\]

\[
j^*=\arg\min_{j:X^{(j)}>0}\left|\log X^{(j)}-\log X_t^e\right|.
\]
ただし \(X\le0\) の場合のログは未定義である。**非正投資を任意の小さな正数に置き換えない**。全均衡候補について \(X>0\) の妥当性をチェックし、正の候補がなければ明示的な数値例外として記録する。零・負の候補を排除すること自体も厳密な実数値モデルに追加された判定なので、サイレントに除外しない。

**論文が実際に採用する選択規則は Section 3.4 の ES1 である。** Section 2.4 の \(\Phi(\Phi(\cdot))\) 型の fictitious tâtonnement は ES2 の best-response dynamics であり、Benchmark Table 2 には置き換え不可。

タイの際は \(|\log K^{(j)}-\log K_{t+1}^*|\) の最小値を二次規則に使い、その後は最小の K を選ぶなど決定的規則を設定する（著者の規則ではなく実装の選択）。

### 7.6 厳密な企業別資本の更新

選択された政策を使って
\[
k_{i,t+1}=d_i\lambda_i^{\kappa_i(K_{t+1})},
\quad X_t=N^{-1}\sum_i(k_{i,t+1}-d_i).
\]
\[
K_{t+1}^{firm}=\left[N^{-1}\sum_i(a_{i,t+1}^{1/\alpha}k_{i,t+1})^\rho\right]^{1/\rho}.
\]
固定点 \(K_{t+1}\) と再集計値が一致することを丸め誤差の範囲で確認する。これはメイン処理に含める invariant であり、別の検証作業ではない。

### 7.7 内生投資ショック

\[
\varepsilon_t=\log X_t-\log X_t^e.
\]
Appendix A の線形近似:
\[
\widetilde K_{t+1}^{lin}
=E\widetilde K_{t+1}+\delta\varepsilon_t
=\xi_s\widetilde K_t+\delta\varepsilon_t.
\]

**重要**: \(K_{t+1}^{firm}\) は式 (13) に基づく実現集計資本、\(K_{t+1}^{lin}\) は一次近似であり厳密には同一でない。元の論文はこの二つの逐次接続のコードレベルの仕様を与えていない。実装では企業資本と整合する \(K_{t+1}^{firm}\) を次期期首の集計資本とし、一次近似の差 \(\widetilde K_{t+1}^{lin}-\log(K_{t+1}^{firm}/\bar K)\) を記録する。これを恣意的に相殺する shock correction は作らない。

### 7.8 マクロ変数を計算

論文の log-linear 条件 (27),(28),(30),(33) により:
\[
\widetilde\mu_t=h\log(K_{t+1}^{firm}/\bar K)
\]
を saddle-path condition として適用する。これは実現有限企業資本を saddle path に写す数値実装上の解釈である。

\[
\widetilde X_t=\log(X_t/\bar X),
\]
\[
\widetilde L_t=\frac{a_x\widetilde X_t-\alpha\widetilde K_t-(a_c/\sigma)\widetilde\mu_t}{A_1},
\]
\[
\widetilde C_t=\frac{-\widetilde\mu_t+(\sigma-1)(c_L/a_c)\widetilde L_t}{\sigma},
\]
\[
\widetilde Y_t=\alpha\widetilde K_t+(1-\alpha)\widetilde L_t,
\]
\[
\widetilde w_t=\widetilde C_t+(\zeta-1+c_L/a_c)\widetilde L_t.
\]

**金利の重要な注意**: Appendix A の Euler 方程式 (32) は、実現した限界効用の対数差に関係するが、(7) は厳密には \(R_t^{-1}=\beta U_{C,t}/U_{C,t-1}\) として記載されている。期待価格の計算に使う情報時点に従って評価し、必要なら実現 \(R_t\) と予測 \(E R_t\) を別の系列として保存する。Benchmark Table 2 の主要評価は Y,C,X,L,K である。

\[
Y_t=\bar Y e^{\widetilde Y_t},\quad C_t=\bar C e^{\widetilde C_t},\quad L_t=\bar L e^{\widetilde L_t}
\]
などを保存する。これらは **対数線形近似から復元した水準**であり、有限変動の非線形財市場の厳密解と混同しない。

### 7.9 状態遷移

\[
(\mathbf k_t,\mathbf a_t,K_t,\widetilde\mu_{t-1})
\longleftarrow
(\mathbf k_{t+1},\mathbf a_{t+1},K_{t+1}^{firm},\widetilde\mu_t)
\]

を更新し、次期へ進む。\(\lambda_i\) は不変。実現ショック・投資・5系列は毎期保存。企業全社×全期間をメモリ上に保持せず、企業の現在状態とマクロ時系列だけ保持する。

## 8. 複数均衡と fictitious tâtonnement の関係

Section 2.4 の fictitious tâtonnement は期待資本を初期点とした best-response dynamics。これは Section 2 の **ES2** へ収束する仕組みとして論じられる。Section 3.4 の Benchmark replication は **ES1**、すなわち実現集計投資の期待からの乖離最小。両者を混同して結果を出さない。

Section 2 の分析に現れる Poisson branching distribution や tail exponent は理論的説明に重要だが、Table 2 を計算するために毎期 Poisson branching process を別途 simulation するわけではない。企業の非連続投資と集計固定点から endogenous avalanche が生成される。

## 9. Monte Carlo と出力

各試行に対し 1,100 期生成し、0–99期を破棄、100–1099期の 1,000 観測を用いる。計算するモーメント:

\[
\operatorname{SD}(\log Y),\operatorname{SD}(\log C),\operatorname{SD}(\log X),\operatorname{SD}(\log L),\operatorname{SD}(\log K),
\]
\[
\operatorname{Corr}(\log Y,\log C),\operatorname{Corr}(\log Y,\log X),\operatorname{Corr}(\log Y,\log L),\operatorname{Corr}(\log Y,\log K).
\]

Table 2 の Benchmark は順に 2.28%, 0.86%, 6.60%, 3.06%, 0.26%, 0.822, 0.979, 0.795, −0.018。各試行でモーメントを求めてから10試行の平均と試行間 standard error（sample SD / √10）を計算する。追加的な HP filtering は論文が指定していないため、初期仕様に含めない。

保存する出力:
- 10試行の各マクロ時系列（time, logY, logC, logX, logL, logK, wage, rate, epsilon）
- 10試行別の9モーメント、平均と standard error
- Table 2 比較表
- 乱数 seed、Julia version、Project.toml / Manifest.toml、数値許容誤差、ES1 選択候補数、例外記録
- 実行可能なエントリーポイントと必要最低限の使用方法

**完了条件**: Benchmark の11万社・1,100期×10試行が実行され、Table 2 比較の結果ファイルが生成されている。どの結果も実行を確認するまで未達成と表記する。

## 10. Julia-native 実装規約

1. Julia の公式サイト・マニュアルを基準に設計する。Python/R/MATLAB等のコードを Julia に機械翻訳することは禁止。
2. モデルの実装に使う言語は Julia のみ。新しい依存ライブラリを追加する前に公式リポジトリとドキュメントの API、サポートバージョンを確認する。
3. ループを関数で包み、グローバル変数での型不安定な計算を避ける。数値は具体型の Float64、整数は Int 等を原則とし、安易な \`Any\` 配列を使わない。
4. 11万社を対象とする期ごとの計算では、作業ベクトルを再利用し、インプレース更新、\`sortperm!\` や \`sort!\` の適切な使用、型安定なスカラー計算、連続アクセスを優先する。ランダムな配列確保や closure の暗黙キャプチャはプロファイラで必要性を判断する。
5. 2N 個の閾値の処理は \(O(N\log N)\)、候補探索は \(O(N)\)、作業記憶は \(O(N)\) を基準とする。各候補で N 社を再計算しない。
6. \`B \ D\` と \`eigen\` に LinearAlgebra を使う。切断正規分布に Distributions.jl、定数 A の一次元数値積分に QuadGK.jl を使用候補とする。乱数再現性を重視するなら StableRNGs.jl を検討する。Project.toml に明示的に記録する。
7. ソート順が違うだけで ES1 の選択が変化しないよう、同一閾値のイベント一括処理、決定的 tie-breaking、有限値・正値の必要最低限の guard を実装する。
8. 安定な差分・対数のため Julia Base の \`expm1\`/\`log1p\` を利用する。コンピューティング環境は Julia 1.x、実際に使った小バージョンを記録する。
9. 推測による外生ショック、恣意的な smoothing、正数への強制置換、ES2 への無断変更、元論文と異なる固定点ルールを採用しない。
10. 計算時間の短縮は benchmarking/profiling の実測に基づき、モデル式を変更しない。

## 11. Julia 公式資料・パッケージ一次資料

- Julia language: https://julialang.org/
- Performance tips: https://docs.julialang.org/en/v1/manual/performance-tips/
- LinearAlgebra: https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/
- Random: https://docs.julialang.org/en/v1/stdlib/Random/
- Sorting: https://docs.julialang.org/en/v1/base/sort/
- Julia Pkg: https://pkgdocs.julialang.org/v1/
- Distributions.jl: https://juliastats.org/Distributions.jl/stable/
- QuadGK.jl: https://juliamath.github.io/QuadGK.jl/stable/
- StableRNGs.jl: https://github.com/JuliaRandom/StableRNGs.jl
- Original article: https://doi.org/10.3982/TE1611

**留意**: 本文は実装のための数値仕様であり、著者の未公開オリジナルコードを複製したものではない。原論文の明示的指定と、再現に必要な未指定事項を区別した。

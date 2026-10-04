#  -- P0-10 连续 transport 指标面板 与 P0-11 联合模型 sanity 复核

生成脚本：`scripts/analysis/transport_summary_joint.py`；seed = 20260911；全部数值为本地真跑复算（frozen // 代数，患者级 S0/标签来自 `data/_workspace`，与  包内副本逐位一致），并与冻结档案交叉验证。

## P0-10 连续 transport 指标面板（修稿包点名表）

within 同管线参照 =  `within_cohort_same_pipeline` Efron 0.632 AUC：IMvigor210 0.567、GSE165252 0.444、GSE78220 0.464（本脚本已用 600-rep OOB 复算逐位复现）。

校准约定（）：logit 尺度 logistic 重校准 logit(P(y=1)) = a + b·logit(p̂)；a = calibration-in-the-large，b = calibration slope。注意  档案中的 calib_slope 为概率尺度回归，与本表不可直接比较。


### 目标队列 IMvigor210

| 模型 | n (ev) | AUC | ΔAUC vs within (95% CI)* | calib 截距 a | calib 斜率 b | Brier | log-loss |
|---|---|---|---|---|---|---|---|
| primary frozen ( F1) | source | source / n-a | n-a | n-a | n-a | n-a | n-a |
| M0 | 298 (68) | 0.500 | -0.067 [-0.166, 0.012] | -0.71 | null | 0.198 | 0.587 |
| M2 | 298 (68) | 0.525 | -0.042 [-0.172, 0.064] | 4.25 | 10.70 | 0.197 | 0.586 |
| M3 | 298 (68) | 0.498† | -0.069 [-0.219, 0.042] | 2.39 | 2.97 | 0.177 | 0.540 |

primary frozen 行：IMvigor210 为训练队列（source），无外部自测，记 source / n-a。其余行 = LOCO（训练 = GSE165252+GSE78220）。


### 目标队列 GSE165252

| 模型 | n (ev) | AUC | ΔAUC vs within (95% CI)* | calib 截距 a | calib 斜率 b | Brier | log-loss |
|---|---|---|---|---|---|---|---|
| primary_frozen_F1 | 71 (23) | 0.501 | 0.057 [-0.231, 0.157] | -0.53 | 0.16 | 0.232 | 0.672 |
| M0 | 71 (23) | 0.500 | 0.056 [-0.172, 0.059] | 0.36 | null | 0.224 | 0.643 |
| M2 | 71 (23) | 0.516 | 0.072 [-0.215, 0.165] | -0.41 | 0.29 | 0.225 | 0.647 |
| M3 | 71 (23) | 0.449† | 0.005 [-0.255, 0.056] | -3.58 | -3.74 | 0.227 | 0.650 |


### 目标队列 GSE78220

| 模型 | n (ev) | AUC | ΔAUC vs within (95% CI)* | calib 截距 a | calib 斜率 b | Brier | log-loss |
|---|---|---|---|---|---|---|---|
| primary_frozen_F1 | 25 (13) | 0.410 | -0.054 [-0.494, 0.142] | -0.56 | -0.52 | 0.347 | 0.956 |
| M0 | 25 (13) | 0.500 | 0.036 [-0.319, 0.096] | 1.20 | null | 0.324 | 0.864 |
| M2 | 25 (13) | 0.449 | -0.015 [-0.454, 0.179] | -0.88 | -0.86 | 0.332 | 0.892 |
| M3 | 25 (13) | 0.451† | -0.013 [-0.461, 0.150] | 0.10 | -1.65 | 0.281 | 0.776 |

† M3（监督式目标域重校准，100 次 50/50 重复）为重复均值；使用目标标签，仅作重校准上限对照。

*ΔAUC = 外部 AUC − within 0.632。CI 方法：外部侧 2000 次患者 bootstrap（冻结预测；M0 常数 0.5 → 外部侧退化；M3 用 100 个 repeat AUC 的 repeat 级重抽样），within 侧为 Efron-0.632 估计量的 2000×600 嵌套 bootstrap，配对差值分布 2.5/97.5 百分位，seed 20260911。敏感性（within 固定为档案点值）CI 见 JSON `delta_auc.*.ci95_within_held_fixed`。


### ΔAUC 关键数字（外部 − within）

| 队列 | primary frozen | M0 | M2 | M3 |
|---|---|---|---|---|
| IMvigor210 | n-a | -0.067 [-0.166, +0.012] | -0.042 [-0.172, +0.064] | -0.069 [-0.219, +0.042] |
| GSE165252 | +0.057 [-0.231, +0.157] | +0.056 [-0.172, +0.059] | +0.072 [-0.215, +0.165] | +0.005 [-0.255, +0.056] |
| GSE78220 | -0.054 [-0.494, +0.142] | +0.036 [-0.319, +0.096] | -0.015 [-0.454, +0.179] | -0.013 [-0.461, +0.150] |

读法：所有 ΔAUC 95% CI 均跨 0（任何模型、任何队列都无对外部优于 within 的证据）；GSE78220 的 primary frozen ΔAUC 为负（−0.054），外部 transport 未获支持，与主文结论一致。


### 汇总 AUC（k=3，descriptive only）

-  DL 随机效应：pooled AUC = 0.498 [0.365, 0.631]，I²(Q) = 0.752，Q 检验 P = 0.0179。仅作描述性汇总，异质性显著，不支持单一可运输 AUC 的解读。


### P0-10 交叉验证

- 复算 vs 冻结档案逐值核对：50/50 项通过（容差 1e-9；含 F1 AUC vs 、M0/M2/M3 全指标 vs 、within 0.632 复算、M3 冻结函数双跑互验）。**P0-10 判定：PASS**。


---

## P0-11 GSE186143 联合模型 sanity 复核（逐项）

- 总判定：**PASS_WITH_NOTES**（25 折联合 vs 独立 CV 全量重跑 1.1s；四格与 Fisher 从 `geo_data/GSE186143_bulk_clinical.json` 原始文件独立重建）。

| # | 复核项 | 冻结值 | 复算值 | 判定 |
|---|---|---|---|---|
| 1 | 2x2 cells (raw file ==  ==  == 10/17/15/5) | {'DCB/severe': 10, 'DCB/no_severe': 17, 'NDB/no_severe': 15, 'NDB/severe': 5} | {'DCB/severe': 10, 'DCB/no_severe': 17, 'NDB/no_severe': 15, 'NDB/severe': 5} | PASS |
| 2 | marginals 27/20/15/32 | {'response_DCB': 27, 'response_NDB': 20, 'severe_irAE_yes': 15, 'severe_irAE_no': 32} | {'response_DCB': 27, 'response_NDB': 20, 'severe_irAE_yes': 15, 'severe_irAE_no': 32} | PASS |
| 3 | Fisher OR = 1.7647 | 1.76471 | 1.76471 | PASS |
| 4 | Fisher exact two-sided P = 0.5293 (hand hypergeometric) | 0.529332 | 0.529332 | PASS |
| 5 | Fisher exact two-sided P (scipy cross-check) | 0.529332 | 0.529332 | PASS |
| 6 | 25-fold per-fold delta log-loss reproduces (fold rule verified) | identical per-fold values (max |diff| = 0) | max |diff| = 0.00e+00 | PASS |
| 7 | pooled joint log-loss 1.3452 | 1.34518 | 1.34518 | PASS |
| 8 | pooled independent log-loss 1.3560 | 1.35597 | 1.35597 | PASS |
| 9 | joint better in 12/25 folds | 12 | 12 | PASS |
| 10 | sklearn replication (1.3549/1.3450, 12/25) | [1.3548948008620516, 1.3449780247480114, 12] | [1.3548948008620516, 1.3449780247480114, 12] | PASS |
| 11 | mean fold-level delta log-loss +0.000479 | 0.000478888 | 0.000478888 | PASS |
| 12 | sign-permutation P = 0.986 (verbatim  code path, frozen deltas) | 0.98626 | 0.98626 | PASS |
| 13 | exact two-sided sign test P = 1.0 | 1 | 1 | PASS |
| 14 | sign-permutation P on the RE-RUN deltas | 0.98626 | 0.98626 | PASS |
| 15 | fold sizes pattern 10/10/9/9/9 per repeat (round-robin on cells 10/17/15/5) | 10,10,9,9,9 x5 | [10, 10, 9, 9, 9] | PASS |
| 16 | identical folds for both model families | shared fold splits in the frozen script | per-fold delta reproduction implies shared splits | PASS |
| 17 | smallest cell NDB/severe n=5 (n=47 joint-evaluable) | 5 of 47 | 5 of 47 | PASS |
| 18 | 4-class multinomial ~9 effective free params vs 6 for the independent pair | 9 vs 6 | 9 vs 6 (derivation verified) | PASS |
| 19 | smallest TRAINING cell n=4 (5 minus 1 held out) | 4 | 4 | PASS |
| 20 | events-per-parameter far below 10 for every joint cell | < 10 | 5.22222 | PASS |
| 21 | claim 'cell proportions carry 95% CIs wider than +/-0.2' | claimed in the frozen identifiability_note | max CI half-width = 0.153 (full widths 0.196-0.288) | FAIL |
| 22 | frozen wording stays within stress-test / descriptive bounds (no joint advantage claim, no mechanism conclusion) | verdict + identifiability_note +  interpretation hedged | all three frozen texts re-read verbatim | PASS |

**问题清单**：

- identifiability_note sentence 'cell proportions carry 95% CIs wider than +/-0.2' is not numerically supported (max half-width 0.153; full widths 0.196-0.288). Word-level fix recommended; conclusion unchanged.


### ED9（joint fold occupancy / identifiability）数值

- 联合可评估 n = 47（NE 6 例剔除、未填补）；四格 DCB/severe 10、DCB/no 17、NDB/no 15、NDB/severe 5；边际 DCB 27 / NDB 20、severe 15 / no 32。
- E–T 关联（描述性）：OR = 1.7647，Fisher 双侧 P = 0.5293（scipy 交叉验证 0.5293）——不支持 efficacy–toxicity 耦合结论。
- 折占用（每 repeat 5 折，类 10/17/15/5 round-robin；DCB/severe+DCB/no+NDB/no+NDB/severe）：2+4+3+1 (n=10)、2+4+3+1 (n=10)、2+3+3+1 (n=9)、2+3+3+1 (n=9)、2+3+3+1 (n=9)；每折最小类占用 = 1，全部 25 折最小训练格 n = 4；两类模型共享同一折划分（逐折复现证明）。
- 逐格比例 95% Clopper-Pearson CI：DCB/severe 10/47 = 0.213 [0.107, 0.357]; DCB/no_severe 17/47 = 0.362 [0.227, 0.515]; NDB/no_severe 15/47 = 0.319 [0.191, 0.471]; NDB/severe 5/47 = 0.106 [0.035, 0.231]。
- 联合 vs 独立：pooled log-loss 1.3452 vs 1.3560，joint 赢 12/25 折，sign-permutation P = 0.9863（重跑折上同值 0.9863）。
- 可辨识度：有效参数 9（联合）vs 6（独立），events-per-parameter = 47/9 = 5.2 < 10；最小格 n=5、最小训练格 n=4。逐折占用矩阵见 JSON `p0_11.ed9.fold_occupancy_by_repeat`。


### 结论

- P0-10：连续面板已产出（表 1–3），ΔAUC 逐队列全部跨 0；全部 50 项档案交叉验证通过。
- P0-11：PASS_WITH_NOTES——四格/Fisher/log-loss/12-of-25 折/permutation P/折构造全部复现；唯一发现是 identifiability_note 中 “95% CI wider than ±0.2”一句的字面数值不成立（最大半宽 0.15），属措辞层修正，不影响“n=47 联合结局不可辨识、只作 stress-test/假设生成”的结论边界。

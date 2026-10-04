# — P0-3 阈值敏感性与 P0-6 零值处理/CLR 敏感性

seed = 20260911；bootstrap = 2000；M3/M6 再校准重复 = 100；全部数值由 `scripts/analysis/threshold_zero_sensitivity.py` 实算（冻结档案交叉验证），无硬编码。

## P0-3 校准斜率容差敏感性（判据②不依赖 0.10 单点值）

规则（ 冻结）：① 外部 AUC bootstrap 95% CI 有利方向排除 0.5；② |logit 尺度校准 slope − 1| ≤ tol。判定矩阵对 tol ∈ {0.05, 0.10, 0.20} 全枚举。

复算与冻结  交叉验证：max |Δslope| = 0.00e+00，max |ΔAUC| = 0.00e+00（逐字复用  代数 + 冻结 S0）。

### 表 1（连续效应）：logit 尺度校准 slope 点值与 bootstrap 95% CI

| 模型 | 队列 | slope（ 冻结点值） | slope（本处复算） | bootstrap 95% CI | n_valid |
|---|---|---|---|---|---|
| M0 | IMvigor210 | N/A | N/A | N/A | 0 |
| M0 | GSE165252 | N/A | N/A | N/A | 0 |
| M0 | GSE78220 | N/A | N/A | N/A | 0 |
| M1 | IMvigor210 | 22.5611 | 22.5611 | [-1.518, 46.158] | 2000 |
| M1 | GSE165252 | 0.2112 | 0.2112 | [-1.999, 2.402] | 2000 |
| M1 | GSE78220 | -0.5513 | -0.5513 | [-6.627, 4.265] | 2000 |
| M2 | IMvigor210 | 10.7039 | 10.7039 | [-10.353, 31.606] | 2000 |
| M2 | GSE165252 | 0.2909 | 0.2909 | [-1.729, 2.212] | 2000 |
| M2 | GSE78220 | -0.8590 | -0.8590 | [-4.667, 2.147] | 2000 |
| M3 | IMvigor210 | 2.9652 | 2.9652 | [-3.629, 4.809] | 2000 |
| M3 | GSE165252 | -3.7390 | -3.7390 | [-4.474, 5.242] | 2000 |
| M3 | GSE78220 | -1.6479 | -1.6479 | [-5.924, 7.962] | 2000 |
| M4 | IMvigor210 | 6.5805 | 6.5805 | [-16.976, 29.327] | 2000 |
| M4 | GSE165252 | 0.2922 | 0.2922 | [-1.628, 2.185] | 2000 |
| M4 | GSE78220 | -0.8635 | -0.8635 | [-5.076, 1.996] | 2000 |
| M5a | IMvigor210 | 0.1509 | 0.1509 | [-13.084, 12.162] | 2000 |
| M5a | GSE165252 | 0.2904 | 0.2904 | [-1.984, 2.474] | 2000 |
| M5a | GSE78220 | -0.8932 | -0.8932 | [-6.363, 2.868] | 2000 |
| M5b | IMvigor210 | -0.1383 | -0.1383 | [-13.201, 12.251] | 2000 |
| M5b | GSE165252 | 0.2920 | 0.2920 | [-1.943, 2.538] | 2000 |
| M5b | GSE78220 | -0.8938 | -0.8938 | [-6.121, 2.898] | 2000 |
| M6 | IMvigor210 | -0.8344 | -0.8344 | [-4.028, 4.571] | 2000 |
| M6 | GSE165252 | 0.4638 | 0.4638 | [-4.002, 5.707] | 2000 |
| M6 | GSE78220 | -0.1960 | -0.1960 | [-6.197, 7.835] | 2000 |

bootstrap 方法：M1/M2/M4/M5a/M5b = 固定预测下对目标队列患者重采样 (y,p) 并重算 logit 尺度 recalibration slope；M3/M6 = 患者重采样后完整重跑 100 次分层 50/50 再校准流程取每 draw 的 eval-half 均值 slope；M0 = 常数预测，slope 不可估。

### 表 2（判定矩阵 A）：判据② |slope−1| ≤ tol 的 PASS/FAIL

| 模型 | 队列 | tol=0.05 | tol=0.1 | tol=0.2 |
|---|---|---|---|---|
| M0 | IMvigor210 | N/A | N/A | N/A |
| M0 | GSE165252 | N/A | N/A | N/A |
| M0 | GSE78220 | N/A | N/A | N/A |
| M1 | IMvigor210 | FAIL | FAIL | FAIL |
| M1 | GSE165252 | FAIL | FAIL | FAIL |
| M1 | GSE78220 | FAIL | FAIL | FAIL |
| M2 | IMvigor210 | FAIL | FAIL | FAIL |
| M2 | GSE165252 | FAIL | FAIL | FAIL |
| M2 | GSE78220 | FAIL | FAIL | FAIL |
| M3 | IMvigor210 | FAIL | FAIL | FAIL |
| M3 | GSE165252 | FAIL | FAIL | FAIL |
| M3 | GSE78220 | FAIL | FAIL | FAIL |
| M4 | IMvigor210 | FAIL | FAIL | FAIL |
| M4 | GSE165252 | FAIL | FAIL | FAIL |
| M4 | GSE78220 | FAIL | FAIL | FAIL |
| M5a | IMvigor210 | FAIL | FAIL | FAIL |
| M5a | GSE165252 | FAIL | FAIL | FAIL |
| M5a | GSE78220 | FAIL | FAIL | FAIL |
| M5b | IMvigor210 | FAIL | FAIL | FAIL |
| M5b | GSE165252 | FAIL | FAIL | FAIL |
| M5b | GSE78220 | FAIL | FAIL | FAIL |
| M6 | IMvigor210 | FAIL | FAIL | FAIL |
| M6 | GSE165252 | FAIL | FAIL | FAIL |
| M6 | GSE78220 | FAIL | FAIL | FAIL |

### 表 3（判定矩阵 B）：整体运输判定（①且②）

| 模型 | 队列 | tol=0.05 | tol=0.1 | tol=0.2 |
|---|---|---|---|---|
| M0 | IMvigor210 | FAIL | FAIL | FAIL |
| M0 | GSE165252 | FAIL | FAIL | FAIL |
| M0 | GSE78220 | FAIL | FAIL | FAIL |
| M1 | IMvigor210 | FAIL | FAIL | FAIL |
| M1 | GSE165252 | FAIL | FAIL | FAIL |
| M1 | GSE78220 | FAIL | FAIL | FAIL |
| M2 | IMvigor210 | FAIL | FAIL | FAIL |
| M2 | GSE165252 | FAIL | FAIL | FAIL |
| M2 | GSE78220 | FAIL | FAIL | FAIL |
| M3 | IMvigor210 | FAIL | FAIL | FAIL |
| M3 | GSE165252 | FAIL | FAIL | FAIL |
| M3 | GSE78220 | FAIL | FAIL | FAIL |
| M4 | IMvigor210 | FAIL | FAIL | FAIL |
| M4 | GSE165252 | FAIL | FAIL | FAIL |
| M4 | GSE78220 | FAIL | FAIL | FAIL |
| M5a | IMvigor210 | FAIL | FAIL | FAIL |
| M5a | GSE165252 | FAIL | FAIL | FAIL |
| M5a | GSE78220 | FAIL | FAIL | FAIL |
| M5b | IMvigor210 | FAIL | FAIL | FAIL |
| M5b | GSE165252 | FAIL | FAIL | FAIL |
| M5b | GSE78220 | FAIL | FAIL | FAIL |
| M6 | IMvigor210 | FAIL | FAIL | FAIL |
| M6 | GSE165252 | FAIL | FAIL | FAIL |
| M6 | GSE78220 | FAIL | FAIL | FAIL |

判据①（AUC CI 有利方向排除 0.5）：全部 24 格均为 FAIL（tolerance 无关， 冻结值）。

**判定随容差翻转的格子数：判据② 0 格；整体运输判定 0 格。**

**主结论是否依赖 0.10 容差：否。** 在 0.05/0.10/0.20 三档容差下，8 模型 × 3 队列共 24 格的判据②与整体判定均无任何翻转；且无任何格在任何容差下被许可运输。所有 slope 点值与 1 的距离最小为 0.536，远超 0.20；结论对容差选择完全不敏感。

### 补充：slope CI 相对容差边界的稳健性（连续视角）

| 模型 | 队列 | tol=0.05 | tol=0.1 | tol=0.2 |
|---|---|---|---|---|
| M0 | IMvigor210 | not_estimable | not_estimable | not_estimable |
| M0 | GSE165252 | not_estimable | not_estimable | not_estimable |
| M0 | GSE78220 | not_estimable | not_estimable | not_estimable |
| M1 | IMvigor210 | inconclusive | inconclusive | inconclusive |
| M1 | GSE165252 | inconclusive | inconclusive | inconclusive |
| M1 | GSE78220 | inconclusive | inconclusive | inconclusive |
| M2 | IMvigor210 | inconclusive | inconclusive | inconclusive |
| M2 | GSE165252 | inconclusive | inconclusive | inconclusive |
| M2 | GSE78220 | inconclusive | inconclusive | inconclusive |
| M3 | IMvigor210 | inconclusive | inconclusive | inconclusive |
| M3 | GSE165252 | inconclusive | inconclusive | inconclusive |
| M3 | GSE78220 | inconclusive | inconclusive | inconclusive |
| M4 | IMvigor210 | inconclusive | inconclusive | inconclusive |
| M4 | GSE165252 | inconclusive | inconclusive | inconclusive |
| M4 | GSE78220 | inconclusive | inconclusive | inconclusive |
| M5a | IMvigor210 | inconclusive | inconclusive | inconclusive |
| M5a | GSE165252 | inconclusive | inconclusive | inconclusive |
| M5a | GSE78220 | inconclusive | inconclusive | inconclusive |
| M5b | IMvigor210 | inconclusive | inconclusive | inconclusive |
| M5b | GSE165252 | inconclusive | inconclusive | inconclusive |
| M5b | GSE78220 | inconclusive | inconclusive | inconclusive |
| M6 | IMvigor210 | inconclusive | inconclusive | inconclusive |
| M6 | GSE165252 | inconclusive | inconclusive | inconclusive |
| M6 | GSE78220 | inconclusive | inconclusive | inconclusive |

（fail_robust = slope CI 完全落在 [1±tol] 之外；pass_robust = 完全落在界内；inconclusive = CI 跨越边界；not_estimable = 常数预测。）

## P0-6 零值处理 / CLR 敏感性（70 未测量基因伪影结论的稳健性）

锚定基因集 =  逐字公共交集 228 基因；其中 70 个在锚定训练集（IMvigor210 标记样本 n=298，即 PCA 实际拟合集，与冻结  一致）中全零（未测量），158 个已测量。补充核查：在全部 348 样本矩阵上全零基因为 67 个——有 3 个基因仅在 50 个未标记样本中出现非零值，不参与锚定 PCA，不影响伪影判定。

### 表 4：三方案对比

| 方案 | 基因数 | 70 基因 PC1 loading 质量占比 | GSE165252 AUC [95% CI] | GSE78220 AUC [95% CI] | 伪影判定 | 运输判定 |
|---|---|---|---|---|---|---|
| S1 原始（log1p/零保持） | 228 | 0.8866 | 0.514 [0.363, 0.664] | 0.455 [0.215, 0.694] | present | 未恢复 |
| S2 multiplicative replacement（0.65×非零最小值） | 228 | 0.5312 | 0.535 [0.384, 0.682] | 0.481 [0.253, 0.731] | present | 未恢复 |
| S3 measurable-only（158 基因） | 158 | 0.0000 | 0.473 [0.324, 0.634] | 0.641 [0.407, 0.865] | removed_by_construction | 未恢复 |

交叉验证  冻结值：S1 mass70 |Δ| = 1.11e-16（冻结 0.886611 vs 复算 0.886611）；S1 外部 AUC |Δ| = GSE165252 0.00e+00, GSE78220 0.00e+00；S3 外部 AUC |Δ| = GSE165252 0.00e+00, GSE78220 0.00e+00。

S2 替换统计（对原始值、CPM 之前执行）：IMvigor210：25715 个零 cells（121 个散发零基因）被替换，70 个全零基因按预设保持为零；GSE165252：417 个零 cells（59 个散发零基因）被替换，0 个全零基因按预设保持为零；GSE78220：195 个零 cells（26 个散发零基因）被替换，0 个全零基因按预设保持为零。

### 稳健性结论

**结论 1（伪影是否存在）：稳健。** 原始流程下 70 未测量基因占锚定 PC1 loading 质量比例 = 0.8866；multiplicative replacement 下 = 0.5312。由于这 70 个基因在锚定训练集（IMvigor210 标记样本）中不存在任何非零值，替换规则按预设无法为其代入值（无非零最小值可乘），保持为零后经 CLR 仍成为 library-size 代理列——伪影由**结构性未测量**驱动，而非零值处理（log1p/伪计数 vs 乘法替换）这一选择驱动；measurable-only 按构造移除该块。

**结论 2（剔除伪影是否恢复运输）：稳健。** 三方案下两个外部队列的 AUC bootstrap 95% CI 均未在有利方向排除 0.5（original_log1p：未恢复；multiplicative_replacement_0p65：未恢复；measured_only_158：未恢复）。即：无论是否剔除伪影、无论零值如何处理，S0 的外部运输失败均持续存在。

### ED5（zero-handling sensitivity）数值

全部数值见 `threshold_zero_sensitivity.json` 的 `ed5_numbers` 块：三方案的 外部 AUC [95% CI]（GSE165252/GSE78220）、70 基因 PC1 loading 质量占比、基因数、零值替换统计与 PC1 解释方差比。

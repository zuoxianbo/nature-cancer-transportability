#  P0-2 主阳性对照（恢复基因身份）— Cristescu 18 基因 T 细胞炎性 GEP

日期：2026-09-24 ｜ seed 20260911 ｜ 脚本：scripts/analysis/positive_control_gep18.py ｜ 冻结 JSON：data/positive_control_gep18.json

**背景：** 原 P0-2 阳性对照（positive_control.py）因 IMvigor210 exmat 基因行身份被上游匿名化破坏而作废。imvigor_gene_recovery.py 已恢复身份（exmat 第 j 行 == fData 按 symbol 大写字母序排序后的第 j 个基因；性别锚定门 PASS、镜像 31085/31085 行数值恒等匹配、上游溯源逐字节一致）。本脚本按原设计重跑：**唯一改动是 IMvigor210 加载器换为恢复钥匙选行**，指标代数、判据、seed 全部原样复用（/ 逐字导入）。

## 1. 恢复钥匙交叉核验

- 脚本内重建的 symbol（大写）字母序钥匙与项目资产 `exmat_row_key.csv`（31085 行）对 6 个签名基因逐一比对：**全部一致**。

## 2. 三队列 6 基因可得性（恢复钥匙下）

| gene | IMvigor210（恢复钥匙，exmat 行号 1-based） | GSE165252（Ensembl） | GSE78220（symbol） |
|---|---|---|---|
| CD8A | YES（行 3133） | YES（ENSG00000153563） | YES |
| CMKLR1 | YES（行 3790） | YES（ENSG00000174600） | YES |
| CXCR6 | YES（行 4433） | YES（ENSG00000172215） | YES |
| HLA-DQA1 | YES（行 8344） | YES（ENSG00000196735） | YES |
| HLA-DRB1 | YES（行 8350） | YES（ENSG00000196126） | YES |
| HLA-E | YES（行 8355） | YES（ENSG00000204592） | YES |
| CCL5 | YES（行 2980） | YES（ENSG00000271503） | YES |
| CXCL9 | YES（行 4426） | YES（ENSG00000138755） | YES |
| CXCL10 | YES（行 4414） | YES（ENSG00000169245） | YES |
| IDO1 | YES（行 8734） | YES（ENSG00000131203） | YES |
| PSMB10 | YES（行 23713） | YES（ENSG00000205220） | YES |
| STAT1 | YES（行 27397） | YES（ENSG00000115415） | YES |
| TIGIT | YES（行 28173） | YES（ENSG00000181847） | YES |
| CD27 | YES（行 3077） | YES（ENSG00000139193） | YES |
| CD274 | YES（行 3079） | YES（ENSG00000120217） | YES |
| CD276 | YES（行 3080） | YES（ENSG00000103855） | YES |
| LAG3 | YES（行 9979） | YES（ENSG00000089692） | YES |
| PDCD1LG2 | YES（行 22453） | YES（ENSG00000197646） | YES |

进入签名的基因：CD8A, CMKLR1, CXCR6, HLA-DQA1, HLA-DRB1, HLA-E, CCL5, CXCL9, CXCL10, IDO1, PSMB10, STAT1, TIGIT, CD27, CD274, CD276, LAG3, PDCD1LG2。缺失：无。
## 3. IMvigor210 锚定（n=298 标注样本）与 within 表现

- 冻结模型：logit(P) = -1.2392 + 0.3712 × signature（ridge λ=1/n=1/298，截距不惩罚）。
- within 表观 AUC = **0.5798**；Efron 0.632 AUC = **0.5728**（n_boot=600，seed 20260911；事件 68/298）。

| 队列（within，常规队列内 z，描述性） | 表观 AUC | mean OOB | Efron 0.632 |
|---|---|---|---|
| GSE165252 | 0.5525 | 0.4804 | **0.5070** |
| GSE78220 | 0.5449 | 0.4337 | **0.4746** |

## 4. 冻结运输：逐队列指标与两判据判定

### primary_logcpm_z_frozen_all6

| target | AUC [95% BCa CI] | slope b | CITL a | Brier | log-loss | 主判据 C1∧C2 | ±0.05 | ±0.10 | ±0.20 |
|---|---|---|---|---|---|---|---|---|---|
| GSE165252 | **0.5462** [0.4034,0.6966] | 0.517 | -0.407 | 0.2201 | 0.6306 | FAIL fail | fail | fail |
| GSE78220 | **0.5705** [0.3179,0.7885] | 0.696 | 1.070 | 0.3384 | 0.9290 | FAIL fail | fail | fail |

### sensitivity_raw_z_frozen

| target | AUC [95% BCa CI] | slope b | CITL a | Brier | log-loss | 主判据 C1∧C2 | ±0.05 | ±0.10 | ±0.20 |
|---|---|---|---|---|---|---|---|---|---|
| GSE165252 | **0.4946** [0.3569,0.6447] | -0.023 | -0.674 | 0.4542 | 2.3411 | FAIL fail | fail | fail |
| GSE78220 | **0.5833** [0.3002,0.7941] | -0.370 | -0.424 | 0.3508 | 0.9408 | FAIL fail | fail | fail |

### sensitivity_logcpm_z_defined_genes

| target | AUC [95% BCa CI] | slope b | CITL a | Brier | log-loss | 主判据 C1∧C2 | ±0.05 | ±0.10 | ±0.20 |
|---|---|---|---|---|---|---|---|---|---|
| GSE165252 | **0.5462** [0.4034,0.6966] | 0.517 | -0.407 | 0.2201 | 0.6306 | FAIL fail | fail | fail |
| GSE78220 | **0.5705** [0.3179,0.7885] | 0.696 | 1.070 | 0.3384 | 0.9290 | FAIL fail | fail | fail |

判据：C1 = AUC 95% CI 下界 > 0.5；C2 = |slope − 1| ≤ 0.10（logit 尺度）；容差档 t 下 C1 放宽为下界 > 0.5 − t、C2 放宽为 |slope − 1| ≤ t。

## 5. 与新 S0 基线对比（同一冻结运输设计）

## 6. 结论：框架是否 trivial rejector？

- IMvigor210-anchored 18-gene T-cell-inflamed GEP positive control (recovered gene identity): GSE165252: AUC 0.546 [0.403,0.697], slope 0.517, main-rule FAIL; GSE78220: AUC 0.571 [0.318,0.788], slope 0.696, main-rule FAIL. Any cohort passes main rule: False; at tiers pm0.05/pm0.10/pm0.20: False/False/False.
- **判定：本数据集内无法证明框架非 trivial rejector**——阳性对照亦未在任何容差档通过。结合 within 表现（见 §3）区分两种解释：(a) 签名在这些小样本外部队列本身信号不足；(b) 框架过严。需要更大基因有效队列才能裁决。

以上全部数字为运行时实算；每患者签名分与预测概率存于 JSON（per_patient_target/per_patient_anchor）。无任何手工数值。

*冻结产出：`data/positive_control_gep18.json`；脚本 `scripts/analysis/positive_control_gep18.py`；seed 20260911；主阳性对照（6 基因 IFN-γ 签名见 data/rerun_positive_control.json，敏感性分析）。*`positive_control.json`。*
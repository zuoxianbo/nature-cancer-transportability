#  P0-2 阳性对照重跑（IMvigor210 基因身份已恢复）— IFN-γ 6 基因签名

日期：2026-09-24 ｜ seed 20260911 ｜ 脚本：scripts/analysis/rerun_positive_control.py ｜ 冻结 JSON：data/rerun_positive_control.json

**背景：** 原 P0-2 阳性对照（positive_control.py）因 IMvigor210 exmat 基因行身份被上游匿名化破坏而作废。imvigor_gene_recovery.py 已恢复身份（exmat 第 j 行 == fData 按 symbol 大写字母序排序后的第 j 个基因；性别锚定门 PASS、镜像 31085/31085 行数值恒等匹配、上游溯源逐字节一致）。本脚本按原设计重跑：**唯一改动是 IMvigor210 加载器换为恢复钥匙选行**，指标代数、判据、seed 全部原样复用（/ 逐字导入）。

## 1. 恢复钥匙交叉核验

- 脚本内重建的 symbol（大写）字母序钥匙与项目资产 `exmat_row_key.csv`（31085 行）对 6 个签名基因逐一比对：**全部一致**。

## 2. 三队列 6 基因可得性（恢复钥匙下）

| gene | IMvigor210（恢复钥匙，exmat 行号 1-based） | GSE165252（Ensembl） | GSE78220（symbol） |
|---|---|---|---|
| IDO1 | YES（行 8734） | YES（ENSG00000131203） | YES |
| CXCL9 | YES（行 4426） | YES（ENSG00000138755） | YES |
| CXCL10 | YES（行 4414） | YES（ENSG00000169245） | YES |
| CCL5 | YES（行 2980） | YES（ENSG00000271503） | YES |
| HLA-DRA | YES（行 8349） | YES（ENSG00000204287） | YES |
| STAT1 | YES（行 27397） | YES（ENSG00000115415） | YES |

进入签名的基因：IDO1, CXCL9, CXCL10, CCL5, HLA-DRA, STAT1。缺失：无。
## 3. IMvigor210 锚定（n=298 标注样本）与 within 表现

- 冻结模型：logit(P) = -1.2604 + 0.4459 × signature（ridge λ=1/n=1/298，截距不惩罚）。
- within 表观 AUC = **0.6065**；Efron 0.632 AUC = **0.6050**（n_boot=600，seed 20260911；事件 68/298）。

| 队列（within，常规队列内 z，描述性） | 表观 AUC | mean OOB | Efron 0.632 |
|---|---|---|---|
| GSE165252 | 0.5598 | 0.4805 | **0.5097** |
| GSE78220 | 0.4936 | 0.3939 | **0.4306** |

## 4. 冻结运输：逐队列指标与两判据判定

### primary_logcpm_z_frozen_all6

| target | AUC [95% BCa CI] | slope b | CITL a | Brier | log-loss | 主判据 C1∧C2 | ±0.05 | ±0.10 | ±0.20 |
|---|---|---|---|---|---|---|---|---|---|
| GSE165252 | **0.5607** [0.4220,0.7119] | 0.489 | -0.044 | 0.2336 | 0.6738 | FAIL fail | fail | fail |
| GSE78220 | **0.4936** [0.2467,0.7214] | 0.250 | 0.588 | 0.4046 | 1.1773 | FAIL fail | fail | fail |

### sensitivity_raw_z_frozen

| target | AUC [95% BCa CI] | slope b | CITL a | Brier | log-loss | 主判据 C1∧C2 | ±0.05 | ±0.10 | ±0.20 |
|---|---|---|---|---|---|---|---|---|---|
| GSE165252 | **0.5417** [0.4043,0.6918] | 0.025 | -0.701 | 0.2346 | 0.6731 | FAIL fail | fail | fail |
| GSE78220 | **0.3590** [0.1466,0.6137] | -48.670 | -74.846 | 0.3680 | 0.9966 | FAIL fail | fail | fail |

### sensitivity_logcpm_z_defined_genes

| target | AUC [95% BCa CI] | slope b | CITL a | Brier | log-loss | 主判据 C1∧C2 | ±0.05 | ±0.10 | ±0.20 |
|---|---|---|---|---|---|---|---|---|---|
| GSE165252 | **0.5607** [0.4220,0.7119] | 0.489 | -0.044 | 0.2336 | 0.6738 | FAIL fail | fail | fail |
| GSE78220 | **0.4936** [0.2467,0.7214] | 0.250 | 0.588 | 0.4046 | 1.1773 | FAIL fail | fail | fail |

判据：C1 = AUC 95% CI 下界 > 0.5；C2 = |slope − 1| ≤ 0.10（logit 尺度）；容差档 t 下 C1 放宽为下界 > 0.5 − t、C2 放宽为 |slope − 1| ≤ t。

## 5. 与新 S0 基线对比（同一冻结运输设计）

| transport | IFN-γ 阳性对照 AUC [CI]（主变体） | 新 S0 AUC [CI] | Δ AUC |
|---|---|---|---|
| IMvigor210 → GSE165252 | **0.5607** [0.4220,0.7119] | 0.5036 [0.3540,0.6585] | +0.0571 |
| IMvigor210 → GSE78220 | **0.4936** [0.2467,0.7214] | 0.6603 [0.4038,0.8617] | -0.1667 |

- within IMvigor210：IFN-γ 0.632 AUC **0.6050** vs 新 S0 0.6571。
- GSE165252：slope 0.489 vs S0 -0.527；Brier 0.2336 vs 0.2627；log-loss 0.6738 vs 0.8142。
- GSE78220：slope 0.250 vs S0 5.678；Brier 0.4046 vs 0.2280；log-loss 1.1773 vs 0.6480。

## 6. 结论：框架是否 trivial rejector？

- IMvigor210-anchored IFN-gamma 6-gene positive control (recovered gene identity): GSE165252: AUC 0.561 [0.422,0.712], slope 0.489, main-rule FAIL; GSE78220: AUC 0.494 [0.247,0.721], slope 0.250, main-rule FAIL. Any cohort passes main rule: False; at tiers pm0.05/pm0.10/pm0.20: False/False/False.
- **判定：本数据集内无法证明框架非 trivial rejector**——阳性对照亦未在任何容差档通过。结合 within 表现（见 §3）区分两种解释：(a) 签名在这些小样本外部队列本身信号不足；(b) 框架过严。需要更大基因有效队列才能裁决。

以上全部数字为运行时实算；每患者签名分与预测概率存于 JSON（per_patient_target/per_patient_anchor）。无任何手工数值。

*冻结产出：`data/rerun_positive_control.json`；脚本 `scripts/analysis/rerun_positive_control.py`；seed 20260911；取代已作废的 `positive_control.json`。*
#  消融网格诚实重规与重跑（恢复基因身份下的受控消融）

日期：2026-09-24 ｜ seed 20260911 ｜ 脚本：scripts/analysis/rerun_ablation.py ｜ 冻结 JSON：data/rerun_ablation.json

## 1. 重规说明（为什么旧  的 26 配置全部作废）

- 基因身份已恢复（_imvigor_gene_recovery）：exmat 行 j == fData 按 symbol（大写归一）字母序排序后的第 j 行；IMvigor210 加载行索引来自恢复钥匙，而非旧的 fData 位置对齐。
- 新锚定面板 = 322 权威 symbol 的跨队列交集 **320 基因**（IMvigor210 可测 320/322，缺 IGHG1, PRKN；GSE165252 322/322；GSE78220 320/322）。
- 面板内 IMvigor210 全零基因数 = **0**（旧 ：70/228），零值基因 PC1 loading 质量旧 88.7% → 新 **0.0**。前提消失，5 个 measured-only 伪影配置删除（removal_reason 见 JSON）。
- λ 规则重规（为逐位复现冻结恢复锚定值）：外部运输 λ=1/n_train=1/298；within 0.632 λ=1/N_total=1/394。 曾对两者统一用 1/N_total；ridge 层扰动则将两个 λ 同设为扰动值。
- 度量引擎逐字复用 ：外部 AUC + percentile bootstrap 95% CI（2000 次重抽外部队列、固定预测、seed 20260911）；within Efron 0.632（600 次、seed 20260911、三队列）；另按恢复脚本逐字附 BCa CI 以保持锚定连续性。

## 2. 锚定自检（baseline 逐位复现冻结恢复值）

| 量 | 冻结期望 | 本次计算 | abs diff | PASS |
|---|---|---|---|---|
| ext_GSE165252_auc | 0.503623188406 | 0.503623188406 | 0.00e+00 | **PASS** |
| ext_GSE165252_ci_bca_lo | 0.354020660429 | 0.354020660429 | 0.00e+00 | **PASS** |
| ext_GSE165252_ci_bca_hi | 0.658472153596 | 0.658472153596 | 0.00e+00 | **PASS** |
| ext_GSE78220_auc | 0.660256410256 | 0.660256410256 | 0.00e+00 | **PASS** |
| ext_GSE78220_ci_bca_lo | 0.403845496435 | 0.403845496435 | 0.00e+00 | **PASS** |
| ext_GSE78220_ci_bca_hi | 0.861658322838 | 0.861658322838 | 0.00e+00 | **PASS** |
| within_imvigor210_auc632 | 0.657074671027 | 0.657074671027 | 0.00e+00 | **PASS** |
| lam_ext | 0.003355704698 | 0.003355704698 | 0.00e+00 | **PASS** |
| lam_within | 0.002538071066 | 0.002538071066 | 0.00e+00 | **PASS** |
| **总门** | | | | **PASS（逐位一致）** |

## 3. 新网格总表（21 配置，全部真跑）

外部 CI 为主引擎 percentile bootstrap；括号内为 BCa（锚定连续性）。within 为 Efron 0.632。

| 配置 | 层 | n_genes | norm | k | λ_ext | GSE165252 AUC [pctCI] (bcaCI) | 排除0.5 | GSE78220 AUC [pctCI] (bcaCI) | 排除0.5 | within IMV 0.632 | 退化 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| baseline | baseline | 320 | clr | 2 | 0.003356 | 0.5036 [0.3519,0.6549] (0.3540,0.6585) | pct:否/bca:否 | 0.6603 [0.4206,0.8718] (0.4038,0.8617) | pct:否/bca:否 | 0.6571 | 否 |
| gene_top100 | gene_set | 100 | clr | 2 | 0.003356 | 0.5353 [0.3897,0.6794] (0.3868,0.6772) | pct:否/bca:否 | 0.5128 [0.2600,0.7598] (0.2598,0.7583) | pct:否/bca:否 | 0.4560 | 否 |
| gene_top50 | gene_set | 50 | clr | 2 | 0.003356 | 0.4656 [0.3288,0.6103] (0.3297,0.6111) | pct:否/bca:否 | 0.5962 [0.3442,0.8162] (0.3269,0.8000) | pct:否/bca:否 | 0.5547 | 否 |
| gene_rand100_s0 | gene_set | 100 | clr | 2 | 0.003356 | 0.5380 [0.3861,0.6880] (0.3874,0.6892) | pct:否/bca:否 | 0.5256 [0.2867,0.7564] (0.2756,0.7500) | pct:否/bca:否 | 0.6085 | 否 |
| gene_rand100_s1 | gene_set | 100 | clr | 2 | 0.003356 | 0.5734 [0.4191,0.7270] (0.4131,0.7253) | pct:否/bca:否 | 0.5705 [0.3269,0.7885] (0.3312,0.7922) | pct:否/bca:否 | 0.6235 | 否 |
| gene_rand100_s2 | gene_set | 100 | clr | 2 | 0.003356 | 0.4864 [0.3352,0.6401] (0.3387,0.6421) | pct:否/bca:否 | 0.6154 [0.3766,0.8401] (0.3636,0.8269) | pct:否/bca:否 | 0.6407 | 否 |
| gene_rand100_s3 | gene_set | 100 | clr | 2 | 0.003356 | 0.5489 [0.4015,0.6974] (0.4015,0.6974) | pct:否/bca:否 | 0.5897 [0.3442,0.8142] (0.3398,0.8060) | pct:否/bca:否 | 0.5802 | 否 |
| gene_rand100_s4 | gene_set | 100 | clr | 2 | 0.003356 | 0.4810 [0.3321,0.6365] (0.3322,0.6367) | pct:否/bca:否 | 0.6154 [0.3733,0.8333] (0.3643,0.8266) | pct:否/bca:否 | 0.6685 | 否 |
| gene_rand50_s0 | gene_set | 50 | clr | 2 | 0.003356 | 0.4674 [0.3260,0.6059] (0.3261,0.6060) | pct:否/bca:否 | 0.6474 [0.4117,0.8590] (0.3902,0.8442) | pct:否/bca:否 | 0.6533 | 否 |
| gene_rand50_s1 | gene_set | 50 | clr | 2 | 0.003356 | 0.4955 [0.3488,0.6447] (0.3514,0.6485) | pct:否/bca:否 | 0.6218 [0.3603,0.8462] (0.3325,0.8333) | pct:否/bca:否 | 0.5842 | 否 |
| gene_rand50_s2 | gene_set | 50 | clr | 2 | 0.003356 | 0.4176 [0.2735,0.5707] (0.2730,0.5701) | pct:否/bca:否 | 0.6538 [0.4122,0.8676] (0.3880,0.8508) | pct:否/bca:否 | 0.6167 | 否 |
| gene_rand50_s3 | gene_set | 50 | clr | 2 | 0.003356 | 0.6024 [0.4565,0.7334] (0.4559,0.7322) | pct:否/bca:否 | 0.5128 [0.2600,0.7500] (0.2540,0.7468) | pct:否/bca:否 | 0.5614 | 否 |
| gene_rand50_s4 | gene_set | 50 | clr | 2 | 0.003356 | 0.5652 [0.4226,0.7114] (0.4158,0.6996) | pct:否/bca:否 | 0.5385 [0.2941,0.7885] (0.3007,0.7924) | pct:否/bca:否 | 0.5045 | 否 |
| norm_log2raw | normalization | 320 | log2raw | 2 | 0.003356 | 0.4755 [0.3330,0.6303] (0.3277,0.6263) | pct:否/bca:否 | 0.6538 [0.4167,0.8681] (0.4067,0.8533) | pct:否/bca:否 | 0.5029 | 否 |
| norm_within_z | normalization | 320 | within_z | 2 | 0.003356 | 0.5091 [0.3539,0.6618] (0.3534,0.6614) | pct:否/bca:否 | 0.6538 [0.4132,0.8637] (0.3898,0.8534) | pct:否/bca:否 | 0.6571 | 否 |
| pc_k1 | pc_k | 320 | clr | 1 | 0.003356 | 0.5000 [0.3603,0.6489] (0.3558,0.6443) | pct:否/bca:否 | 0.6218 [0.3782,0.8667] (0.3595,0.8458) | pct:否/bca:否 | 0.5183 | 否 |
| pc_k3 | pc_k | 320 | clr | 3 | 0.003356 | 0.4846 [0.3359,0.6321] (0.3369,0.6345) | pct:否/bca:否 | 0.6795 [0.4416,0.8846] (0.4167,0.8718) | pct:否/bca:否 | 0.6509 | 否 |
| pc_k5 | pc_k | 320 | clr | 5 | 0.003356 | 0.4638 [0.3276,0.6121] (0.3275,0.6110) | pct:否/bca:否 | 0.5833 [0.3421,0.8015] (0.3301,0.7982) | pct:否/bca:否 | 0.6806 | 否 |
| ridge_lam0.01 | ridge_lambda | 320 | clr | 2 | 0.01 | 0.5036 [0.3519,0.6549] (0.3540,0.6585) | pct:否/bca:否 | 0.6603 [0.4206,0.8718] (0.4038,0.8617) | pct:否/bca:否 | 0.6571 | 否 |
| ridge_lam0.1 | ridge_lambda | 320 | clr | 2 | 0.1 | 0.5036 [0.3519,0.6549] (0.3540,0.6585) | pct:否/bca:否 | 0.6603 [0.4206,0.8718] (0.4038,0.8617) | pct:否/bca:否 | 0.6571 | 否 |
| ridge_lam1 | ridge_lambda | 320 | clr | 2 | 1 | 0.5036 [0.3519,0.6549] (0.3540,0.6585) | pct:否/bca:否 | 0.6603 [0.4206,0.8718] (0.4038,0.8617) | pct:否/bca:否 | 0.6571 | 否 |

退化配置数：0（无）；每配置子集内全零基因数均为 0（n_unmeasured_in_subset 全 0）。

## 4. 逐层外部 AUC 方差分解（对齐旧  口径）

每层含共享 baseline；统计基于非退化配置；均值为两外部队列方差平均。

| 层 | 队列 | n | min | max | range | IQR | variance |
|---|---|---|---|---|---|---|---|
| gene_set | GSE165252 | 13 | 0.4176 | 0.6024 | 0.1848 | 0.0679 | 0.002694 |
| gene_set | GSE78220 | 13 | 0.5128 | 0.6603 | 0.1474 | 0.0833 | 0.002811 |
| normalization | GSE165252 | 3 | 0.4755 | 0.5091 | 0.0335 | 0.0168 | 0.000324 |
| normalization | GSE78220 | 3 | 0.6538 | 0.6603 | 0.0064 | 0.0032 | 0.000014 |
| pc_k | GSE165252 | 4 | 0.4638 | 0.5036 | 0.0399 | 0.0215 | 0.000329 |
| pc_k | GSE78220 | 4 | 0.5833 | 0.6795 | 0.0962 | 0.0529 | 0.001818 |
| ridge_lambda | GSE165252 | 4 | 0.5036 | 0.5036 | 0.0000 | 0.0000 | 0.000000 |
| ridge_lambda | GSE78220 | 4 | 0.6603 | 0.6603 | 0.0000 | 0.0000 | 0.000000 |

层方差排序（主 4 层，均值跨两外部队列）：**gene_set > pc_k > normalization > ridge_lambda**；最大方差层 = **gene_set**（旧 ：gene_set）。

| 层 | 旧  均值方差 | 新  均值方差 | 倍数变化 |
|---|---|---|---|
| gene_set | 0.008383 | 0.002752 | 0.33x |
| normalization | 0.007544 | 0.000169 | 0.02x |
| pc_k | 0.003125 | 0.001074 | 0.34x |
| ridge_lambda | 0.000000 | 0.000000 | 两代均恰为 0 |

## 5. 恢复运输判定（任何配置 CI 有利方向排除 0.5？）

- **GSE165252**：percentile CI 下界>0.5 的配置：无；BCa CI 下界>0.5 的配置：无；（不利方向，CI 上界<0.5：percentile 无 / BCa 无）
- **GSE78220**：percentile CI 下界>0.5 的配置：无；BCa CI 下界>0.5 的配置：无；（不利方向，CI 上界<0.5：percentile 无 / BCa 无）

**判定：无任何非退化配置的外部 AUC bootstrap 95% CI 在有利方向排除 0.5（percentile 与 BCa 两法一致）——恢复身份后运输失败依旧系统性存在。**

## 6. 与旧  结论对照（存活 / 翻转 / 撤回）

### 撤回（前提已消失）

1. **“70 未测量基因伪影”结构诊断**：旧  报告 70/228 基因在 IMvigor210 全零、承载 88.7% PC1 loading 质量、top-50 |PC1 loading| 子集完全退化；恢复身份后 320 基因面板全零基因数 = 0，零基因 loading 质量 = 0.0，整条机制叙事撤回。
2. **5 个 measured-only 伪影配置**（gene_measured_only / gene_top100_measured / gene_top50_measured / norm_log2raw_measured / norm_within_z_measured）：removal_reason = "70-gene artefact retracted: zero-gene count = 0 under recovered identity"。
3. **旧 228 基因“交集”**：由假 Entrez 选择产生（Supplementary_Table_1 列 1 为行号非 Entrez），被 320 基因真实交集取代。
4. **gene_top50 退化结论**：旧退化配置 gene_top50；新网格退化配置 无——top50 现为真实可判别配置。

### 翻转（数值方向性改变）

| 量 | 旧 （错误身份） | 新 （恢复身份） | 解读 |
|---|---|---|---|
| baseline 外部 AUC GSE78220 | 0.4551 [0.2153,0.6945] | 0.6603 [0.4206,0.8718] | **方向翻转**（<0.5 → >0.5，点估计越过无效线；CI 仍含 0.5，n=25） |
| baseline within 0.632 IMvigor210 | 0.5669 | 0.6571 | 锚定内部判别大幅改善（+0.0902） |
| baseline 外部 AUC GSE165252 | 0.5136 [0.3627,0.6644] | 0.5036 [0.3519,0.6549] | 数值微降但结论不变：仍≈0.5、CI 横跨 0.5 |

### 存活（重跑后仍成立）

- “无任何配置外部 CI 排除 0.5 —— 运输失败是系统性的，而非特定预处理选择的伪影”这一  核心结论在恢复身份、撤回伪影配置后的 21 配置网格中**存活**（percentile 与 BCa 两法一致）。
- “外部 AUC 跨配置方差最大层 = gene_set”**存活**（旧 gene_set → 新 gene_set）。
- “退化配置不构成信息”结论以空集形式存活：新网格 0 个退化配置。

## 7. 汇总（计算所得 key finding）

Across 21 non-degenerate ablation configurations spanning gene set, normalization, PC dimension and ridge lambda under the recovered gene identity (320-gene panel, 0 all-zero genes), external AUC ranged 0.418-0.602 (GSE165252) and 0.513-0.679 (GSE78220). NO non-degenerate configuration yielded an external AUC whose bootstrap 95% CI excluded 0.5 under either the percentile or the BCa method -- the transport failure remains SYSTEMATIC under the recovered identity, not an artefact of any specific preprocessing choice tested. The largest cross-configuration external-AUC variance was attributable to the 'gene_set' layer (ranking: gene_set > pc_k > normalization > ridge_lambda). The  '70 unmeasured genes carrying 88.7% of PC1 loading mass' structural artefact is RETRACTED: under the recovered identity the anchor panel has 0 all-zero genes and the PC1 loading mass on (nonexistent) zero-genes is exactly 0; the 5 measured-only artefact configurations were removed (removal_reason recorded).

## 8. 复现性

- Python：managed venv（/Users/zuoxianbo/.workbuddy/binaries/python/envs/default/bin/python）；seed 20260911 全程。
- 队列样本（标注 n / 事件 n）：IMvigor210 298/68；GSE165252 71/23；GSE78220 25/13；N_total=394。
- baseline 锚定自检逐位通过后方可产出本报告（否则脚本中止，见 anchor_self_check）。
- 旧  冻结对照源：/Users/zuoxianbo/Desktop/SCI论文/Sandy_Nature/outputs/analysis_base/data/alignment_ablation.json。

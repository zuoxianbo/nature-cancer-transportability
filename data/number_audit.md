#  投稿前数字审计报告

- **审计对象**：`reviews_clean_manuscript_extract.md`（ GPT 清洁版手稿，57 行提取稿）
- **权威冻结档案**：`data/` 下  /  /  / _meta_sensitivity / _eff_meta / _tox_meta /  系列 JSON
- **审计日期**：2026-09-24
- **审计性质**：投稿门槛审计——手稿中每一个数字型断言（AUC/CI/slope/Brier/log-loss/n/占比/OR/P/基因数/配置数）逐项与冻结 JSON 比对；不复查算术。
- **位置标注**：L<n> = 提取稿行号（L5=Abstract，L13=Results 第 1 段，L15=Results 第 2 段，L17=Results 第 3 段，L19=Results 第 4 段，L21=Results 第 5 段，L32=Methods 数据集段，L36=Methods 表示段，L40=Methods 敏感性段）。
- **判定定义**：MATCH=与冻结值一致；ROUNDING_OK=冻结值按手稿精度正确舍入（或由冻结值正确导出的区间/跨度）；MISMATCH=与冻结值不一致（含错误舍入）；NOT_IN_FROZEN=冻结档案中不存在该数字的来源。

---

## 一、审计表

### Abstract（L5）

| 手稿位置 | 断言原文 | 手稿值 | 冻结来源（文件:键） | 冻结值 | 判定 |
|---|---|---|---|---|---|
| L5 | "A 322-gene baseline immune-state representation" | 322 | /external_validation.json: representation | "322-gene CLR -> PCA k=2" | MATCH |
| L5 | "cohort-specific AUCs of 0.41–0.58" | 0.41–0.58 | :F1–F4 全部 9 个 cohort-specific AUC | min 0.4103（F1→GSE78220），max 0.5796（F2→IMvigor210） | ROUNDING_OK（注 A1） |

### Results 第 1 段（L13）

| 手稿位置 | 断言原文 | 手稿值 | 冻结来源（文件:键） | 冻结值 | 判定 |
|---|---|---|---|---|---|
| L13 | "comprised 322 genes" | 322 | : representation | 322 | MATCH |
| L13 | "projected through a 228-gene cross-cohort intersection" | 228 | /alignment_ablation.json: unmeasured_gene_artifact.n_common_genes | 228 | MATCH |
| L13 | "158 genes were measurably expressed in the anchor" | 158 | : unmeasured_gene_artifact.n_measured_in_imvigor210 | 158 | MATCH |
| L13 | "whereas 70 were unmeasured" | 70 | : unmeasured_gene_artifact.n_unmeasured_in_imvigor210 | 70 | MATCH |
| L13 | "the anchor-to-target AUCs were **0.444** in GSE165252" | 0.444 | 0.444 = /competitor_benchmark.json: within_cohort_same_pipeline.per_cohort.GSE165252.efron_0632_auc（**within 值**）；真外部冻结值 = : folds.F1_train_IMvigor210.test_GSE165252.auc | within 0.4441；**外部 0.5009 [0.3456, 0.6481]** | **MISMATCH（严重，注 A2）** |
| L13 | "and **0.464** in GSE78220" | 0.464 | 0.464 = : within_cohort_same_pipeline.per_cohort.GSE78220.efron_0632_auc（**within 值**）；真外部冻结值 = : folds.F1_train_IMvigor210.test_GSE78220.auc | within 0.4640；**外部 0.4103 [0.1948, 0.6603]** | **MISMATCH（严重，注 A2）** |
| L13 | "spanned approximately 0.41–0.58" | 0.41–0.58 | : F1–F4（含  true_external_pairwise_transports 复核一致） | 0.4103–0.5796 | ROUNDING_OK（注 A1） |
| L13 | "A pooled AUC of 0.498 (95% CI 0.365–0.631)" | 0.498 [0.365, 0.631] | meta_sensitivity.json: efficacy_dlmeta.pooled_auc / ci_low / ci_high | 0.498 [0.365, 0.631] | MATCH |
| L13 | "with I²=0.75" | 0.75 | : efficacy_dlmeta.i2_qbased | 0.7516 | ROUNDING_OK（注 A3） |
| L13 | "only three cohorts contributed to that summary" | 3 | : efficacy_dlmeta（k=3，无独立键，由 cohorts 推出） | k=3 | MATCH |

### Results 第 2 段（L15）

| 手稿位置 | 断言原文 | 手稿值 | 冻结来源（文件:键） | 冻结值 | 判定 |
|---|---|---|---|---|---|
| L15 | "compressed near the source prevalence, spanning approximately 0.371–0.379" | 0.371–0.379 | /supporting_analyses.json: calibration_curves.results.IMvigor210.M2.bins[].mean_predicted | 0.3706–0.3792（5 bin 均值预测） | ROUNDING_OK |
| L15 | "the target event rate was approximately 0.228" | 0.228 | : calibration_curves.results.IMvigor210.M2.overall_observed_rate | 0.22819（68/298） | ROUNDING_OK |
| L15 | "representative estimates of 10.70" | 10.70 | : ...IMvigor210.M2.calib_slope | 10.7039 | ROUNDING_OK |
| L15 | "0.29" | 0.29 | : ...GSE165252.M2.calib_slope | 0.2909 | ROUNDING_OK |
| L15 | "and −0.86" | −0.86 | : ...GSE78220.M2.calib_slope | −0.8590 | ROUNDING_OK |
| L15 | "improved log-loss from 0.586 to 0.540" | 0.586→0.540 | : IMvigor210 M2.log_loss → M3.repeat_means.log_loss | 0.58613→0.53980 | ROUNDING_OK |
| L15 | "and Brier score from **0.198** to 0.177" | 0.198→0.177 | : IMvigor210 M2.brier → M3.repeat_means.brier | **0.19747**→0.17717 | **MISMATCH（轻微，注 A4：0.198 应为 0.197）** |
| L15 | "improved log-loss from 0.892 to 0.776 in GSE78220" | 0.892→0.776 | : GSE78220 M2.log_loss → M3.repeat_means.log_loss | 0.89237→0.77612 | ROUNDING_OK |
| L15 | "calibration-slope confidence interval spanning −26.2 to 2.9" | [−26.2, 2.9] | : GSE165252.M3.repeat_means.calib_slope_ci95 | [−26.2288, 2.9201] | MATCH（**存在于冻结档案**，注 A5） |

### Results 第 3 段（L17）

| 手稿位置 | 断言原文 | 手稿值 | 冻结来源（文件:键） | 冻结值 | 判定 |
|---|---|---|---|---|---|
| L17 | "25 non-degenerate representation and preprocessing configurations" | 25 | : per_config（共 26 配置，degenerate_configs=[gene_top50]） | 26 总数 − 1 退化 = 25 | MATCH |
| L17 | "External AUC ranged from 0.333 to 0.575 in GSE165252" | 0.333–0.575 | : per_config 非退化 external.GSE165252.auc 最小/最大 | 0.3333（gene_top100）–0.5752（gene_rand50_s1） | ROUNDING_OK |
| L17 | "and from 0.301 to 0.724 in GSE78220" | 0.301–0.724 | : per_config 非退化 external.GSE78220.auc 最小/最大 | 0.3013（pc_k3）–0.7244（gene_top50_measured） | ROUNDING_OK |
| L17 | "the 70 genes unmeasured in the anchor" | 70 | : unmeasured_gene_artifact | 70 | MATCH |
| L17 | "accounted for 88.7% of the anchor PC1 loading mass" | 88.7% | : pc1_loading_squared_mass_on_unmeasured | 0.88661 | ROUNDING_OK |
| L17 | "Restricting ... to the 158 measured genes" | 158 | : unmeasured_gene_artifact.n_measured_in_imvigor210 | 158 | MATCH |
| L17 | "increased the observed AUC in GSE78220 to 0.583–0.724" | 0.583–0.724 | : artifact_removal 层 5 配置 external.GSE78220.auc 范围 | 0.5833（norm_within_z_measured）–0.7244（gene_top50_measured） | ROUNDING_OK |
| L17 | "no confidence interval excluded 0.5" | — | : artifact_removal 各配置 GSE78220.ci_excludes_0p5 | 全部 false | MATCH |

### Results 第 4 段（L19）

| 手稿位置 | 断言原文 | 手稿值 | 冻结来源（文件:键） | 冻结值 | 判定 |
|---|---|---|---|---|---|
| L19 | "TCR diversity alone produced an AUC of 0.570" | 0.570 | /model_competition.json: models.M1_tcr_diversity.out_of_sample.auc（= /joint_stress_test.json: model_competition_extract.M1...auc） | 0.57037 | ROUNDING_OK |
| L19 | "compared with 0.544 for the combined ... model" | 0.544 | : models.M4_S0_plus_features.out_of_sample.auc | 0.54444 | ROUNDING_OK |
| L19 | "log-loss was 0.696 versus 0.767" | 0.696 / 0.767 | : M1 / M4 out_of_sample.log_loss | 0.69571 / 0.76670 | ROUNDING_OK |
| L19 | "Brier score **0.250** versus 0.271" | 0.250 / 0.271 | : M1 / M4 out_of_sample.brier | **0.24950** / 0.27106 | **MISMATCH（轻微，注 A6：0.250 应为 0.249）** |

### Results 第 5 段（L21）

| 手稿位置 | 断言原文 | 手稿值 | 冻结来源（文件:键） | 冻结值 | 判定 |
|---|---|---|---|---|---|
| L21 | "10 patients with durable clinical benefit and severe toxicity" | 10 | joint_stress_test.json: two_by_two.cells."DCB/severe" | 10 | MATCH |
| L21 | "17 with durable benefit without severe toxicity" | 17 | _joint: cells."DCB/no_severe" | 17 | MATCH |
| L21 | "5 without durable benefit with severe toxicity" | 5 | _joint: cells."NDB/severe" | 5 | MATCH |
| L21 | "15 with neither outcome" | 15 | _joint: cells."NDB/no_severe" | 15 | MATCH |
| L21 | "Fisher's exact test gave OR=1.76" | 1.76 | _joint: association_E_T.odds_ratio_DCB_vs_severe | 1.76471 | ROUNDING_OK |
| L21 | "and P=0.53" | 0.53 | _joint: association_E_T.fisher_exact_two_sided_p | 0.52933 | ROUNDING_OK |
| L21 | "log-loss was 1.345 versus 1.356 for the comparator" | 1.345 / 1.356 | _joint: joint_vs_independent.log_loss_joint / log_loss_independent | 1.34518 / 1.35597 | ROUNDING_OK |
| L21 | "the model won 12 of 25 folds" | 12/25 | _joint: n_folds_joint_better / n_folds_total | 12 / 25 | MATCH |
| L21 | "the sign-permutation P value was 0.986" | 0.986 | supporting_analyses.json: joint_formal_test.results.sign_permutation_p | 0.98626 | ROUNDING_OK |

### Methods 数据集段（L32）及队列 n

| 手稿位置 | 断言原文 | 手稿值 | 冻结来源（文件:键） | 冻结值 | 判定 |
|---|---|---|---|---|---|
| L32 | "IMvigor210 (n=298 labelled of 348)" | 298 / 348 | eff_meta.json: cohorts[0].n_labeled / n_total | 298 / 348 | MATCH |
| L32 | "GSE165252 (n=71)" | 71 | _eff_meta: cohorts[1].n_labeled | 71 | MATCH |
| L32 | "GSE78220 (n=25)" | 25 | _eff_meta: cohorts[2].n_labeled | 25 | MATCH |
| L32 | "GSE166181 (n=20)" | 20 | /gse166181.json: n_patients_baseline_T0 | 20 | MATCH |
| L32 | "GSE186143 (n=47 for efficacy; 53 for toxicity)" | 47 / 53 | /gse186143.json: n_with_response_label / n_with_toxicity_label | 47 / 53 | MATCH |
| L32 | "GSE189125 (n=13)" | 13 | _joint: toxicity_cohort_recheck.GSE189125.cohort.n_patients | 13 | MATCH |
| L32 | "GSE216329 (n=24 patients; 48 samples)" | 24 / 48 | 24 = tox_meta.json: output.per_cohort[0].n；**48 = 项目冻结档案中无**（仅 `_reftmp/gds_summ.json` GEO 系列记录 n_samples=48） | 24 冻结；48 仅 GEO 记录 | **NOT_IN_FROZEN（GEO 可佐证，注 A7）** |
| L32 | "GSE253720 (n=12)" | 12 | _joint: toxicity_cohort_recheck.GSE253720.auc_toxicity_CPI_colitis_vs_HC.n_patients | 12 | MATCH |

### Methods 其余数字（L36/L40）与 Discussion（L23）

| 手稿位置 | 断言原文 | 手稿值 | 冻结来源（文件:键） | 冻结值 | 判定 |
|---|---|---|---|---|---|
| L36 | "158 measurably expressed genes from the 228-gene intersection and 70 unmeasured" | 158/228/70 | : unmeasured_gene_artifact | 158/228/70 | MATCH |
| L36 | "the primary analysis using k=2 components" | 2 | : representation（PCA k=2）；: baseline.k | 2 | MATCH |
| L40 | "25 non-degenerate configurations varying ... dimension (k=1,2,3,5)" | 25；k∈{1,2,3,5} | : grid.pc_k（baseline(k=2), pc_k1, pc_k3, pc_k5）+ degenerate_configs | 25 非退化；k=1,2,3,5 | MATCH |
| L23 | "a frozen 322-gene representation" | 322 | : representation | 322 | MATCH |

---

## 二、汇总

| 判定 | 条数 | 明细 |
|---|---|---|
| **MATCH** | 29 | 322/228/158/70（多处）、0.498 [0.365,0.631]、k=3、slope CI [−26.2,2.9]、25 配置、no-CI-excludes-0.5、10/17/5/15、12/25、全部队列 n（除 48）、k=2、k=1,2,3,5 |
| **ROUNDING_OK** | 21 | 0.41–0.58（Abstract 与 Results 两处）、I²=0.75、0.371–0.379、0.228、10.70/0.29/−0.86、0.586→0.540、0.892→0.776、0.333–0.575、0.301–0.724、88.7%、0.583–0.724、0.570/0.544、0.696/0.767、OR 1.76、P 0.53、1.345/1.356、P 0.986 |
| **MISMATCH** | **4** | 0.444（L13，严重）、0.464（L13，严重）、0.198（L15，轻微舍入）、0.250（L19，轻微舍入） |
| **NOT_IN_FROZEN** | **1** | GSE216329 "48 samples"（L32；GEO 记录可佐证但项目冻结档案无） |

### MISMATCH 逐条修正建议

| # | 位置 | 手稿值 | 冻结依据 | 修正建议值 |
|---|---|---|---|---|
| 1 | L13 | 0.444（GSE165252 anchor-to-target） |  F1 真外部 = 0.5009 [0.3456, 0.6481]；0.444 实为  within Efron 0.632 值 0.4441 | **0.501（95% CI 0.346–0.648）**；若采用  消融基线管线则为 0.514 [0.363, 0.664]（须与管线选择一致，见注 A2） |
| 2 | L13 | 0.464（GSE78220 anchor-to-target） |  F1 真外部 = 0.4103 [0.1948, 0.6603]；0.464 实为  within Efron 0.632 值 0.4640 | **0.410（95% CI 0.195–0.660）**；若采用  消融基线管线则为 0.455 [0.215, 0.694] |
| 3 | L15 | Brier 0.198 |  冻结 0.19747（正确三位舍入为 0.197；0.198 系 0.1975→0.198 的双重舍入产物） | **0.197** |
| 4 | L19 | Brier 0.250（TCR） |  冻结 M1 brier 0.24950（正确舍入 0.249；注意  中 M2_activated_cd4_memory 的 brier=0.2502，疑为错行转录或双重舍入） | **0.249** |

### NOT_IN_FROZEN 列表

| # | 位置 | 数字 | 状态 | 处置建议 |
|---|---|---|---|---|
| 1 | L32 | GSE216329 "48 samples" | 项目冻结 JSON（/_tox_meta/）只记录 n=24 患者；48 仅见于 `_reftmp/gds_summ.json` 的 GEO 系列记录（n_samples=48），非项目冻结档案 | 数字本身与 GEO 一致，**无需改数**；建议在分析 manifest 或 Data Availability 中固化 GEO n_samples 出处，或将该冻结参考记录纳入正式档案 |

### 注释

- **A1（0.41–0.58 的口径）**：该跨度 =  F1–F4 全部 9 个 cohort-specific AUC 的 min/max（0.4103–0.5796）。上限 0.58 来自源互换敏感性转运（GSE165252→IMvigor210，0.5796），**不是**主估计量（IMvigor210 冻结→外部）的值；主外部（F1）仅 0.410/0.501。Results 正文已用 "broader pre-specified cross-cohort sensitivity transports" 限定，与冻结值可对上；Abstract 的 "cohort-specific AUCs of 0.41–0.58" 建议保留但确保与 Results 的限定语呼应（若 Abstract 只指主外部值，则应为 0.41–0.50）。
- **A2（已知硬伤，已确认）**：手稿 Results 第 1 段 "anchor-to-target AUCs were 0.444 in GSE165252 and 0.464 in GSE78220" 中的 0.444/0.464 与  `within_cohort_same_pipeline` 的 GSE165252/GSE78220 Efron 0.632 within 值（0.4441/0.4640）逐位吻合，属 "within 值被当外部值" 的错位。按手稿自己定义的主估计量（S0 在 IMvigor210 冻结、目标队列不重估），对应冻结值是  F1 的 **0.501 / 0.410**。注意修正后方向一升一降（GSE165252 由 0.444→0.501 升过 0.5；GSE78220 由 0.464→0.410 降），"external discrimination was not stably reproduced" 的结论不变（两 CI 均含 0.5）。若作者意在引用  消融基线（锚 PCA 旋转 + lambda=1/N_total 管线），则冻结值为 0.514/0.455——但那不是  语义下的主外部验证，二者不可混用。
- **A3（I²）**： `efficacy_dlmeta.i2_qbased` = 0.7516（Q-based，主）；_eff_meta 中另有 tau-based i2 = 0.8548（ 已明示 0.85 是 tau²-型 I² 且非首选）。手稿 0.75 与  主值一致，正确。
- **A4（0.198）**：/ 冻结 IMvigor210 M2 Brier = 0.19747，三位舍入应为 0.197。0.198 只能由先舍到 0.1975 再舍到 0.198 的双重舍入产生。属轻微但必须改（Nature Cancer 级别的数字一致性要求）。
- **A5（−26.2 到 2.9）**：任务书怀疑此 CI 可能无冻结来源——**已确认存在**：supporting_analyses.json GSE165252 M3 `calib_slope_ci95` = [−26.2288, 2.9201]，且语境（recalibration 不稳定）与手稿用法一致。
- **A6（0.250）**： 冻结 M1（TCR diversity）Brier = 0.24950，正确舍入 0.249。同表中 M2_activated_cd4_memory 的 Brier 恰为 0.2502——0.250 疑为错行转录（取了 M2 行的值）或双重舍入；无论哪种都必须改为 0.249。
- **A7（48 samples）**：GSE216329 患者 n=24 在 _tox_meta 冻结；"48 samples" 仅能溯源到 GEO 系列记录（`_reftmp/gds_summ.json` → n_samples=48）。数字与 GEO 一致，无需改值，但应补一个正式可引用的出处。

---

## 三、修正清单（定稿时逐项执行）

1. **【必须改·严重】L13（Results 第 1 段）**："the anchor-to-target AUCs were 0.444 in GSE165252 and 0.464 in GSE78220"
   → 改为 "the anchor-to-target AUCs were **0.501 (95% CI 0.346–0.648)** in GSE165252 and **0.410 (95% CI 0.195–0.660)** in GSE78220"（ F1 冻结真外部值）。
   备选（如整段改挂  消融基线管线并明示管线差异）：0.514 (0.363–0.664) / 0.455 (0.215–0.694)。二选一，不得保留 0.444/0.464。
2. **【必须改·轻微】L15**："Brier score from 0.198 to 0.177" → "Brier score from **0.197** to 0.177"。
3. **【必须改·轻微】L19**："Brier score 0.250 versus 0.271" → "Brier score **0.249** versus 0.271"。
4. **【建议核对口径】L5 Abstract / L13**："0.41–0.58" 依赖含源互换敏感性的 F1–F4 跨度（冻结 min 0.4103 / max 0.5796，舍入正确）。保留可以，但确保 Abstract 与 Results 的 "broader pre-specified ... sensitivity transports" 限定语一致；若 Abstract 意图只报主外部值，应改为 0.41–0.50。
5. **【建议补出处】L32**：GSE216329 "48 samples" 在分析 manifest / Data Availability 中固化 GEO n_samples 出处（或纳入冻结参考档案）。数字本身无需改。
6. **【无需改·已确认存在】L15**：GSE165252 校准 slope CI −26.2 到 2.9（ 冻结 [−26.2288, 2.9201]）——保留。

### 审计覆盖与范围说明

- 已逐段扫描全文 57 行提取稿的全部数字型断言（含 Abstract、Results 5 段、Discussion、Methods、数据/代码可用性声明）；参考文献中的卷期页码属书目信息，不在数据审计范围。
- 所有 MATCH/ROUNDING_OK 条目的冻结键已在审计表列明，可逐键复核。
- 未发现除上述 4 处 MISMATCH 与 1 处 NOT_IN_FROZEN 之外的数字问题；特别地，pooled 0.498 [0.365, 0.631] I²=0.75、校准压缩 0.371–0.379 vs 0.228、slopes 10.70/0.29/−0.86、recalibration 0.586→0.540 与 0.892→0.776、25 配置范围 0.333–0.575 / 0.301–0.724、158 measured → 0.583–0.724、88.7% PC1 loading、GSE186143 10/17/5/15、OR 1.76 / P 0.53、1.345 vs 1.356、12/25 folds、P=0.986、TCR 0.570 vs 0.544（log-loss 0.696/0.767、Brier 0.271）及全部队列 n 均与冻结档案一致（或在正确舍入范围内）。

# — Clinical-only 竞争基线（P0-4 投稿门槛）

**脚本**: `scripts/analysis/clinical_only.py` → **数据**: `data/clinical_only.json`
**seed** 20260911；指标代数（`fit_logistic` / `metrics_block` / 2000-rep bootstrap BCa AUC CI）自 `competitor_benchmark.py` **逐字导入**，与  完全一致。
**患者对齐**：三队列标签向量在运行时重建并**断言等于**冻结的 `_eff_{IMvigor210,GSE165252,GSE78220}_y.npy`（298 / 71 / 25，事件 68 / 23 / 13）——与冻结 S0 /  基准**完全同一批 held-out 患者**。

---

## 1. 变量可得性矩阵（基于源文件实际解析，非假设）

| 变量 | IMvigor210 (n=298) | GSE165252 (n=71) | GSE78220 (n=25) | 可调和性 |
|---|---|---|---|---|
| Sex | ✓ pData `Sex`（M/F，0 缺失） | ✗ GEO 未存 | ✓ `gender`（25/25） | **IMvigor+GSE78220**（同构念同编码） |
| 患者年龄 | ✗ **pData 无患者年龄**（`Sample.age` 是采样时间分类，非年龄） | ✗ | ✓ `age (yrs)` 数值 | **否**——训练队列即缺失，无法进入任何训练模型 |
| 分期 | ✓ `TCGA.Subtype` I–IV（0 缺失） | ✗ | ✓ AJCC M 类别 M0/M1a/b/c | **仅粗糙二值**（构念不同：TNM stage group vs M category，已声明为 pre-declared crude 调和） |
| ECOG PS | ✓ 0/1/2（0 缺失） | ✗ | ✗ | 否（单队列） |
| TMB | ✓ FMOne（**64/298 缺失**） | ✗ | ✗（仅个别基因突变，无 TMB） | 否——IMvigor-only 敏感性 |
| 转移状态 | ✓ `Met.Disease.Status`（27/298 缺失） | ✗ | 部分（并入 M 类别，构念不同） | 否——IMvigor-only 敏感性 |
| 既往治疗 | ✓ platinum / BCG | ✗ | ✓ prior MAPKi | 否（药物/瘤种不同，非共同构念） |
| 吸烟 | ✓ | ✗ | ✗ | 否（单队列） |

**结论**：
- **跨三队列共同患者协变量集 = 空集**（如实记录，P0-4 步骤 4 触发）。
- **最大可调和集（IMvigor210 训练 → GSE78220 运输）= {sex, 进展期/转移期二值}**；GSE165252 的 GEO 存档**不含任何患者级协变量**（特征行仅 tissue / treatment / response——已实际解析验证）。
- 患者年龄无法进入任何训练模型（训练队列 IMvigor210 无此变量）。

**调和规则**：`sex_male`（M=1）；`stage_adv`：IMvigor210 `TCGA.Subtype=="IV"`→1，否则 0；GSE78220 M1a/b/c→1，M0→0。C_main/C_sex 的特征**零缺失**（已验证）；C_rich 的缺失处理见 §4。

## 2. Clinical-only 模型（IMvigor210 n=298 训练，λ=1/298，截距不惩罚，冻结）

| 模型 | 特征 | 冻结系数 |
|---|---|---|
| C_main | sex_male, stage_adv | 截距 −1.619；sex 0.463；stage 0.148 |
| C_sex | sex_male | 截距 −1.591；sex 0.463 |

### 外部指标（与  同一批 held-out 患者、同一指标集）

| 队列 | 模型 | AUC [95% BCa CI] | logit slope | calib-in-large (a) | Brier | log-loss |
|---|---|---|---|---|---|---|
| GSE78220 (n=25, 13 事件) | C_main | 0.487 [0.305, 0.698] | 0.011 | +0.093 | 0.3292 | 0.8805 |
| GSE78220 | C_sex | 0.513 [0.325, 0.700] | 0.254 | +0.404 | 0.3396 | 0.9105 |
| GSE165252 (n=71, 23 事件) | C_main | 0.500（常数，退化） | 不可估 | 0.483 | 0.2282 | 0.6538 |
| GSE165252 | C_sex | 0.500（常数，退化） | 不可估 | 0.483 | 0.2282 | 0.6538 |

**GSE165252 退化说明（数据约束，非建模偷懒）**：该队列 GEO 元数据无任何患者协变量，冻结临床模型无法产生个体化预测，退化为协变量边缘常数预测（=训练集平均预测概率 0.2282 = 训练患病率；log-loss 与 M0 prevalence-only 完全一致，AUC 按构造为 0.5，slope 不可估）。

### 双准则裁定（ 规则，同口径）

临床模型在**全部 2 个目标队列均未获运输许可**（C_main 与 C_sex：AUC CI 下界 ≤0.5 且 |slope−1|>0.1），与 M0–M6 全体竞争者的  裁定一致。

## 3. 与 M2（冻结 S0）/ M0 对比（同为 IMvigor210 训练、同一批患者；M2/M0 数值取自  `true_external_pairwise_transports`）

| 目标队列 | C_main AUC | M2 (S0) AUC | ΔAUC (C−M2) | C log-loss | M2 log-loss | M0 log-loss | C Brier | M2 Brier | M0 Brier |
|---|---|---|---|---|---|---|---|---|---|
| GSE165252 | 0.500 | 0.501 [0.346, 0.648]（slope 0.162） | −0.001 | 0.6538 | 0.6717 | **0.6538** | 0.2282 | 0.2320 | **0.2282** |
| GSE78220 | 0.487 [0.305, 0.698] | 0.410 [0.195, 0.660]（slope −0.517） | +0.077 | **0.8805** | 0.9564 | 0.8927 | **0.3292** | 0.3471 | 0.3348 |

**对比结论**：
1. **无任何模型在任一外部队列通过双准则**。临床基线与冻结 S0 同样无法运输。
2. GSE78220 上 C_main 的 AUC（0.487）与概率精度（log-loss 0.881 / Brier 0.329）**略优于** M2（0.410 / 0.956 / 0.347）且略优于 M0（0.893 / 0.335），但差异远小于该队列最小可检测 AUC（：0.768，80% 检验效能）——不构成可检出的优势，CI 跨 0.5。
3. GSE165252 上一切临床模型退化 ≡ prevalence-only；M2 相对 M0/C 在 log-loss/Brier 上更差（ 已有结论）。
4. **总体结论（P0-4）**：建立真正的 clinical-only 基线后，"prevalence-only 是唯一简单基线"的  借口被推翻，但推翻的结论方向不变——**可调和临床变量（sex、粗糙分期）在外部队列同样不携带超越患病率的判别信息**；S0 与 clinical-only 谁也不具备运输资格。

## 4. IMvigor210-only 敏感性（C_rich，不运输）

全本地协变量模型（sex、分期、ECOG、log10(TMB+1)（中位数填补 + 缺失指示，64/298 缺失）、转移状态 one-hot + 缺失指示（27/298 缺失）、platinum、BCG、吸烟），同一 ridge 代数，within-cohort Efron 0.632（600 次 bootstrap 重拟合，seed 20260911）：

| 模型 | apparent AUC | OOB 均值 | **0.632 AUC** |
|---|---|---|---|
| C_rich（14 特征） | 0.839 | 0.766 | **0.793** |
| C_main | 0.546 | 0.507 | 0.521 |
| C_sex | 0.537 | 0.522 | 0.528 |
| 参考：S0（ within-cohort 同管线） | 0.584 | 0.557 | 0.567 |

主导系数：log10 TMB +3.60、ECOG1 −1.17、met_liver −1.12（正向 TMB 效应与文献一致）。**解读**：IMvigor210 队列内，富临床模型（0.793）大幅超过 S0（0.567）与可调和临床模型（0.521）——但其信息载体（TMB、ECOG、转移状态）在两个 GEO 队列均不可得，**富临床模型不可运输**。这把  的结论精确化：不是"临床变量无用"，而是"**可运输的变量集（sex、粗糙分期）无用，有用的变量（TMB 等）不可运输**"。

## 5. Case-mix 描述性对比（P0-9 descriptive 部分）

| 共同变量 | IMvigor210 (n=298) | GSE165252 (n=71) | GSE78220 (n=25) |
|---|---|---|---|
| 应答率 | 0.228 | 0.324 | 0.520 |
| 男性比例 | 78.2% | 不可得 | 68.0% |
| 患者年龄 | 不可得（无此变量） | 不可得 | 均值 60.3（SD 12.2，27–84） |
| 分期分布 | I:107 / II:75 / III:60 / IV:56 | 不可得 | M1c:20 / M1b:2 / M1a:2 / M0:1 |
| 进展/转移比例 | 18.8%（stage IV） | 不可得 | 96.0%（M1） |

**Case-mix shift 要点**：
- **患病率漂移**：0.228 → 0.324 → 0.520（GSE78220 应答率是 IMvigor210 的 2.3 倍）；IMvigor210 训练的临床模型在 GSE78220 的预测均值 ≈0.23 而实际应答率 0.52，系统性低估应答（重校准 slope≈0.01、intercept≈logit(0.52)，即重校准后退化为患病率常数——模型信号在目标队列消失）。
- **分期构成极端漂移**：进展期比例 SMD = **−2.50**（18.8% vs 96.0%）——IMvigor210 混杂大量 non-stage-IV，GSE78220 几乎全为 M1；两队列的"分期"不可直接比较，印证 §1 的 crude 调和警示。
- 性别构成中度漂移：SMD = 0.23（78.2% vs 68.0% 男性）。
- GSE165252 在所有协变量维度上**不可比**（无可得变量），仅患病率可量化。

## 6. Limitation 语句草稿（供手稿引用）

> Clinical-only baseline limitation: the patient-level covariates deposited with the three efficacy cohorts are not mutually harmonizable. GSE165252 (PERFECT) deposits no patient-level covariate of any kind (only tissue, treatment arm and response), the IMvigor210 pData contains no patient age (its "Sample.age" column is a sample-timing category), and the staging constructs differ (IMvigor210 TCGA stage group vs GSE78220 AJCC M category). The covariate set common to all three cohorts is therefore empty, and the maximal set harmonizable between the training cohort (IMvigor210) and any external cohort is {sex, a crude advanced/metastatic-stage binary}, transportable only to GSE78220. On GSE165252 every clinical-only model necessarily degenerates to a constant (covariate-marginal) prediction with AUC exactly 0.5, so for that cohort prevalence-only is the strongest clinical baseline the deposited data permit. These constraints are properties of the deposited data, not modelling choices; a richer clinical baseline would require patient-level covariate tables from the trial sponsors/publications that are not part of the GEO deposits.

## 7. 复现

- 脚本：`outputs/analysis_base/scripts/analysis/clinical_only.py`
- 数据：`outputs/analysis_base/data/clinical_only.json`（可得性矩阵、调和规则、逐队列指标、0.632、 对比、case-mix、limitation 全量入库）
- Python：`/Users/zuoxianbo/.workbuddy/binaries/python/envs/default/bin/python`；seed 20260911；bootstrap 2000（AUC BCa CI）/ 600（Efron 0.632）
- 标签重建断言：三队列重建 y 与冻结 `_eff_*_y.npy` 逐元素相等（脚本内 `assert`，运行通过）

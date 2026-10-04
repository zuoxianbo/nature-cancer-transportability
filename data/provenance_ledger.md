#  Provenance Ledger — P0-7 cohort/source-lineage audit

**日期**: 2026-09-24 · **范围**: GSE186143, GSE166181, GSE216329, GSE253720 ·
**机器可读版**: `data/provenance_ledger.json` · **冻结证据**: `data/evidence/`（六个 family.soft.gz 官方记录）

**核验通道说明（如实记录）**：审计时 NCBI 主站（www.ncbi.nlm.nih.gov 的 acc.cgi 页面与 eutils）被 NCBI 反滥用机制拦截（reCAPTCHA / Access Denied），**未标注为 UNVERIFIED_WEB**，因为官方 GEO 记录改由 `ftp.ncbi.nlm.nih.gov` 的 family.soft.gz 全文取得（四个目标 series + GSE186144 超系列 + GSE189125 子系列全部 200 OK，内容等同 acc.cgi 页面的权威记录）；关联文献经 Europe PMC REST API（元数据 + 全文 XML）与出版商记录核验。唯一未能取得的全文是 Lozano Nat Med 2022 的正文（Europe PMC fullTextXML 对 PMC8866214 返回 500，出版商站点被浏览器验证拦截），相关两处字段标 UNVERIFIED。

**判定统计**：CONFIRMED 27 项 · CORRECTED 10 项 · UNVERIFIED 3 项（含 1 项 "not claimed" 确认项计入 CONFIRMED）。

---

## 1. GSE186143（重中之重）

**GEO 记录**（`evidence/GSE186143_family.soft.gz`）：
- Title: *"T cell characteristics associated with toxicity to immune checkpoint blockade in patients with melanoma [bulk RNA-seq]"*
- Overall design（原文）: *"bulk RNA sequencing (bulk RNA-seq) of peripheral blood mononuclear cell (PBMC) samples obtained at baseline from 61 patients with metastatic melanoma and five healthy controls. The 61 patient samples consisted of two independently profiled bulk cohorts (bulk cohort 1, n = 26; bulk cohort 2, n = 27) and eight patient samples used for technical experiments."*
- 原始数据：*"raw data will be submitted to dbGaP due to patient privacy concerns"*
- 平台：GPL16791（Illumina HiSeq 2500），TruSeq RNA Exome，2×150 bp
- 公开处理矩阵：TPM（`GSE186143_bulk_TPM_samples.txt.gz`，56,809 基因 × 66 样本）
- 结构：GSE186144 超系列 = GSE186143（bulk）+ GSE189125（scRNA）+ GSE189126（TCR/BCR）；BioProject PRJNA772729
- **关联文献**：Lozano AX, Chaudhuri AA, Nene A, et al. *Nat Med* 2022;28(2). doi:10.1038/s41591-021-01623-z；PMID 35027754；PMC8866214

**逐项判定**：

| 字段 | 包内声称 | 核实值 | 判定 |
|---|---|---|---|
| **癌种** | melanoma（Table S15）； 曾记 "cancer type not locally verifiable" | **转移性黑色素瘤**——GEO 标题/摘要/逐样本 `diagnosis: Melanoma`（61/61）+ Nat Med 摘要。** 遗留问题就此了结：答案是 melanoma，证据链完整。** | **CONFIRMED** |
| 组织 | PBMC/血液 | 外周血 PBMC（66/66 样本 `tissue: Blood`） | CONFIRMED |
| **患者数 n** | Table S15: 47；README: 53； 联合分析: 47 | **系列共 61 例黑色素瘤患者 + 5 例健康对照**（bulk cohort 1 = 26 + bulk cohort 2 = 27 + 8 个技术实验样本）；本地临床注释恰覆盖 53（=26+27）；47 = 53 − 6 NE。三个数字都真实但指代不同集合，手稿必须写明链条 | **CORRECTED** |
| 表达矩阵样本 | 未声明 | 66 列 = 61 患者 + 5 健康对照（H-1~H-5）；分析须显式剔除 HC 与技术样本列 | CORRECTED |
| **疗效终点** | response DCB/NDB/NE（27/20/6） | **GEO 元数据（含超系列）无任何 response 字段**；DCB/NDB 标签源自 **Nat Med 2022 Supplementary Table 3**（同源注释文件 `GSE189125_clinical_annotation.json` 明示 source）。手稿引用 E 标签必须引论文补充材料而非 GEO | **CORRECTED**（来源归属） |
| **毒性终点** | severe_irAE（二分类），grade 0–4 | GEO 逐样本 `highest irae grade`；severe = grade ≥ 3；本地 JSON 与 GEO **逐人零不一致**（53/53 匹配；severe 17/53 = grade 3×14 + grade 4×3） | **CONFIRMED** |
| **ICI 方案** | "ICI" | **anti-PD-1 单药（29/53）或 anti-PD-1 + anti-CTLA-4 联合（24/53）**（GEO `therapy type` + Nat Med 摘要） | **CORRECTED** |
| 基线采样 | 治疗前基线 | 基线 PBMC（GEO overall design "obtained at baseline"；Nat Med "pretreatment blood"） | CONFIRMED |
| 治疗线数 | 未声称 | GEO 无此字段；论文正文被拦截未能核验 | UNVERIFIED |
| 模态 | bulk(TPM) | PBMC bulk RNA-seq；公开矩阵为 TPM；原始 counts 在 dbGaP 受控 | CONFIRMED |
| 平台 | 未声明 | GPL16791 / HiSeq 2500 / TruSeq RNA Exome | CONFIRMED（新增细节） |
| 本地临床 JSON 一致性 | 53 例 ground truth | 53/53 与 GEO 匹配；DCB 27/NDB 20/NE 6；severe 17； 2×2（n=47）= 10/17/5/15（DCB&sev / DCB&非sev / NDB&sev / NDB&非sev），NE 含 2 severe + 4 非 severe——与冻结数字完全一致 | CONFIRMED |
| 论文-队列映射 | 未涉及 | Nat Med 摘要称 3 队列 n=27/26/18，与 GEO bulk 26/27 的对应（尤其 n=18）未能从公开元数据完全对齐（论文全文未取得） | UNVERIFIED |

---

## 2. GSE166181

**GEO 记录**：Title *"scRNA-Seq of peripheral CD8 T cells from patients with metastatic melanoma taking ICI therapy."* · Public on Feb 05 2021 · 平台 GPL21697（NextSeq 550 + Bio-Rad ddSEQ/SureCell WTA 3'）·
**关联文献**：De Biasi S, Gibellini L, Lo Tartaro D, et al. *Nat Commun* 2021;12:1669. doi:10.1038/s41467-021-21928-4；PMID 33723257；PMC7961017（全文已核验）。

| 字段 | 包内声称 | 核实值 | 判定 |
|---|---|---|---|
| 癌种 | melanoma | 转移性黑色素瘤 | CONFIRMED |
| 细胞/组织 | scRNA、外周 | 外周血 FACS 分选 CD3⁺CD8⁺ T 细胞（纯度 96–99%），ddSEQ，~4 测序重复/样本 | CONFIRMED |
| **n** | 20 | **20 例**（11 R / 9 NR；本地元数据重解析 20 名 Mela* 受试者、患者层面 11/9，与 GEO 一致）。注意：**论文总队列 28 例（17 R/11 NR，含流式）**，GEO 为 scRNA 子集；本地矩阵 51,701 细胞 vs 论文 QC 后 56,142（过滤不同，引用时以 GEO 存档数为准） | CONFIRMED（附注） |
| 样本数 | 未声明 | 183 个 GSM 文库 = 20 患者 × 3 时点 × ~3–4 重复 | CONFIRMED |
| **疗效终点** | efficacy | **RECIST、6 个月随访**：R = CR/PR/SD/混合反应 >6 个月无进展；NR = 疾病进展（PD）——论文 Methods 原文 | CONFIRMED（新增细节） |
| 毒性终点 | 无 | 无 irAE 注释 | CONFIRMED |
| **ICI 方案** | "ICI" | **nivolumab 或 pembrolizumab（anti-PD-1 单药；标准治疗）** | **CORRECTED** |
| 基线采样 | 基线使用 | T0 = 治疗前，T1 = 第 1 周期后，T2 = 第 2 周期后；基线分析须取 T0 子集 | CONFIRMED |
| 模态 | scRNA | scRNA（分选 CD8 T）；UMI 计数矩阵（原始+标准化） | CONFIRMED |
| **元数据列错位** | 未涉及 | **GEO 存档 Metadata.tsv 表头与数据错位一列**（header 的 `subject` 列实为 T0/T1/T2；真实受试者 ID 在 `sample` 列/replicate 前缀）。按表头直读会静默错赋受试者 | **CORRECTED**（预警） |

---

## 3. GSE216329

**GEO 记录**：Title *"Single cell RNA sequencing reveals distinct T cell populations in immune-related adverse events of checkpoint inhibitors"* · Public on Dec 14 2022 · 平台 GPL24676（NovaSeq 6000；10x Chromium V(D)J  5'；Cell Ranger .1.1）·
**关联文献**：Bukhari S, Henick BS, Winchester RJ, et al. *Cell Rep Med* 2023;4(1):100868. doi:10.1016/j.xcrm.2022.100868；PMID 36513074；PMC9873824。

| 字段 | 包内声称 | 核实值 | 判定 |
|---|---|---|---|
| **癌种** | mixed solid tumours | **以 NSCLC 为主：18/24 为 NSCLC**（出版物 + 独立综述 PMC11232920："Approximately 2 thirds of the patients had a diagnosis of NSCLC (18/24)"）；GEO 仅写 "cancer patients" | **CORRECTED** |
| 细胞/组织 | PBMC scRNA | 外周血 RosetteSep 富集未触碰 CD3 T 细胞（非分选） | CONFIRMED |
| **n** | 24（患者）/48（样本） | **24 患者 × 2 时点 = 48 样本**；本地文件清点完全一致（P1–P24 × T/TN）；**后缀语义：_TN = 基线（Baseline），_T = 治疗中（On treatment）** | CONFIRMED |
| **毒性终点** | toxicity（12 维 S0 AUC 0.804 的二分类标签） | **irAE vs 无 irAE（15/24 患者，12 个月随访）**；器官特异性（肺炎/关节炎/甲状腺炎/肾炎/结肠炎/神经毒性）；**分级为 grade 2–3（中-重度），非 grade ≥ 3 "severe irAE"**。不得与 GSE186143 的 grade ≥ 3 队列无说明地并置 | **CORRECTED** |
| 疗效终点 | 无 | 无 response 注释 | CONFIRMED |
| ICI 方案 | anti-PD-1 | anti-PD-1（22/24 单药，出版物） | CONFIRMED |
| 基线采样 | 基线使用 | 每患者一份治疗前基线（_TN）+ 一份治疗中（_T）；基线分析须取 24 个 _TN | CONFIRMED |
| **模态** | **bulk（Table S15 "Type: bulk"）** | **scRNA（10x Chromium V(D)J ）**；本地即 48 套 matrix.mtx/barcodes/genes 三元组。除非明确说明为 pseudo-bulk 聚合，"bulk" 标注错误 | **CORRECTED** |
| 平台 | 未声明 | GPL24676 / NovaSeq 6000 / Chromium V(D)J  / Cell Ranger .1.1 | CONFIRMED |

---

## 4. GSE253720

**GEO 记录**：Title *"Reverse translational discovery reveals dysregulation of CD4+ and CD8+ resident memory T, myeloid, and stromal cells in steroid-refractory, checkpoint inhibitor colitis [biopsy scRNA_TCR-seq]"* · Public on Mar 20 2024 · 平台 GPL24676（NovaSeq 6000；Chromium 5' ；6 个 pooled 文库经 demuxlet/freemuxlet 去卷积）·
**结构**：本 series 为活检 scRNA 子系列；出版物数据声明将全部数据存于超系列 **GSE253723** ·
**关联文献**：He JY, Kim YJ, Mennillo E, et al. *J Immunother Cancer* 2024;12(4):e008628. doi:10.1136/jitc-2023-008628；PMID 38642938；PMC11033653（全文已核验）。

| 字段 | 包内声称 | 核实值 | 判定 |
|---|---|---|---|
| 癌种 | ICI-colitis（混合原发瘤） | CPI-colitis 患者原发瘤混杂（h5ad `cancer`）：黑色素瘤 3、RCC 2、肺鳞癌 1、前列腺癌 1、HNSCC 1、Merkel 细胞癌 1 | CONFIRMED |
| 组织 | biopsy | 结肠内镜冷活检（每人 6 块）；配对血液 CITE-seq 在超系列中、不在本 h5ad | CONFIRMED |
| **n** | **12** | **活检数据集 14 人：9 CPI-colitis（HS7–HS15）+ 3 健康对照（HS1–HS3）+ 2 溃疡性结肠炎（HS4–HS5）**；出版物共入组 15 人（1 名 UC 的 HS6 活检无活细胞、仅入血液数据）。**12 与任何公开计数都不符**（最接近 9+3 HC，疑似把病例-对照比较误记为 n）。Table S15 须改为 9/3/2（h5ad 内 14 人） | **CORRECTED** |
| **毒性终点** | toxicity； 已标健康对照混杂 | 类固醇经治/难治性 CPI-colitis（**组织学+内镜确诊的病例定义**），对照组为 HC 与 UC——**是活动性结肠炎生物学的病例-对照研究，不是治疗前预测性生物标志物设计**。健康对照混杂确认并量化：3 HC + 2 UC 共存于同一 h5ad；活检时免疫抑制：激素 3 人、TNF 阻断 3 人、无 3 人 | **CORRECTED**（设计定性） |
| ICI 方案 | anti-PD-1 / combo | 6/9 anti-PD-1 单药 + 3/9 anti-PD-1 + anti-CTLA-4 联合（h5ad `cpi` + 论文原文） | CONFIRMED |
| **基线采样** | 未明确（毒性臂语境隐含基线） | **非治疗前基线**：活检取自活动性、以类固醇经治为主的 CPI-colitis 时点（8/9 类固醇经治，多人已用 TNF 阻断/vedolizumab）。**不得将该队列为"基线预测"证据；如引用须明标为 on-toxicity/结局组织证据** | **CORRECTED**（关键） |
| 模态 | scRNA(h5ad, biopsy) | 多重化活检 scRNA + TCR-seq（50,845 细胞；h5ad 另含 TCR 克隆型与 ADT 字段；CITE-seq/CyTOF 在超系列） | CONFIRMED |
| 平台 | 未声明 | GPL24676 / NovaSeq 6000 / Chromium 5'  / 6 pooled 文库 | CONFIRMED |
| 健康对照混杂 |  已标 | 确认（见毒性终点行） | CONFIRMED |

---

## 5. 修正清单（包内需要改动的地方）

1. **Table S15 GSE186143 行 n**：47 → 写明完整链条："61 patients in series (26 bulk-cohort-1 + 27 bulk-cohort-2 + 8 technical) + 5 healthy controls; 53 with clinical annotation; 47 evaluable (excl. 6 NE)"。
2. **Table S15 GSE186143 Agent**："ICI" → "anti-PD-1 (29/53) 或 anti-PD-1+anti-CTLA-4 (24/53)"。
3. **GSE186143 E 标签来源**：所有引用 DCB/NDB 的地方注明"labels from Nat Med 2022 Supplementary Table 3（不在 GEO 元数据中）"。
4. **Table S15 GSE166181 Agent**："ICI" → "nivolumab/pembrolizumab (anti-PD-1 mono)"。
5. **GSE166181 元数据预警**：任何复用 Metadata.tsv 的代码须注明表头错位（subject 列实为时点）。
6. **Table S15 GSE216329 Cancer**："mixed solid tumours" → "predominantly NSCLC (18/24)"。
7. **Table S15 GSE216329 Type**："bulk" → "scRNA (10x; pseudo-bulk aggregation if used must be stated)"。
8. **GSE216329 毒性终点描述**：正文/S3 表述为 "any irAE, grade 2–3, organ-specific, 12-month follow-up, 15/24"——非 severe (grade ≥ 3) irAE；与 GSE186143 的 grade ≥ 3 终点并置时须注明定义差异。
9. **Table S15 GSE253720 n**：12 → "9 CPI-colitis + 3 HC + 2 UC (14 donors in biopsy h5ad)"。
10. **GSE253720 采样时点**：明确"活动性类固醇经治/难治 CPI-colitis 时点的结肠活检（非治疗前基线）；病例-对照设计（vs HC/UC）"—— 的健康对照混杂警示升格为设计层面限定。

---

## 6. Supplementary Table 1（cohort/source-lineage ledger）可直接用草稿

**Table S1. Cohort and source-lineage ledger for the four audited public cohorts.**

| Cohort (GEO) | Subseries of | Cancer | Tissue / cell | Patients (donors) | Samples | ICI regimen | Sampling time | Endpoints (as defined by source) | Modality / platform | Source publication |
|---|---|---|---|---|---|---|---|---|---|---|
| GSE186143 | GSE186144 (with GSE189125 scRNA, GSE189126 TCR/BCR); BioProject PRJNA772729 | Metastatic melanoma | Peripheral blood (PBMC) | 61 in series (26 bulk-cohort-1 + 27 bulk-cohort-2 + 8 technical) + 5 healthy controls; 53 with clinical annotation; 47 evaluable for joint E×T | 66 (TPM matrix) | anti-PD-1 mono (29/53) or anti-PD-1 + anti-CTLA-4 (24/53) | Pretreatment baseline | T: highest irAE grade 0–4 (GEO); severe = grade ≥ 3 (17/53). E: DCB/NDB/NE (27/20/6) from Nat Med Suppl. Table 3 (not in GEO metadata) | Bulk RNA-seq (TruSeq RNA Exome); GPL16791 / HiSeq 2500; TPM processed matrix public, raw counts dbGaP-controlled | Lozano AX et al. Nat Med 2022;28(2). PMID 35027754 |
| GSE166181 | — (standalone) | Metastatic melanoma | Peripheral blood; FACS-sorted CD3⁺CD8⁺ T cells | 20 (11 responders / 9 non-responders); parent study 28 (17/11) | 183 libraries (20 × T0/T1/T2 × ~4 replicates); 51,701 cells | Nivolumab or pembrolizumab (anti-PD-1 mono) | T0 before therapy; T1 after 1st, T2 after 2nd cycle | E: RECIST at 6-month follow-up (R = CR/PR/SD/MR > 6 mo no progression; NR = PD). No toxicity annotation | scRNA-seq (Bio-Rad ddSEQ SureCell WTA 3′); GPL21697 / NextSeq 550 | De Biasi S et al. Nat Commun 2021;12:1669. PMID 33723257 |
| GSE216329 | — (BioProject PRJNA892990) | Predominantly NSCLC (18/24; other solid tumours) | Peripheral blood; enriched CD3 T cells | 24 (15 irAE / 9 no-irAE) | 48 (24 baseline _TN + 24 on-treatment _T) | anti-PD-1 (22/24 mono) | Baseline + on-treatment (paired) | T: any irAE, grade 2–3, organ-specific (pneumonitis/arthritis/thyroiditis/nephritis/colitis/neurotoxicity), 12-month follow-up — not grade ≥ 3 severe irAE. No efficacy annotation | scRNA-seq (10x Chromium V(D)J ; Cell Ranger 2.1.1); GPL24676 / NovaSeq 6000 | Bukhari S et al. Cell Rep Med 2023;4(1):100868. PMID 36513074 |
| GSE253720 | GSE253723 (biopsy scRNA subseries) | CPI-colitis with mixed primary tumours (melanoma 3, RCC 2, SCC-lung 1, prostate 1, HNSCC 1, Merkel 1) | Colon biopsy (endoscopic) | 14 donors in biopsy h5ad: 9 CPI-colitis + 3 healthy controls + 2 UC (study enrolled 15; 1 UC biopsy failed) | 6 pooled libraries, demuxlet-deconvoluted; 50,845 cells | anti-PD-1 mono (6/9) or anti-PD-1 + anti-CTLA-4 (3/9) | At active, steroid-experienced/refractory CPI colitis — NOT pretreatment baseline | T: histoendoscopically confirmed steroid-experienced CPI colitis (case-control vs HC and UC). No efficacy annotation; not a predictive-baseline design | Multiplexed biopsy scRNA + TCR-seq (10x Chromium 5′ ); GPL24676 / NovaSeq 6000 | He JY et al. J Immunother Cancer 2024;12(4):e008628. PMID 38642938 |

---

## 7. 残留 UNVERIFIED 项（3 项，均不阻碍投稿，建议 camera-ready 前补核）

1. **GSE186143 治疗线数**：GEO 无字段；Nat Med 正文/补充表未能取得（站点拦截）。
2. **GSE186143 论文-队列映射**：Nat Med 摘要 "3 cohorts (n=27, 26, 18)" 与 GEO bulk 26/27 的精确对应（尤其 n=18 与 bulk cohort 2 n=27 的关系）待对论文 Supplementary 核验。
3. **（结构性）NCBI acc.cgi 页面直核**：本次以 ftp.ncbi.nlm.nih.gov 的 family.soft.gz 官方记录替代（内容权威等同）；如需页面截图存档，需在网络恢复后补做。

**审计自检**：所有 verified_value 均可溯源到 `data/evidence/` 的冻结 SOFT 文件、`geo_data/` 本地数据或本审计的脚本化重解析（jq/h5py/awk，命令见会话记录）；GSE186143 临床标签与 GEO 的逐人交叉验证为零不一致。无任何值凭记忆转录。

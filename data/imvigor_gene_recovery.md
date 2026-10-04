#  IMvigor210 基因行身份恢复与锚定臂重建（Go/No-Go 级数据取证）

日期：2026-09-24 ｜ seed 20260911 ｜ 脚本：scripts/analysis/imvigor_gene_recovery.py ｜ 冻结 JSON：data/imvigor_gene_recovery.json

## 1. 恢复路径

**成功路径：阶段 1a（本地排序假设检验）**；1b/1c 作为独立确认（1c 已完成）。

恢复的钥匙（key）：**exmat 第 j 行 == fData 按 symbol（大写归一）字母序排序后的第 j 个基因**。
即上游匿名化只替换了行标签（gene_1..gene_31085）而保留了 `counts(cds)` 的行序，而该行序恰为 symbol 字母序；
项目管线旧假设（exmat 行 i == fData 行 i，即 entrez 字符串序位置对齐）因此证伪。

### 1a 排序假设检验（性别锚定门）

| 假设 | PASS | DDX3Y f/m | UTY f/m | RPS4Y1 f/m | XIST f/m | GAPDH 行和 | ACTB 行和 |
|---|---|---|---|---|---|---|---|
| position_fData (=entrez string sort, the falsified o | fail | 0.9749 | 0.9263 | 0.0000 | 0.9807 | 0 | 36 |
| symbol_alphabetical_upper | **PASS** | 0.0014 | 0.0008 | 0.0014 | 672.8287 | 16533113 | 118077389 |
| symbol_alphabetical_rawcase | **PASS** | 0.0014 | 0.0008 | 0.0014 | 672.8287 | 16533113 | 118077389 |
| entrez_numeric_desc | fail | 1.0842 | 0.4311 | 0.7234 | 1.4041 | 248111 | 1529517 |
| entrez_numeric_asc | fail | 0.6153 | 0.9558 | 1.0066 | 0.8519 | 71377 | 519353 |
| symbol_desc | fail | 0.9752 | 0.6255 | 0.0000 | 0.8703 | 2 | 26273 |
| length_asc | fail | 1.1246 | 0.9712 | 1.2322 | 0.9967 | 5622694 | 0 |
| length_desc | fail | 0.6111 | 0.9887 | 0.9288 | 0.1405 | 235571 | 0 |
| n_exons_asc | fail | 0.4332 | 0.8820 | 0.0000 | 1.3014 | 1910430 | 0 |

### 1a 深度验证（胜出假设 symbol 字母序）

- **Y 连锁基因族**（8 基因全部男性特异）：
  - DDX3Y: 女性均值 11.3 / 男性均值 8082.6（f/m=0.0014）
  - UTY: 女性均值 7.1 / 男性均值 9125.6（f/m=0.0008）
  - RPS4Y1: 女性均值 19.9 / 男性均值 13830.3（f/m=0.0014）
  - KDM5D: 女性均值 3.1 / 男性均值 6654.7（f/m=0.0005）
  - ZFY: 女性均值 8.3 / 男性均值 3682.9（f/m=0.0023）
  - EIF1AY: 女性均值 1.6 / 男性均值 1129.9（f/m=0.0014）
  - TMSB4Y: 女性均值 0.1 / 男性均值 38.2（f/m=0.0014）
  - USP9Y: 女性均值 34.6 / 男性均值 18840.6（f/m=0.0018）
- **顶级男性偏倚行（m≥50 且 f/m<0.05，共 25 行）中 chrY 标记命中 14 个**：DDX3Y, EIF1AY, KDM5D, NLGN4Y, PCDH11Y, PRKY, RPS4Y1, SEMG1, SEMG2, TBL1Y, TXLNGY, USP9Y, UTY, ZFY
- **XIST**：行 30113，f/m = 672.8（女均 430148 vs 男均 639）
- **管家/高表达**：全矩阵行和 Top30 = MALAT1, COL1A2, FN1, COL1A1, AHNAK, ACTB, COL3A1, RNF213, ACTG1, PABPC1, MUC4, SNORA73A…（MALAT1/COL1A2/FN1/ACTB 等，教科书级）
- **免疫浸润相关性（Spearman，对 pData TC/IC 分级）**：
  - CD8A: ρ(TC)=+0.299 (p=1.4e-08)，ρ(IC)=+0.554 (p=2.8e-29)
  - CD3D: ρ(TC)=+0.273 (p=2.4e-07)，ρ(IC)=+0.487 (p=4.6e-22)
  - GZMB: ρ(TC)=+0.326 (p=4.7e-10)，ρ(IC)=+0.542 (p=7.6e-28)
  - PRF1: ρ(TC)=+0.341 (p=7.0e-11)，ρ(IC)=+0.540 (p=1.2e-27)
  - CD274: ρ(TC)=+0.523 (p=9.7e-26)，ρ(IC)=+0.520 (p=2.0e-25)
  - PDCD1: ρ(TC)=+0.264 (p=5.9e-07)，ρ(IC)=+0.539 (p=1.4e-27)
  - CTLA4: ρ(TC)=+0.245 (p=4.0e-06)，ρ(IC)=+0.552 (p=4.5e-29)
  - CXCL9: ρ(TC)=+0.319 (p=1.2e-09)，ρ(IC)=+0.608 (p=1.9e-36)
  - CXCL10: ρ(TC)=+0.389 (p=5.5e-14)，ρ(IC)=+0.569 (p=3.3e-31)
- **Pearson 抽查**：CD8A~GZMB=+0.610，COL1A1~COL1A2=+0.982，ACTB~GAPDH=+0.449，XIST~DDX3Y=-0.369
- 排序键大小写歧义：raw-case 与 upper-case 排序在 145 行上不同（均为含小写 orf 的非面板基因）；**322 面板基因两种排序行号完全一致**（逐一核验，0 受影响）。

### 1a 跨队列身份验证（仅用本地数据）

对 320 个面板基因，比较 IMvigor210（恢复钥匙选行）与外部队列（自带真实标签）的**基因平均 log-CPM 秩相关**：

| 对照 | 恢复钥匙 Spearman ρ | 旧位置对齐 Spearman ρ |
|---|---|---|
| IMvigor210 vs GSE165252（320 基因） | **+0.813** (p=1.3e-76) | -0.170 |
| IMvigor210 vs GSE78220（320 基因） | **+0.619** (p=3.1e-35) | -0.066 |

同一批 symbol 在两个队列的平均表达量高度一致（ρ≈0.81/0.62），而旧位置对齐下≈0——行身份恢复正确性的又一独立证据。

### 1b 溯源

- 本地 exmat.csv / fData.csv / pData.csv 与上游 `mimifp/tfm_mUC`（2_data/1_IMvigor210/）**逐字节一致**（fData/pData SHA256 比对通过；exmat 尺寸 34,121,302 字节一致）。
- 上游 `1_download_data.R`（已存档 geo_data/imvigor210_labeled/upstream_1_download_data.R）显示：exmat 以 `counts(cds)` 行序写出、首列为 gene symbol；censored 版仅把行标签替换为 gene_N 占位符。
- 带标签镜像：`SiYangming/IMvigor210CoreBiologies`（官方 research-pub.gene.com 包的第三方镜像）data/cds.RData。

### 1c 数值恒等匹配（镜像 counts(cds) vs 本地 exmat，逐行 348 样本全值相等）

- **匹配 31085/31085 行（100%，0 行未匹配）**；值重复歧义行 5627（全零/等值行，取与 1a 钥匙一致候选）；唯一值匹配 25458 行中 25454 与钥匙完全一致。
- 与 1a symbol(upper) 钥匙一致：**31081/31085（99.987%）**；与 raw-case 排序一致：30936 行。
- 4 行不一致均为包版本间注记差异（C4BPA/C4BPB/C4B_2 家族 entrez 互换、CSNK1E 去重保留行），**均非 322 面板基因**。
- 排序键裁决：**upper**（大写归一 symbol 字母序）。

## 2. 验证门

| 门 | 结果 |
|---|---|
| 性别锚定 | PASS |
| 样本数 = 348 | True |
| 样本 ID 与 pData 一致（集合/顺序） | True / True |
| binaryResponse 计数（CR/PR=68，SD/PD=230，NA=50） | {'CR/PR': 68, 'SD/PD': 230, 'NA': 50} |
| 1c 数值恒等匹配 | True |
| 322 面板可测量基因数 | 新 320/322（旧假 entrez 选择 229；跨队列交集 320） |
| **总门** | **PASS** |

## 3. 主数字新旧对照（诚实呈现，数字变就是变）

对照基线：旧 = 冻结  S0（错误对齐+假 entrez 选择）经同一套度量代码复算（外部 AUC 与冻结  值逐位一致：GSE165252 0.5009==0.5009，GSE78220 0.4103==0.4103）。

| 量 | 旧值（错误对齐） | 新值（真实标签） | 变化 |
|---|---|---|---|
| 面板基因集（IMv 可测量/跨队列交集） | 229（假 entrez）→交集 228 | 320/322 →交集 320 | 选择修正（缺 IGHG1, PRKN） |
| 面板内全零基因数 | 70 | 0 | 消除 |
| 零值基因占 PC1 loading 平方质量 | 88.7% | 0.0% | 消除 |
| PC1/PC2 解释方差比 | 34.6% / 2.9% | 26.5% / 7.3% | 结构改变 |
| IMvigor210 within 表观 AUC | 0.5836 | 0.6659 | +0.0822 |
| IMvigor210 within Efron 0.632 AUC [95% CI] | 0.5669 [冻结 ] | 0.6571 [0.5857, 0.7284] | +0.0902 |
| GSE165252 外部 AUC | 0.5009 [0.3456, 0.6481] | 0.5036 [0.3540, 0.6585] | +0.0027 |
| GSE165252 校准 slope | 0.980 | -0.527 | -1.507 |
| GSE165252 Brier | 0.2320 | 0.2627 | +0.0307 |
| GSE165252 log-loss | 0.6717 | 0.8142 | +0.1425 |
| GSE78220 外部 AUC | 0.4103 [0.1948, 0.6603] | 0.6603 [0.4038, 0.8617] | +0.2500 |
| GSE78220 校准 slope | -2.481 | 5.678 | +8.159 |
| GSE78220 Brier | 0.3471 | 0.2280 | -0.1191 |
| GSE78220 log-loss | 0.9564 | 0.6480 | -0.3084 |

要点（不挑好的说）：

- GSE165252 外部 AUC 基本没动（0.5009→0.5036），CI 仍横跨 0.5；校准 slope 由 0.98 恶化为 -0.53。
- GSE78220 外部 AUC 大幅改善（0.4103→0.6603），但 CI 仍宽（n=25）；slope 由 -2.48 变为 5.68。
- 锚定队列内部判别改善（0.632 AUC 0.5669→0.6571），且 70 基因全零伪影与 88.7% loading 质量被彻底消除——旧结果中的“未测量基因驱动 PC1”结构性缺陷不复存在。

## 4. 对投稿的影响评估

一句话：此前冻结的主数字（外部 AUC 0.501/0.410、70 零值基因、88.7% PC1 伪影）建立在双重错位的基因行上，
按真实标签重建后锚定臂内部与 GSE78220 运输显著改善、GSE165252 运输仍为无效（AUC≈0.50，CI 横跨 0.5），
因此论文必须以本  数字替换全部受影响结论后方可投稿；任何基于旧 S0 的图、表与敏感性分析都需重算。

## 5. 溯源与哈希

- exmat.csv SHA256：`fc51849b9e923b6a0ac8341c443d83d35d9fb2643d19b15554b3e444f74b5445`
- fData.csv SHA256：`d0f781b5209748ed8475918c3db055db72d83f2838a9ce6636cbfe674b0f196b`（== 上游 fData_IMvigor210.csv `d0f781b5209748ed8475918c3db055db72d83f2838a9ce6636cbfe674b0f196b`）
- pData.csv SHA256：`33aafccabc6ea95fa95aa8049ace5261c7e5089ecde28811e877f0f5c010f0de`（== 上游 `33aafccabc6ea95fa95aa8049ace5261c7e5089ecde28811e877f0f5c010f0de`）
- 镜像 cds.RData（https://raw.githubusercontent.com/SiYangming/IMvigor210CoreBiologies/master/data/cds.RData）SHA256：`b7f1cfb648cc6125e48ba9ff5094e4ed547ad4758cf68fca5117767ebf65df2e`（仅用于 1c 恒等核对，不进入分析数值）
- 上游提取与匿名化脚本存档：geo_data/imvigor210_labeled/upstream_1_download_data.R、upstream_2_preprocess_data.R、upstream_README.md

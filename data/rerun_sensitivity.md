#  重算敏感性报告（P0-3 / P0-10 新数字版）+ P0-12 manifest

日期：2026-09-24 ｜ seed 20260911 ｜ 脚本：scripts/analysis/rerun_sensitivity_manifest.py ｜ 冻结 JSON：data/rerun_sensitivity.json

基线：IMvigor210 基因行身份已恢复（symbol 字母序钥匙，_imvigor_gene_recovery 总门 PASS）；primary frozen = 320 基因 CLR→PCA k=2 锚定臂（298 labelled），ridge λ=1/298 运输。全部数字由本脚本从原始矩阵重算，并与恢复报告冻结值逐位交叉校验（tol 1e-9）。

## 1. 重建一致性（先验证再分析）

| 队列 | AUC（重算） | AUC（冻结） | slope（重算/冻结） | Brier | log-loss | 全部一致 |
|---|---|---|---|---|---|---|
| GSE165252 | 0.5036 | 0.5036 | -0.527 / -0.527 | 0.2627 | 0.8142 | **是** |
| GSE78220 | 0.6603 | 0.6603 | 5.678 / 5.678 | 0.2280 | 0.6480 | **是** |

## 2. P0-3 阈值敏感性（新数字版）

覆盖：primary frozen 两个方向（GSE165252 / GSE78220）。M0-M6 基准模型行 **待补**：rerun_benchmark.json 尚未落盘（manifest 中记 pending）。

判定规则：C1 = 外部 AUC BCa 95% CI（2000，seed 20260911）向上排除 0.5（与容差无关）；C2 = |slope−1| ≤ tol（slope 用恢复报告口径：以预测概率为特征的重标定斜率，即引文基线 −0.527 / 5.678；logit 口径同时报告）。

| 队列 | AUC [BCa CI] | C1 | slope(概率口径) [boot CI] | slope(logit口径) [boot CI] | |slope−1| | tol=0.05 | 0.10 | 0.20 |
|---|---|---|---|---|---|---|---|---|
| GSE165252 | 0.5036 [0.3540, 0.6585] | 不过 | -0.527 [-9.20, 7.99] | -0.06 [-0.90, 0.85] | 1.527 | 不过 | 不过 | 不过 |
| GSE78220 | 0.6603 [0.4038, 0.8617] | 不过 | 5.678 [-1.88, 18.43] | 1.35 [-0.43, 4.52] | 4.678 | 不过 | 不过 | 不过 |

- **判定翻转格数（overall，概率口径）：0**（翻转队列：无）；C2 单项翻转 0 格。
- **主结论是否依赖 0.10 容差：否** —— 任何容差下均无任何格允许运输（any_transport_permitted = False）：两个方向的 |slope−1|（1.53 / 4.68）远超最大容差 0.20，且 C1 在两个方向均不过（AUC CI 含 0.5）。
- 两口径（概率特征 / logit 特征）slope 结论一致；slope bootstrap 95% CI 在两个方向均极宽，任何容差下均为 inconclusive（CI 同时跨越 [1−tol, 1+tol] 边界）。

## 3. P0-10 ΔAUC 连续面板（新数字版）

within 侧（任务规格：恢复钥匙 + 同一加载逻辑 + CLR + PCA k=2 + λ=1/n + 600 bootstrap，Efron 0.632 optimism 形式，与 IMvigor210 基线 0.6571 同一估计器）：

| 队列 | within 0.632 [CI] | λ | n（事件） | nested bootstrap |
|---|---|---|---|---|
| IMvigor210 | 0.6571 [0.5857, 0.7284] | 0.002538 | 298（68） | 未重跑（冻结 600-boot CI） |
| GSE165252 | 0.4096 [0.3002, 0.5191] | 0.014085 | 71（23） | 2000 抽样，CI [0.4046, 0.6706] |
| GSE78220 | 0.4627 [0.2738, 0.6515] | 0.040000 | 25（13） | 2000 抽样，CI [0.3553, 0.8280] |

primary frozen 行（M0 / M2 / M3 **待补**——rerun_benchmark.json 未落盘，已在 manifest 记 pending，落盘后同表追加）：

| 目标队列 | 外部 AUC [BCa CI] | ΔAUC（外−内） | ΔAUC 95% CI（配对 bootstrap） | logit 截距 | logit slope | 概率口径 slope | Brier | log-loss |
|---|---|---|---|---|---|---|---|---|
| GSE165252 | 0.5036 [0.3540, 0.6585] | 0.0940 | [-0.2234, 0.1700] | -0.87 | -0.06 | -0.527 | 0.2627 | 0.8142 |
| GSE78220 | 0.6603 [0.4038, 0.8617] | 0.1976 | [-0.2509, 0.3951] | -0.04 | 1.35 | 5.678 | 0.2280 | 0.6480 |

解读（不挑好的说）：

- GSE165252：ΔAUC = 0.0940 [-0.2234, 0.1700]，CI 均含 0——同管线 within 参考下，外部表现与队列内表现无统计可分辨的差异。
- GSE78220：ΔAUC = 0.1976 [-0.2509, 0.3951]，CI 均含 0——同管线 within 参考下，外部表现与队列内表现无统计可分辨的差异。

### Pooled AUC（DL 随机效应，k=3）

- **新 320 基因口径（主规格：IMvigor within 0.632 + 两外部队列运输 AUC）：pooled AUC = 0.615 [0.515, 0.715]，I²(Q) = 0.387，Q = 3.26（df=2，P = 0.195），tau² = 0.00326**
- 次规格（三队列全部 within 0.632，镜像旧  设计）：pooled = 0.517 [0.335, 0.699]，I²(Q) = 0.868
- 历史  descriptive **0.498 [0.365, 0.631] 保留在案**，但注明为旧 229 基因（错误位置对齐）口径的历史值，已由上列新口径合成取代；异质性 caveat 不变（I² 高，单一可运输 AUC 解释仍不可取）。

## 4. P0-12 manifest 摘要

- 冻结映射：data/analysis_manifest.json —— 10 个  分析条目（-01 … -09，含 phase1a 取证组与 cds 镜像提取子条目），每条含脚本/输入/输出/seed/bootstrap 次数/关键参数/状态（done / pending / superseded / split）。
- IMvigor210 三文件 SHA256 与恢复报告逐一断言相等：exmat `fc51849b…`、fData `d0f781b5…`、pData `33aafcca…`；上游溯源（mimifp/tfm_mUC censored、SiYangming/IMvigor210CoreBiologies 镜像）已固化。
- 伪影撤回文档：data/artefact_retraction.md（ED10/provenance 教训表）。

## 5. 文件

- 脚本：scripts/analysis/rerun_sensitivity_manifest.py
- 数据：data/rerun_sensitivity.json；data/analysis_manifest.json；within 嵌套 bootstrap 缓存 data/__within_boot_newkey_cache.npz
- 报告：data/rerun_sensitivity.md；撤回文档 data/artefact_retraction.md
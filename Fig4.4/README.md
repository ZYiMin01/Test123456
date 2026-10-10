# Figure 4.4 进展模式与死亡方式关联分析

## 运行
```bash
# 将 raw.xlsx / baseline.xlsx / cod.xlsx 放到 data/ 下（患者级数据，不入库）
Rscript R/fig4_4_analysis.R                 # 默认 B=1000 bootstrap
N_BOOT=200 N_CORES=1 Rscript R/fig4_4_analysis.R   # Windows 自动单核
```
R 包：survival, riskRegression, prodlim, coxphf, ggplot2, patchwork, openxlsx, readxl, gridExtra, scales。
脚本为 UTF-8 编码（含中文列名），在 R 中请 `source("R/fig4_4_analysis.R", encoding = "UTF-8")`。

## 输出（output/）
| 文件 | 内容 |
|---|---|
| `Fig4.4A_adjusted_CIF_matrix.pdf` | 主图 A：调整后 12 个月 CIF 矩阵 + 死因构成条形图 |
| `Fig4.4B_csHR_matrix_refIHNL.pdf` | 主图 B：IH-NL 参照 HR/csHR 矩阵 + 整体 Wald P + Lunn–McNeil P |
| `Fig4.4_Supplementary.pdf` | sHR、ARD、PMR、未调整 AJ、6/24 个月 CIF、Firth、排除无法判定死因、纳入标记病例 |
| `Fig4.4_Full_Report.pdf` | 汇总报告（21 页）：主图 + 所有表格 + 补充图 |
| `Fig4.4_Tables.xlsx` | 所有估计值（CIF/ARD/Q/PMR/HR，含 CI、事件数）和检验 P 值、验证、QC |
| `png/` | 所有图的 PNG |

`Methods_Fig4.4.md`：英文 Statistical analysis 段落。`output_run.log`：完整运行日志。

## 采用的默认设定（交接文档第 6 节待确认问题）
1. 模式归类：癌栓升级 = 肝内 TG；新发血管侵犯 = 肝内 NL（交接文档规则）。
2. 死因合并：IH / EH / Systemic HCC-related 各为一类；Non-HCC = Liver-related + Other Causes + Treatment-related；
   Undetermined = Unclassified HCC-related + Unknown + 死亡但死因缺失（共 16 例，作为竞争事件保留，敏感性分析中排除）。
3. 协变量：TTP + 年龄、性别、HBV、ECOG、ALBI、AFP≥400、PIVKA≥1000、肿瘤数>3、肿瘤大小、Vp 分级、HVTT、EHS、肝硬化；strata(Group)。
4. 排除：cod `Delete` 标记的 11 例进展患者（"可纳入"标记除外）和 1 例 PFS > OS 的数据错误；纳入标记病例作为敏感性分析。
5. 主时间点 12 个月（进展后中位随访 18.7 个月；12 个月时各模式风险人数 ≥ 7）。

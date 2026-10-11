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
| `Fig4.4A_main_CIF12_plus_composition.pdf` | 正文主图 A：12 个月调整后 CIF 矩阵 + 死因构成 |
| `Fig4.4A1_adjusted_CIF12_matrix.pdf` | A1：调整后 12 个月 CIF 矩阵（9 模式 × 4 死因 + 全因） |
| `Fig4.4A2_ARD12_vs_population_average.pdf` | A2：ARD vs 人群标准化平均 |
| `Fig4.4A3_cause_composition_6_12_24m.pdf` | A3：6/12/24 个月死因构成 |
| `Fig4.4A4_PMR12_enrichment.pdf` | A4：PMR 富集 |
| `Fig4.4B1_csHR_matrix_refSysTG.pdf` | B1：Sys-TG 参照 HR/csHR 矩阵 + 整体 Wald P + Lunn–McNeil P |
| `Fig4.4B1_sHR_matrix_refSysTG.pdf` | B1：Sys-TG 参照 sHR 矩阵 |
| `Fig4.4B2_LunnMcNeil_tests.pdf` | B2：Lunn–McNeil 检验 + 后续检验 |
| `Fig4.4_Supplementary.pdf` | 补充图（未调整 AJ、6/24 个月 CIF/ARD/PMR、Firth、敏感性分析） |
| `Fig4.4_Full_Report.pdf` | 汇总报告：全部图 + 表格 |
| `Fig4.4_Tables.xlsx` | 所有估计值和检验 P 值 |
| `png/` | 所有图的 PNG |

`Fig4.4_结果说明.docx`：各结果的分析方法、目的与结果阐释（Word）。`Methods_Fig4.4.md`：英文 Statistical analysis 段落。`output_run.log`：运行日志。
4.4B 参照组由环境变量 `REF_HR` 控制（默认 `Sys-TG`）。

## 采用的默认设定（交接文档第 6 节待确认问题）
1. 模式归类：按《解剖分类 × 肿瘤模式分类》文档。TTU（癌栓升级）并入 IHG，nMVI（新发血管侵犯）并入 NIH；Sys-TG = 仅 STG，Sys-NL = NIH + NEH，其余组合均归为 Sys-TG+NL。
   与早期规则相比，有 22 例改判（STG+TTU 14 例，EHG+TTU 8 例：原 Sys-TG → 现 Sys-TG+NL）；早期规则结果见补充图 figS_v0 和 xlsx 的 `Sens_earlier_rule_HR`，逐组合对照见 `Pattern_classification`。
2. 死因合并：IH / EH / Systemic HCC-related 各为一类；Non-HCC = Liver-related + Other Causes + Treatment-related；
   Undetermined = Unclassified HCC-related + Unknown（分析队列中共 8 例，作为竞争事件保留，敏感性分析中排除；全队列 16 例死因缺失者均不在进展队列中）。
3. 协变量：TTP + 年龄、性别、HBV、ECOG、ALBI、AFP≥400、PIVKA≥1000、肿瘤数>3、肿瘤大小、Vp 分级、HVTT、EHS、肝硬化；strata(Group)。
4. 排除：cod `Delete` 标记的 11 例进展患者（"可纳入"标记除外）和 1 例 PFS > OS 的数据错误；纳入标记病例作为敏感性分析。
5. 主时间点 12 个月（进展后中位随访 18.7 个月；12 个月时各模式风险人数 ≥ 7）。

# zhangsiyu Skill

生物信息学分析 Skill 集合，共 18 个模块。维护者：Zhang Siyu。

[ORCID 0009-0005-5383-541X](https://orcid.org/0009-0005-5383-541X) 仅用于署名标识，不代表账号认证或科学背书，也不是运行依赖。

## Skill 目录

- [atac-seq-analysis](skills/atac-seq-analysis/SKILL.md)：ATAC-seq 与 scATAC-seq。
- [biofigure-polish](skills/biofigure-polish/SKILL.md)：生信图表整理与美化。
- [bulkrna-analysis](skills/bulkrna-analysis/SKILL.md)：Bulk RNA-seq。
- [crispr-screen-analysis](skills/crispr-screen-analysis/SKILL.md)：CRISPR 筛选。
- [drug-target-repurposing](skills/drug-target-repurposing/SKILL.md)：药物靶点与老药新用候选排序。
- [enrichment-analysis](skills/enrichment-analysis/SKILL.md)：富集分析。
- [gbd-burden-analysis](skills/gbd-burden-analysis/SKILL.md)：GBD 负担质控、EAPC 与变化分解。
- [gwas-causal-target-analysis](skills/gwas-causal-target-analysis/SKILL.md)：GWAS、QTL 与因果靶点分析指导及质控。
- [immune-repertoire-analysis](skills/immune-repertoire-analysis/SKILL.md)：免疫组库。
- [large-scale-singlecell-engineering](skills/large-scale-singlecell-engineering/SKILL.md)：大规模单细胞资源规划与工程化。
- [microbiome-analysis](skills/microbiome-analysis/SKILL.md)：微生物组。
- [network-pathway-analysis](skills/network-pathway-analysis/SKILL.md)：网络与通路分析。
- [proteomics-metabolomics-analysis](skills/proteomics-metabolomics-analysis/SKILL.md)：蛋白组与代谢组。
- [public-db-survival](skills/public-db-survival/SKILL.md)：公共队列表达与生存分析。
- [public-omics-meta-analysis](skills/public-omics-meta-analysis/SKILL.md)：公共多队列整理与 Meta 分析。
- [singlecell-analysis](skills/singlecell-analysis/SKILL.md)：单细胞转录组。
- [spatial-transcriptomics-analysis](skills/spatial-transcriptomics-analysis/SKILL.md)：空间转录组。
- [variant-cnv-signature-analysis](skills/variant-cnv-signature-analysis/SKILL.md)：突变、拷贝数与突变特征。

## 使用方式与边界

下载仓库后，保留所需 `skills/<名称>/` 的完整目录结构，交给支持 Agent Skills 的工具读取其中的 `SKILL.md`。依赖和输入要求以每个模块的文档为准。

这些模块包含工作流指引、输入质控、分析模板和部分可执行脚本，不是已在所有真实队列上验证的一键全流程软件。大规模单细胞分析仍需适配实际服务器、数据与软件环境。统计关联、遗传证据和药物候选排序不能替代因果验证、实验验证或临床判断。

GBD 核心检查可在仓库根目录运行：

```sh
Rscript --vanilla skills/gbd-burden-analysis/tests/test_core.R
```

本次发布进行了文本文件安全筛查和 Python 语法检查；未宣称所有模块完成真实数据端到端测试。`tests/fixtures` 中的小型示例仅供测试，不是研究证据。GBD 具体测试范围见其说明。

## 发布范围与许可

包含 Skill 说明、脚本、参考文档、配置与小型测试样例。不包含本地分析输出、真实研究数据、课程原稿、直播课件、缓存或凭证。

既有 GBD 模块的 MIT 许可保留，见 [LICENSE](LICENSE)。新增模块保留原有声明；未提供独立许可的模块不因公开上传而自动获得 MIT 许可。使用者应核实其来源和再分发权限。仓库根目录的 [CITATION.cff](CITATION.cff) 目前仅对应 GBD 核心工具，不代表其他模块的原创归属。

`release-manifest.json` 记录本次发布文件的 SHA-256 校验值。

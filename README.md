# zhangsiyu Skill

Zhang Siyu 的科研分析 Skills。

作者标识：[ORCID 0009-0005-5383-541X](https://orcid.org/0009-0005-5383-541X)

目前包含 [gbd-burden-analysis](skills/gbd-burden-analysis/SKILL.md)：GBD 输入质控、分层 EAPC 趋势估计、人口增长/年龄结构/年龄别率变化分解，以及其他 GBD 分析的方法指导。

## 快速开始

安装 R 后，在仓库根目录执行：

```sh
Rscript --vanilla skills/gbd-burden-analysis/scripts/check_dependencies.R core
Rscript --vanilla skills/gbd-burden-analysis/tests/test_core.R
Rscript --vanilla skills/gbd-burden-analysis/scripts/run_analysis.R eapc your-rates.csv new-results
```

将 `skills/gbd-burden-analysis` 完整目录复制到你的 Agent Skills 目录后重新加载任务；也可以直接让代理读取其 `SKILL.md`。输入要求见 [数据格式](skills/gbd-burden-analysis/references/input-schema.md)。

核心工具仅需 base R，可离线运行，无需 GitHub/ORCID 登录。地图、Excel 和 BAPC 需要额外软件及项目级验证。发布内容不含真实 GBD 数据和课程原稿。测试使用合成数据；通过测试不代表所有统计模型适用于所有研究问题。

代码与维护文档采用 MIT 许可，见 [LICENSE](LICENSE)。引用信息见 [CITATION.cff](CITATION.cff)。

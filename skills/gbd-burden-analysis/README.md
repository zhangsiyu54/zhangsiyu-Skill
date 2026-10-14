# GBD Burden Analysis Skill

面向 Codex 等兼容 Agent Skills 的 GBD 分析工具。可直接运行：输入质控、分层 EAPC、人口增长/年龄结构/年龄别率变化三因素分解。核心计算仅需 R，无第三方 R 包依赖。

作者：Zhang Siyu · [ORCID: 0009-0005-5383-541X](https://orcid.org/0009-0005-5383-541X)。此链接用于署名；不代表 ORCID 认证、背书或账号登录绑定。GitHub/ORCID 都不是运行依赖。

## 使用

在本目录执行：

```sh
Rscript scripts/check_dependencies.R core
Rscript scripts/validate_gbd_input.R your-gbd.csv
Rscript scripts/run_analysis.R eapc your-rates.csv new-eapc-output
Rscript scripts/run_analysis.R decompose your-age-table.csv new-decomposition-output 100000
Rscript tests/test_core.R
```

输出目录必须尚不存在；每次输出包含结果、输入与脚本 MD5、R 会话信息。命令支持包含空格的路径，请加引号。输入字段与合成示例见 [输入格式](references/input-schema.md)。直接在 Codex 中指定本目录的 `SKILL.md` 也可使用；安装时将完整目录复制到个人 Skills 目录并重新加载任务。

## 能力边界

EAPC 区间是回归模型的置信区间，不是 GBD 原始不确定性区间。三因素分解使用六种替换顺序平均的 Shapley 贡献，仅输出点估计。地图、SDI、风险归因摘要、APC/BAPC 有方法指导，但仍需针对研究数据编写并验证分析代码。

发布版不包含课程原稿、真实数据库导出、人口数据或视频。课程原稿可保留在作者本地，核心运行不依赖它们。全部测试数据均为合成数据。

# Top References

- Colaprico A et al. TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. Nucleic Acids Research. 2016. https://doi.org/10.1093/nar/gkv1507
- Liu J et al. An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. Cell. 2018. https://doi.org/10.1016/j.cell.2018.02.052
- Therneau TM, Grambsch PM. Modeling Survival Data: Extending the Cox Model. Springer. 2000. https://doi.org/10.1007/978-1-4757-3294-8
- Kassambara A, Kosinski M, Biecek P. survminer: Drawing Survival Curves using ggplot2. R package. https://rpkgs.datanovia.com/survminer/
- Cerami E et al. The cBio Cancer Genomics Portal: An Open Platform for Exploring Multidimensional Cancer Genomics Data. Cancer Discovery. 2012. https://doi.org/10.1158/2159-8290.CD-12-0095
- Goldman MJ et al. Visualizing and interpreting cancer genomics data via the Xena platform. Nature Biotechnology. 2020. https://doi.org/10.1038/s41587-020-0546-8
- Cox DR. Regression Models and Life-Tables. Journal of the Royal Statistical Society Series B. 1972. https://doi.org/10.1111/j.2517-6161.1972.tb00899.x (foundational Cox proportional hazards model paper)
- Grambsch PM, Therneau TM. Proportional hazards tests and diagnostics based on weighted residuals. Biometrika. 1994. https://doi.org/10.1093/biomet/81.3.515 (Schoenfeld-residual-based test for the proportional hazards assumption; implemented as `cox.zph` in R)
- Hothorn T, Lausen B. On the exact distribution of maximally selected rank statistics. Computational Statistics & Data Analysis. 2003. https://doi.org/10.1016/S0167-9473(02)00225-6 (statistical basis for data-driven "optimal" expression cutoffs; must be reported with correction-adjusted p values, see `references/pitfalls.md`)
- Budczies J et al. Cutoff Finder: A Comprehensive and Straightforward Web Application Enabling Rapid Biomarker Cutoff Optimization. PLOS ONE. 2012. https://doi.org/10.1371/journal.pone.0051862
- Fine JP, Gray RJ. A Proportional Hazards Model for the Subdistribution of a Competing Risk. Journal of the American Statistical Association. 1999. Vol 94, 496-509. (subdistribution/Fine-Gray model for competing risks; contrast with cause-specific Cox)
- Peduzzi P, Concato J, Kemper E, Holford TR, Feinstein AR. A simulation study of the number of events per variable in logistic regression analysis. Journal of Clinical Epidemiology. 1996. https://doi.org/10.1016/S0895-4356(96)00236-3 (origin of the ~10 events-per-variable heuristic for Cox/logistic models)
- Vittinghoff E, McCulloch CE. Relaxing the Rule of Ten Events per Variable in Logistic and Cox Regression. American Journal of Epidemiology. 2007. https://doi.org/10.1093/aje/kwk052 (re-examines and relaxes the EPV=10 rule under some conditions)
- Suissa S. Immortal Time Bias in Pharmacoepidemiology. American Journal of Epidemiology. 2008. https://doi.org/10.1093/aje/kwm324
- McShane LM, Altman DG, Sauerbrei W, Taube SE, Gion M, Clark GM. REporting recommendations for tumour MARKer prognostic studies (REMARK). British Journal of Cancer. 2005. https://doi.org/10.1038/sj.bjc.6602678 (reporting standard for prognostic biomarker studies; motivates transparent methods text in `references/deliverables.md`)
- Tang Z, Kang B, Li C, Chen T, Zhang Z. GEPIA2: an enhanced web server for large-scale expression profiling and interactive analysis. Nucleic Acids Research. 2019. https://doi.org/10.1093/nar/gkz430
- Lánczky A, Győrffy B. Web-Based Survival Analysis Tool Tailored for Medical Research (KMplot): Development and Implementation. Journal of Medical Internet Research. 2021. https://doi.org/10.2196/27633

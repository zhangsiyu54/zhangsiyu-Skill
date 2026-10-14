# Paper Practice

Follow MiXCR/TRUST4/scRepertoire/immunarch-style repertoire papers: clonotype definition and chain handling determine the biological meaning.

## Reference Sources

- MiXCR GitHub: https://github.com/milaboratory/mixcr
- MiXCR Nature Methods: https://www.nature.com/articles/nmeth.3364
- TRUST4 GitHub: https://github.com/liulab-dfci/TRUST4
- TRUST4 Nature Methods: https://www.nature.com/articles/s41592-020-0957-5
- scRepertoire GitHub: https://github.com/BorchLab/scRepertoire
- immunarch: https://immunarch.com/

## Paper-style Flow

1. Confirm chain type, clonotype definition, CDR3 columns, V/J genes, sample IDs, and count/frequency columns.
2. Summarize clonotype abundance, expansion, diversity, and sample overlap.
3. Add V/J usage, CDR3 length, alluvial clone tracking, or cell type coupling when fields exist.
4. State that expansion does not prove antigen specificity unless supported by antigen data.
5. Document clonotype definition and filtering thresholds.

## Visualization Pattern

- Expansion plots should group rare clones and label only dominant clones.
- V/J usage should use horizontal bars when labels are long.
- Alluvial/circos-style plots should be limited to top clonotypes.

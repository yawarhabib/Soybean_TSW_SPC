# SoybeanGWAS

Genome-wide association study and genomic prediction of thousand-seed weight (TSW) and seed protein content (SPC)
in the Khabarovsk soybean diversity panel, Russian Far East. Field trials 2017-2024.
Code for a manuscript in preparation. Scripts are added as each step of the analysis is finalised.

## Layout

```
Paper_TSW_SPC/
  scripts/     analysis scripts, numbered in pipeline order
  data/        raw/ inputs and processed/ files read by later steps (not in the repository yet)
  results/     figures and tables written by the scripts (not in the repository yet)
```

Scripts use `here::here()` paths from the repository root, marked by the `.here` file, so the
repository runs from any location. They are written for interactive use in RStudio with the
repository root as the project folder.

## Pipeline

| Step | Script | Content | Status |
|------|--------|---------|--------|
| 01 | `01_Manuscript_Analysis_Phenotypes.R` | descriptive statistics; outlier screen (Bernal-Vasquez et al. 2016); genotype BLUEs with season as a fixed effect; entry-mean and Cullis heritability; bivariate model for the genetic correlation between TSW and SPC; residual diagnostics | done |
| 02 | | genotype quality control | to come |
| 03 | | population structure and linkage disequilibrium | to come |
| 04 | | GWAS | to come |
| 05 | | post-GWAS | to come |
| 06 | | genomic prediction | to come |

## Requirements

R 4.6.1 with dplyr, ggplot2, tidyr, lme4, emmeans, sommer, patchwork and here.

## Data

The phenotype and genotype data are not part of the repository yet. They will be added, or
deposited with a DOI, when the manuscript is submitted.

## Author

Yawar Habib, Center for Bio and Medical Technologies, Skolkovo Institute of Science and Technology, Moscow.

Code is released under the MIT licence, see `LICENSE`.

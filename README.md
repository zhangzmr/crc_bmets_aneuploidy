# CRC brain metastases analysis

R scripts + **demo data** for the manuscript:

**Aneuploidy profiling among different types of brain metastases reveals a colon cancer-specific chromosomal signature**

Aparajita Khan, Mengrui Zhang, Susan M. Grimes, Billy T. Lau, Anuja Sathe, Summer Han, Hanlee P. Ji  

Stanford University, Stanford, CA

## Layout

```
code/GitHub/
├── README.md
├── .gitignore
├── data/                 # demo ASCETS arm calls, clinical, MAF (small)
└── scripts/
    ├── 01_cancer_type_arm_CNV_CRC.R
    ├── 01b–01d_*.R
    ├── 01e_plot_cancer_type_heatmaps_from_csv.R
    ├── 02_metastatic_site_arm_CNV.R
    ├── 03_crc_arm_cosegregation.R
    └── 04_mutation_coBarplot_brain_vs_primary.R
```

No external path config: each script sets `REPO_ROOT` to this folder and reads `data/`.

## Demo data

Synthetic, small cohorts under `data/` (same relative paths as the full project):

- `data/ASCETS/ascets-master/...` — arm AMP/DEL/NEUTRAL calls (CRC + LC/MELAN/BRCA)
- `data/msk_met_2021_clinical_data.tsv`
- `data/msk_met_2021/data_mutations.maf` + clinical sample file

Enough to smoke-test pipelines end-to-end. Replace with real ASCETS / MSK-MET files for manuscript figures.

## Figure → script map

| Panel (approx.) | Script |
|---|---|
| Fig 2 arm AMP/DEL vs other cancer BMets | `01_cancer_type_arm_CNV_CRC.R` |
| Fig 2 / Supp arm co-occurrence | `03_crc_arm_cosegregation.R` |
| Fig 4 CRC Brain vs other sites | `02_metastatic_site_arm_CNV.R` |
| Fig 5 mutations Brain vs Primary | `04_mutation_coBarplot_brain_vs_primary.R` |

## Requirements

```r
install.packages(c(
  "tidyverse", "dplyr", "tidyr", "stringr", "ggplot2", "scales",
  "data.table", "RColorBrewer", "readxl", "writexl", "circlize"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "maftools"))
```

## Run

```bash
cd code/GitHub
Rscript scripts/01_cancer_type_arm_CNV_CRC.R
Rscript scripts/02_metastatic_site_arm_CNV.R
Rscript scripts/03_crc_arm_cosegregation.R
Rscript scripts/04_mutation_coBarplot_brain_vs_primary.R
```

Scripts create an `outputs/` directory when run (not committed).

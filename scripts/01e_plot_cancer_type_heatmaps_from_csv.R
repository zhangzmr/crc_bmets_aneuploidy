rm(list=ls())


# Self-contained repo root (directory containing scripts/ and data/)
.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", grep("^--file=", .args, value = TRUE))
.script_dir <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
REPO_ROOT <- normalizePath(file.path(.script_dir, ".."), mustWork = TRUE)
setwd(REPO_ROOT)
OUT_ROOT <- file.path(REPO_ROOT, "outputs")
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)
DIR_CANCER_TYPE_ARM <- file.path(OUT_ROOT, "cancer_type_arm")
DIR_METSITE_ARM <- file.path(OUT_ROOT, "metastatic_site_arm")
DIR_COSEG <- file.path(OUT_ROOT, "crc_arm_cosegregation")
DIR_MUTATION <- file.path(OUT_ROOT, "mutation_figure")

library(ComplexHeatmap)
library(ggplot2)
library(readxl)
library(circlize)
library(RColorBrewer)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Requires prop_* CSVs from scripts/01_cancer_type_arm_CNV_CRC.R
.csv_dir <- DIR_CANCER_TYPE_ARM
if (!file.exists(file.path(.csv_dir, "prop_AMP_CRC_LC.csv"))) {
  stop("Missing prop_AMP_CRC_LC.csv in outputs/cancer_type_arm/. Run 01_cancer_type_arm_CNV_CRC.R first.")
}
dir.create(DIR_CANCER_TYPE_ARM, showWarnings = FALSE, recursive = TRUE)

######################################################
## amplification
LC <- read.csv(file.path(.csv_dir, "prop_AMP_CRC_LC.csv"))
LC <- LC[-1, ]
rownames(LC) <- LC[, 1]
LC <- LC[, -1]

MELAN <- read.csv(file.path(.csv_dir, "prop_AMP_CRC_MELAN.csv"))
MELAN <- MELAN[-1, ]
rownames(MELAN) <- MELAN[, 1]
MELAN <- MELAN[, -1]

BRCA <- read.csv(file.path(.csv_dir, "prop_AMP_CRC_BRCA.csv"))
BRCA <- BRCA[-1, ]
rownames(BRCA) <- BRCA[, 1]
BRCA <- BRCA[, -1]

ALL <- read.csv(file.path(.csv_dir, "prop_AMP_CRC_ALL.csv"))
ALL <- ALL[-1, ]
rownames(ALL) <- ALL[, 1]
ALL <- ALL[, -1]

mat1 <- data.frame(LC = LC$Prop_CRC - LC$Prop_LC, 
                   MELAN = MELAN$Prop_CRC - MELAN$Prop_MELAN, 
                   BRCA = BRCA$Prop_CRC - BRCA$Prop_BRCA,
                   ALL = ALL$Prop_CRC - ALL$Prop_ALL, 
                   row.names = rownames(LC))
colnames(mat1) <- c("CRC-LC", "CRC-MELAN", "CRC-BRCA", "CRC-ALL")
# mat1[is.na(mat1)] <- 0

mat2 <- data.frame(LC = p.adjust(LC$pval, method = "fdr"),
                   MELAN = p.adjust(MELAN$pval, method = "fdr"),
                   BRCA = p.adjust(BRCA$pval, method = "fdr"),
                   ALL = p.adjust(ALL$pval, method = "fdr"),
                   row.names = rownames(LC))
mat2[is.na(mat2)] <- 1

mat1 <- t(as.matrix(mat1))
mat2 <- t(as.matrix(mat2))

colnames(mat1) <- gsub("Arm_", "", colnames(mat1))
colnames(mat2) <- gsub("Arm_", "", colnames(mat2))


cell_fun <- function(j, i, x, y, w, h, fill) {
  if(mat2[i, j] < 0.01) {
    grid.text("**", x, y)
  } else if(mat2[i, j] < 0.05) {
    grid.text("*", x, y)
  } else if(mat2[i, j] < 0.001) {
    grid.text("***", x, y)
  }
}

ann <- data.frame(Comparisons = rownames(mat1))
ann$Comparisons <- factor(ann$Comparisons, levels = c("CRC-LC", "CRC-MELAN", "CRC-BRCA", "CRC-ALL"))
colors <- col_vector[11:14]
names(colors) = ann$Comparisons

colours <- list('Comparisons' = colors)
colAnn <- HeatmapAnnotation(df = ann, which = "row", col = colours, annotation_name_gp= gpar(fontsize = 15), show_annotation_name = FALSE,
                            annotation_legend_param = list(labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize = 15)))

max(mat1)
min(mat1)
col_fun = colorRamp2(c(-0.4, 0, 0.72), c("green", "white", "red"))
gra <- Heatmap(mat1, cell_fun = cell_fun, col = col_fun, name = "PropChange", right_annotation = colAnn,
               cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE)
gra

pdf(file.path(DIR_CANCER_TYPE_ARM, "cancer_type_ArmAmplification_heatmap.pdf"), height = 3, width = 12)
print(gra)
dev.off()


######################################################
## Deletion
LC <- read.csv(file.path(.csv_dir, "prop_DEL_CRC_LC.csv"))
LC <- LC[-1, ]
rownames(LC) <- LC[, 1]
LC <- LC[, -1]

MELAN <- read.csv(file.path(.csv_dir, "prop_DEL_CRC_MELAN.csv"))
MELAN <- MELAN[-1, ]
rownames(MELAN) <- MELAN[, 1]
MELAN <- MELAN[, -1]

BRCA <- read.csv(file.path(.csv_dir, "prop_DEL_CRC_BRCA.csv"))
BRCA <- BRCA[-1, ]
rownames(BRCA) <- BRCA[, 1]
BRCA <- BRCA[, -1]

ALL <- read.csv(file.path(.csv_dir, "prop_DEL_CRC_ALL.csv"))
ALL <- ALL[-1, ]
rownames(ALL) <- ALL[, 1]
ALL <- ALL[, -1]

mat1 <- data.frame(LC = LC$Prop_CRC - LC$Prop_LC, 
                   MELAN = MELAN$Prop_CRC - MELAN$Prop_MELAN, 
                   BRCA = BRCA$Prop_CRC - BRCA$Prop_BRCA,
                   ALL = ALL$Prop_CRC - ALL$Prop_ALL, 
                   row.names = rownames(LC))
colnames(mat1) <- c("CRC-LC", "CRC-MELAN", "CRC-BRCA", "CRC-ALL")
# mat1[is.na(mat1)] <- 0

mat2 <- data.frame(LC = p.adjust(LC$pval, method = "fdr"),
                   MELAN = p.adjust(MELAN$pval, method = "fdr"),
                   BRCA = p.adjust(BRCA$pval, method = "fdr"),
                   ALL = p.adjust(ALL$pval, method = "fdr"),
                   row.names = rownames(LC))
mat2[is.na(mat2)] <- 1

mat1 <- t(as.matrix(mat1))
mat2 <- t(as.matrix(mat2))

colnames(mat1) <- gsub("Arm_", "", colnames(mat1))
colnames(mat2) <- gsub("Arm_", "", colnames(mat2))


cell_fun <- function(j, i, x, y, w, h, fill) {
  if(mat2[i, j] < 0.01) {
    grid.text("**", x, y)
  } else if(mat2[i, j] < 0.05) {
    grid.text("*", x, y)
  } else if(mat2[i, j] < 0.001) {
    grid.text("***", x, y)
  }
}

ann <- data.frame(Comparisons = rownames(mat1))
ann$Comparisons <- factor(ann$Comparisons, levels = c("CRC-LC", "CRC-MELAN", "CRC-BRCA", "CRC-ALL"))
colors <- col_vector[11:14]
names(colors) = ann$Comparisons

colours <- list('Comparisons' = colors)
colAnn <- HeatmapAnnotation(df = ann, which = "row", col = colours, annotation_name_gp= gpar(fontsize = 15), show_annotation_name = FALSE,
                            annotation_legend_param = list(labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize = 15)))

max(mat1)
min(mat1)

col_fun = colorRamp2(c(-0.5, 0, 0.5), c("green", "white", "#00008B"))
gra <- Heatmap(mat1, cell_fun = cell_fun, col = col_fun, name = "PropChange", right_annotation = colAnn,
               cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE)
gra

pdf(file.path(DIR_CANCER_TYPE_ARM, "cancer_type_ArmDeletion_heatmap.pdf"), height = 3, width = 12)
print(gra)
dev.off()

















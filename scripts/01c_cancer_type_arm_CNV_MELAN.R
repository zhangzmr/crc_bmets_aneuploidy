# MELAN vs Other Cancer Types: Arm-Level CNV Proportion Analysis + Heatmaps
# Reference cancer type: Melanoma (labelled 'MELAN')
#
# Cohorts (all brain metastasis samples):
#   MELAN = Melanoma                    (MSK-MET, n=63)
#   CRC   = Colorectal Cancer           (Integrated CRC dataset, n=70)
#   BRCA  = Breast Cancer               (MSK-MET, n=41)
#   LC    = Non-Small Cell Lung Cancer  (MSK-MET, n=196)
#
# Comparisons (prop.test, one-sided: MELAN proportion > comparator):
#   AMP: MELAN vs CRC | MELAN vs BRCA | MELAN vs LC | MELAN vs ALL (CRC+BRCA+LC)
#   DEL: MELAN vs CRC | MELAN vs BRCA | MELAN vs LC | MELAN vs ALL (CRC+BRCA+LC)
#
# Outputs (saved to cancer_type_arm/):
#   - prop_AMP_MELAN_*.csv / prop_DEL_MELAN_*.csv
#   - MELAN_arm_proptest_all.csv
#   - MELAN_arm_proptest_significant.csv
#   - MELAN_cancer_type_ArmAmplification_heatmap.pdf
#   - MELAN_cancer_type_ArmDeletion_heatmap.pdf

# ============================================================
# Libraries
# ============================================================


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

library(readxl)
library(writexl)
library(stringr)
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

dir.create(DIR_CANCER_TYPE_ARM, showWarnings = FALSE, recursive = TRUE)

REF <- "MELAN"

# ============================================================
# Load clinical data
# ============================================================
msk_clinical <- read.table("data/msk_met_2021_clinical_data.tsv", header = TRUE, sep = '\t', check.names = FALSE)
colnames(msk_clinical) <- str_replace_all(colnames(msk_clinical), ': ', '_')
colnames(msk_clinical) <- str_replace_all(colnames(msk_clinical), "\\(", "")
colnames(msk_clinical) <- str_replace_all(colnames(msk_clinical), "\\)", "")
colnames(msk_clinical) <- str_replace_all(colnames(msk_clinical), "\\/", "")
colnames(msk_clinical) <- str_replace_all(colnames(msk_clinical), ' ', '_')

primary_samples <- msk_clinical$Patient_ID[intersect(which(msk_clinical$Metastatic_patient == FALSE),
                                                      which(is.na(msk_clinical$Metastatic_Site)))]
Met_site_CNSBrain <- subset(msk_clinical, Metastatic_Site == 'CNS/Brain')
Met_site_Other    <- subset(msk_clinical, !is.na(Metastatic_Site) & !(Patient_ID %in%
                              msk_clinical$Patient_ID[msk_clinical$Distant_Mets_CNSBrain == 'Yes']))
primary           <- subset(msk_clinical, Patient_ID %in% primary_samples)

msk_clinical <- rbind(primary, Met_site_CNSBrain, Met_site_Other)
met_unspecified_sample <- subset(msk_clinical, Metastatic_Site == 'Unspecified')$Patient_ID
msk_clinical <- msk_clinical[!msk_clinical$Patient_ID %in% met_unspecified_sample, ]

msk_clinical$BrainMet <- as.factor(
  ifelse(msk_clinical$Patient_ID %in% primary_samples, 'Primary',
  ifelse(msk_clinical$Metastatic_Site == 'CNS/Brain', 'BrainMet', 'OtherMets')))

cat("Clinical data loaded:", dim(msk_clinical), "\n")

# ============================================================
# Load arm-level CNV data for each cancer type
# ============================================================

load_arm_cnv <- function(fname, sample_meta) {
  arm_cnv <- read.table(fname, header = TRUE)
  colnames(arm_cnv)[1] <- 'Sample_ID'
  colnames(arm_cnv) <- str_replace(colnames(arm_cnv), 'X', 'Arm_')
  nc <- ncol(arm_cnv)
  colnames(arm_cnv)[c(nc-1, nc)] <- c('Arm_Xp', 'Arm_Xq')

  colorder <- 'Sample_ID'
  for (i in 1:22)
    colorder <- c(colorder, paste0('Arm_', i, 'p'), paste0('Arm_', i, 'q'))
  colorder <- c(colorder, 'Arm_Xp', 'Arm_Xq')
  colorder <- intersect(colorder, colnames(arm_cnv))
  colorder <- setdiff(colorder, c('Arm_13p', 'Arm_14p', 'Arm_15p', 'Arm_21p', 'Arm_22p'))

  arm_cnv_ordered <- arm_cnv[, colorder]
  rownames(arm_cnv_ordered) <- arm_cnv_ordered$Sample_ID
  merge(arm_cnv_ordered, sample_meta)
}

### arm_cnv_mCRC_BrM  (Integrated CRC data)
SampleMeta_CRC <- read.csv('data/ASCETS/ascets-master/IntegratedData/Sample_Metadata.csv')
SampleMeta_CRC$X <- NULL

arm_cnv_CRC <- read.table("data/ASCETS/ascets-master/IntegratedData/IntegratedData_arm_level_calls.txt",
                           header = TRUE)
colnames(arm_cnv_CRC)[1] <- 'Sample_ID'
colnames(arm_cnv_CRC) <- str_replace(colnames(arm_cnv_CRC), 'X', 'Arm_')
nc <- ncol(arm_cnv_CRC); colnames(arm_cnv_CRC)[(nc-3):nc] <- c('Arm_Xp', 'Arm_Xq', 'Arm_Yp', 'Arm_Yq')
colorder_crc <- 'Sample_ID'
for (i in 1:22)
  colorder_crc <- c(colorder_crc, paste0('Arm_', i, 'p'), paste0('Arm_', i, 'q'))
colorder_crc <- c(colorder_crc, 'Arm_Xp', 'Arm_Xq', 'Arm_Yp', 'Arm_Yq')
colorder_crc <- setdiff(colorder_crc, c('Arm_13p'))
arm_cnv_CRC <- arm_cnv_CRC[, colorder_crc]
rownames(arm_cnv_CRC) <- arm_cnv_CRC$Sample_ID
arm_cnv_CRC <- merge(arm_cnv_CRC, SampleMeta_CRC)

JiQC_reject <- c('P10033_20283A', 'P10030_20281A')
arm_cnv_CRC <- subset(arm_cnv_CRC, !(Sample_ID %in% JiQC_reject) & Metastatic_Site == 'Brain')
arm_cnv_CRC$CancerType <- 'CRC'
arm_cnv_CRC$Subtype    <- 'CRC'
cols_to_drop <- c('Arm_14p', 'Arm_15p', 'Arm_21p', 'Arm_22p', 'Arm_Yp', 'Arm_Yq')
arm_cnv_CRC <- arm_cnv_CRC %>% select(-any_of(cols_to_drop))
arm_cnv_mCRC_BrM <- arm_cnv_CRC
cat("CRC brain mets:", nrow(arm_cnv_mCRC_BrM), "samples\n")

### arm_cnv_mLC_BrM
SampleMeta_LC <- subset(msk_clinical, Cancer_Type == 'Non-Small Cell Lung Cancer' & BrainMet == "BrainMet")
SampleMeta_LC$Dataset        <- 'MSK'
SampleMeta_LC$CancerType     <- 'LC'
SampleMeta_LC$Subtype        <- SampleMeta_LC$Subtype_Abbreviation
SampleMeta_LC$SampleType     <- SampleMeta_LC$Sample_Type
SampleMeta_LC$Metastatic_Site <- ifelse(SampleMeta_LC$Metastatic_Site == 'CNS/Brain', 'Brain',
                                         SampleMeta_LC$Metastatic_Site)
SampleMeta_LC <- SampleMeta_LC[, c('Sample_ID', 'Patient_ID', 'SampleType',
                                    'Metastatic_Site', 'Dataset', 'CancerType', 'Subtype')]

arm_cnv_mLC_BrM <- load_arm_cnv("data/ASCETS/ascets-master/MSK_NSCLC_Data/Output/LC_arm_level_calls.txt",
                                  SampleMeta_LC)
cat("LC brain mets:", nrow(arm_cnv_mLC_BrM), "samples\n")

### arm_cnv_mMELAN_BrM
SampleMeta_MELAN <- subset(msk_clinical, Cancer_Type == 'Melanoma' & BrainMet == "BrainMet")
SampleMeta_MELAN$Dataset        <- 'MSK'
SampleMeta_MELAN$CancerType     <- 'MELAN'
SampleMeta_MELAN$Subtype        <- SampleMeta_MELAN$Subtype_Abbreviation
SampleMeta_MELAN$SampleType     <- SampleMeta_MELAN$Sample_Type
SampleMeta_MELAN$Metastatic_Site <- ifelse(SampleMeta_MELAN$Metastatic_Site == 'CNS/Brain', 'Brain',
                                            SampleMeta_MELAN$Metastatic_Site)
SampleMeta_MELAN <- SampleMeta_MELAN[, c('Sample_ID', 'Patient_ID', 'SampleType',
                                          'Metastatic_Site', 'Dataset', 'CancerType', 'Subtype')]

arm_cnv_mMELAN_BrM <- load_arm_cnv("data/ASCETS/ascets-master/MSK_MELAN_Data/Output/MELAN_arm_level_calls.txt",
                                     SampleMeta_MELAN)
cat("MELAN brain mets:", nrow(arm_cnv_mMELAN_BrM), "samples\n")

### arm_cnv_mBRCA_BrM
SampleMeta_BRCA <- subset(msk_clinical, Cancer_Type == 'Breast Cancer' & BrainMet == "BrainMet")
SampleMeta_BRCA$Dataset        <- 'MSK'
SampleMeta_BRCA$CancerType     <- 'BRCA'
SampleMeta_BRCA$Subtype        <- SampleMeta_BRCA$Subtype_Abbreviation
SampleMeta_BRCA$SampleType     <- SampleMeta_BRCA$Sample_Type
SampleMeta_BRCA$Metastatic_Site <- ifelse(SampleMeta_BRCA$Metastatic_Site == 'CNS/Brain', 'Brain',
                                           SampleMeta_BRCA$Metastatic_Site)
SampleMeta_BRCA <- SampleMeta_BRCA[, c('Sample_ID', 'Patient_ID', 'SampleType',
                                        'Metastatic_Site', 'Dataset', 'CancerType', 'Subtype')]

arm_cnv_mBRCA_BrM <- load_arm_cnv("data/ASCETS/ascets-master/MSK_BRCA_Data/Output/BRCA_arm_level_calls.txt",
                                    SampleMeta_BRCA)
cat("BRCA brain mets:", nrow(arm_cnv_mBRCA_BrM), "samples\n")

### Combine all four
arm_cnv_DF_saved <- rbind(arm_cnv_mCRC_BrM,
                           arm_cnv_mLC_BrM,
                           arm_cnv_mMELAN_BrM,
                           arm_cnv_mBRCA_BrM)
cat("Combined data:", nrow(arm_cnv_DF_saved), "samples,", length(unique(arm_cnv_DF_saved$Patient_ID)), "patients\n")
cat("Cancer types:", unique(arm_cnv_DF_saved$CancerType), "\n")

arm_cols <- colnames(arm_cnv_DF_saved)
arm_cols <- arm_cols[grepl("^Arm_", arm_cols)]

# ============================================================
# Helper: run proportion test for one CNV type (AMP or DEL)
# ============================================================
run_prop_test <- function(arm_cnv_DF, sites_interest, cnv_type = "AMP") {
  met_n <- c(
    length(which(arm_cnv_DF$CancerType == sites_interest[1])),
    length(which(arm_cnv_DF$CancerType == sites_interest[2]))
  )
  prop_mat <- matrix(NA, length(arm_cols) + 1, length(sites_interest) + 4)
  rownames(prop_mat) <- c('Total_n', arm_cols)
  colnames(prop_mat) <- c(sites_interest, paste0('Prop_', sites_interest), 'statistic', 'pval')
  prop_mat[1, ] <- c(met_n, rep(NA, 4))
  ntest   <- length(arm_cols)
  confInt <- rep('', nrow(prop_mat))

  for (j in seq_along(arm_cols)) {
    for (k in 1:2) {
      site_sub    <- subset(arm_cnv_DF, CancerType == sites_interest[k])
      aneu_count  <- length(which(site_sub[, arm_cols[j]] == cnv_type))
      prop_mat[arm_cols[j], sites_interest[k]]               <- aneu_count
      prop_mat[arm_cols[j], paste0('Prop_', sites_interest[k])] <- round(aneu_count / met_n[k], 3)
    }
    if (sum(prop_mat[arm_cols[j], c(1, 2)], na.rm = TRUE) > 0) {
      propTest <- prop.test(x = as.numeric(prop_mat[arm_cols[j], 1:2]),
                            n = as.numeric(prop_mat[1, 1:2]),
                            alternative = "greater")
      prop_mat[arm_cols[j], 'statistic'] <- round(propTest$statistic, 3)
      prop_mat[arm_cols[j], 'pval']      <- propTest$p.value
      confInt[j + 1] <- paste0('[ ', round(propTest$conf.int[1], 3),
                                ' - ', round(propTest$conf.int[2], 3), ' ]')
    }
  }

  prop_df           <- as.data.frame(prop_mat)
  prop_df$confInt   <- confInt
  p_th_05  <- 0.05 / ntest
  p_th_01  <- 0.01 / ntest
  p_th_001 <- 0.001 / ntest
  prop_df$padj_sign <- ''
  prop_df$padj_sign <- ifelse(prop_df$pval < p_th_05,  '*',  prop_df$padj_sign)
  prop_df$padj_sign <- ifelse(prop_df$pval < p_th_01,  '**', prop_df$padj_sign)
  prop_df$padj_sign <- ifelse(prop_df$pval < p_th_001, '***', prop_df$padj_sign)
  prop_df$padj_thr  <- p_th_001
  prop_df
}

# ============================================================
# Statistical tests: MELAN as reference vs CRC, BRCA, LC, ALL
# ============================================================

# --- AMP: MELAN vs CRC ---
sites_interest <- c('MELAN', 'CRC')
arm_cnv_DF     <- subset(arm_cnv_DF_saved, CancerType %in% sites_interest)
prop_AMP_MELAN_CRC <- run_prop_test(arm_cnv_DF, sites_interest, "AMP")
write.csv(prop_AMP_MELAN_CRC, file.path(DIR_CANCER_TYPE_ARM, "prop_AMP_MELAN_CRC.csv"), row.names = TRUE)
cat("Saved: prop_AMP_MELAN_CRC.csv\n")

# --- AMP: MELAN vs BRCA ---
sites_interest  <- c('MELAN', 'BRCA')
arm_cnv_DF      <- subset(arm_cnv_DF_saved, CancerType %in% sites_interest)
prop_AMP_MELAN_BRCA <- run_prop_test(arm_cnv_DF, sites_interest, "AMP")
write.csv(prop_AMP_MELAN_BRCA, file.path(DIR_CANCER_TYPE_ARM, "prop_AMP_MELAN_BRCA.csv"), row.names = TRUE)
cat("Saved: prop_AMP_MELAN_BRCA.csv\n")

# --- AMP: MELAN vs LC (NSCLC) ---
sites_interest  <- c('MELAN', 'LC')
arm_cnv_DF      <- subset(arm_cnv_DF_saved, CancerType %in% sites_interest)
prop_AMP_MELAN_LC <- run_prop_test(arm_cnv_DF, sites_interest, "AMP")
write.csv(prop_AMP_MELAN_LC, file.path(DIR_CANCER_TYPE_ARM, "prop_AMP_MELAN_LC.csv"), row.names = TRUE)
cat("Saved: prop_AMP_MELAN_LC.csv\n")

# --- AMP: MELAN vs ALL (CRC + BRCA + LC combined) ---
arm_cnv_DF_temp <- arm_cnv_DF_saved
arm_cnv_DF_temp$CancerType <- ifelse(arm_cnv_DF_temp$CancerType == 'MELAN', 'MELAN', 'ALL')
sites_interest  <- c('MELAN', 'ALL')
prop_AMP_MELAN_ALL <- run_prop_test(arm_cnv_DF_temp, sites_interest, "AMP")
write.csv(prop_AMP_MELAN_ALL, file.path(DIR_CANCER_TYPE_ARM, "prop_AMP_MELAN_ALL.csv"), row.names = TRUE)
cat("Saved: prop_AMP_MELAN_ALL.csv\n")

# --- DEL: MELAN vs CRC ---
sites_interest <- c('MELAN', 'CRC')
arm_cnv_DF     <- subset(arm_cnv_DF_saved, CancerType %in% sites_interest)
prop_DEL_MELAN_CRC <- run_prop_test(arm_cnv_DF, sites_interest, "DEL")
write.csv(prop_DEL_MELAN_CRC, file.path(DIR_CANCER_TYPE_ARM, "prop_DEL_MELAN_CRC.csv"), row.names = TRUE)
cat("Saved: prop_DEL_MELAN_CRC.csv\n")

# --- DEL: MELAN vs BRCA ---
sites_interest  <- c('MELAN', 'BRCA')
arm_cnv_DF      <- subset(arm_cnv_DF_saved, CancerType %in% sites_interest)
prop_DEL_MELAN_BRCA <- run_prop_test(arm_cnv_DF, sites_interest, "DEL")
write.csv(prop_DEL_MELAN_BRCA, file.path(DIR_CANCER_TYPE_ARM, "prop_DEL_MELAN_BRCA.csv"), row.names = TRUE)
cat("Saved: prop_DEL_MELAN_BRCA.csv\n")

# --- DEL: MELAN vs LC ---
sites_interest  <- c('MELAN', 'LC')
arm_cnv_DF      <- subset(arm_cnv_DF_saved, CancerType %in% sites_interest)
prop_DEL_MELAN_LC <- run_prop_test(arm_cnv_DF, sites_interest, "DEL")
write.csv(prop_DEL_MELAN_LC, file.path(DIR_CANCER_TYPE_ARM, "prop_DEL_MELAN_LC.csv"), row.names = TRUE)
cat("Saved: prop_DEL_MELAN_LC.csv\n")

# --- DEL: MELAN vs ALL ---
arm_cnv_DF_temp <- arm_cnv_DF_saved
arm_cnv_DF_temp$CancerType <- ifelse(arm_cnv_DF_temp$CancerType == 'MELAN', 'MELAN', 'ALL')
sites_interest  <- c('MELAN', 'ALL')
prop_DEL_MELAN_ALL <- run_prop_test(arm_cnv_DF_temp, sites_interest, "DEL")
write.csv(prop_DEL_MELAN_ALL, file.path(DIR_CANCER_TYPE_ARM, "prop_DEL_MELAN_ALL.csv"), row.names = TRUE)
cat("Saved: prop_DEL_MELAN_ALL.csv\n")

# ============================================================
# Helper: build heatmap matrices from a list of prop tables
# ============================================================
make_heatmap_matrices <- function(prop_list, comp_labels, ref = REF) {
  arms <- rownames(prop_list[[1]])

  prop_col <- paste0("Prop_", ref)
  mat1 <- do.call(data.frame, lapply(seq_along(prop_list), function(i) {
    df  <- prop_list[[i]]
    comparator_col <- setdiff(grep("^Prop_", colnames(df), value = TRUE), prop_col)
    as.numeric(df[arms, prop_col]) - as.numeric(df[arms, comparator_col])
  }))
  rownames(mat1) <- arms
  colnames(mat1) <- comp_labels

  mat2 <- do.call(data.frame, lapply(prop_list, function(df) {
    p.adjust(as.numeric(df[arms, "pval"]), method = "fdr")
  }))
  rownames(mat2) <- arms
  colnames(mat2) <- comp_labels
  mat2[is.na(mat2)] <- 1

  list(mat1 = t(as.matrix(mat1)), mat2 = t(as.matrix(mat2)))
}

strip_arm <- function(mat) {
  colnames(mat) <- gsub("Arm_", "", colnames(mat))
  mat
}

cell_fun_factory <- function(mat2_env) {
  function(j, i, x, y, w, h, fill) {
    v <- mat2_env[i, j]
    if (!is.na(v)) {
      if      (v < 0.001) grid.text("***", x, y)
      else if (v < 0.01)  grid.text("**",  x, y)
      else if (v < 0.05)  grid.text("*",   x, y)
    }
  }
}

make_row_annotation <- function(comp_labels, col_indices) {
  ann <- data.frame(Comparisons = factor(comp_labels, levels = comp_labels))
  colors <- col_vector[col_indices]
  names(colors) <- comp_labels
  HeatmapAnnotation(
    df  = ann, which = "row",
    col = list(Comparisons = colors),
    annotation_name_gp = gpar(fontsize = 15),
    show_annotation_name = FALSE,
    annotation_legend_param = list(labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize = 15))
  )
}

# ============================================================
# Heatmap: Amplifications
# ============================================================
comp_labels <- c("MELAN-CRC", "MELAN-BRCA", "MELAN-LC", "MELAN-ALL")

amp_list <- list(
  subset(prop_AMP_MELAN_CRC,  rownames(prop_AMP_MELAN_CRC)  != 'Total_n'),
  subset(prop_AMP_MELAN_BRCA, rownames(prop_AMP_MELAN_BRCA) != 'Total_n'),
  subset(prop_AMP_MELAN_LC,   rownames(prop_AMP_MELAN_LC)   != 'Total_n'),
  subset(prop_AMP_MELAN_ALL,  rownames(prop_AMP_MELAN_ALL)  != 'Total_n')
)

mats_amp <- make_heatmap_matrices(amp_list, comp_labels, ref = REF)
mat1_amp <- strip_arm(mats_amp$mat1)
mat2_amp <- strip_arm(mats_amp$mat2)

col_fun_amp <- colorRamp2(c(-0.4, 0, 0.72), c("green", "white", "red"))
colAnn_amp  <- make_row_annotation(comp_labels, col_indices = 15:18)

gra_amp <- Heatmap(mat1_amp,
                   cell_fun          = cell_fun_factory(mat2_amp),
                   col               = col_fun_amp,
                   name              = "PropChange",
                   right_annotation  = colAnn_amp,
                   cluster_rows      = FALSE,
                   cluster_columns   = FALSE,
                   show_row_names    = FALSE,
                   column_title      = "Chromosome Arm",
                   row_title         = "Comparison")
print(gra_amp)

pdf(file.path(DIR_CANCER_TYPE_ARM, "MELAN_cancer_type_ArmAmplification_heatmap.pdf"), height = 3, width = 12)
print(gra_amp)
dev.off()
cat("Saved: MELAN_cancer_type_ArmAmplification_heatmap.pdf\n")

# ============================================================
# Heatmap: Deletions
# ============================================================
del_list <- list(
  subset(prop_DEL_MELAN_CRC,  rownames(prop_DEL_MELAN_CRC)  != 'Total_n'),
  subset(prop_DEL_MELAN_BRCA, rownames(prop_DEL_MELAN_BRCA) != 'Total_n'),
  subset(prop_DEL_MELAN_LC,   rownames(prop_DEL_MELAN_LC)   != 'Total_n'),
  subset(prop_DEL_MELAN_ALL,  rownames(prop_DEL_MELAN_ALL)  != 'Total_n')
)

mats_del <- make_heatmap_matrices(del_list, comp_labels, ref = REF)
mat1_del <- strip_arm(mats_del$mat1)
mat2_del <- strip_arm(mats_del$mat2)

col_fun_del <- colorRamp2(c(-0.5, 0, 0.5), c("green", "white", "#00008B"))
colAnn_del  <- make_row_annotation(comp_labels, col_indices = 15:18)

gra_del <- Heatmap(mat1_del,
                   cell_fun          = cell_fun_factory(mat2_del),
                   col               = col_fun_del,
                   name              = "PropChange",
                   right_annotation  = colAnn_del,
                   cluster_rows      = FALSE,
                   cluster_columns   = FALSE,
                   show_row_names    = FALSE,
                   column_title      = "Chromosome Arm",
                   row_title         = "Comparison")
print(gra_del)

pdf(file.path(DIR_CANCER_TYPE_ARM, "MELAN_cancer_type_ArmDeletion_heatmap.pdf"), height = 3, width = 12)
print(gra_del)
dev.off()
cat("Saved: MELAN_cancer_type_ArmDeletion_heatmap.pdf\n")

# ============================================================
# Summary table: significant arms across all comparisons
# ============================================================
make_summary <- function(prop_df, comparison, cnv_type, ref = REF) {
  df <- subset(prop_df, rownames(prop_df) != 'Total_n')
  ref_col       <- paste0("Prop_", ref)
  comp_prop_col <- setdiff(grep("^Prop_", colnames(df), value = TRUE), ref_col)
  comparator    <- sub("Prop_", "", comp_prop_col)

  data.frame(
    Arm              = gsub("Arm_", "", rownames(df)),
    CNV_Type         = cnv_type,
    Comparison       = comparison,
    N_Ref            = as.numeric(df[, ref]),
    N_Comparator     = as.numeric(df[, comparator]),
    Prop_Ref         = as.numeric(df[, ref_col]),
    Prop_Comparator  = as.numeric(df[, comp_prop_col]),
    Prop_Diff        = round(as.numeric(df[, ref_col]) - as.numeric(df[, comp_prop_col]), 3),
    Pval             = as.numeric(df$pval),
    Padj_FDR         = round(p.adjust(as.numeric(df$pval), method = "fdr"), 5),
    Significance     = df$padj_sign,
    stringsAsFactors = FALSE
  )
}

summary_all <- rbind(
  make_summary(prop_AMP_MELAN_CRC,  "MELAN-CRC",  "AMP"),
  make_summary(prop_AMP_MELAN_BRCA, "MELAN-BRCA", "AMP"),
  make_summary(prop_AMP_MELAN_LC,   "MELAN-LC",   "AMP"),
  make_summary(prop_AMP_MELAN_ALL,  "MELAN-ALL",  "AMP"),
  make_summary(prop_DEL_MELAN_CRC,  "MELAN-CRC",  "DEL"),
  make_summary(prop_DEL_MELAN_BRCA, "MELAN-BRCA", "DEL"),
  make_summary(prop_DEL_MELAN_LC,   "MELAN-LC",   "DEL"),
  make_summary(prop_DEL_MELAN_ALL,  "MELAN-ALL",  "DEL")
)

summary_sig <- subset(summary_all, Significance != '')
summary_sig <- summary_sig[order(summary_sig$CNV_Type, summary_sig$Comparison, summary_sig$Pval), ]
rownames(summary_sig) <- NULL

write.csv(summary_all, file.path(DIR_CANCER_TYPE_ARM, "MELAN_arm_proptest_all.csv"),        row.names = FALSE)
write.csv(summary_sig, file.path(DIR_CANCER_TYPE_ARM, "MELAN_arm_proptest_significant.csv"), row.names = FALSE)

cat("\nSummary table (all arms):", nrow(summary_all), "rows\n")
cat("Significant arms (Bonferroni):", nrow(summary_sig), "rows\n")
print(summary_sig[, c("Arm","CNV_Type","Comparison","Prop_Ref","Prop_Comparator","Prop_Diff","Pval","Significance")])
cat("Saved: MELAN_arm_proptest_all.csv and MELAN_arm_proptest_significant.csv\n")

cat("\nDone! All CSVs written to cancer_type_arm/ and heatmap PDFs saved.\n")

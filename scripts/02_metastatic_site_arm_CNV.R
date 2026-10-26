# Metastatic-site arm CNV: Brain vs Primary / Liver / Lung / Lymph / ALL (LLL)
# Within each cancer type (same design as CRC ArmAmplification/Deletion heatmaps)
#
# Cancer types:
#   CRC   = Integrated CRC (MSK + Curtis + Ji)
#   LC    = MSK Non-Small Cell Lung Cancer
#   MELAN = MSK Melanoma
#   BRCA  = MSK Breast Cancer
#
# Comparisons (prop.test, one-sided: Brain proportion > comparator):
#   AMP & DEL for each of:
#     Brain vs Primary | Brain vs Liver | Brain vs Lung | Brain vs Lymph | Brain vs ALL
#   ALL = Liver + Lung + Lymph (LLL_Mets)
#
# Outputs → metastatic_site_arm/  (metsite_ prefix to distinguish from cancer_type_arm/)
#   - {CANCER}_metsite_sample_counts.csv
#   - prop_AMP_metsite_{CANCER}_{Primary,Liver,Lung,Lymph,ALL}.csv
#   - prop_DEL_metsite_{CANCER}_{Primary,Liver,Lung,Lymph,ALL}.csv
#   - {CANCER}_metsite_arm_proptest_all.csv / _significant.csv
#   - {CANCER}_metsite_ArmAmplification_heatmap.pdf
#   - {CANCER}_metsite_ArmDeletion_heatmap.pdf



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

suppressPackageStartupMessages({
  library(stringr)
  library(tidyverse)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

outDir <- DIR_METSITE_ARM
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Clinical (for MSK LC / MELAN / BRCA site labels)
# ============================================================
msk_clinical <- read.table("data/msk_met_2021_clinical_data.tsv",
                           header = TRUE, sep = '\t', check.names = FALSE)
colnames(msk_clinical) <- str_replace_all(colnames(msk_clinical), ': ', '_')
colnames(msk_clinical) <- str_replace_all(colnames(msk_clinical), "\\(", "")
colnames(msk_clinical) <- str_replace_all(colnames(msk_clinical), "\\)", "")
colnames(msk_clinical) <- str_replace_all(colnames(msk_clinical), "\\/", "")
colnames(msk_clinical) <- str_replace_all(colnames(msk_clinical), ' ', '_')

primary_samples <- msk_clinical$Patient_ID[
  intersect(which(msk_clinical$Metastatic_patient == FALSE),
            which(is.na(msk_clinical$Metastatic_Site)))
]

msk_clinical$SiteLabel <- ifelse(
  msk_clinical$Patient_ID %in% primary_samples, 'Primary',
  ifelse(msk_clinical$Metastatic_Site == 'CNS/Brain', 'Brain',
         as.character(msk_clinical$Metastatic_Site))
)

# ============================================================
# Helpers to load arm CNV
# ============================================================
load_msk_arm <- function(fname, cancer_type) {
  arm_cnv <- read.table(fname, header = TRUE)
  colnames(arm_cnv)[1] <- 'Sample_ID'
  colnames(arm_cnv) <- str_replace(colnames(arm_cnv), 'X', 'Arm_')
  nc <- ncol(arm_cnv)
  colnames(arm_cnv)[c(nc - 1, nc)] <- c('Arm_Xp', 'Arm_Xq')

  colorder <- 'Sample_ID'
  for (i in 1:22)
    colorder <- c(colorder, paste0('Arm_', i, 'p'), paste0('Arm_', i, 'q'))
  colorder <- c(colorder, 'Arm_Xp', 'Arm_Xq')
  colorder <- intersect(colorder, colnames(arm_cnv))
  colorder <- setdiff(colorder, c('Arm_13p', 'Arm_14p', 'Arm_15p', 'Arm_21p', 'Arm_22p'))
  arm_cnv <- arm_cnv[, colorder]

  meta <- subset(msk_clinical, Cancer_Type == cancer_type)
  meta <- meta[, c('Sample_ID', 'Patient_ID', 'SiteLabel')]
  colnames(meta)[colnames(meta) == 'SiteLabel'] <- 'Metastatic_Site'
  merge(arm_cnv, meta, by = 'Sample_ID')
}

load_crc_integrated <- function() {
  SampleMeta_CRC <- read.csv('data/ASCETS/ascets-master/IntegratedData/Sample_Metadata.csv')
  SampleMeta_CRC$X <- NULL

  arm_cnv_CRC <- read.table(
    "data/ASCETS/ascets-master/IntegratedData/IntegratedData_arm_level_calls.txt",
    header = TRUE
  )
  colnames(arm_cnv_CRC)[1] <- 'Sample_ID'
  colnames(arm_cnv_CRC) <- str_replace(colnames(arm_cnv_CRC), 'X', 'Arm_')
  nc <- ncol(arm_cnv_CRC); colnames(arm_cnv_CRC)[(nc-3):nc] <- c('Arm_Xp', 'Arm_Xq', 'Arm_Yp', 'Arm_Yq')

  colorder_crc <- 'Sample_ID'
  for (i in 1:22)
    colorder_crc <- c(colorder_crc, paste0('Arm_', i, 'p'), paste0('Arm_', i, 'q'))
  colorder_crc <- c(colorder_crc, 'Arm_Xp', 'Arm_Xq', 'Arm_Yp', 'Arm_Yq')
  colorder_crc <- setdiff(colorder_crc, c('Arm_13p'))
  arm_cnv_CRC <- arm_cnv_CRC[, colorder_crc]
  arm_cnv_CRC <- merge(arm_cnv_CRC, SampleMeta_CRC, by = 'Sample_ID')

  JiQC_reject <- c('P10033_20283A', 'P10030_20281A')
  arm_cnv_CRC <- subset(arm_cnv_CRC, !(Sample_ID %in% JiQC_reject))
  cols_to_drop <- c('Arm_14p', 'Arm_15p', 'Arm_21p', 'Arm_22p', 'Arm_Yp', 'Arm_Yq')
  arm_cnv_CRC <- arm_cnv_CRC %>% select(-any_of(cols_to_drop))
  arm_cnv_CRC
}

# ============================================================
# Proportion test: Brain vs one comparator
# ============================================================
run_site_prop_test <- function(arm_cnv_DF, arm_cols, comparator, cnv_type = "AMP") {
  # comparator is Primary / Liver / Lung / Lymph / ALL
  df <- arm_cnv_DF
  if (comparator == "ALL") {
    df$SiteCmp <- ifelse(df$Metastatic_Site == "Brain", "Brain",
                         ifelse(df$Metastatic_Site %in% c("Liver", "Lung", "Lymph"), "ALL", NA))
    df <- subset(df, !is.na(SiteCmp))
    sites_interest <- c("Brain", "ALL")
    site_col <- "SiteCmp"
  } else {
    df <- subset(df, Metastatic_Site %in% c("Brain", comparator))
    sites_interest <- c("Brain", comparator)
    site_col <- "Metastatic_Site"
  }

  met_n <- c(
    length(which(df[[site_col]] == sites_interest[1])),
    length(which(df[[site_col]] == sites_interest[2]))
  )

  prop_mat <- matrix(NA, length(arm_cols) + 1, length(sites_interest) + 4)
  rownames(prop_mat) <- c("Total_n", arm_cols)
  colnames(prop_mat) <- c(sites_interest, paste0("Prop_", sites_interest), "statistic", "pval")
  prop_mat[1, ] <- c(met_n, rep(NA, 4))
  ntest <- length(arm_cols)
  confInt <- rep("", nrow(prop_mat))

  for (j in seq_along(arm_cols)) {
    for (k in 1:2) {
      site_sub <- df[df[[site_col]] == sites_interest[k], , drop = FALSE]
      aneu_count <- length(which(site_sub[[arm_cols[j]]] == cnv_type))
      prop_mat[arm_cols[j], sites_interest[k]] <- aneu_count
      prop_mat[arm_cols[j], paste0("Prop_", sites_interest[k])] <-
        if (met_n[k] > 0) round(aneu_count / met_n[k], 3) else NA
    }
    if (sum(prop_mat[arm_cols[j], c(1, 2)], na.rm = TRUE) > 0 && all(met_n > 0)) {
      propTest <- prop.test(
        x = as.numeric(prop_mat[arm_cols[j], 1:2]),
        n = as.numeric(prop_mat[1, 1:2]),
        alternative = "greater"
      )
      prop_mat[arm_cols[j], "statistic"] <- round(propTest$statistic, 3)
      prop_mat[arm_cols[j], "pval"] <- propTest$p.value
      confInt[j + 1] <- paste0(
        "[ ", round(propTest$conf.int[1], 3), " - ",
        round(propTest$conf.int[2], 3), " ]"
      )
    }
  }

  prop_df <- as.data.frame(prop_mat)
  prop_df$confInt <- confInt
  p_th_05  <- 0.05 / ntest
  p_th_01  <- 0.01 / ntest
  p_th_001 <- 0.001 / ntest
  prop_df$padj_sign <- ""
  prop_df$padj_sign <- ifelse(prop_df$pval < p_th_05,  "*",   prop_df$padj_sign)
  prop_df$padj_sign <- ifelse(prop_df$pval < p_th_01,  "**",  prop_df$padj_sign)
  prop_df$padj_sign <- ifelse(prop_df$pval < p_th_001, "***", prop_df$padj_sign)
  prop_df$padj_thr <- p_th_001
  prop_df
}

# ============================================================
# Heatmap helpers
# ============================================================
make_heatmap_matrices <- function(prop_list, comp_labels) {
  arms <- rownames(prop_list[[1]])
  mat1 <- do.call(data.frame, lapply(seq_along(prop_list), function(i) {
    df <- prop_list[[i]]
    as.numeric(df[arms, "Prop_Brain"]) - as.numeric(df[arms, grep("^Prop_", colnames(df), value = TRUE)[2]])
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
    df = ann, which = "row",
    col = list(Comparisons = colors),
    annotation_name_gp = gpar(fontsize = 15),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      labels_gp = gpar(fontsize = 15),
      title_gp = gpar(fontsize = 15)
    )
  )
}

make_summary <- function(prop_df, comparison, cnv_type) {
  df <- subset(prop_df, rownames(prop_df) != "Total_n")
  prop_cols <- grep("^Prop_", colnames(df), value = TRUE)
  ref_col <- "Prop_Brain"
  comp_prop_col <- setdiff(prop_cols, ref_col)
  comparator <- sub("Prop_", "", comp_prop_col)
  count_cols <- setdiff(colnames(df)[1:2], c("statistic", "pval"))
  # first two data cols are Brain + comparator counts
  data.frame(
    Arm = gsub("Arm_", "", rownames(df)),
    CNV_Type = cnv_type,
    Comparison = comparison,
    N_Brain = as.numeric(df[, "Brain"]),
    N_Comparator = as.numeric(df[, comparator]),
    Prop_Brain = as.numeric(df[, ref_col]),
    Prop_Comparator = as.numeric(df[, comp_prop_col]),
    Prop_Diff = round(as.numeric(df[, ref_col]) - as.numeric(df[, comp_prop_col]), 3),
    Pval = as.numeric(df$pval),
    Padj_FDR = round(p.adjust(as.numeric(df$pval), method = "fdr"), 5),
    Significance = df$padj_sign,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Run one cancer type end-to-end
# ============================================================
run_cancer <- function(cancer_code, arm_cnv_DF, col_indices) {
  cat("\n====================", cancer_code, "====================\n")

  keep_sites <- c("Brain", "Primary", "Liver", "Lung", "Lymph")
  arm_cnv_DF <- subset(arm_cnv_DF, Metastatic_Site %in% keep_sites)

  counts <- as.data.frame(table(arm_cnv_DF$Metastatic_Site), stringsAsFactors = FALSE)
  colnames(counts) <- c("Metastatic_Site", "N")
  counts$Cancer <- cancer_code
  # ensure all sites present
  for (s in keep_sites) {
    if (!s %in% counts$Metastatic_Site)
      counts <- rbind(counts, data.frame(Metastatic_Site = s, N = 0, Cancer = cancer_code))
  }
  counts <- counts[match(keep_sites, counts$Metastatic_Site), ]
  lll_n <- sum(counts$N[counts$Metastatic_Site %in% c("Liver", "Lung", "Lymph")])
  counts <- rbind(counts, data.frame(Metastatic_Site = "ALL_LLL", N = lll_n, Cancer = cancer_code))
  write.csv(counts, file.path(outDir, paste0(cancer_code, "_metsite_sample_counts.csv")), row.names = FALSE)
  print(counts)

  arm_cols <- colnames(arm_cnv_DF)[grepl("^Arm_", colnames(arm_cnv_DF))]
  comparators <- c("Primary", "Liver", "Lung", "Lymph", "ALL")
  heatmap_labels <- c("Brain-Primary", "Brain-Liver", "Brain-Lung", "Brain-Lymph", "Brain-ALL")

  amp_list <- list()
  del_list <- list()
  summary_rows <- list()

  for (comp in comparators) {
    prop_amp <- run_site_prop_test(arm_cnv_DF, arm_cols, comp, "AMP")
    prop_del <- run_site_prop_test(arm_cnv_DF, arm_cols, comp, "DEL")
    write.csv(prop_amp, file.path(outDir, paste0("prop_AMP_metsite_", cancer_code, "_", comp, ".csv")),
              row.names = TRUE)
    write.csv(prop_del, file.path(outDir, paste0("prop_DEL_metsite_", cancer_code, "_", comp, ".csv")),
              row.names = TRUE)
    cat("Saved prop_AMP/DEL_metsite_", cancer_code, "_", comp, ".csv\n", sep = "")

    amp_list[[comp]] <- subset(prop_amp, rownames(prop_amp) != "Total_n")
    del_list[[comp]] <- subset(prop_del, rownames(prop_del) != "Total_n")
    label <- paste0("Brain-", ifelse(comp == "ALL", "ALL", comp))
    summary_rows[[paste0("AMP_", comp)]] <- make_summary(prop_amp, label, "AMP")
    summary_rows[[paste0("DEL_", comp)]] <- make_summary(prop_del, label, "DEL")
  }

  # AMP heatmap
  mats_amp <- make_heatmap_matrices(amp_list[comparators], heatmap_labels)
  mat1_amp <- strip_arm(mats_amp$mat1)
  mat2_amp <- strip_arm(mats_amp$mat2)
  # dynamic color scale based on data
  amp_rng <- range(mat1_amp, na.rm = TRUE)
  col_fun_amp <- colorRamp2(
    c(min(amp_rng[1], -0.1), 0, max(amp_rng[2], 0.1)),
    c("green", "white", "red")
  )
  gra_amp <- Heatmap(
    mat1_amp,
    cell_fun = cell_fun_factory(mat2_amp),
    col = col_fun_amp,
    name = "PropChange",
    right_annotation = make_row_annotation(heatmap_labels, col_indices),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    column_title = paste0(cancer_code, " — Chromosome Arm (AMP)"),
    row_title = "Comparison"
  )
  pdf(file.path(outDir, paste0(cancer_code, "_metsite_ArmAmplification_heatmap.pdf")),
      height = 3, width = 12)
  print(gra_amp)
  dev.off()
  cat("Saved:", cancer_code, "_metsite_ArmAmplification_heatmap.pdf\n")

  # DEL heatmap
  mats_del <- make_heatmap_matrices(del_list[comparators], heatmap_labels)
  mat1_del <- strip_arm(mats_del$mat1)
  mat2_del <- strip_arm(mats_del$mat2)
  del_rng <- range(mat1_del, na.rm = TRUE)
  col_fun_del <- colorRamp2(
    c(min(del_rng[1], -0.1), 0, max(del_rng[2], 0.1)),
    c("green", "white", "#00008B")
  )
  gra_del <- Heatmap(
    mat1_del,
    cell_fun = cell_fun_factory(mat2_del),
    col = col_fun_del,
    name = "PropChange",
    right_annotation = make_row_annotation(heatmap_labels, col_indices),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    column_title = paste0(cancer_code, " — Chromosome Arm (DEL)"),
    row_title = "Comparison"
  )
  pdf(file.path(outDir, paste0(cancer_code, "_metsite_ArmDeletion_heatmap.pdf")),
      height = 3, width = 12)
  print(gra_del)
  dev.off()
  cat("Saved:", cancer_code, "_metsite_ArmDeletion_heatmap.pdf\n")

  summary_all <- do.call(rbind, summary_rows)
  rownames(summary_all) <- NULL
  summary_sig <- subset(summary_all, Significance != "")
  summary_sig <- summary_sig[order(summary_sig$CNV_Type, summary_sig$Comparison, summary_sig$Pval), ]
  rownames(summary_sig) <- NULL

  write.csv(summary_all, file.path(outDir, paste0(cancer_code, "_metsite_arm_proptest_all.csv")),
            row.names = FALSE)
  write.csv(summary_sig, file.path(outDir, paste0(cancer_code, "_metsite_arm_proptest_significant.csv")),
            row.names = FALSE)
  cat("Significant arms (Bonferroni):", nrow(summary_sig), "\n")
  if (nrow(summary_sig) > 0) {
    print(summary_sig[, c("Arm", "CNV_Type", "Comparison", "Prop_Brain",
                          "Prop_Comparator", "Prop_Diff", "Pval", "Significance")])
  }
  invisible(summary_all)
}

# ============================================================
# Load data and run
# ============================================================
cat("Loading CRC (Integrated)...\n")
crc_df <- load_crc_integrated()
cat("CRC samples:", nrow(crc_df), "\n")

cat("Loading LC...\n")
lc_df <- load_msk_arm(
  "data/ASCETS/ascets-master/MSK_NSCLC_Data/Output/LC_arm_level_calls.txt",
  "Non-Small Cell Lung Cancer"
)
cat("LC samples:", nrow(lc_df), "\n")

cat("Loading MELAN...\n")
melan_df <- load_msk_arm(
  "data/ASCETS/ascets-master/MSK_MELAN_Data/Output/MELAN_arm_level_calls.txt",
  "Melanoma"
)
cat("MELAN samples:", nrow(melan_df), "\n")

cat("Loading BRCA...\n")
brca_df <- load_msk_arm(
  "data/ASCETS/ascets-master/MSK_BRCA_Data/Output/BRCA_arm_level_calls.txt",
  "Breast Cancer"
)
cat("BRCA samples:", nrow(brca_df), "\n")

run_cancer("CRC",   crc_df,   col_indices = 6:10)
run_cancer("LC",    lc_df,    col_indices = 11:15)
run_cancer("MELAN", melan_df, col_indices = 16:20)
run_cancer("BRCA",  brca_df,  col_indices = 21:25)

cat("\nDone! All outputs written to", outDir, "/\n")

#!/usr/bin/env Rscript
# Rebuild maftools coBarplot: MSK mCRC Brain mets vs Primary
# Manuscript: Figure 5 (mutation comparison)
#
# Legend note: stock maftools::coBarplot allocates a short bottom panel and often
# clips Variant_Classification labels. This script uses a local wrapper with more
# legend space.

suppressPackageStartupMessages({
  library(maftools)
  library(dplyr)
  library(data.table)
})

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", grep("^--file=", .args, value = TRUE))
if (!length(.file)) stop("Run with Rscript so --file= is set.")
.script_dir <- dirname(normalizePath(.file))
REPO_ROOT <- normalizePath(file.path(.script_dir, ".."), mustWork = TRUE)
setwd(REPO_ROOT)

OUT_ROOT <- file.path(REPO_ROOT, "outputs")
DIR_MUTATION <- file.path(OUT_ROOT, "mutation_figure")
dir.create(DIR_MUTATION, showWarnings = FALSE, recursive = TRUE)
out_dir <- DIR_MUTATION

MSK_MUT_MAF <- file.path("data", "msk_met_2021", "data_mutations.maf")
MSK_CLINICAL_SAMPLE <- file.path("data", "msk_met_2021", "data_clinical_sample_renamed.txt")
SAMPLE_META <- file.path("data", "ASCETS", "ascets-master", "IntegratedData", "Sample_Metadata.csv")

stopifnot(file.exists(MSK_MUT_MAF), file.exists(MSK_CLINICAL_SAMPLE), file.exists(SAMPLE_META))

message("MAF:      ", MSK_MUT_MAF)
message("Clinical: ", MSK_CLINICAL_SAMPLE)
message("Meta:     ", SAMPLE_META)


# ---------------------------------------------------------------------------
# Local coBarplot with readable legend (based on maftools::coBarplot)
# ---------------------------------------------------------------------------
get_col_df <- getFromNamespace("get_col_df", "maftools")
get_vcColors <- getFromNamespace("get_vcColors", "maftools")

coBarplot_legible <- function(
  m1, m2, genes = NULL, orderBy = NULL,
  m1Name = "m1", m2Name = "m2",
  colors = NULL, normalize = TRUE, yLims = NULL,
  borderCol = "gray", titleSize = 1, geneSize = 0.8,
  showPct = TRUE, pctSize = 0.7, axisSize = 0.8,
  showLegend = TRUE, legendTxtSize = 0.85, geneMar = 5,
  legendPanelHeight = 2.2
) {
  if (is.null(genes)) {
    m1.genes <- getGeneSummary(m1)[1:5, Hugo_Symbol]
    m2.genes <- getGeneSummary(m2)[1:5, Hugo_Symbol]
    genes <- rev(unique(c(m1.genes, m2.genes)))
  }
  if (is.null(colors)) colors <- get_vcColors()

  m1.gs <- get_col_df(m = m1, g = genes)
  m1.ss <- as.numeric(m1@summary[ID %in% "Samples", summary])
  gs1.load <- getGeneSummary(x = m1)[Hugo_Symbol %in% genes, .(Hugo_Symbol, AlteredSamples)]
  gs1.load[, `:=`(AlteredSamples, round(AlteredSamples / m1.ss, digits = 2) * 100)]
  data.table::setDF(x = gs1.load, rownames = gs1.load$Hugo_Symbol)
  gs1.load <- gs1.load[, "AlteredSamples", drop = FALSE]
  m1.missing <- genes[!genes %in% rownames(gs1.load)]
  m1.missing <- m1.missing[!is.na(m1.missing)]
  if (length(m1.missing) > 0) {
    add1 <- data.frame(AlteredSamples = rep(0, length(m1.missing)))
    rownames(add1) <- m1.missing
    gs1.load <- rbind(gs1.load, add1)
  }

  m2.gs <- get_col_df(m = m2, g = genes)
  if (ncol(m1.gs) == 0L && ncol(m2.gs) == 0L) stop("Cannot plot genes without mutation!")
  if (ncol(m1.gs) == 0L) { m1.gs <- m2.gs; m1.gs[] <- 0 }
  if (ncol(m2.gs) == 0L) { m2.gs <- m1.gs; m2.gs[] <- 0 }

  m2.ss <- as.numeric(m2@summary[ID %in% "Samples", summary])
  gs2.load <- getGeneSummary(x = m2)[Hugo_Symbol %in% genes, .(Hugo_Symbol, AlteredSamples)]
  gs2.load[, `:=`(AlteredSamples, round(AlteredSamples / m2.ss, digits = 2) * 100)]
  data.table::setDF(x = gs2.load, rownames = gs2.load$Hugo_Symbol)
  gs2.load <- gs2.load[, "AlteredSamples", drop = FALSE]
  m2.missing <- genes[!genes %in% rownames(gs2.load)]
  m2.missing <- m2.missing[!is.na(m2.missing)]
  if (length(m2.missing) > 0) {
    add2 <- data.frame(AlteredSamples = rep(0, length(m2.missing)))
    rownames(add2) <- m2.missing
    gs2.load <- rbind(gs2.load, add2)
  }

  if (!is.null(orderBy)) {
    orderBy <- match.arg(orderBy, c("m1", "m2"))
    if (orderBy == "m1") {
      genes <- rev(rownames(gs1.load[order(gs1.load$AlteredSamples, decreasing = TRUE), , drop = FALSE]))
    } else {
      genes <- rev(rownames(gs2.load[order(gs2.load$AlteredSamples, decreasing = TRUE), , drop = FALSE]))
    }
  }

  m1.gs <- t(m1.gs[genes, , drop = FALSE])
  gs1.load <- gs1.load[genes, , drop = FALSE]
  m2.gs <- t(m2.gs[genes, , drop = FALSE])
  gs2.load <- gs2.load[genes, , drop = FALSE]

  if (normalize) {
    m1.gs <- (m1.gs / m1.ss) * 100
    m2.gs <- (m2.gs / m2.ss) * 100
    xat <- if (!is.null(yLims)) pretty(c(-yLims[1], 0, yLims[2])) else seq(-100, 100, 20)
  } else if (!is.null(yLims)) {
    xat <- pretty(c(-yLims[1], 0, yLims[2]))
  } else {
    m12max <- max(
      max(apply(m1.gs, 2, function(x) cumsum(x))),
      max(apply(m2.gs, 2, function(x) cumsum(x)))
    )
    xat <- pretty(c(-m12max, 0, m12max))
  }

  lo <- matrix(c(1, 1, 2, 2), nrow = 2, byrow = TRUE)
  if (showLegend) {
    graphics::layout(mat = lo, heights = c(4, legendPanelHeight))
  } else {
    graphics::layout(mat = lo, heights = c(20, 1))
  }

  par(mar = c(3, geneMar, 1.5, 2))
  b1 <- barplot(
    height = -m1.gs, horiz = TRUE, las = TRUE,
    col = colors[rownames(m1.gs)], xlim = range(xat),
    axes = FALSE, names.arg = rep(NA, ncol(m1.gs)), border = borderCol
  )
  if (showPct) {
    text(
      x = -colSums(m1.gs), y = b1,
      labels = paste0(gs1.load$AlteredSamples, "%"),
      adj = 1.1, cex = pctSize, xpd = TRUE
    )
  }
  barplot(
    height = m2.gs, horiz = TRUE, las = TRUE,
    col = colors[rownames(m2.gs)], add = TRUE, axes = FALSE,
    names.arg = rep(NA, ncol(m2.gs)), border = borderCol
  )
  if (showPct) {
    text(
      x = colSums(m2.gs), y = b1,
      labels = paste0(gs2.load$AlteredSamples, "%"),
      adj = -0.1, cex = pctSize, xpd = TRUE
    )
  }
  axis(
    side = 1, at = xat,
    labels = if (normalize) paste0(abs(xat), "%") else abs(xat),
    cex.axis = axisSize
  )
  mtext(text = colnames(m1.gs), side = 2, font = 3, las = 2, at = b1, line = 0.75, cex = geneSize)
  mtext(
    text = if (normalize) "Percent of cases" else "Number of cases",
    side = 1, line = 2.2, cex = axisSize * 1.15
  )
  title(main = paste0(m1Name, " [N=", m1.ss, "]"), cex.main = titleSize, font = 2, adj = 0)
  title(main = paste0(m2Name, " [N=", m2.ss, "]"), cex.main = titleSize, font = 2, adj = 1)

  par(mar = c(0.2, 0.8, 0.2, 0.5), xpd = TRUE)
  plot(NULL, ylab = "", xlab = "", xlim = 0:1, ylim = 0:1, axes = FALSE)
  if (showLegend) {
    vcs <- unique(c(rownames(m1.gs), rownames(m2.gs)))
    col <- colors[vcs]
    label_map <- c(
      Frame_Shift_Del = "Frame shift del",
      Frame_Shift_Ins = "Frame shift ins",
      In_Frame_Del = "In-frame del",
      In_Frame_Ins = "In-frame ins",
      Missense_Mutation = "Missense",
      Nonsense_Mutation = "Nonsense",
      Nonstop_Mutation = "Nonstop",
      Splice_Site = "Splice site",
      Translation_Start_Site = "Translation start",
      Multi_Hit = "Multi-hit"
    )
    leg_labels <- ifelse(
      names(col) %in% names(label_map),
      label_map[names(col)],
      gsub("_", " ", names(col))
    )
    legend(
      "topleft", legend = leg_labels, col = as.character(col), bty = "n",
      border = NA, xpd = TRUE, pch = 15,
      cex = legendTxtSize, y.intersp = 1.35, x.intersp = 1.1,
      pt.cex = 1.6 * legendTxtSize,
      text.font = 1,
      ncol = 2
    )
  }
  invisible(list(m1 = m1.gs, m2 = m2.gs))
}

# ---------------------------------------------------------------------------
# Load MAF + sample metadata; subset Brain vs Primary (MSK CRC)
# ---------------------------------------------------------------------------
SampleMeta <- read.csv(SAMPLE_META, stringsAsFactors = FALSE)

MSKBrM <- subset(SampleMeta, Metastatic_Site %in% c("Brain", "CNS/Brain") & Dataset == "MSK")
MSKPrimary <- subset(SampleMeta, Metastatic_Site == "Primary" & Dataset == "MSK")
message("MSK Brain n=", nrow(MSKBrM), "; MSK Primary n=", nrow(MSKPrimary))

message("Reading MAF...")
maf_data <- read.maf(maf = MSK_MUT_MAF, clinicalData = MSK_CLINICAL_SAMPLE)

mut_BrM <- subsetMaf(maf = maf_data, tsb = MSKBrM$Sample_ID)
mut_BrM@clinical.data$SampleType <- "mCRC-Brain"

mut_Primary <- subsetMaf(maf = maf_data, tsb = MSKPrimary$Sample_ID)
mut_Primary@clinical.data$SampleType <- "Primary"

BrM_Primary <- mafCompare(
  m1 = mut_BrM,
  m2 = mut_Primary,
  m1Name = "mCRC_Brain",
  m2Name = "Primary",
  minMut = 5
)

genes <- BrM_Primary$results$Hugo_Symbol
genes <- genes[!is.na(genes) & nzchar(genes)]
genes <- head(genes, 20)
if (length(genes) < 1) stop("No genes available from mafCompare.")
genes_rev <- rev(genes)

write.csv(
  BrM_Primary$results,
  file.path(out_dir, "mafCompare_mCRC_Brain_vs_Primary.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(rank = seq_along(genes), Hugo_Symbol = genes),
  file.path(out_dir, "coBarplot_genes_top20.csv"),
  row.names = FALSE
)

pdf_file <- file.path(out_dir, "compare_coBar_mCRC_Brain_vs_Primary.pdf")
pdf(file = pdf_file, width = 6, height = 9)
coBarplot_legible(
  m1 = mut_BrM,
  m2 = mut_Primary,
  m1Name = "mCRC_Brain",
  m2Name = "Primary",
  genes = genes_rev,
  geneSize = 0.8,
  showPct = TRUE,
  pctSize = 0.85,
  axisSize = 1.0,
  showLegend = TRUE,
  legendTxtSize = 1.35,
  geneMar = 5.2,
  titleSize = 1.0,
  legendPanelHeight = 2.6
)
invisible(dev.off())
message("Wrote: ", pdf_file)
message("Done.")

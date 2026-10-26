# CRC chromosome-arm co-segregation
# Rank-order pairs / sets of AMP (gain) and DEL (loss) that co-occur
# in the same sample, by metastatic site, with CRC Brain vs others.
#
# Sites (Integrated CRC: MSK + Curtis + Ji):
#   Brain, Primary, Liver, Lung, Lymph, ALL (Liver + Lung + Lymph)
#
# Also compares the top CRC-Brain pairs in LC / MELAN / BRCA brain mets.
#
# Ranking: co-occurrence frequency (n_both / n_samples), then n_both, Jaccard.
# Extra stats (not used for primary rank): lift, phi, Fisher one-sided p (enrichment).
# same_chromosome = TRUE marks p+q of the same chromosome (often whole-chr aneuploidy).
#
# Outputs -> crc_arm_cosegregation/



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
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(scales)
})

set.seed(20260911)

outDir <- DIR_COSEG
figDir <- file.path(outDir, "figures")
tabDir <- file.path(outDir, "tables")
dir.create(figDir, showWarnings = FALSE, recursive = TRUE)
dir.create(tabDir, showWarnings = FALSE, recursive = TRUE)

site_levels <- c("Brain", "Primary", "Liver", "Lung", "Lymph", "ALL")
site_cols <- c(
  Brain = "#C0392B", Primary = "#7F8C8D", Liver = "#16A085",
  Lung = "#2980B9", Lymph = "#8E44AD", ALL = "#F39C12"
)
amp_col <- "#D73027"
del_col <- "#313695"
type_cols <- c(AMP = amp_col, DEL = del_col, Mixed = "#AF7AC5")

theme_pub <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    axis.text = element_text(color = "black")
  )

# ============================================================
# Data loaders (same filters as metastatic_site_CNV_analysis.R)
# ============================================================
msk_clinical <- read.table(
  "data/msk_met_2021_clinical_data.tsv",
  header = TRUE, sep = "\t", check.names = FALSE
)
colnames(msk_clinical) <- str_replace_all(colnames(msk_clinical), ": ", "_")
colnames(msk_clinical) <- str_replace_all(colnames(msk_clinical), "[\\(\\)]", "")
colnames(msk_clinical) <- str_replace_all(colnames(msk_clinical), "/", "")
colnames(msk_clinical) <- str_replace_all(colnames(msk_clinical), " ", "_")

primary_samples <- msk_clinical$Patient_ID[
  intersect(
    which(msk_clinical$Metastatic_patient == FALSE),
    which(is.na(msk_clinical$Metastatic_Site))
  )
]
msk_clinical$SiteLabel <- ifelse(
  msk_clinical$Patient_ID %in% primary_samples, "Primary",
  ifelse(msk_clinical$Metastatic_Site == "CNS/Brain", "Brain",
         as.character(msk_clinical$Metastatic_Site))
)

load_msk_arm <- function(fname, cancer_type) {
  arm_cnv <- read.table(fname, header = TRUE)
  colnames(arm_cnv)[1] <- "Sample_ID"
  colnames(arm_cnv) <- str_replace(colnames(arm_cnv), "X", "Arm_")
  nc <- ncol(arm_cnv)
  colnames(arm_cnv)[c(nc - 1, nc)] <- c("Arm_Xp", "Arm_Xq")

  colorder <- "Sample_ID"
  for (i in 1:22)
    colorder <- c(colorder, paste0("Arm_", i, "p"), paste0("Arm_", i, "q"))
  colorder <- c(colorder, "Arm_Xp", "Arm_Xq")
  colorder <- intersect(colorder, colnames(arm_cnv))
  colorder <- setdiff(colorder, c("Arm_13p", "Arm_14p", "Arm_15p", "Arm_21p", "Arm_22p"))
  arm_cnv <- arm_cnv[, colorder]

  meta <- subset(msk_clinical, Cancer_Type == cancer_type)
  meta <- meta[, c("Sample_ID", "Patient_ID", "SiteLabel")]
  colnames(meta)[colnames(meta) == "SiteLabel"] <- "Metastatic_Site"
  merge(arm_cnv, meta, by = "Sample_ID")
}

load_crc_integrated <- function() {
  SampleMeta_CRC <- read.csv("data/ASCETS/ascets-master/IntegratedData/Sample_Metadata.csv")
  SampleMeta_CRC$X <- NULL
  arm_cnv_CRC <- read.table(
    "data/ASCETS/ascets-master/IntegratedData/IntegratedData_arm_level_calls.txt",
    header = TRUE
  )
  colnames(arm_cnv_CRC)[1] <- "Sample_ID"
  colnames(arm_cnv_CRC) <- str_replace(colnames(arm_cnv_CRC), "X", "Arm_")
  nc <- ncol(arm_cnv_CRC); colnames(arm_cnv_CRC)[(nc-3):nc] <- c("Arm_Xp", "Arm_Xq", "Arm_Yp", "Arm_Yq")

  colorder_crc <- "Sample_ID"
  for (i in 1:22)
    colorder_crc <- c(colorder_crc, paste0("Arm_", i, "p"), paste0("Arm_", i, "q"))
  colorder_crc <- c(colorder_crc, "Arm_Xp", "Arm_Xq", "Arm_Yp", "Arm_Yq")
  colorder_crc <- setdiff(colorder_crc, "Arm_13p")
  arm_cnv_CRC <- arm_cnv_CRC[, colorder_crc]
  arm_cnv_CRC <- merge(arm_cnv_CRC, SampleMeta_CRC, by = "Sample_ID")

  JiQC_reject <- c("P10033_20283A", "P10030_20281A")
  arm_cnv_CRC <- subset(arm_cnv_CRC, !(Sample_ID %in% JiQC_reject))
  drop_cols <- c("Arm_14p", "Arm_15p", "Arm_21p", "Arm_22p", "Arm_Yp", "Arm_Yq")
  arm_cnv_CRC <- arm_cnv_CRC[, setdiff(colnames(arm_cnv_CRC), drop_cols)]
  arm_cnv_CRC
}

pretty_arm <- function(x) gsub("Arm_", "", x)
arm_chr <- function(x) gsub("[pq]$", "", pretty_arm(x))

# ============================================================
# Event matrices and pairwise / triple stats
# ============================================================
binary_event <- function(df, arm_cols, cnv_type) {
  mat <- as.matrix(df[, arm_cols, drop = FALSE])
  storage.mode(mat) <- "character"
  out <- matrix(as.integer(mat == cnv_type), nrow = nrow(mat), ncol = ncol(mat))
  out[is.na(out)] <- 0L
  colnames(out) <- arm_cols
  rownames(out) <- df$Sample_ID
  out
}

fisher_greater <- function(n11, n10, n01, n00) {
  m <- matrix(c(n11, n01, n10, n00), nrow = 2)
  if (any(is.na(m)) || any(m < 0) || sum(m) == 0) return(NA_real_)
  tryCatch(fisher.test(m, alternative = "greater")$p.value, error = function(e) NA_real_)
}

pairs_from_mats <- function(mat_a, mat_b, type, n, min_both = 0L) {
  arms_a <- colnames(mat_a)
  arms_b <- colnames(mat_b)
  co <- crossprod(mat_a, mat_b)
  n_a <- colSums(mat_a)
  n_b <- colSums(mat_b)
  same_type <- identical(arms_a, arms_b) && type != "Mixed"

  rows <- list()
  k <- 0L
  for (i in seq_along(arms_a)) {
    j_start <- if (same_type) i + 1L else 1L
    if (j_start > length(arms_b)) next
    for (j in j_start:length(arms_b)) {
      if (!same_type && arms_a[i] == arms_b[j]) next
      n11 <- as.integer(co[i, j])
      if (n11 < min_both) next
      n10 <- as.integer(n_a[i] - n11)
      n01 <- as.integer(n_b[j] - n11)
      n00 <- as.integer(n - n11 - n10 - n01)
      freq_a <- n_a[i] / n
      freq_b <- n_b[j] / n
      freq_both <- n11 / n
      expected <- freq_a * freq_b
      union_n <- n11 + n10 + n01
              denom_phi <- sqrt(as.numeric(n11 + n10) * as.numeric(n01 + n00) *
                                 as.numeric(n11 + n01) * as.numeric(n10 + n00))
      ev_a <- if (type == "DEL") "DEL" else "AMP"
      ev_b <- if (type == "AMP") "AMP" else "DEL"
      if (type == "Mixed") {
        ev_a <- "AMP"
        ev_b <- "DEL"
      }
      k <- k + 1L
      rows[[k]] <- data.frame(
        Pair = paste0(pretty_arm(arms_a[i]), " ", ev_a, " + ", pretty_arm(arms_b[j]), " ", ev_b),
        Event1 = paste(pretty_arm(arms_a[i]), ev_a),
        Event2 = paste(pretty_arm(arms_b[j]), ev_b),
        Arm1 = pretty_arm(arms_a[i]),
        Arm2 = pretty_arm(arms_b[j]),
        CNV1 = ev_a,
        CNV2 = ev_b,
        PairType = type,
        same_chromosome = arm_chr(arms_a[i]) == arm_chr(arms_b[j]),
        N_samples = n,
        N_event1 = as.integer(n_a[i]),
        N_event2 = as.integer(n_b[j]),
        N_both = n11,
        Freq_event1 = round(freq_a, 4),
        Freq_event2 = round(freq_b, 4),
        Freq_both = round(freq_both, 4),
        Jaccard = if (union_n > 0) round(n11 / union_n, 4) else NA_real_,
        Lift = if (expected > 0) round(freq_both / expected, 3) else NA_real_,
        Phi = if (denom_phi > 0) round((n11 * n00 - n10 * n01) / denom_phi, 4) else NA_real_,
        Fisher_p = fisher_greater(n11, n10, n01, n00),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) {
    return(NULL)
  }
  out <- bind_rows(rows)
  out$Fisher_FDR <- p.adjust(out$Fisher_p, method = "BH")
  out
}

rank_pairs <- function(df) {
  df %>%
    arrange(desc(Freq_both), desc(N_both), desc(Jaccard), Pair) %>%
    mutate(Rank = row_number(), .before = 1)
}

site_pairs <- function(df, arm_cols, site_name) {
  n <- nrow(df)
  amp <- binary_event(df, arm_cols, "AMP")
  del <- binary_event(df, arm_cols, "DEL")
  pairs <- bind_rows(
    pairs_from_mats(amp, amp, "AMP", n),
    pairs_from_mats(del, del, "DEL", n),
    pairs_from_mats(amp, del, "Mixed", n)
  )
  pairs$Site <- site_name
  rank_pairs(pairs)
}

triples_from_events <- function(df, event_df, site_name, min_freq = 0.20, max_events = 14L, min_both = 3L) {
  keep <- event_df %>%
    filter(Freq >= min_freq) %>%
    arrange(desc(Freq), desc(N)) %>%
    slice_head(n = max_events)
  if (nrow(keep) < 3) return(NULL)

  bin <- matrix(0L, nrow = nrow(df), ncol = nrow(keep))
  colnames(bin) <- keep$Event
  for (i in seq_len(nrow(keep))) {
    bin[, i] <- as.integer(df[[keep$Arm_col[i]]] == keep$CNV[i])
  }
  bin[is.na(bin)] <- 0L
  n <- nrow(df)
  combos <- utils::combn(seq_len(nrow(keep)), 3)
  rows <- vector("list", ncol(combos))
  for (j in seq_len(ncol(combos))) {
    idx <- combos[, j]
    both <- as.integer(sum(bin[, idx[1]] * bin[, idx[2]] * bin[, idx[3]]))
    if (both < min_both) {
      rows[[j]] <- NULL
      next
    }
    rows[[j]] <- data.frame(
      Set = paste(keep$Event[idx], collapse = " + "),
      Event1 = keep$Event[idx[1]],
      Event2 = keep$Event[idx[2]],
      Event3 = keep$Event[idx[3]],
      N_samples = n,
      N_set = both,
      Freq_set = round(both / n, 4),
      Freq_e1 = keep$Freq[idx[1]],
      Freq_e2 = keep$Freq[idx[2]],
      Freq_e3 = keep$Freq[idx[3]],
      Lift = {
        expc <- keep$Freq[idx[1]] * keep$Freq[idx[2]] * keep$Freq[idx[3]]
        if (expc > 0) round((both / n) / expc, 3) else NA_real_
      },
      Site = site_name,
      stringsAsFactors = FALSE
    )
  }
  out <- bind_rows(rows)
  if (is.null(out) || nrow(out) == 0) return(NULL)
  out %>%
    arrange(desc(Freq_set), desc(N_set), desc(Lift)) %>%
    mutate(Rank = row_number(), .before = 1)
}

single_events <- function(df, arm_cols) {
  n <- nrow(df)
  rows <- lapply(arm_cols, function(ac) {
    v <- df[[ac]]
    data.frame(
      Arm = pretty_arm(ac),
      Arm_col = ac,
      CNV = c("AMP", "DEL"),
      Event = paste(pretty_arm(ac), c("AMP", "DEL")),
      N = c(sum(v == "AMP", na.rm = TRUE), sum(v == "DEL", na.rm = TRUE)),
      Freq = c(
        sum(v == "AMP", na.rm = TRUE) / n,
        sum(v == "DEL", na.rm = TRUE) / n
      ),
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows) %>%
    arrange(desc(Freq), desc(N), Event) %>%
    mutate(Rank = row_number(), .before = 1, N_samples = n)
}

onco_matrix <- function(df, arm_cols) {
  mat <- t(as.matrix(df[, arm_cols, drop = FALSE]))
  storage.mode(mat) <- "character"
  mat[!(mat %in% c("AMP", "DEL"))] <- ""
  mat[is.na(mat)] <- ""
  rownames(mat) <- pretty_arm(arm_cols)
  colnames(mat) <- df$Sample_ID
  mat
}

alter_fun <- list(
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w * 0.88, h * 0.88, gp = gpar(fill = amp_col, col = NA))
  },
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w * 0.88, h * 0.88, gp = gpar(fill = del_col, col = NA))
  }
)

# ============================================================
# Load CRC and restrict to analysis sites
# ============================================================
cat("Loading Integrated CRC...\n")
crc_all <- load_crc_integrated()
arm_cols <- grep("^Arm_", colnames(crc_all), value = TRUE)
cat("Arm columns:", length(arm_cols), "\n")

keep_sites <- c("Brain", "Primary", "Liver", "Lung", "Lymph")
crc <- subset(crc_all, Metastatic_Site %in% keep_sites)
crc$Metastatic_Site <- factor(crc$Metastatic_Site, levels = keep_sites)

counts <- as.data.frame(table(crc$Metastatic_Site), stringsAsFactors = FALSE)
colnames(counts) <- c("Metastatic_Site", "N")
counts$Metastatic_Site <- as.character(counts$Metastatic_Site)
counts <- rbind(
  counts,
  data.frame(
    Metastatic_Site = "ALL",
    N = sum(counts$N[counts$Metastatic_Site %in% c("Liver", "Lung", "Lymph")])
  )
)
counts <- counts[match(site_levels, counts$Metastatic_Site), ]
write.csv(counts, file.path(tabDir, "CRC_sample_counts.csv"), row.names = FALSE)
print(counts)

site_dfs <- list(
  Brain = subset(crc, Metastatic_Site == "Brain"),
  Primary = subset(crc, Metastatic_Site == "Primary"),
  Liver = subset(crc, Metastatic_Site == "Liver"),
  Lung = subset(crc, Metastatic_Site == "Lung"),
  Lymph = subset(crc, Metastatic_Site == "Lymph"),
  ALL = subset(crc, Metastatic_Site %in% c("Liver", "Lung", "Lymph"))
)

# ============================================================
# Ranked pair / triple tables
# ============================================================
cat("Computing pairwise co-occurrence by site...\n")
pair_list <- lapply(site_levels, function(s) {
  cat("  pairs:", s, "n =", nrow(site_dfs[[s]]), "\n")
  site_pairs(site_dfs[[s]], arm_cols, s)
})
names(pair_list) <- site_levels
pairs_all <- bind_rows(pair_list)
write.csv(pairs_all, file.path(tabDir, "CRC_pair_cosegregation_all.csv"), row.names = FALSE)

top50 <- bind_rows(lapply(pair_list, function(d) {
  d %>% slice_head(n = 50)
}))
write.csv(top50, file.path(tabDir, "CRC_pair_cosegregation_top50_by_site.csv"), row.names = FALSE)

top50_inter <- bind_rows(lapply(pair_list, function(d) {
  d %>% filter(!same_chromosome) %>% mutate(Rank = row_number()) %>% slice_head(n = 50)
}))
write.csv(top50_inter, file.path(tabDir, "CRC_pair_cosegregation_top50_interchromosomal.csv"), row.names = FALSE)

cat("Computing single-event frequencies and triples...\n")
event_list <- lapply(site_levels, function(s) {
  ev <- single_events(site_dfs[[s]], arm_cols)
  ev$Site <- s
  ev
})
names(event_list) <- site_levels
events_all <- bind_rows(event_list)
write.csv(events_all, file.path(tabDir, "CRC_single_event_frequency_by_site.csv"), row.names = FALSE)

triple_list <- lapply(site_levels, function(s) {
  cat("  triples:", s, "\n")
  triples_from_events(site_dfs[[s]], event_list[[s]], s)
})
names(triple_list) <- site_levels
triples_all <- bind_rows(triple_list)
write.csv(triples_all, file.path(tabDir, "CRC_triple_cosegregation_by_site.csv"), row.names = FALSE)

# Wide frequency table for the top Brain pairs
brain_top_n <- 30L
brain_pairs <- pair_list$Brain
brain_top <- brain_pairs %>% slice_head(n = brain_top_n)

wide_freq <- pairs_all %>%
  filter(Pair %in% brain_top$Pair) %>%
  select(Pair, PairType, same_chromosome, Site, Freq_both, N_both, Jaccard, Fisher_FDR) %>%
  pivot_wider(
    names_from = Site,
    values_from = c(Freq_both, N_both, Jaccard, Fisher_FDR),
    names_sep = "_"
  )
wide_freq <- wide_freq[match(brain_top$Pair, wide_freq$Pair), ]
wide_freq$Rank_Brain <- seq_len(nrow(wide_freq))
wide_freq <- wide_freq[, c("Rank_Brain", setdiff(colnames(wide_freq), "Rank_Brain"))]
write.csv(wide_freq, file.path(tabDir, "CRC_Brain_top_pairs_across_sites.csv"), row.names = FALSE)

# Brain vs other CRC sites: pair-frequency enrichment
enrich_rows <- list()
for (comp in c("Primary", "Liver", "Lung", "Lymph", "ALL")) {
  a <- pair_list$Brain %>%
    select(Pair, PairType, same_chromosome, N_samples, N_both, Freq_both)
  b <- pair_list[[comp]] %>%
    select(Pair, N_samples, N_both, Freq_both)
  colnames(a) <- c("Pair", "PairType", "same_chromosome", "N_Brain", "Nboth_Brain", "Freq_Brain")
  colnames(b) <- c("Pair", "N_Comp", "Nboth_Comp", "Freq_Comp")
  m <- merge(a, b, by = "Pair", all.x = TRUE)
  m$N_Comp[is.na(m$N_Comp)] <- nrow(site_dfs[[comp]])
  m$Nboth_Comp[is.na(m$Nboth_Comp)] <- 0
  m$Freq_Comp[is.na(m$Freq_Comp)] <- 0
  m$Freq_diff <- round(m$Freq_Brain - m$Freq_Comp, 4)
  m$Pval <- mapply(function(x1, n1, x2, n2) {
    if (n1 < 2 || n2 < 2 || (x1 + x2) == 0) return(NA_real_)
    tryCatch(
      prop.test(c(x1, x2), c(n1, n2), alternative = "greater")$p.value,
      error = function(e) NA_real_
    )
  }, m$Nboth_Brain, m$N_Brain, m$Nboth_Comp, m$N_Comp)
  m$FDR <- p.adjust(m$Pval, method = "BH")
  m$Comparator <- comp
  enrich_rows[[comp]] <- m
}
enrich_all <- bind_rows(enrich_rows) %>%
  arrange(Comparator, desc(Freq_diff), Pval)
write.csv(enrich_all, file.path(tabDir, "CRC_Brain_pair_enrichment_vs_sites.csv"), row.names = FALSE)

# ============================================================
# Other-cancer brain mets: frequency of the same CRC-Brain top pairs
# ============================================================
cat("Loading other-cancer brain mets for pair comparison...\n")
lc_brain <- subset(
  load_msk_arm("data/ASCETS/ascets-master/MSK_NSCLC_Data/Output/LC_arm_level_calls.txt",
               "Non-Small Cell Lung Cancer"),
  Metastatic_Site == "Brain"
)
melan_brain <- subset(
  load_msk_arm("data/ASCETS/ascets-master/MSK_MELAN_Data/Output/MELAN_arm_level_calls.txt",
               "Melanoma"),
  Metastatic_Site == "Brain"
)
brca_brain <- subset(
  load_msk_arm("data/ASCETS/ascets-master/MSK_BRCA_Data/Output/BRCA_arm_level_calls.txt",
               "Breast Cancer"),
  Metastatic_Site == "Brain"
)

# Same Fisher FDR procedure as CRC sites / CRC_Brain_top_pairs_across_sites:
# full AMP/DEL/Mixed pair scan per cohort, BH within each PairType (pairs_from_mats),
# then look up the top CRC-Brain pairs.
cancer_brain <- list(
  CRC = site_dfs$Brain,
  LC = lc_brain,
  MELAN = melan_brain,
  BRCA = brca_brain
)
cancer_levels <- names(cancer_brain)

cat("Computing full pairwise co-occurrence in other-cancer brain mets (for FDR)...\n")
cancer_pair_list <- list(
  CRC = pair_list$Brain  # already genome-wide within PairType; same FDR as Freq_both_Brain table
)
for (cn in c("LC", "MELAN", "BRCA")) {
  dfc <- cancer_brain[[cn]]
  ac <- grep("^Arm_", colnames(dfc), value = TRUE)
  cat("  pairs:", cn, "Brain n =", nrow(dfc), "\n")
  cancer_pair_list[[cn]] <- site_pairs(dfc, ac, cn)
}

cross_long <- lapply(cancer_levels, function(cn) {
  pl <- cancer_pair_list[[cn]]
  hit <- pl[match(brain_top$Pair, pl$Pair), ]
  data.frame(
    Pair = brain_top$Pair,
    Cancer = cn,
    Freq_both = hit$Freq_both,
    N_both = hit$N_both,
    Jaccard = hit$Jaccard,
    Fisher_FDR = hit$Fisher_FDR,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

# Pairs absent from a cohort scan (should be rare with min_both=0): fill zeros / NA FDR
cross_long$Freq_both[is.na(cross_long$Freq_both)] <- 0
cross_long$N_both[is.na(cross_long$N_both)] <- 0L

cross_cancer <- cross_long %>%
  pivot_wider(
    names_from = Cancer,
    values_from = c(Freq_both, N_both, Jaccard, Fisher_FDR),
    names_sep = "_"
  )
cross_cancer <- merge(
  brain_top[, c("Pair", "PairType", "same_chromosome")],
  cross_cancer,
  by = "Pair",
  all.x = TRUE,
  sort = FALSE
)
cross_cancer <- cross_cancer[match(brain_top$Pair, cross_cancer$Pair), ]
cross_cancer$Rank_CRC <- seq_len(nrow(cross_cancer))
metric_order <- unlist(lapply(
  c("Freq_both", "N_both", "Jaccard", "Fisher_FDR"),
  function(m) paste0(m, "_", cancer_levels)
))
cross_cancer <- cross_cancer[, c(
  "Rank_CRC", "Pair", "PairType", "same_chromosome", metric_order
)]
write.csv(
  cross_cancer,
  file.path(tabDir, "CRC_BrM_vs_otherCancers_top_pairs.csv"),
  row.names = FALSE
)
cat("Brain-met n: CRC", nrow(site_dfs$Brain),
    "LC", nrow(lc_brain), "MELAN", nrow(melan_brain), "BRCA", nrow(brca_brain), "\n")
# Sanity: CRC FDR in this table should match Fisher_FDR_Brain in across-sites table
cat("CRC Fisher_FDR vs Brain across-sites (max abs diff):",
    max(abs(cross_cancer$Fisher_FDR_CRC - wide_freq$Fisher_FDR_Brain), na.rm = TRUE), "\n")

# ============================================================
# Core multi-event signatures in CRC Brain
# ============================================================
core_n <- 8L
core_events <- event_list$Brain %>% slice_head(n = core_n)
brain_df <- site_dfs$Brain
core_bin <- matrix(0L, nrow = nrow(brain_df), ncol = core_n)
colnames(core_bin) <- core_events$Event
for (i in seq_len(core_n)) {
  core_bin[, i] <- as.integer(brain_df[[core_events$Arm_col[i]]] == core_events$CNV[i])
}
core_bin[is.na(core_bin)] <- 0L
pattern <- apply(core_bin, 1, paste, collapse = "")
pat_tab <- as.data.frame(sort(table(pattern), decreasing = TRUE), stringsAsFactors = FALSE)
colnames(pat_tab) <- c("Pattern", "N")
pat_tab$Freq <- round(pat_tab$N / nrow(brain_df), 4)
pat_tab$Rank <- seq_len(nrow(pat_tab))
bit_mat <- do.call(rbind, strsplit(as.character(pat_tab$Pattern), ""))
bit_mat <- apply(bit_mat, 2, as.integer)
colnames(bit_mat) <- core_events$Event
sig_tbl <- cbind(pat_tab[, c("Rank", "N", "Freq")], as.data.frame(bit_mat))
write.csv(sig_tbl, file.path(tabDir, "CRC_Brain_core_signatures.csv"), row.names = FALSE)

# ============================================================
# Figures
# ============================================================
save_ht <- function(ht, path, width, height) {
  pdf(path, width = width, height = height)
  draw(ht, merge_legend = TRUE)
  dev.off()
  cat("Wrote", path, "\n")
}

# --- Oncoprints (sample x arm): default column order groups co-segregating samples
cat("Drawing oncoprints...\n")
pdf(file.path(figDir, "CRC_oncoprint_by_site.pdf"), width = 14, height = 8.5)
for (s in site_levels) {
  mat <- onco_matrix(site_dfs[[s]], arm_cols)
  ht <- oncoPrint(
    mat,
    alter_fun = alter_fun,
    col = c(AMP = amp_col, DEL = del_col),
    show_column_names = FALSE,
    show_pct = TRUE,
    pct_gp = gpar(fontsize = 7),
    row_names_gp = gpar(fontsize = 8),
    column_title = paste0("CRC ", s, "  (n = ", ncol(mat), ")"),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_param = list(title = "Alteration"),
    top_annotation = HeatmapAnnotation(
      Site = rep(s, ncol(mat)),
      col = list(Site = site_cols),
      show_annotation_name = FALSE,
      simple_anno_size = unit(2.5, "mm")
    )
  )
  draw(ht, newpage = TRUE)
}
dev.off()
cat("Wrote", file.path(figDir, "CRC_oncoprint_by_site.pdf"), "\n")

mat_brain <- onco_matrix(site_dfs$Brain, arm_cols)
ht_brain <- oncoPrint(
  mat_brain,
  alter_fun = alter_fun,
  col = c(AMP = amp_col, DEL = del_col),
  show_column_names = FALSE,
  show_pct = TRUE,
  pct_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 9),
  column_title = paste0("CRC brain metastases - arm AMP / DEL (n = ", ncol(mat_brain), ")"),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(title = "Alteration")
)
save_ht(ht_brain, file.path(figDir, "CRC_oncoprint_Brain.pdf"), 13, 8.5)

# --- Rank heatmap: top Brain pairs x CRC sites
cat("Drawing rank / comparison heatmaps...\n")
rank_mat <- as.matrix(wide_freq[, paste0("Freq_both_", site_levels)])
rownames(rank_mat) <- wide_freq$Pair
colnames(rank_mat) <- site_levels
rank_mat[is.na(rank_mat)] <- 0
pair_anno <- data.frame(
  Type = factor(wide_freq$PairType, levels = c("AMP", "DEL", "Mixed")),
  Chr = factor(ifelse(wide_freq$same_chromosome, "same chr", "inter-chr"),
               levels = c("same chr", "inter-chr"))
)
ha_pair <- rowAnnotation(
  Type = pair_anno$Type,
  Chr = pair_anno$Chr,
  col = list(
    Type = type_cols,
    Chr = c("same chr" = "#BDC3C7", "inter-chr" = "#2C3E50")
  ),
  annotation_name_gp = gpar(fontsize = 9)
)
ht_rank <- Heatmap(
  rank_mat,
  name = "Co-occur\nfreq",
  col = colorRamp2(c(0, 0.3, 0.65), c("#F7F7F7", "#FDBB84", "#B30000")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  right_annotation = ha_pair,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 45,
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", rank_mat[i, j]), x, y, gp = gpar(fontsize = 7))
  },
  column_title = "Top 30 CRC Brain pairs - co-occurrence frequency by site",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  width = unit(6, "cm")
)
save_ht(ht_rank, file.path(figDir, "CRC_top_pairs_rank_heatmap.pdf"), 11, 10)

# CRC BrM vs other cancer BrMs
cc_mat <- as.matrix(cross_cancer[, paste0("Freq_both_", cancer_levels)])
rownames(cc_mat) <- cross_cancer$Pair
colnames(cc_mat) <- c("CRC Brain", "LC Brain", "MELAN Brain", "BRCA Brain")
cc_mat[is.na(cc_mat)] <- 0
ha_cc <- rowAnnotation(
  Type = factor(cross_cancer$PairType, levels = c("AMP", "DEL", "Mixed")),
  Chr = factor(ifelse(cross_cancer$same_chromosome, "same chr", "inter-chr"),
               levels = c("same chr", "inter-chr")),
  col = list(
    Type = type_cols,
    Chr = c("same chr" = "#BDC3C7", "inter-chr" = "#2C3E50")
  )
)
n_lab <- paste0(
  colnames(cc_mat), "\n(n=",
  c(nrow(site_dfs$Brain), nrow(lc_brain), nrow(melan_brain), nrow(brca_brain)),
  ")"
)
colnames(cc_mat) <- n_lab
ht_cc <- Heatmap(
  cc_mat,
  name = "Co-occur\nfreq",
  col = colorRamp2(c(0, 0.3, 0.65), c("#F7F7F7", "#FDBB84", "#B30000")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  right_annotation = ha_cc,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", cc_mat[i, j]), x, y, gp = gpar(fontsize = 7))
  },
  column_title = "Top CRC Brain pairs in brain mets of other cancer types",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  width = unit(5.5, "cm")
)
save_ht(ht_cc, file.path(figDir, "CRC_BrM_vs_otherCancers_pair_heatmap.pdf"), 11, 10)

# --- Ranked lollipop plot
cat("Drawing ranked lollipop plot...\n")
lolli <- pair_list$Brain %>%
  filter(!same_chromosome) %>%
  slice_head(n = 20)
p_lolli <- ggplot(lolli, aes(x = Freq_both, y = reorder(Pair, Freq_both), color = PairType)) +
  geom_segment(aes(x = 0, xend = Freq_both, yend = Pair), linewidth = 0.7) +
  geom_point(aes(size = N_both)) +
  scale_color_manual(values = type_cols) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_size_continuous(range = c(2.5, 6), name = "n both") +
  labs(
    title = "CRC Brain - top inter-chromosomal co-occurring pairs",
    subtitle = "Ranked by co-occurrence frequency; point size = number of samples with both events",
    x = "Co-occurrence frequency", y = NULL, color = "Pair type"
  ) +
  theme_pub +
  theme(legend.position = "right")
ggsave(file.path(figDir, "CRC_Brain_interchrom_pairs_lollipop.pdf"), p_lolli,
       width = 10, height = 7)

# --- UpSet of top Brain events (multi-event sets)
cat("Drawing UpSet and core-signature plots...\n")
top_ev <- event_list$Brain %>% slice_head(n = 8)
lt <- list()
for (i in seq_len(nrow(top_ev))) {
  hit <- brain_df$Sample_ID[brain_df[[top_ev$Arm_col[i]]] == top_ev$CNV[i]]
  lt[[top_ev$Event[i]]] <- as.character(hit)
}
comb <- make_comb_mat(lt)
pdf(file.path(figDir, "CRC_Brain_UpSet_top_events.pdf"), width = 12, height = 6.5)
ht_up <- UpSet(
  comb,
  comb_col = "#2C3E50",
  bg_col = c("#F4F6F7", "white"),
  top_annotation = upset_top_annotation(comb, add_numbers = TRUE, annotation_name_rot = 90),
  right_annotation = upset_right_annotation(comb, add_numbers = TRUE),
  row_names_gp = gpar(fontsize = 10)
)
draw(ht_up)
grid.text(
  "CRC Brain - intersections among the 8 most frequent arm events",
  x = 0.5, y = 0.98, gp = gpar(fontsize = 13, fontface = "bold")
)
dev.off()
cat("Wrote", file.path(figDir, "CRC_Brain_UpSet_top_events.pdf"), "\n")

# Core signature heatmap (unique binary patterns among top 8 events)
sig_show <- sig_tbl %>% slice_head(n = min(20, nrow(sig_tbl)))
sig_mat <- as.matrix(sig_show[, core_events$Event])
storage.mode(sig_mat) <- "numeric"
rownames(sig_mat) <- paste0("S", sig_show$Rank, "  n=", sig_show$N, "  (", percent(sig_show$Freq, 0.1), ")")
ha_sig <- rowAnnotation(
  n = anno_barplot(sig_show$N, gp = gpar(fill = "#2C3E50"), border = FALSE, width = unit(2.2, "cm")),
  annotation_name_gp = gpar(fontsize = 9)
)
col_sig <- c("0" = "#F4F6F7", "1" = "#1A5276")
# color by event type on columns
ha_col <- HeatmapAnnotation(
  CNV = core_events$CNV,
  col = list(CNV = type_cols[c("AMP", "DEL")]),
  annotation_name_gp = gpar(fontsize = 9)
)
ht_sig <- Heatmap(
  sig_mat,
  name = "Present",
  col = colorRamp2(c(0, 1), c("#F4F6F7", "#1A5276")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  right_annotation = ha_sig,
  top_annotation = ha_col,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  column_title = "CRC Brain - most frequent combinations of the top 8 arm events",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  heatmap_legend_param = list(at = c(0, 1), labels = c("absent", "present"))
)
save_ht(ht_sig, file.path(figDir, "CRC_Brain_core_signature_heatmap.pdf"), 10, 7.5)

# ============================================================
# Console rank-order summary
# ============================================================
cat("\n========== RANK-ORDER SUMMARY ==========\n")
for (s in site_levels) {
  cat("\n---", s, "n =", nrow(site_dfs[[s]]), "---\n")
  cat("Top 8 pairs (all):\n")
  print(pair_list[[s]] %>%
          slice_head(n = 8) %>%
          select(Rank, Pair, PairType, Freq_both, N_both, Jaccard, same_chromosome),
        row.names = FALSE)
  cat("Top 8 inter-chromosomal pairs:\n")
  print(pair_list[[s]] %>%
          filter(!same_chromosome) %>%
          mutate(Rank = row_number()) %>%
          slice_head(n = 8) %>%
          select(Rank, Pair, PairType, Freq_both, N_both, Jaccard),
        row.names = FALSE)
  if (!is.null(triple_list[[s]])) {
    cat("Top 5 triples (events with site freq >= 20%):\n")
    print(triple_list[[s]] %>%
            slice_head(n = 5) %>%
            select(Rank, Set, Freq_set, N_set, Lift),
          row.names = FALSE)
  }
}

cat("\nDone. Tables:", tabDir, "\nFigures:", figDir, "\n")

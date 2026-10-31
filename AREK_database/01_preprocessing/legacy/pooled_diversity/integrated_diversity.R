#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  need <- c("vegan", "ggplot2")
  missing <- need[!vapply(need, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing R packages: ", paste(missing, collapse=", "),
      "\nInstall with: conda install -n alcohol -c conda-forge r-vegan r-ggplot2 -y"
    )
  }
  library(vegan)
  library(ggplot2)
})

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2)stop("usage: Rscript integrated_diversity.R OUTPUT_BASE EXISTING_INTEGRATION_ROOT")
if(Sys.getenv("AREK_POOLED_SCOPE")!="verified_same_host_site_pipeline")stop("Set AREK_POOLED_SCOPE=verified_same_host_site_pipeline only after checking source study/host/site/pipeline conditions; no automatic merging")

BASE_DEFAULT <- "."
INTEGRATION_DEFAULT <- ""

base_dir <- if (length(args) >= 1) args[[1]] else BASE_DEFAULT
integration_root <- if (length(args) >= 2) args[[2]] else INTEGRATION_DEFAULT

out_dir <- file.path(base_dir, "results")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

set.seed(20261005)

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep="")
}

pick_existing <- function(paths) {
  for (p in paths) if (file.exists(p) && file.info(p)$size > 0) return(p)
  stop("None of the expected input files exists:\n", paste(paths, collapse="\n"))
}

unfiltered_file <- pick_existing(c(
  file.path(integration_root, "FINAL", "GUT_STOOL_all_datasets_genus_relative_unfiltered.tsv"),
  file.path(integration_root, "02_merged_cohorts", "all_ready", "merged_relative_unfiltered.tsv")
))

filtered_file <- pick_existing(c(
  file.path(integration_root, "FINAL", "GUT_STOOL_all_datasets_genus_relative_filtered.tsv"),
  file.path(integration_root, "02_merged_cohorts", "all_ready", "merged_relative_filtered_tss.tsv")
))

metadata_file <- pick_existing(c(
  file.path(integration_root, "FINAL", "GUT_STOOL_all_datasets_metadata.tsv"),
  file.path(integration_root, "02_merged_cohorts", "all_ready", "metadata.tsv")
))

std_root <- file.path(integration_root, "01_standardized_by_dataset")

msg("Unfiltered matrix: ", unfiltered_file)
msg("Filtered matrix:   ", filtered_file)
msg("Metadata:          ", metadata_file)
msg("Output:            ", out_dir)

read_matrix <- function(path) {
  x <- read.delim(
    path, row.names=1, check.names=FALSE, quote="", comment.char="",
    stringsAsFactors=FALSE
  )
  x[] <- lapply(x, function(z) as.numeric(as.character(z)))
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if(any(!is.finite(x))||any(x<0))stop("matrix contains missing/nonfinite/negative values; no silent zero filling")
  x
}

rel_unf <- read_matrix(unfiltered_file)
rel_flt <- read_matrix(filtered_file)
meta_raw <- read.delim(
  metadata_file, check.names=FALSE, stringsAsFactors=FALSE,
  quote="", comment.char=""
)

# ============================================================
# Robust metadata-column discovery
# Different integration versions may call the same fields
# Run_accession / run_accession / sample_id / SampleID, etc.
# Resolve them from actual overlap with matrix sample names rather
# than requiring one hard-coded header.
# ============================================================

clean_chr <- function(x) trimws(as.character(x))

sample_names <- intersect(colnames(rel_unf), colnames(rel_flt))
if (length(sample_names) < 4) {
  stop("Too few shared sample columns between unfiltered and filtered matrices: ",
       length(sample_names))
}

sample_candidates <- c(
  "Run_accession", "run_accession", "Run", "run",
  "sample_id", "SampleID", "sampleid", "Sample_ID", "sample",
  "sample_name", "SampleName", "id", "ID", "accession"
)

sample_scores <- sapply(names(meta_raw), function(nm) {
  vals <- unique(clean_chr(meta_raw[[nm]]))
  sum(vals %in% sample_names)
})

# Prefer a conventional ID column when ties occur.
candidate_rank <- match(names(meta_raw), sample_candidates)
candidate_rank[is.na(candidate_rank)] <- length(sample_candidates) + 1000
ord <- order(-sample_scores, candidate_rank)
sample_col <- names(meta_raw)[ord[1]]
best_overlap <- sample_scores[ord[1]]

if (best_overlap < 4) {
  top <- head(data.frame(
    column=names(meta_raw)[order(-sample_scores)],
    overlap=sort(sample_scores, decreasing=TRUE),
    stringsAsFactors=FALSE
  ), 10)
  stop(
    "Could not identify the metadata sample-ID column. ",
    "Best overlap with matrix sample names was ", best_overlap, ".\n",
    "Metadata columns: ", paste(names(meta_raw), collapse=", "), "\n",
    "Top overlap candidates:\n",
    paste(capture.output(print(top, row.names=FALSE)), collapse="\n")
  )
}

meta <- meta_raw
meta$Run_accession <- clean_chr(meta[[sample_col]])

normalize_group <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y <- gsub("[ -]+", "_", y)
  y[y %in% c("NONALCOHOL", "NON_ALCOHOLIC", "NO_ALCOHOL")] <- "NON_ALCOHOL"
  y
}

group_candidates <- c(
  "Alcohol", "analysis_group", "alcohol_group", "exposure_group",
  "Group", "group", "phenotype", "Phenotype"
)

group_score <- sapply(names(meta), function(nm) {
  y <- normalize_group(meta[[nm]])
  sum(y %in% c("ALCOHOL", "NON_ALCOHOL"))
})
group_rank <- match(names(meta), group_candidates)
group_rank[is.na(group_rank)] <- length(group_candidates) + 1000
gord <- order(-group_score, group_rank)
group_col <- names(meta)[gord[1]]

if (group_score[gord[1]] < 4) {
  stop(
    "Could not identify ALCOHOL/NON_ALCOHOL metadata column.\n",
    "Metadata columns: ", paste(names(meta_raw), collapse=", ")
  )
}
meta$Alcohol <- normalize_group(meta[[group_col]])
for(nm in intersect(c("subject_id","SubjectID","participant_id","participant_ID"),names(meta))){
 v<-meta[[nm]];v<-v[!is.na(v)&v!=""]
 if(anyDuplicated(v))stop("repeated subjects in pooled legacy branch; use corrected canonical metadata/block route or explicit study-specific specimen selection")
}

dataset_candidates <- c(
  "Dataset", "study_accession", "Study_accession", "study",
  "Study", "project_accession", "Project", "project"
)

dataset_col <- NULL
for (nm in dataset_candidates) {
  if (nm %in% names(meta)) {
    dataset_col <- nm
    break
  }
}
if (is.null(dataset_col)) {
  # Fallback: choose a column that most often looks like a study/project accession.
  acc_scores <- sapply(names(meta), function(nm) {
    y <- clean_chr(meta[[nm]])
    sum(grepl("^(PRJ|CRA|ERP|SRP|DRP|SAMN|SAMEA)", y))
  })
  if (max(acc_scores) > 0) {
    dataset_col <- names(meta)[which.max(acc_scores)]
  }
}

if (is.null(dataset_col)) {
  warning("No dataset/study column detected; Dataset will be set to UNKNOWN. ",
          "Dataset-adjusted models will not be informative.")
  meta$Dataset <- "UNKNOWN"
  dataset_col <- "<derived:UNKNOWN>"
} else {
  meta$Dataset <- clean_chr(meta[[dataset_col]])
}

# Save exactly which columns were resolved, so the run is auditable.
mapping <- data.frame(
  canonical=c("Run_accession", "Alcohol", "Dataset"),
  source_column=c(sample_col, group_col, dataset_col),
  stringsAsFactors=FALSE
)
write.table(mapping, file.path(out_dir, "metadata_column_mapping.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

msg("Metadata mapping: Run_accession <- ", sample_col,
    "; Alcohol <- ", group_col,
    "; Dataset <- ", dataset_col,
    "; sample-ID overlap=", best_overlap)

# Restrict to the two requested exposure groups.
meta <- meta[meta$Alcohol %in% c("ALCOHOL", "NON_ALCOHOL"), , drop=FALSE]

# Remove empty sample IDs; retain one row per sample. If duplicated IDs exist,
# write them for inspection before using the first occurrence.
meta <- meta[nzchar(meta$Run_accession), , drop=FALSE]
if (anyDuplicated(meta$Run_accession)) {
  dup <- meta[duplicated(meta$Run_accession) |
              duplicated(meta$Run_accession, fromLast=TRUE), , drop=FALSE]
  write.table(dup, file.path(out_dir, "duplicated_metadata_sample_ids.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  warning("Duplicated sample IDs found in metadata; first occurrence retained. ",
          "See duplicated_metadata_sample_ids.tsv")
  meta <- meta[!duplicated(meta$Run_accession), , drop=FALSE]
}

common <- Reduce(intersect, list(meta$Run_accession, colnames(rel_unf), colnames(rel_flt)))
if (length(common) < 4) {
  stop("Too few shared samples between matrices and metadata after automatic ID mapping: ",
       length(common),
       "\nResolved metadata sample column: ", sample_col)
}

meta <- meta[match(common, meta$Run_accession), , drop=FALSE]
rel_unf <- rel_unf[, common, drop=FALSE]
rel_flt <- rel_flt[, common, drop=FALSE]

meta$Alcohol <- factor(meta$Alcohol, levels=c("NON_ALCOHOL", "ALCOHOL"))
meta$Dataset <- factor(meta$Dataset)

# Re-normalize defensively.
renorm <- function(x) {
  s <- colSums(x)
  good <- s > 0
  x[, good] <- sweep(x[, good, drop=FALSE], 2, s[good], "/")
  x
}
rel_unf <- renorm(rel_unf)
rel_flt <- renorm(rel_flt)

group_counts <- as.data.frame(table(meta$Alcohol), stringsAsFactors=FALSE)
names(group_counts) <- c("Group", "N")
write.table(group_counts, file.path(out_dir, "group_counts.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

dataset_group <- as.data.frame.matrix(table(meta$Dataset, meta$Alcohol))
dataset_group$Dataset <- rownames(dataset_group)
dataset_group <- dataset_group[, c("Dataset", setdiff(names(dataset_group), "Dataset")), drop=FALSE]
write.table(dataset_group, file.path(out_dir, "dataset_group_counts.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

msg("Samples used: ", nrow(meta),
    " (ALCOHOL=", sum(meta$Alcohol=="ALCOHOL"),
    ", NON_ALCOHOL=", sum(meta$Alcohol=="NON_ALCOHOL"), ")")

# ============================================================
# Alpha diversity
# Biomni-compatible set:
# Observed richness, Shannon, Simpson, Chao1
#
# Observed/Shannon/Simpson use the UNFILTERED genus relative table.
# Chao1 is reconstructed from per-dataset standardized abundance tables
# only when the abundance is count/count-like; relative-only samples are NA.
# ============================================================

community_rel <- t(rel_unf)  # samples x genera

alpha <- data.frame(
  Run_accession = rownames(community_rel),
  Dataset = as.character(meta$Dataset),
  Alcohol = as.character(meta$Alcohol),
  Observed = rowSums(community_rel > 0),
  Shannon = diversity(community_rel, index="shannon"),
  Simpson = diversity(community_rel, index="simpson"),
  Chao1 = NA_real_,
  Chao1_source = "NO_COUNT_LIKE_MATRIX",
  stringsAsFactors=FALSE
)

# Rebuild count/count-like sample profiles from standardized per-dataset outputs.
# This avoids inventing pseudo-counts from relative abundance.
for (ds in if(Sys.getenv("AREK_INTEGRATION_INTEGER_COUNTS")=="1") unique(alpha$Dataset) else character()) {
  f <- file.path(std_root, ds, "genus_abundance_unfiltered.tsv")
  if (!file.exists(f) || file.info(f)$size <= 0) next

  ab <- tryCatch(read_matrix(f), error=function(e) NULL)
  if (is.null(ab)) next

  ss <- intersect(alpha$Run_accession[alpha$Dataset == ds], colnames(ab))
  if (length(ss) == 0) next

  for (sid in ss) {
    v <- ab[, sid]
    v[!is.finite(v) | v < 0] <- 0
    total <- sum(v)
    if (total <= 0) next

    # Reject obvious relative-fraction / relative-percent tables.
    relative_like <- (total >= 0.85 && total <= 1.15) ||
                     (total >= 85 && total <= 115)

    if (relative_like) next

    # Count-like abundance (including redistributed read estimates such as Bracken).
    # Round only for Chao1's singleton/doubleton calculation.
    if(any(v != floor(v))) next
    vr <- v # repaired: never fabricate singleton/doubleton counts by rounding
    vr[vr < 0] <- 0
    if (sum(vr) <= 0) next

    er <- tryCatch(estimateR(matrix(vr, nrow=1)), error=function(e) NULL)
    if (is.null(er) || !"S.chao1" %in% rownames(er)) next

    ii <- match(sid, alpha$Run_accession)
    alpha$Chao1[ii] <- as.numeric(er["S.chao1", 1])
    alpha$Chao1_source[ii] <- "COUNT_OR_COUNT_LIKE"
  }
}

write.table(alpha, file.path(out_dir, "alpha_diversity_values.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# Group summaries.
metrics <- c("Observed", "Shannon", "Simpson", "Chao1")
summaries <- list()
tests <- list()

for (metric in metrics) {
  d <- alpha[is.finite(alpha[[metric]]), , drop=FALSE]
  if (nrow(d) == 0) next

  for (g in c("NON_ALCOHOL", "ALCOHOL")) {
    v <- d[[metric]][d$Alcohol == g]
    if (length(v) == 0) next
    summaries[[length(summaries)+1]] <- data.frame(
      Metric=metric, Group=g, N=length(v),
      Mean=mean(v), SD=sd(v), Median=median(v),
      Q1=unname(quantile(v, 0.25)), Q3=unname(quantile(v, 0.75))
    )
  }

  d$Alcohol <- factor(d$Alcohol, levels=c("NON_ALCOHOL", "ALCOHOL"))
  d$Dataset <- factor(d$Dataset)

  w <- tryCatch(wilcox.test(d[[metric]] ~ d$Alcohol, exact=FALSE),
                error=function(e) NULL)

  p_w <- if (is.null(w)) NA_real_ else w$p.value

  # Dataset-adjusted sensitivity test; coefficient is ALCOHOL vs NON_ALCOHOL.
  p_adj <- NA_real_
  beta_adj <- NA_real_
  fit <- tryCatch(lm(d[[metric]] ~ Alcohol + Dataset, data=d), error=function(e) NULL)
  if (!is.null(fit)) {
    cf <- summary(fit)$coefficients
    if ("AlcoholALCOHOL" %in% rownames(cf)) {
      beta_adj <- cf["AlcoholALCOHOL", "Estimate"]
      p_adj <- cf["AlcoholALCOHOL", "Pr(>|t|)"]
    }
  }

  med_a <- median(d[[metric]][d$Alcohol=="ALCOHOL"], na.rm=TRUE)
  med_n <- median(d[[metric]][d$Alcohol=="NON_ALCOHOL"], na.rm=TRUE)

  tests[[length(tests)+1]] <- data.frame(
    Metric=metric,
    N=nrow(d),
    N_ALCOHOL=sum(d$Alcohol=="ALCOHOL"),
    N_NON_ALCOHOL=sum(d$Alcohol=="NON_ALCOHOL"),
    Median_ALCOHOL=med_a,
    Median_NON_ALCOHOL=med_n,
    Median_difference_ALCOHOL_minus_NON=med_a-med_n,
    Wilcoxon_p=p_w,
    Dataset_adjusted_beta_ALCOHOL=beta_adj,
    Dataset_adjusted_p=p_adj
  )
}

alpha_summary <- do.call(rbind, summaries)
alpha_tests <- do.call(rbind, tests)
if (!is.null(alpha_tests) && nrow(alpha_tests) > 0) {
  alpha_tests$Wilcoxon_FDR <- p.adjust(alpha_tests$Wilcoxon_p, method="BH")
  alpha_tests$Dataset_adjusted_FDR <- p.adjust(alpha_tests$Dataset_adjusted_p, method="BH")
}

write.table(alpha_summary, file.path(out_dir, "alpha_diversity_group_summary.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(alpha_tests, file.path(out_dir, "alpha_diversity_tests.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# Alpha plot.
long_alpha <- do.call(rbind, lapply(metrics, function(m) {
  data.frame(
    Run_accession=alpha$Run_accession,
    Dataset=alpha$Dataset,
    Alcohol=alpha$Alcohol,
    Metric=m,
    Value=alpha[[m]],
    stringsAsFactors=FALSE
  )
}))
long_alpha <- long_alpha[is.finite(long_alpha$Value), , drop=FALSE]
long_alpha$Alcohol <- factor(long_alpha$Alcohol, levels=c("NON_ALCOHOL", "ALCOHOL"))

p_alpha <- ggplot(long_alpha, aes(x=Alcohol, y=Value, fill=Alcohol)) +
  geom_boxplot(outlier.shape=NA, width=0.65) +
  facet_wrap(~Metric, scales="free_y", ncol=2) +
  labs(x=NULL, y=NULL, title="Alpha diversity: ALCOHOL vs NON_ALCOHOL") +
  theme_bw(base_size=12) +
  theme(legend.position="none")

ggsave(file.path(out_dir, "alpha_diversity_boxplots.png"),
       p_alpha, width=9, height=7, dpi=300)
ggsave(file.path(out_dir, "alpha_diversity_boxplots.pdf"),
       p_alpha, width=9, height=7)

# ============================================================
# Beta diversity
# Biomni-compatible:
#   Bray-Curtis
#   Aitchison (CLR + Euclidean)
# Tests:
#   PERMANOVA, pooled + Dataset-adjusted
#   PERMDISP
# ============================================================

comm <- t(rel_flt)  # samples x genera
comm <- comm[meta$Run_accession, , drop=FALSE]

# Remove all-zero taxa/samples defensively.
comm <- comm[, colSums(comm) > 0, drop=FALSE]
keep_s <- rowSums(comm) > 0
comm <- comm[keep_s, , drop=FALSE]
meta_beta <- meta[keep_s, , drop=FALSE]
rownames(meta_beta) <- meta_beta$Run_accession

meta_beta$Alcohol <- factor(meta_beta$Alcohol, levels=c("NON_ALCOHOL", "ALCOHOL"))
meta_beta$Dataset <- factor(meta_beta$Dataset)

run_adonis <- function(dist_obj, label) {
  rows <- list()

  pooled <- tryCatch(
    adonis2(dist_obj ~ Alcohol, data=meta_beta, permutations=999, by="margin"),
    error=function(e) e
  )
  if (inherits(pooled, "error")) {
    rows[[length(rows)+1]] <- data.frame(
      Distance=label, Model="Alcohol_only", Term="ERROR",
      Df=NA, R2=NA, F=NA, P=NA, Note=conditionMessage(pooled)
    )
  } else {
    rr <- as.data.frame(pooled)
    rr$Term <- rownames(rr)
    a <- rr[rr$Term=="Alcohol", , drop=FALSE]
    rows[[length(rows)+1]] <- data.frame(
      Distance=label, Model="Alcohol_only", Term="Alcohol",
      Df=a$Df, R2=a$R2, F=a$F, P=a$`Pr(>F)`, Note=""
    )
  }

  adjusted <- tryCatch(
    adonis2(dist_obj ~ Dataset + Alcohol, data=meta_beta, permutations=999, by="margin"),
    error=function(e) e
  )
  if (inherits(adjusted, "error")) {
    rows[[length(rows)+1]] <- data.frame(
      Distance=label, Model="Dataset_plus_Alcohol", Term="ERROR",
      Df=NA, R2=NA, F=NA, P=NA, Note=conditionMessage(adjusted)
    )
  } else {
    rr <- as.data.frame(adjusted)
    rr$Term <- rownames(rr)
    for (term in intersect(c("Dataset","Alcohol"), rr$Term)) {
      a <- rr[rr$Term==term, , drop=FALSE]
      rows[[length(rows)+1]] <- data.frame(
        Distance=label, Model="Dataset_plus_Alcohol", Term=term,
        Df=a$Df, R2=a$R2, F=a$F, P=a$`Pr(>F)`, Note=""
      )
    }
  }
  do.call(rbind, rows)
}

run_dispersion <- function(dist_obj, label) {
  bd <- tryCatch(betadisper(dist_obj, meta_beta$Alcohol), error=function(e) e)
  if (inherits(bd, "error")) {
    return(data.frame(Distance=label, F=NA, P=NA, Note=conditionMessage(bd)))
  }
  pt <- tryCatch(permutest(bd, permutations=999), error=function(e) e)
  if (inherits(pt, "error")) {
    return(data.frame(Distance=label, F=NA, P=NA, Note=conditionMessage(pt)))
  }
  tab <- pt$tab
  data.frame(
    Distance=label,
    F=as.numeric(tab[1, "F"]),
    P=as.numeric(tab[1, "Pr(>F)"]),
    Note=""
  )
}

# Bray-Curtis
bray <- vegdist(comm, method="bray")
bray_test <- run_adonis(bray, "Bray-Curtis")
bray_disp <- run_dispersion(bray, "Bray-Curtis")

# Bray PCoA
bray_ord <- wcmdscale(bray, k=2, eig=TRUE, add="lingoes")
bray_pts <- as.data.frame(bray_ord$points[, 1:2, drop=FALSE])
names(bray_pts) <- c("Axis1", "Axis2")
bray_pts$Run_accession <- rownames(bray_pts)
bray_pts <- merge(bray_pts, meta_beta[, c("Run_accession","Alcohol","Dataset")],
                  by="Run_accession", sort=FALSE)

pos_eig <- bray_ord$eig[bray_ord$eig > 0]
bray_var <- if (length(pos_eig) >= 2) 100 * pos_eig[1:2] / sum(pos_eig) else c(NA, NA)

p_bray <- ggplot(bray_pts, aes(x=Axis1, y=Axis2, color=Alcohol)) +
  geom_point(alpha=0.70, size=1.8) +
  stat_ellipse(level=0.95, linewidth=0.7, na.rm=TRUE) +
  labs(
    title="Bray-Curtis PCoA: ALCOHOL vs NON_ALCOHOL",
    x=sprintf("PCoA1 (%.1f%%)", bray_var[1]),
    y=sprintf("PCoA2 (%.1f%%)", bray_var[2])
  ) +
  theme_bw(base_size=12)

ggsave(file.path(out_dir, "beta_bray_pcoa.png"),
       p_bray, width=7.5, height=6, dpi=300)
ggsave(file.path(out_dir, "beta_bray_pcoa.pdf"),
       p_bray, width=7.5, height=6)

write.table(bray_pts, file.path(out_dir, "beta_bray_pcoa_coordinates.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# Aitchison: pseudocount + CLR + Euclidean
pseudo <- 1e-6
comm_pc <- comm + pseudo
comm_pc <- comm_pc / rowSums(comm_pc)
clr <- log(comm_pc) - rowMeans(log(comm_pc))
aitch <- dist(clr, method="euclidean")

aitch_test <- run_adonis(aitch, "Aitchison")
aitch_disp <- run_dispersion(aitch, "Aitchison")

pca <- prcomp(clr, center=TRUE, scale.=FALSE)
aitch_pts <- data.frame(
  Run_accession=rownames(pca$x),
  PC1=pca$x[,1],
  PC2=pca$x[,2],
  stringsAsFactors=FALSE
)
aitch_pts <- merge(aitch_pts, meta_beta[, c("Run_accession","Alcohol","Dataset")],
                   by="Run_accession", sort=FALSE)
pvar <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

p_aitch <- ggplot(aitch_pts, aes(x=PC1, y=PC2, color=Alcohol)) +
  geom_point(alpha=0.70, size=1.8) +
  stat_ellipse(level=0.95, linewidth=0.7, na.rm=TRUE) +
  labs(
    title="Aitchison PCA: ALCOHOL vs NON_ALCOHOL",
    x=sprintf("PC1 (%.1f%%)", pvar[1]),
    y=sprintf("PC2 (%.1f%%)", pvar[2])
  ) +
  theme_bw(base_size=12)

ggsave(file.path(out_dir, "beta_aitchison_pca.png"),
       p_aitch, width=7.5, height=6, dpi=300)
ggsave(file.path(out_dir, "beta_aitchison_pca.pdf"),
       p_aitch, width=7.5, height=6)

write.table(aitch_pts, file.path(out_dir, "beta_aitchison_pca_coordinates.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

beta_tests <- rbind(bray_test, aitch_test)
beta_disp <- rbind(bray_disp, aitch_disp)

write.table(beta_tests, file.path(out_dir, "beta_PERMANOVA_tests.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(beta_disp, file.path(out_dir, "beta_PERMDISP_tests.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# ============================================================
# Concise text summary
# ============================================================

summary_file <- file.path(out_dir, "DIVERSITY_SUMMARY.txt")
sink(summary_file)
cat("GUT_STOOL ALCOHOL vs NON_ALCOHOL diversity analysis\n")
cat("===================================================\n\n")
cat("Input unfiltered:", unfiltered_file, "\n")
cat("Input filtered:  ", filtered_file, "\n")
cat("Metadata:        ", metadata_file, "\n\n")
cat("Samples:", nrow(meta), "\n")
cat("ALCOHOL:", sum(meta$Alcohol=="ALCOHOL"), "\n")
cat("NON_ALCOHOL:", sum(meta$Alcohol=="NON_ALCOHOL"), "\n")
cat("Datasets:", length(unique(meta$Dataset)), "\n\n")

cat("Alpha diversity metrics:\n")
cat("  Observed richness\n")
cat("  Shannon\n")
cat("  Simpson (1-D)\n")
cat("  Chao1 (only where count/count-like standardized abundance is available)\n\n")
print(alpha_tests, row.names=FALSE)

cat("\nBeta diversity:\n")
cat("  Bray-Curtis\n")
cat("  Aitchison (CLR Euclidean)\n\n")
cat("PERMANOVA:\n")
print(beta_tests, row.names=FALSE)
cat("\nPERMDISP:\n")
print(beta_disp, row.names=FALSE)

cat("\nInterpretation rule:\n")
cat("  PERMANOVA tests centroid/compositional separation.\n")
cat("  PERMDISP tests group dispersion; a significant PERMDISP means a PERMANOVA result\n")
cat("  should be interpreted cautiously because dispersion differs between groups.\n")
cat("  Dataset_plus_Alcohol is the preferred pooled sensitivity model for the merged cohort.\n")
sink()

msg("Done. Results written to: ", out_dir)
msg("Summary: ", summary_file)

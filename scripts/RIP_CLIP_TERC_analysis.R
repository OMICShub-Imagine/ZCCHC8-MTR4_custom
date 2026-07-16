#!/usr/bin/env Rscript
# =============================================================================
# Locus-level TERC/Terc signal and top targets from public RIP-seq / iCLIP
# =============================================================================
#
# Purpose
#   Quantify publicly available ZCCHC8 RIP-seq (mouse) and RBM7 iCLIP (human)
#   signal at the TERC/Terc locus, and export top enriched/bound genes.
#
# Public data
#   - ZCCHC8 RIP-seq (mouse ES cells): GEO GSE127790
#   - RBM7 iCLIP (human HeLa):         GEO GSE63791
#
# Annotation
#   - Mouse: GENCODE vM9 (mm9 / GRCm38 matching GSE127790 bigWigs)
#   - Human: GENCODE hg19 gene annotation matching GSE63791 bigWigs
#
# Outputs (written to OUTDIR)
#   - RIP_Terc_signal_summary.tsv
#   - RIP_Terc_log2_enrichment.txt
#   - RIP_top100_enriched_genes_mouse.tsv
#   - CLIP_TERC_signal_summary_minus.tsv
#   - CLIP_top200_RBM7_targets.tsv
#   - sessionInfo.txt
#
# Usage
#   Rscript scripts/RIP_CLIP_TERC_analysis.R
#   # or from R: source("scripts/RIP_CLIP_TERC_analysis.R")
#
# =============================================================================

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(data.table)
})

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
# Resolve repository root from this script location when possible.
.args <- commandArgs(trailingOnly = FALSE)
.file_arg <- grep("^--file=", .args, value = TRUE)
ROOT <- if (length(.file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", .file_arg))))
} else if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
           rstudioapi::isAvailable() &&
           !is.null(rstudioapi::getSourceEditorContext())) {
  dirname(dirname(rstudioapi::getSourceEditorContext()$path))
} else {
  # Fallback: assume current working directory is the repository root.
  normalizePath(".", winslash = "/", mustWork = TRUE)
}

# Input directories / files -----------------------------------------------------
GSE127790_DIR <- file.path(ROOT, "data", "GSE127790")
GSE63791_DIR  <- file.path(ROOT, "data", "GSE63791")

GTF_MOUSE <- file.path(ROOT, "data", "annotation", "gencode.vM9.annotation.gtf")
GTF_HUMAN <- file.path(ROOT, "data", "annotation", "gencode.hg19.annotation.gtf")

BW_RIP1 <- file.path(GSE127790_DIR, "GSM3638745_RIP_Zcchc8_1.bw")
BW_RIP2 <- file.path(GSE127790_DIR, "GSM3638746_RIP_Zcchc8_2.bw")
BW_INP  <- file.path(GSE127790_DIR, "GSM3638747_RIP_Zcchc8_input.bw")

RIP_GENE_ENRICHMENT <- file.path(
  GSE127790_DIR, "mouse_Zcchc8_RIP_gene_enrichment.tsv"
)

BW_CLIP_MINUS <- file.path(
  GSE63791_DIR,
  c(
    "GSM1557476_NNNGGCG_R7_minus.bw",
    "GSM1557477_NNNTTGC_R7-2_minus.bw",
    "GSM1557478_NNNCCGG_R740_minus.bw",
    "GSM1557479_NNNGGTC_R740-2_minus.bw"
  )
)

RBM7_GENE_SCORES <- file.path(
  GSE63791_DIR, "RBM7_iCLIP_gene_scores_HGNC.tsv"
)

# Analysis parameters ---------------------------------------------------------
PSEUDOCOUNT <- 1e-6
TOP_RIP_GENES  <- 100L
TOP_CLIP_GENES <- 200L

# Output directory ------------------------------------------------------------
OUTDIR <- file.path(ROOT, "outputs")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  GTF_MOUSE, GTF_HUMAN,
  BW_RIP1, BW_RIP2, BW_INP,
  RIP_GENE_ENRICHMENT,
  BW_CLIP_MINUS,
  RBM7_GENE_SCORES
)
missing <- required_files[!file.exists(required_files)]
if (length(missing)) {
  stop(
    "Missing required input file(s).\n",
    "See data/README.md for download instructions.\n",
    paste0("  - ", missing, collapse = "\n")
  )
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
bw_summary <- function(bw, gr) {
  x <- import(bw, which = gr, as = "NumericList")
  v <- unlist(x, use.names = FALSE)
  v <- v[is.finite(v)]
  if (!length(v)) {
    return(data.table(
      track = basename(bw),
      mean = NA_real_, max = NA_real_, sum = NA_real_, n = 0L
    ))
  }
  data.table(
    track = basename(bw),
    mean = mean(v),
    max = max(v),
    sum = sum(v),
    n = length(v)
  )
}

get_gene_locus <- function(gtf_path, gene_symbols) {
  gtf <- import(gtf_path)
  genes <- gtf[gtf$type == "gene"]
  hit <- genes[genes$gene_name %in% gene_symbols]
  if (length(hit) == 0) {
    stop(
      "Gene not found in ", basename(gtf_path),
      ". Looked for: ", paste(gene_symbols, collapse = ", ")
    )
  }
  range(hit)
}

detect_symbol_column <- function(dt) {
  cand <- c(
    "hgnc_symbol", "HGNC_symbol", "hgnc", "HGNC",
    "symbol", "SYMBOL", "gene", "Gene", "gene_name", "GeneName"
  )
  hit <- intersect(cand, names(dt))
  if (!length(hit)) {
    # Fallback: first column that looks gene-like and is character/factor
    char_cols <- names(dt)[vapply(dt, function(x) {
      is.character(x) || is.factor(x)
    }, logical(1))]
    if (!length(char_cols)) {
      stop("Could not detect a gene-symbol column in RBM7 score table.")
    }
    return(char_cols[[1]])
  }
  hit[[1]]
}

message("Working directory / repo root: ", ROOT)
message("Writing outputs to: ", OUTDIR)

# =============================================================================
# 1) Mouse ZCCHC8 RIP-seq signal at Terc
# =============================================================================
message("1) Mouse RIP-seq: Terc locus signal (GSE127790)...")

terc_mm <- get_gene_locus(GTF_MOUSE, c("Terc", "TERC"))
message("  Terc locus: ", as.character(terc_mm))

terc_rip_stats <- rbindlist(list(
  bw_summary(BW_RIP1, terc_mm),
  bw_summary(BW_RIP2, terc_mm),
  bw_summary(BW_INP,  terc_mm)
))

rip_mean <- mean(
  terc_rip_stats[track %in% basename(c(BW_RIP1, BW_RIP2))]$mean,
  na.rm = TRUE
)
inp_mean <- terc_rip_stats[track == basename(BW_INP)]$mean
terc_log2_enrich <- log2(rip_mean + PSEUDOCOUNT) - log2(inp_mean + PSEUDOCOUNT)

terc_rip_stats[, sample_type := fcase(
  track == basename(BW_INP), "Input",
  track %in% basename(c(BW_RIP1, BW_RIP2)), "RIP",
  default = "other"
)]

fwrite(
  terc_rip_stats,
  file.path(OUTDIR, "RIP_Terc_signal_summary.tsv"),
  sep = "\t"
)
writeLines(
  sprintf("Terc_log2_enrichment_RIP_vs_Input\t%0.4f", terc_log2_enrich),
  file.path(OUTDIR, "RIP_Terc_log2_enrichment.txt")
)
message(sprintf("  Terc log2(RIP/Input) enrichment: %.4f", terc_log2_enrich))

# =============================================================================
# 2) Top mouse RIP-enriched genes (precomputed gene-level table)
# =============================================================================
message("2) Mouse RIP-seq: top enriched genes...")

sig_gene <- fread(RIP_GENE_ENRICHMENT)
needed <- c("log2_enrich", "RIP_mean")
missing_cols <- setdiff(needed, names(sig_gene))
if (length(missing_cols)) {
  stop(
    "Expected columns in ", basename(RIP_GENE_ENRICHMENT), ": ",
    paste(needed, collapse = ", "),
    ". Missing: ", paste(missing_cols, collapse = ", ")
  )
}

# Keep finite enrichments with detectable RIP signal (avoids Input == 0 artefacts)
sig_gene2 <- sig_gene[
  is.finite(log2_enrich) & is.finite(RIP_mean) & RIP_mean > 0
]
top_rip <- sig_gene2[order(-log2_enrich)][seq_len(min(TOP_RIP_GENES, .N))]
fwrite(
  top_rip,
  file.path(OUTDIR, "RIP_top100_enriched_genes_mouse.tsv"),
  sep = "\t"
)
message("  Wrote top ", nrow(top_rip), " RIP-enriched genes.")

# =============================================================================
# 3) Human RBM7 iCLIP signal at TERC (minus strand; TERC-encoded strand)
# =============================================================================
message("3) Human RBM7 iCLIP: TERC locus signal (GSE63791)...")

terc_hs <- get_gene_locus(GTF_HUMAN, c("TERC", "Terc"))
message("  TERC locus: ", as.character(terc_hs))

clip_terc_minus <- rbindlist(
  lapply(BW_CLIP_MINUS, bw_summary, gr = terc_hs)
)
fwrite(
  clip_terc_minus,
  file.path(OUTDIR, "CLIP_TERC_signal_summary_minus.tsv"),
  sep = "\t"
)
message(
  "  Mean minus-strand CLIP signal at TERC: ",
  sprintf("%.2f", mean(clip_terc_minus$mean, na.rm = TRUE))
)

# =============================================================================
# 4) Top RBM7 iCLIP targets (control mean binding strength)
# =============================================================================
message("4) Human RBM7 iCLIP: top targets by control mean...")

rbm7 <- fread(RBM7_GENE_SCORES)
sym_col <- detect_symbol_column(rbm7)

if (!("ctrl_mean" %in% names(rbm7))) {
  stop(
    "Expected column 'ctrl_mean' in ", basename(RBM7_GENE_SCORES),
    ". Available columns: ", paste(names(rbm7), collapse = ", ")
  )
}

setorder(rbm7, -ctrl_mean)
keep_cols <- c(sym_col, "ctrl_mean")
if ("kd_mean" %in% names(rbm7)) keep_cols <- c(keep_cols, "kd_mean")

top_clip <- rbm7[seq_len(min(TOP_CLIP_GENES, .N)), ..keep_cols]
setnames(top_clip, sym_col, "gene")
top_clip <- top_clip[!is.na(gene) & gene != ""]

fwrite(
  top_clip,
  file.path(OUTDIR, "CLIP_top200_RBM7_targets.tsv"),
  sep = "\t"
)
message("  Wrote top ", nrow(top_clip), " RBM7 iCLIP targets.")

# =============================================================================
# Session info
# =============================================================================
sink(file.path(OUTDIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("Done. Outputs in: ", OUTDIR)

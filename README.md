# Reanalysis of public ZCCHC8 RIP-seq and RBM7 iCLIP datasets

This repository reproduces the reanalysis described in Proper ZCCHC8-MTR4 Interaction Within the Nuclear Exosome Targeting (NEXT) Complex
is Critical for TERC Biogenesis and Telomere Homeostasis (de Tocqueville et al.)

## Requirements

R >= 4.3

Packages

- rtracklayer
- GenomicRanges
- data.table

## Input

- GSE127790
- GSE63791

Download processed bigWig files from GEO.

## Run

Rscript scripts/reanalyze_RIP_CLIP.R

## Output

- RIP_Terc_signal_summary.tsv
- RIP_top100_enriched_genes_mouse.tsv
- CLIP_TERC_signal_summary_minus.tsv
- CLIP_top200_RBM7_targets.tsv

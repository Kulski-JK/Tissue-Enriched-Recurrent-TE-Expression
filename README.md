# Tissue-Enriched and Recurrent Expression of Repetitive Elements Across 32 Human Tissues

This repository is being developed to support the manuscript **“Tissue-Enriched and Recurrent Expression of Repetitive Elements Across 32 Human Tissues.”**

It contains R scripts and documentation for analysing significant transposable-element (TE) expression results across 32 human tissues using a one-tissue-versus-all-others framework.

## Analysis overview

The analysis uses tissue-level result files named `TEresSig*.csv`. For each tissue, associations are classified as:

- **Tissue-enriched:** adjusted (p)-value (`padj`) < 0.05 and `log2FoldChange` > 0.
- **Tissue-depleted:** adjusted (p)-value (`padj`) < 0.05 and `log2FoldChange` < 0.
- **Cross-tissue recurrent enriched TE:** a harmonised enriched TE identity detected in at least two distinct tissues.

The recurrence analysis counts distinct tissues, not simply the number of rows. The original SalmonTE TE names are retained alongside a separate harmonised identity field. Harmonisation is a naming convention for analysis and does not establish that two targets are biologically identical.

## Data and directory layout

Place the 32 input CSV files in `data/32files/`. See [data/README.md](data/README.md) for input requirements and data provenance notes.

```text
.
├── README.md
├── data/
│   ├── README.md
│   └── 32files/       # Place the 32 TEresSig*.csv input files here
├── scripts/           # Reproducible R analysis scripts
└── outputs/           # Generated result tables (created by scripts)
```

The input result CSVs are not included in this repository at this stage. Do not assume that the files in another repository are the same TE inputs; confirm their location and permissions before copying or redistributing them.

## Software

The main analysis script uses R and the `dplyr` and `stringr` packages. Install these packages once in R if they are not already available:

```r
install.packages(c("dplyr", "stringr"))
```

Run scripts from the repository root so relative paths such as `data/32files` resolve correctly.

## Outputs

The analysis is designed to produce tables for the combined significant results, tissue-enriched associations, tissue-depleted associations, cross-tissue recurrence, within-tissue duplicate/alias identities, and a summary of counts.

## Scope and interpretation

These analyses describe differential TE expression in the supplied tissue-level results. “Tissue-enriched” is used rather than “tissue-specific”; positive and negative log2 fold changes indicate enrichment and depletion, respectively, under the one-versus-all comparison.

## Authors and status

Analysis scripts associated with work by **Sulev Koks and Jerzy Kulski**. Repository contents and documentation are being assembled; scripts and input provenance should be checked before treating this repository as a complete release.

No software licence has been added yet. Licensing and attribution will be clarified before a formal public release.

# Input data

## Expected input files

The analysis expects 32 CSV files in this directory:

`TEresSig*.csv`

Each file represents one tissue and should contain these columns:

- `TE`
- `baseMean`
- `log2FoldChange`
- `lfcSE`
- `stat`
- `pvalue`
- `padj`

The script derives the tissue label from each filename and adds it as a `Tissue` column.

## Provenance and access

The inputs are not currently committed to this repository. Before adding them, confirm their source, authorship, redistribution permissions, and whether the files are already available in the associated tissue-enriched gene-expression repository or another project location. Do not copy data merely because a related repository exists.

Place approved input files in `data/32files/`. Keep generated output tables in `outputs/` rather than mixing them with the input files.

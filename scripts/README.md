# Analysis scripts

## Primary analysis

[`analyze_all_significant_TEs_and_recurrence.R`](analyze_all_significant_TEs_and_recurrence.R) combines the 32 tissue-level `TEresSig*.csv` files, applies the adjusted-(p)-value and log2-fold-change criteria, retains the original SalmonTE names, creates a separate harmonised identity field, and summarises cross-tissue recurrence.

Run it from the repository root in R or RStudio. It expects the input files in `data/32files/` and writes its six CSV output tables to `outputs/`.

Required packages: `dplyr` and `stringr`. Install them once, if needed, using:

```r
install.packages(c("dplyr", "stringr"))
```

The script checks for 32 input files and stops if required columns are missing. Review the output and methods against the final manuscript before citing the repository as a formal release.

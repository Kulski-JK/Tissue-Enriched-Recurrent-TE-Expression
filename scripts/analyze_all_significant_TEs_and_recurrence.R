############################################################
# ANALYSIS: ALL SIGNIFICANT TEs
# Authors: Sulev Koks and Jerzy Kulski
# Working analysis script adapted for the public repository
# Date of original working script: 17.07.2026
#Uses ALL significant TE results (not just top-6)
# from the 32 TEresSig*.csv files.
#
# Main analysis:
#   Tissue-enriched TEs = padj < 0.05 & log2FoldChange > 0
#
# Secondary analysis:
#   Tissue-depleted TEs = padj < 0.05 & log2FoldChange < 0
#
# Outputs:
#   1. All_Significant_TE_Master.csv
#   2. All_Significant_Tissue_Enriched_TE.csv
#   3. All_Significant_Tissue_Depleted_TE.csv
#   4. All_Significant_Cross_Tissue_Recurrent_TE.csv
#   5. All_Significant_Within_Tissue_Duplicate_TE.csv
#   6. All_Significant_TE_Summary.csv
#
# Original SalmonTE names are retained.
############################################################


############################
# 1. Load 2 packages
############################

# Install once if needed: install.packages(c("dplyr", "stringr"))

library(dplyr)
library(stringr)

# 1b. Run 2 packages

packageVersion("dplyr")
packageVersion("stringr")

############################
# 2. Set input directory
############################

input_dir <- "data/32files"
output_dir <- "outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Input directory:\n", input_dir, "\n\n")


############################
# 3. Find the 32 TE result files
############################

files <- list.files(
  input_dir,
  pattern = "^TEresSig.*\\.csv$",
  full.names = TRUE
)

cat("Number of TE result files found:", length(files), "\n\n")

if (length(files) != 32) {
  warning(
    "Expected 32 TEresSig*.csv files, but found ",
    length(files)
  )
}


############################
# 4. Read each file
############################

TE_list <- lapply(files, function(f) {
  tissue <- sub(
    "^TEresSig(.*)\\.csv$",
    "\\1",
    basename(f)
  )
  
  x <- read.csv(
    f,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  x$Tissue <- tissue
  return(x)
})

############################
# 5. Check required columns
############################

required_columns <- c(
  "TE",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj",
  "Tissue"
)

for (i in seq_along(TE_list)) {
  
  missing_columns <- setdiff(
    required_columns,
    colnames(TE_list[[i]])
  )
  
  if (length(missing_columns) > 0) {
    
    stop(
      "Missing columns in file: ",
      basename(files[i]),
      "\nMissing: ",
      paste(missing_columns, collapse = ", ")
    )
    
  }
}


############################
# 6. Combine all 32 files
############################

all_TE <- bind_rows(TE_list)

cat(
  "Total rows across all 32 files:",
  nrow(all_TE),
  "\n\n"
)


############################
# 7. Basic validation
############################

cat(
  "Number of tissues represented:",
  n_distinct(all_TE$Tissue),
  "\n"
)

cat(
  "Number of unique SalmonTE TE names:",
  n_distinct(all_TE$TE),
  "\n\n"
)


############################
# 8. Define significance
############################

# Explicitly apply the significance criterion.
#
# Although the files are named TEresSig*, we retain the
# explicit padj < 0.05 filter for reproducibility.

all_TE <- all_TE %>%
  mutate(
    Significant = !is.na(padj) & padj < 0.05
  )

cat(
  "Total significant rows:",
  sum(all_TE$Significant, na.rm = TRUE),
  "\n"
)


############################
# 9. Separate enriched and depleted TEs
############################

TE_enriched <- all_TE %>%
  filter(
    Significant,
    !is.na(log2FoldChange),
    log2FoldChange > 0
  )

TE_depleted <- all_TE %>%
  filter(
    Significant,
    !is.na(log2FoldChange),
    log2FoldChange < 0
  )

cat(
  "Tissue-enriched TE associations:",
  nrow(TE_enriched),
  "\n"
)

cat(
  "Tissue-depleted TE associations:",
  nrow(TE_depleted),
  "\n"
)

cat(
  "Enriched + depleted:",
  nrow(TE_enriched) + nrow(TE_depleted),
  "\n"
)

############################
# 10. Identity harmonisation
############################
#
# IMPORTANT:
#
# The original SalmonTE TE name is NOT changed.
#
# TE_harmonised is a separate analytical field.
#
# For now we use only the two conservative transformations
# already used in Analysis 1:
#
#   LTR18B_dup2  -> LTR18B
#   HERV-K14I    -> HERVK14I
#
# and removal of a trailing "_dup2".
#
# This should be regarded as nomenclature harmonisation,
# NOT proof that the two SalmonTE targets are biologically
# identical.
############################

library(stringr)

harmonise_TE <- function(x) {
  
  x %>%
    str_replace("_dup2$", "") %>%
    str_replace("^HERV-K14I$", "HERVK14I")
}


TE_enriched <- TE_enriched %>%
  mutate(
    TE_harmonised = harmonise_TE(TE)
  )

TE_depleted <- TE_depleted %>%
  mutate(
    TE_harmonised = harmonise_TE(TE)
  )

cat(
  "Unique enriched SalmonTE names:",
  n_distinct(TE_enriched$TE),
  "\n"
)

cat(
  "Unique harmonised enriched TE identities:",
  n_distinct(TE_enriched$TE_harmonised),
  "\n"
)





############################
# 11. Save master tables
############################

write.csv(
  all_TE,
  file.path(
    output_dir,
    "All_Significant_TE_Master.csv"
  ),
  row.names = FALSE
)

write.csv(
  TE_enriched,
  file.path(
    input_dir,
    "All_Significant_Tissue_Enriched_TE.csv"
  ),
  row.names = FALSE
)

write.csv(
  TE_depleted,
  file.path(
    input_dir,
    "All_Significant_Tissue_Depleted_TE.csv"
  ),
  row.names = FALSE
)


############################
# 12. Cross-tissue recurrence
#
# IMPORTANT:
# Count distinct tissues, not simply number of rows.
#
# A TE occurring twice within one tissue should NOT
# automatically be interpreted as cross-tissue recurrence.
############################

cross_tissue <- TE_enriched %>%
  group_by(TE_harmonised) %>%
  summarise(
    
    Number_of_Tissues =
      n_distinct(Tissue),
    
    Tissues =
      paste(
        sort(unique(Tissue)),
        collapse = "; "
      ),
    
    Total_Topological_Positions =
      n(),
    
    Original_SalmonTE_Names =
      paste(
        sort(unique(TE)),
        collapse = "; "
      ),
    
    .groups = "drop"
    
  ) %>%
  
  filter(
    Number_of_Tissues >= 2
  ) %>%
  
  arrange(
    desc(Number_of_Tissues),
    TE_harmonised
  )

cat(
  "Cross-tissue recurrent enriched TE identities:",
  nrow(cross_tissue),
  "\n"
)

head(cross_tissue, 20)

############################
# 13. Save cross-tissue recurrence
############################

write.csv(
  cross_tissue,
  file.path(
    input_dir,
    "All_Significant_Cross_Tissue_Recurrent_TE.csv"
  ),
  row.names = FALSE
)


############################
# 14. Within-tissue duplicate / alias analysis
#
# These are cases where the harmonised identity occurs
# more than once within the SAME tissue.
#
# Example:
#   HERVK
#   HERVK_dup2
#
# Both may become:
#   HERVK
#
# but their original SalmonTE names are retained.
############################

within_tissue_duplicates <- TE_enriched %>%
  
  group_by(
    Tissue,
    TE_harmonised
  ) %>%
  
  summarise(
    
    Number_of_Rows =
      n(),
    
    Original_SalmonTE_Names =
      paste(
        sort(unique(TE)),
        collapse = "; "
      ),
    
    Ranks_by_log2FC =
      paste(
        sort(
          rank(
            -log2FoldChange,
            ties.method = "first"
          )
        ),
        collapse = "; "
      ),
    
    log2FC_values =
      paste(
        round(log2FoldChange, 3),
        collapse = "; "
      ),
    
    padj_values =
      paste(
        signif(padj, 4),
        collapse = "; "
      ),
    
    .groups = "drop"
    
  ) %>%
  
  filter(
    Number_of_Rows >= 2
  ) %>%
  
  arrange(
    Tissue,
    TE_harmonised
  )

cat(
  "Within-tissue duplicate/alias identities:",
  nrow(within_tissue_duplicates),
  "\n"
)

head(within_tissue_duplicates, 20)


############################
# 15. Save within-tissue duplicates
############################

write.csv(
  within_tissue_duplicates,
  file.path(
    input_dir,
    "All_Significant_Within_Tissue_Duplicate_TE.csv"
  ),
  row.names = FALSE
)


############################
# 16. Overall summary
############################

summary_table <- tibble(
  
  Metric = c(
    
    "Number of TE result files",
    "Number of tissues",
    "Total rows in all files",
    "Total significant rows",
    "Total tissue-enriched TE associations",
    "Total tissue-depleted TE associations",
    "Unique SalmonTE TE names in master data",
    "Unique enriched SalmonTE TE names",
    "Unique depleted SalmonTE TE names",
    "Unique harmonised enriched TE identities",
    "Cross-tissue recurrent enriched TE identities",
    "Within-tissue duplicate/alias identities"
    
  ),
  
  Value = c(
    
    length(files),
    
    n_distinct(all_TE$Tissue),
    
    nrow(all_TE),
    
    sum(all_TE$Significant, na.rm = TRUE),
    
    nrow(TE_enriched),
    
    nrow(TE_depleted),
    
    n_distinct(all_TE$TE),
    
    n_distinct(TE_enriched$TE),
    
    n_distinct(TE_depleted$TE),
    
    n_distinct(TE_enriched$TE_harmonised),
    
    nrow(cross_tissue),
    
    nrow(within_tissue_duplicates)
    
  )
  
)

print(summary_table)

############################
# 17. Save summary
############################

write.csv(
  summary_table,
  file.path(
    input_dir,
    "All_Significant_TE_Summary.csv"
  ),
  row.names = FALSE
)
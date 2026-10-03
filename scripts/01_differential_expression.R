## =============================================================
## Bulk RNA-seq Differential Expression: Osteoarthritis vs Normal
## Dataset: GSE114007 (via GREIN)
## =============================================================
## HOW TO USE: Run this script one STEP at a time (highlight the step's
## code and press Ctrl+Enter, or run line-by-line). Each step prints a
## checkpoint so you can see what's happening before moving to the next.
## If a step errors out, paste me the exact error + the printed output
## above it and I'll fix that step.
##
## Memory note: written to be light on RAM (8GB systems) — uses fread(),
## drops objects with rm() once no longer needed, and calls gc() after
## heavy steps.

## ---- STEP 0: Libraries ----
library(data.table)
library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)

## ---- STEP 1: Load the two files ----
counts_raw <- fread("data/GSE114007_GeneLevel_Raw_data.csv")
metadata   <- fread("data/GSE114007_filtered_metadata.csv")

## CHECKPOINT 1 — look at this output before continuing
cat("=== counts_raw: dimensions ===\n")
print(dim(counts_raw))
cat("=== counts_raw: first 4 column names ===\n")
print(colnames(counts_raw)[1:4])
cat("\n=== metadata: column names ===\n")
print(colnames(metadata))
cat("\n=== metadata: first 6 rows ===\n")
print(head(metadata))

## >>> STOP HERE and read the checkpoint output before running Step 2 <<<


## ---- STEP 2: Standardize metadata column names ----
## This step auto-detects common column name variants so it works whether
## GREIN exported "characteristics", "Characteristics", "oa grade", etc.

names(metadata) <- tolower(trimws(names(metadata)))   # lowercase + trim whitespace
print(colnames(metadata))  # confirm after cleaning

## Identify the sample ID column (first column, usually GSM ids)
setnames(metadata, old = names(metadata)[1], new = "sample_id")

## Identify the characteristics column (may be named "characteristics" or similar)
char_col <- grep("characteristic", names(metadata), value = TRUE)[1]
stopifnot(!is.na(char_col))  # stop with error if not found — tell me if this fails
setnames(metadata, old = char_col, new = "characteristics")

cat("Sample ID column renamed to 'sample_id'. Characteristics column found:", char_col, "\n")
print(head(metadata[, .(sample_id, characteristics)]))

## >>> STOP HERE and confirm sample_id + characteristics look right <<<


## ---- STEP 3: Build the group variable (Normal vs OA) ----
metadata[, group := ifelse(grepl("^Normal", characteristics, ignore.case = TRUE), "Normal",
                     ifelse(grepl("^OA", characteristics, ignore.case = TRUE), "OA", NA_character_))]
metadata[, group := factor(group, levels = c("Normal", "OA"))]

cat("Group counts:\n")
print(table(metadata$group, useNA = "always"))

## If you see any NA here, some characteristics values don't start with
## "Normal" or "OA" — paste me the unique values below so I can fix the pattern:
print(unique(metadata$characteristics[is.na(metadata$group)]))

stopifnot(all(!is.na(metadata$group)))  # must pass before continuing

## >>> STOP HERE — group counts should read Normal: 18, OA: 20 <<<


## ---- STEP 4: Prepare the counts matrix ----
## counts_raw structure: col1 = Ensembl ID, col2 = Gene symbol, then one
## column per sample (GSM IDs)
gene_id_col  <- names(counts_raw)[1]
gene_sym_col <- names(counts_raw)[2]

gene_annot <- counts_raw[, c(gene_id_col, gene_sym_col), with = FALSE]
setnames(gene_annot, c("ensembl_id", "gene_symbol"))

count_mat <- as.matrix(counts_raw[, -c(1,2), with = FALSE])
rownames(count_mat) <- counts_raw[[gene_id_col]]

rm(counts_raw); gc()

cat("count_mat dimensions:", dim(count_mat)[1], "genes x", dim(count_mat)[2], "samples\n")
print(colnames(count_mat)[1:5])

## >>> STOP HERE — confirm count_mat has 38 sample columns <<<


## ---- STEP 5: Align sample order between counts and metadata ----
common_samples <- intersect(colnames(count_mat), metadata$sample_id)
cat("Samples matched between counts and metadata:", length(common_samples), "\n")

## If this is NOT 38, the sample_id format differs between the two files
## (e.g. counts use GSM IDs, metadata uses GSM IDs but with extra text).
## Paste me: head(colnames(count_mat)) and head(metadata$sample_id)
print(head(colnames(count_mat)))
print(head(metadata$sample_id))

stopifnot(length(common_samples) == 38)

count_mat <- count_mat[, common_samples]
metadata  <- metadata[match(common_samples, metadata$sample_id)]
stopifnot(all(colnames(count_mat) == metadata$sample_id))

## >>> STOP HERE — must say "38" above before continuing <<<


## ---- STEP 6: Filter low-count genes ----
keep <- rowSums(count_mat >= 10) >= 3
count_mat <- count_mat[keep, ]
gene_annot <- gene_annot[gene_annot$ensembl_id %in% rownames(count_mat), ]
cat("Genes retained after filtering:", nrow(count_mat), "out of", length(keep), "\n")


## ---- STEP 7: Build the DESeq2 dataset ----
## age is included as a covariate (Normal donors ~36.6 yrs vs OA donors ~66.2 yrs)
age_col <- grep("^age$", names(metadata), value = TRUE)[1]
stopifnot(!is.na(age_col))

col_data <- data.frame(
  row.names = metadata$sample_id,
  group = metadata$group,
  age   = as.numeric(metadata[[age_col]])
)

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = col_data,
  design    = ~ age + group
)

rm(count_mat); gc()
cat("DESeqDataSet built successfully.\n")


## ---- STEP 8: Run DESeq2 ----
dds <- DESeq(dds)
cat("DESeq2 model fitting complete.\n")


## ---- STEP 9: Extract results (OA vs Normal) ----
res <- results(dds, contrast = c("group", "OA", "Normal"), alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = "group_OA_vs_Normal", type = "apeglm")

res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("ensembl_id") %>%
  left_join(gene_annot, by = "ensembl_id") %>%
  arrange(padj)

cat("Top 6 DE genes:\n")
print(head(res_df))


## ---- STEP 10: Save results ----
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
fwrite(res_df, "results/tables/DEGs_OA_vs_Normal.csv")

sig_degs <- res_df %>% filter(padj < 0.05, abs(log2FoldChange) >= 1)
fwrite(sig_degs, "results/tables/DEGs_significant_padj0.05_LFC1.csv")
cat("Significant DEGs (padj<0.05, |LFC|>=1):", nrow(sig_degs), "\n")


## ---- STEP 11: QC plots ----
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

vsd <- vst(dds, blind = FALSE)

pca_plot <- plotPCA(vsd, intgroup = c("group")) + theme_minimal() +
  ggtitle("PCA: Normal vs OA (GSE114007)")
ggsave("results/figures/PCA_plot.png", pca_plot, width = 6, height = 5, dpi = 150)

volcano_df <- res_df %>% filter(!is.na(padj))
volcano_plot <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.4, size = 0.8,
             color = ifelse(volcano_df$padj < 0.05 & abs(volcano_df$log2FoldChange) >= 1,
                             "red", "grey60")) +
  theme_minimal() +
  labs(title = "Volcano Plot: OA vs Normal", x = "log2 Fold Change", y = "-log10 adjusted p-value")
ggsave("results/figures/volcano_plot.png", volcano_plot, width = 6, height = 5, dpi = 150)

cat("Plots saved to results/figures/\n")


## ---- STEP 12: Cleanup ----
rm(dds, vsd); gc()
cat("\nAll done. Check results/tables/ and results/figures/\n")

## =============================================================
## GO / KEGG Pathway Enrichment: Osteoarthritis vs Normal DEGs
## Dataset: GSE114007 (via GREIN)
## Input: results/tables/DEGs_significant_padj0.05_LFC1.csv
##        (output of 01_differential_expression.R)
## =============================================================
## HOW TO USE: Run STEP BY STEP, same as script 01. Checkpoints print
## after each step — read them before continuing.
##
## Memory note: clusterProfiler + org.Hs.eg.db load a moderate amount of
## annotation data into memory (a few hundred MB) — fine for 8GB, but
## close other heavy apps/tabs while this runs, same as script 01.

## ---- STEP 0: Libraries ----
## If any of these aren't installed yet, run this first:
## if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
## BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"))

library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

## ---- STEP 1: Load DEG results ----
sig_degs <- fread("results/tables/DEGs_significant_padj0.05_LFC1.csv")
all_degs <- fread("results/tables/DEGs_OA_vs_Normal.csv")   # full tested gene list = background

cat("Significant DEGs loaded:", nrow(sig_degs), "\n")
cat("Background (all tested genes):", nrow(all_degs), "\n")
print(head(sig_degs))

## >>> STOP HERE — confirm sig_degs has 1248 rows <<<


## ---- STEP 2: Split into up- and down-regulated gene sets ----
## Enriching these separately often gives cleaner biological interpretation
## than lumping all DEGs together (e.g. "cartilage breakdown" vs "immune
## activation" pathways tend to split by direction in OA).
up_genes   <- sig_degs %>% filter(log2FoldChange > 0) %>% pull(ensembl_id)
down_genes <- sig_degs %>% filter(log2FoldChange < 0) %>% pull(ensembl_id)
all_sig    <- sig_degs$ensembl_id

cat("Up-regulated in OA:", length(up_genes), "\n")
cat("Down-regulated in OA:", length(down_genes), "\n")

## >>> STOP HERE — the two numbers above should sum to 1248 <<<


## ---- STEP 3: Convert Ensembl IDs to Entrez IDs ----
## GO/KEGG databases key on Entrez IDs, not Ensembl — this mapping step
## is required and will drop a small number of genes with no Entrez match
## (normal and expected, usually <5-10% loss).

convert_ids <- function(ensembl_ids) {
  bitr(ensembl_ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}

up_entrez   <- convert_ids(up_genes)
down_entrez <- convert_ids(down_genes)
bg_entrez   <- convert_ids(all_degs$ensembl_id)

cat("Up-regulated mapped to Entrez:", nrow(up_entrez), "/", length(up_genes), "\n")
cat("Down-regulated mapped to Entrez:", nrow(down_entrez), "/", length(down_genes), "\n")
cat("Background mapped to Entrez:", nrow(bg_entrez), "/", nrow(all_degs), "\n")

## >>> STOP HERE — check the mapping rates look reasonable (>85-90%) <<<


## ---- STEP 4: GO enrichment (Biological Process) — up-regulated genes ----
go_up <- enrichGO(
  gene          = up_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",          # Biological Process; can also try "MF", "CC"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

cat("GO-BP terms enriched (up-regulated):", nrow(as.data.frame(go_up)), "\n")
print(head(as.data.frame(go_up)))


## ---- STEP 5: GO enrichment (Biological Process) — down-regulated genes ----
go_down <- enrichGO(
  gene          = down_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

cat("GO-BP terms enriched (down-regulated):", nrow(as.data.frame(go_down)), "\n")
print(head(as.data.frame(go_down)))

## >>> STOP HERE — both up/down should return at least a handful of terms.
## If one comes back empty, that's biologically possible (small gene list,
## weak signal in that direction) — not necessarily an error. <<<


## ---- STEP 6: KEGG pathway enrichment (up- and down-regulated) ----
kegg_up <- enrichKEGG(
  gene          = up_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

kegg_down <- enrichKEGG(
  gene          = down_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

cat("KEGG pathways enriched (up-regulated):", nrow(as.data.frame(kegg_up)), "\n")
cat("KEGG pathways enriched (down-regulated):", nrow(as.data.frame(kegg_down)), "\n")

## KEGG results store Entrez IDs in the gene column by default; convert to
## readable gene symbols for the saved tables (readable=TRUE isn't
## supported directly in enrichKEGG, so we do it after the fact)
kegg_up   <- setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_down <- setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

## >>> STOP HERE before saving outputs <<<


## ---- STEP 7: Save result tables ----
dir.create("results/tables/enrichment", showWarnings = FALSE, recursive = TRUE)

fwrite(as.data.frame(go_up),   "results/tables/enrichment/GO_BP_upregulated.csv")
fwrite(as.data.frame(go_down), "results/tables/enrichment/GO_BP_downregulated.csv")
fwrite(as.data.frame(kegg_up),   "results/tables/enrichment/KEGG_upregulated.csv")
fwrite(as.data.frame(kegg_down), "results/tables/enrichment/KEGG_downregulated.csv")

cat("Enrichment tables saved to results/tables/enrichment/\n")


## ---- STEP 8: Plots ----
dir.create("results/figures/enrichment", showWarnings = FALSE, recursive = TRUE)

## Dotplots: top 15 terms each, up- and down-regulated, GO-BP
if (nrow(as.data.frame(go_up)) > 0) {
  p1 <- dotplot(go_up, showCategory = 15) + ggtitle("GO-BP: Up-regulated in OA")
  ggsave("results/figures/enrichment/GO_BP_up_dotplot.png", p1, width = 8, height = 7, dpi = 150)
}

if (nrow(as.data.frame(go_down)) > 0) {
  p2 <- dotplot(go_down, showCategory = 15) + ggtitle("GO-BP: Down-regulated in OA")
  ggsave("results/figures/enrichment/GO_BP_down_dotplot.png", p2, width = 8, height = 7, dpi = 150)
}

## KEGG dotplots
if (nrow(as.data.frame(kegg_up)) > 0) {
  p3 <- dotplot(kegg_up, showCategory = 15) + ggtitle("KEGG: Up-regulated in OA")
  ggsave("results/figures/enrichment/KEGG_up_dotplot.png", p3, width = 8, height = 7, dpi = 150)
}

if (nrow(as.data.frame(kegg_down)) > 0) {
  p4 <- dotplot(kegg_down, showCategory = 15) + ggtitle("KEGG: Down-regulated in OA")
  ggsave("results/figures/enrichment/KEGG_down_dotplot.png", p4, width = 8, height = 7, dpi = 150)
}

cat("Enrichment plots saved to results/figures/enrichment/\n")


## ---- STEP 9: Cleanup ----
rm(go_up, go_down, kegg_up, kegg_down); gc()

cat("\nEnrichment analysis complete.\n")
cat("Check results/tables/enrichment/ and results/figures/enrichment/\n")

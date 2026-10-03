# Bulk RNA-seq Analysis: Osteoarthritis vs Normal Cartilage

A reproducible bulk RNA-seq differential expression and pathway enrichment analysis comparing osteoarthritic (OA) and normal human knee cartilage, using public data retrieved via [GREIN](https://www.ilincs.org/apps/grein/).

## Dataset

- **Source**: GEO accession [GSE114007](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114007), retrieved via GREIN
- **Tissue**: Human knee cartilage
- **Samples**: 38 total — 18 Normal, 20 Osteoarthritic (OA)
- **Design note**: Normal and OA donors differ in mean age (~36.6 vs ~66.2 years), so age was included as a covariate in the differential expression model to reduce confounding.
- **Data type**: Gene-level raw counts (28,125 genes), downloaded directly from GREIN

## Workflow

1. **Differential expression analysis** (`scripts/01_differential_expression.R`)
   - Import raw counts and metadata from GREIN exports
   - Filter low-count genes
   - Model: `~ age + group` using DESeq2
   - Shrunk log2 fold changes, significance threshold: padj < 0.05, |log2FC| ≥ 1
2. **Pathway enrichment analysis** (`scripts/02_enrichment_analysis.R`)
   - Significant DEGs split into up- and down-regulated sets
   - Ensembl → Entrez ID mapping
   - GO Biological Process enrichment (clusterProfiler, against full tested-gene background)
   - KEGG pathway enrichment (same background)

## Results

- **1,248 significant DEGs** identified between OA and Normal cartilage (padj < 0.05, |log2FC| ≥ 1)
- PCA shows clear separation of Normal and OA samples along PC1 (50% of variance explained), supporting a strong, consistent transcriptional difference between groups
- GO-BP and KEGG enrichment returned significant terms/pathways for both up- and down-regulated gene sets (see `results/tables/enrichment/`)

| Output | Location |
|---|---|
| Full DE results | `results/tables/DEGs_OA_vs_Normal.csv` |
| Significant DEGs | `results/tables/DEGs_significant_padj0.05_LFC1.csv` |
| GO-BP enrichment (up/down) | `results/tables/enrichment/GO_BP_*.csv` |
| KEGG enrichment (up/down) | `results/tables/enrichment/KEGG_*.csv` |
| PCA plot | `results/figures/PCA_plot.png` |
| Volcano plot | `results/figures/volcano_plot.png` |
| Enrichment dotplots | `results/figures/enrichment/` |

## Repository Structure

```
.
├── data/                  # Raw counts + metadata from GREIN (GSE114007)
├── scripts/
│   ├── 01_differential_expression.R
│   └── 02_enrichment_analysis.R
├── results/
│   ├── figures/           # PCA, volcano, enrichment dotplots
│   └── tables/            # DEG and enrichment result tables
├── LICENSE
└── README.md
```

## Requirements

- R (≥ 4.x)
- Bioconductor packages: `DESeq2`, `apeglm` (or use `normal` shrinkage as a lighter alternative), `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`, `DOSE`
- CRAN packages: `data.table`, `dplyr`, `tibble`, `ggplot2`

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "apeglm", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"))
install.packages(c("data.table", "dplyr", "tibble", "ggplot2"))
```

## How to Reproduce

1. Clone this repository and open `bulk-rnaseq-meta-analysis-osteoarthritis.Rproj` in RStudio
2. Run `scripts/01_differential_expression.R`
3. Run `scripts/02_enrichment_analysis.R`
4. Outputs will populate in `results/figures/` and `results/tables/`

## Acknowledgements

This repository follows the analysis framework taught by [DeepBio Academy](https://github.com/deepbioacademy), applied here to an independently selected GREIN dataset.

## License

MIT — see `LICENSE`.

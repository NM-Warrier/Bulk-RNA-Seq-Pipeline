Bulk RNA-seq Downstream Analysis

A general R pipeline for gene annotation, protein-coding gene selection, count quality checks, normalization, PCA, and pairwise differential expression analysis using DESeq2, with optional batch adjustment.

## Requirements

- R
- Required packages: `DESeq2`, `dplyr`, `ggplot2`
- Optional packages: `biomaRt` for annotation retrieval, `edgeR` for filtering comparison, and `ashr` for log2 fold-change shrinkage

Install the required packages in R:

```r
install.packages(c("BiocManager", "dplyr", "ggplot2"))
BiocManager::install("DESeq2")
```

Optional packages:

```r
BiocManager::install(c("biomaRt", "edgeR"))
install.packages("ashr")
```

## Input files

Provide tab-separated files with headers.

| File | Required contents |
| --- | --- |
| Raw counts | First column: unique gene IDs. Remaining columns: sample IDs containing raw, nonnegative integer counts. |
| Sample metadata | `SampleID`, a phenotype column such as `Condition`, and `Batch` when batch adjustment is enabled. |
| Gene annotation | `gene_id`, `gene_symbol`, and `gene_biotype`, with one row per gene ID. |

Count column names must match metadata sample IDs exactly. Annotation must match the organism, gene identifier format, and annotation release used for counting. Resolve ambiguous annotation mappings before running.

Do not use TPM, FPKM, log-transformed, or previously normalized expression values as count input.

## Usage

1. Open `general_bulk_RNA_downstream_pipeline.R`.
2. Replace the working-directory and input-file placeholders.
3. Edit the comparison table with your phenotype column, numerator, and denominator. Add rows for additional comparisons.
4. Set `use_batch` to `TRUE` or `FALSE` as appropriate.
5. Review the annotation, filtering, and optional analysis settings.
6. Run the script in RStudio or from a terminal:

```bash
Rscript general_bulk_RNA_downstream_pipeline.R
```

To skip annotation, set `annotation_file <- NULL` and `protein_coding_only <- FALSE`. To retain annotated genes of all biotypes, set only `protein_coding_only <- FALSE`.

## Workflow

1. Validate raw counts and annotate gene IDs.
2. Optionally select protein-coding genes.
3. Export duplicate-symbol checks and a descriptive table of counts summed by symbol.
4. Match counts to sample metadata and apply explicitly specified exclusions.
5. Subset each pairwise comparison and validate its design.
6. Estimate DESeq2 size factors and filter low-count genes.
7. Fit the group effect, including batch when enabled and multiple batches are present.
8. Generate PCA and dispersion plots.
9. Export differential expression results, optional shrunken effects, and diagnostics.

Stable gene IDs are retained for statistical testing. The symbol-summed table is an additional export and is not used as DESeq2 input.

## Analysis settings

- Each comparison requires at least two independent biological replicates per group and residual degrees of freedom.
- Positive log2 fold changes indicate higher expression in the numerator group.
- The default prefilter retains genes with a normalized count sum of at least 10 within each comparison; this setting is editable.
- Significance-only tests use adjusted P-value thresholds of 0.05, 0.01, and 0.10.
- Additional tests evaluate effects greater than fold-change ratios of 1.5 and 2 at an adjusted P-value threshold of 0.05.
- DESeq2's default independent filtering and Cook-distance result filtering remain enabled. Automatic replacement of outlier counts is disabled.
- Optional shrinkage reports effect estimates separately from the original Wald test statistics.

Batch is modeled during differential expression analysis. PCA uses variance-stabilized values and retains batch variation. Adapt the design for paired samples, repeated measures, or additional covariates. Completely confounded batch and group effects cannot be separated by this pipeline.

The optional edgeR step compares filtering results; it does not select which statistical package should be used.

## Outputs

Results are saved under `RNAseq_results/`, including:

- Annotated and selected count tables
- Duplicate-symbol reports and symbol-summed counts
- Sample-matching and exclusion records
- Comparison-specific metadata and design matrices
- Size factors, filtering summaries, and normalized counts
- PCA scores, PCA plots, and dispersion plots
- Full, significant, upregulated, and downregulated result tables
- Cook-distance and model-convergence diagnostics
- Fitted DESeq2 objects, a combined DEG summary, and R session information

Use a fresh output directory for each analysis to avoid mixing results from different runs. Comparisons are fitted separately, and adjusted P-values are calculated within each result table rather than across all comparisons.

## Validation

The script has passed an R syntax check. Full execution and scientific validation require suitable input data and the required packages.

## Reference

[DESeq2 documentation and workflow](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
# Bulk-RNA-Seq-Pipeline

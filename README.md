# Thy1-haSyn-Gut-Brain

Code to accompany "Distinct patterns of gene expression changes in the colon and striatum of young mice overexpressing alpha-synuclein support Parkinson Disease as a multi-system process", Journal of Parkinson's Disease

-   RNA-seq processing follows recommendations from Lexogen for QuantSeq
-   Correction for sequencing artifacts and WGCNA code is heavily adapted from code published by the laboratory of Dan Geschwind, much of which can be found at [Geschwind Lab (github.com)](https://github.com/dhglab)

## Abbreviations

-   In the code, the abbreviation ASO for "alpha-synuclein overexpressing" is frequently used rather that Thy1-haSyn
-   dc distal colon
-   str striatum
-   cons consensus

## Data

Processed gene expression data is distributed as R `ExpressionSet` objects (`.rds` files). For each tissue × timepoint grouping, three objects are produced by the scripts in [`1-RNAseqWorkflow/C_QCandNormalization/`](1-RNAseqWorkflow/C_QCandNormalization) and committed to that grouping's `data/` folder (see the [datasets table](#datasets-available) below for direct links). They share the same samples and gene set but differ in how much processing has been applied.

| File | Description | When to use |
|---|---|---|
| `filteredEset.rds` | Raw counts, restricted to genes with **>5 counts in at least 30% of samples**. No normalization, no batch/technical correction. | **This is what to use if you want "raw counts."** It is the count matrix the rest of the pipeline starts from. We do not redistribute the unfiltered counts because low-expression genes are uninformative for QuantSeq and inflate multiple-testing burden in any downstream analysis. |
| `filtNormEset.rds` | `filteredEset` after DESeq2 variance-stabilizing transformation (`estimateSizeFactors` → `estimateDispersions` → `getVarianceStabilizedData`). No technical-covariate regression. | Intermediate object, exposed mainly for reproducibility. Most users should prefer `CorrectedEset.rds` below. |
| `CorrectedEset.rds` | `filtNormEset` after regressing out sequencing-artifact principal components derived from the Picard QC metrics. | **Recommended for downstream analysis.** A substantial amount of technical variance was removed at this step (see the per-dataset script for the specific PCs regressed, e.g. seqPCs 1, 2, 3, and 7 for `dc1m`). Differential expression, WGCNA, heatmaps, and all figures in the paper use this object. |

The PCs chosen for regression were selected per dataset by inspecting the association of each sequencing PC with biological covariates of interest versus its association with technical Picard metrics; PCs dominated by technical signal were removed. The selection is hard-coded in each tissue/timepoint script.

### Datasets available

The same three-object pattern is produced for each of the six tissue × timepoint groupings. Each row's **Data folder** link opens the directory containing that grouping's `filteredEset.rds`, `filtNormEset.rds`, and `CorrectedEset.rds`.

| Group | Tissue | Timepoint(s) | Script | Data folder |
|---|---|---|---|---|
| `dcAll` | distal colon | 1 mo + 3 mo | [dcAll.R](1-RNAseqWorkflow/C_QCandNormalization/dcAll.R) | [`dc/data/`](1-RNAseqWorkflow/C_QCandNormalization/dc/data) |
| `dc1m`  | distal colon | 1 mo         | [dc1m.R](1-RNAseqWorkflow/C_QCandNormalization/dc1m.R)   | [`dc1m/data/`](1-RNAseqWorkflow/C_QCandNormalization/dc1m/data) |
| `dc3m`  | distal colon | 3 mo         | [dc3m.R](1-RNAseqWorkflow/C_QCandNormalization/dc3m.R)   | [`dc3m/data/`](1-RNAseqWorkflow/C_QCandNormalization/dc3m/data) |
| `strAll`| striatum     | 1 mo + 3 mo | [strAll.R](1-RNAseqWorkflow/C_QCandNormalization/strAll.R) | [`str/data/`](1-RNAseqWorkflow/C_QCandNormalization/str/data) |
| `str1m` | striatum     | 1 mo         | [str1m.R](1-RNAseqWorkflow/C_QCandNormalization/str1m.R)  | [`str1m/data/`](1-RNAseqWorkflow/C_QCandNormalization/str1m/data) |
| `str3m` | striatum     | 3 mo         | [str3m.R](1-RNAseqWorkflow/C_QCandNormalization/str3m.R)  | [`str3m/data/`](1-RNAseqWorkflow/C_QCandNormalization/str3m/data) |

> Note: the `dcAll` and `strAll` scripts write to folders named `dc/` and `str/` (not `dcAll/` and `strAll/`) — the links above already point to the right place.

> **If you are looking for data:** for "raw counts," use `filteredEset.rds`. For *any* downstream analysis (differential expression, WGCNA, clustering, visualization, machine learning), please use `CorrectedEset.rds` instead — the uncorrected versions retain a meaningful amount of sequencing-batch signal that will leak into your results.

### Working with `ExpressionSet` objects

These files are R `ExpressionSet` objects from Bioconductor's [Biobase](https://bioconductor.org/packages/release/bioc/html/Biobase.html) package. To load one and pull out the expression matrix, sample metadata, or feature annotation:

```r
# install.packages("BiocManager"); BiocManager::install("Biobase")
library(Biobase)

eset <- readRDS("CorrectedEset.rds")

exprs(eset)     # gene × sample expression matrix
pData(eset)     # sample metadata (phenoData) as a data.frame
fData(eset)     # gene/feature annotation (featureData) as a data.frame
sampleNames(eset)
featureNames(eset)
```

For a full walk-through of the class and its accessors, see the Bioconductor vignette **["An Introduction to Bioconductor's ExpressionSet Class"](https://bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf)** (PDF) and the [Biobase reference manual](https://bioconductor.org/packages/release/bioc/manuals/Biobase/man/Biobase.pdf).

## Contents

*Code executed in the following order*

### 1. RNA-seq Workflow

#### A. Trim FASTQs

1.  bbduk - [A1_bbdukQS.sh](1-RNAseqWorkflow/A_TrimFastqs/A1_bbdukQS.sh)
2.  FastQC - [A2_runFastQC.sh](1-RNAseqWorkflow/A_TrimFastqs/A2runFastQC.sh)

#### B. Align reads with STAR and run Picard tools

1.  generate index - [B1_genomeGenMM22](1-RNAseqWorkflow/B_STARandPicard/B1_genomeGenMM22.sh)
2.  align - [B2_runStar.sh](1-RNAseqWorkflow/B_STARandPicard/B2_runStar.sh)
3.  samtools - [B3_samtoolsAndCountTable.sh](1-RNAseqWorkflow/B_STARandPicard/B3_samtoolsAndCountTable.sh)
4.  Picard tools - [B4_runPicard.sh](1-RNAseqWorkflow/B_STARandPicard/B4_runPicard.sh)
5.  Picard table - [B5_makePicardTable.R](1-RNAseqWorkflow/B_STARandPicard/B5_makePicardTable.R)
6.  Count table - [makeCountTable.R](1-RNAseqWorkflow/B_STARandPicard/B6_makeCountTable.R)


### C. QC and Normalization

1.  outlier detection - [QCfunctions.R](1-RNAseqWorkflow/C_QCandNormalization/QCfunctions.R) is sourced for the subsequent steps
2.  colon, both 1 and 3 months - [dcAll.R](1-RNAseqWorkflow/C_QCandNormalization/dcAll.R)
3.  colon, 1 month - [dc1m.R](1-RNAseqWorkflow/C_QCandNormalization/dc1m.R)
4.  colon, 3 months - [dc3m.R](1-RNAseqWorkflow/C_QCandNormalization/dc3m.R)
5.  striatum, both 1 and 3 months - [strAll.R](1-RNAseqWorkflow/C_QCandNormalization/strAll.R)
6.  striatum, 1 month - [str1m.R](1-RNAseqWorkflow/C_QCandNormalization/str1m.R)
7.  striatum, 3 months - [str3m.R](1-RNAseqWorkflow/C_QCandNormalization/str3m.R)

### 2. Analysis

1.  Differential Expression - [dea.R](2-Analysis/dea.R)

2.  Heatmaps - [expressionHeatmaps.R](2-Analysis/expressionHeatmaps.R)

3.  WGCNA

    1.  Network construction, module characterization
        -   Colon
            -   [a-NetworkConstruction.R](2-Analysis/WGCNA/Colon/a-NetworkConstruction.R)
            -   [b-ModuleCharacterization.R](2-Analysis/WGCNA/Colon/a-ModuleCharacterization.R)
        -   Striatum
            -   [a-NetworkConstruction.R](2-Analysis/WGCNA/Striatum/a-NetworkConstruction.R)
            -   [b-ModuleCharacterization.R](2-Analysis/WGCNA/Striatum/a-ModuleCharacterization.R)
        -   Consensus
            -   [ConsensusWGCNA.R](2-Analysis/WGCNA/cons/ConsensusWGCNA.R)
    2.  [CellTypeEnrichment.R](2-Analysis/WGCNA//CellTypeEnrichment.R)
    3.  [CytoscapeExport.R](2-Analysis/WGCNA/CytoscapeExport.R)
    4.  [ModulePreservation.R](2-Analysis/WGCNA/ModulePreservation.R)
    5.  [RiskGenes.R](2-Analysis/WGCNA/RiskGenes.R)

4.  GSEA - [GSEA-commandSetup.R](2-Analysis/GSEA-commandSetup.R) (writes commands for GSEA command line)

5.  Rank-rank hypergeometric overlap - [RRHO.R](2-Analysis/RRHO.R)

6.  Validation

    -   nCounter - [Nanstring.R](2-Analysis/Nanostring.R)
    -   Comparison with Gries et al - [matchGries.R](2-Analysis/matchGries.R)

7.  FinalProcessedData is a list containing processed corrected data and metadata

### Interactive WGCNA networks are available on NDEx Bio:

Full networks (All genes with module membership at least 0.6).

-   Colon: <https://doi.org/10.18119/N9D31X>

-   Striatum: <https://doi.org/10.18119/N93K77>

Networks of selected genes in selected notable modules:

-   Gut modules:

    -   Gut-Proliferation (cyan) module: <https://doi.org/10.18119/N98C9G>

    -   Gut-Proteolysis (darkmagenta) module: <https://doi.org/10.18119/N94K65>

    -   Gut-Metabolism (blue) module: <https://doi.org/10.18119/N90W48>

    -   Gut-ENS (yellow) module: <https://doi.org/10.18119/N9RC9T>

    -   Gut-Lysosome (greenyellow) module: <https://doi.org/10.18119/N9MK7W>

    -   Gut-Mito-QC (violet) module: <https://doi.org/10.18119/N97C8S>

-   Brain modules:

    -   Brain-Snca (cyan) module: <https://doi.org/10.18119/N9GW37>
    -   Brain-UPR (red) module: <https://doi.org/10.18119/N9C320>

Consensus networks

-   Colon adjacency: <https://doi.org/10.18119/N9V60M>

-   Striatum adjacency: <https://doi.org/10.18119/N9ZW3K>

-   Midnight blue module (colon adjacency): <https://doi.org/10.18119/N9W318>

### Respirometry

Analysis of respirometry data is included as a vignette in [videlock/SeahorseR (github.com)](https://github.com/videlock/SeahorseR)

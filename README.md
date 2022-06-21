# Thy1-haSyn-Gut-Brain

Code to accompany "Distinct patterns of gene expression changes in the colon and striatum of young mice overexpressing alpha-synuclein support Parkinson Disease as a multi-system process"

* RNA-seq processing follows recommendations from Lexogen for QuantSeq
* Correction for sequencing artifacts and WGCNA code is heavily adapted from code published by the laboratory of Dan Geschwind, much of which can be found at [Geschwind Lab (github.com)](https://github.com/dhglab) 

## Abbreviations
* In the code, the abbreviation ASO for "alpha-synuclein overexpressing" is frequently used rather that Thy1-haSyn
* dc distal colon
* str striatum
* cons consensus

## Contents
*Code executed in the following order*

### 1. RNA-seq Workflow 
#### A. Trim FASTQs  
1. bbduk - [A1_bbdukQS.sh](1-RNAseqWorkflow/A_TrimFastqs/A1_bbdukQS.sh)  
2. FastQC - [A2_runFastQC.sh](1-RNAseqWorkflow/A_TrimFastqs/A2_runFastQC.sh)  
#### B. Align reads with STAR and run Picard tools
1. generate index - [B1_genomeGenMM22](1-RNAseqWorkflow/B_STARandPicard/B1_genomeGenMM22.sh)
2. align - [B2_runStar.sh](1-RNAseqWorkflow/B_STARandPicard/B2_runStar.sh)
3. samtools - [B3_samtoolsAndCountTable.sh](1-RNAseqWorkflow/B_STARandPicard/B3_samtoolsAndCountTable.sh)
4. count table - [makeCountTable.R](1-RNAseqWorkflow/B_STARandPicard/makeCountTable.R)
5. Picard tools - [B4_runPicard.sh](1-RNAseqWorkflow/B_STARandPicard/B4_runPicard.sh)
6. Picard table - [B5_makePicardTable.R](1-RNAseqWorkflow/B_STARandPicard/B5_makePicardTable.R)
#### C. QC and Normalization
1. outlier detection - [outlierDetection.R](1-RNAseqWorkflow/C_QCandNormalization/outlierDetection.R)
2. colon, both 1 and 3 months - [dc_1and3m.R](1-RNAseqWorkflow/C_QCandNormalization/dc_1and3m.R)
3. colon, 1 month - [dc_1m.R](1-RNAseqWorkflow/C_QCandNormalization/dc_1m.R)
4. colon, 3 months - [dc_3m.R](1-RNAseqWorkflow/C_QCandNormalization/dc_3m.R)
5. striatum, both 1 and 3 months - [str_1and3m.R](1-RNAseqWorkflow/C_QCandNormalization/str_1and3m.R)
6. striatum, 1 month - [str_1m.R](1-RNAseqWorkflow/C_QCandNormalization/str_1m.R)
7. striatum, 3 months - [str_3m.R](1-RNAseqWorkflow/C_QCandNormalization/str_3m.R)
### 2. Analysis
1. Differential Expression - [dea.R](2-Analysis/dea.R)
2. Heatmaps - [expressionHeatmaps.R](2-Analysis/expressionHeatmaps.R)
3. WGCNA
	1. colon and striatum
          1. a-NetworkConstruction.R
          2. b-ModuleCharacterization.R
    2. Consensus - [ConsensusWGCNA.R](2-Analysis/WGCNA/cons/ConsensusWGCNA.R)
    3. [CellTypeEnrichment.R](2-Analysis/WGCNA//CellTypeEnrichment.R)
    4. [CytoscapeExport.R](2-Analysis/WGCNA/CytoscapeExport.R)
    5. [ModulePreservation.R](2-Analysis/WGCNA/ModulePreservation.R)
    6. [RiskGenes.R](2-Analysis/WGCNA/RiskGenes.R)
4. GSEA
5. Validation

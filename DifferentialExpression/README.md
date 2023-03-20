<img src="https://raw.githubusercontent.com/oicr-gsi/analysisProcedures/DEanalysisProcedure/docs/oicr_logo.png" width=25% height=25%>

# Differential Expression Analysis Procedure  

### Introduction
By sequencing RNA we can estimate gene abundance using tools like RSEM1, the next step in RNA-seq analysis is often to run a Differential Expression analysis to determine whether there are significant differences between the gene expression of groups of interests. This necessitates sequencing of replicates for statistical power. The goal of this procedure is to guide the analyst on running a Differential Expression analysis to convert gene counts to gene fold change with associated statistical significance.  

### Summary
Estimating differential gene expression between defined groups.   
Analysis: Alignment of RNA sequences to reference genome with STAR aligner, followed by quantification of RNA abundance per gene using RSEM, then differential expression is estimated using DESeq2.   

### Inputs  
- Gene counts: RSEM output from pipeline. Example: [TEST_0001_Ov_C_WT_S1.genes.results](inputs/TEST_0001_Ov_C_WT_S1.genes.results)
- Metadata: map between library names and groups or conditions to be compared. Example: [metadata](inputs/TEST_metadata.txt)
- Study design and experimental question, designs can be evaluated at case by case basis. Example:

<br />

|Design|Example|Design formula|Experimental Question|
|----------------|----------------|-------------------------------|-----------------------------|
|Single factor|Study of treatments and control|`design = ~ group`|What is the difference in gene expression between treatment and control?|
|Multiple factors|Study of treatments on different cell types|`design= ~ group` <br>**Convert to a single factor and create contrasts* |What are the differences between treatments within each cell type? What are the differences between cell types? |
|Accounting for effects not of interest|Time series experiment|`design= ~ batch + group`|What is the difference in gene expression between treatments accounting for batch effect?|
### Procedure
 1. Create metadata table  

	group IDs recorded in the File Provenance Report (FPR) may not be sufficient, collaborators will have to provide 	guidance on the groups to be compared.   

 2. Copy Deseq2.R to your working directory and update the paths to input, metadata and output

	`cp /.mounts/labs/gsiprojects/gsi/AnalysisProcedures/analysisProcedures/DifferentialExpression/script/`  
   `Deseq2.R /path/to/workingDir`

 3. Run the analysis and visualize results  
`$ module load deseqanalysis`
`$ Rscript Deseq2.R`

### Deliverables
 - Differentially expressed genes: a table with Ensembl gene IDs and respective log2FoldChange and pvalues

### Resources
RSEM - https://github.com/deweylab/RSEM
A guide to creating design matrices for gene expression experiments - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/
DE_QC - https://github.com/hbctraining/DGE_workshop_salmon/blob/master/lessons/03_DGE_QC_analysis.md
DESeq2 - http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
tximport - http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#RSEM

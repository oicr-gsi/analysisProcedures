#!/usr/bin/env Rscript

library(optparse)

# read user info
option_list = list(
  make_option(c("-i", "--inputDir"), type="character", default=NA, help="path to directory containing input RSEM data", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NA, help="path to metadata file"),
  make_option(c("-o", "--outputDir"),  type="character", default=NA, help="path to where output files will be written", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

# error if missing argument: if not, return an error
if (is.na(opt$i) || is.na(opt$m) || is.na(opt$o)) {
  print("Missing argument")
  print_help(opt_parser)
  stop("", call.=FALSE)
}

suppressMessages(library(DESeq2))
suppressMessages(library(tximport))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))

print("STEP 1: reading inputs")

# read inputs
#inputDir <- "/.mounts/labs/gsiprojects/gsi/AnalysisProcedures/analysisProcedures/DifferentialExpression/inputs/rsem"
inputDir <- opt$i
#outputDir <- "/.mounts/labs/gsiprojects/gsi/AnalysisProcedures/analysisProcedures/DifferentialExpression/outputs/"
outputDir <- opt$o
metadata <- read.csv(opt$m, sep = "\t")
#metadata <- read.csv("/.mounts/labs/gsiprojects/gsi/AnalysisProcedures/analysisProcedures/DifferentialExpression/inputs/TEST_metadata.txt", sep = "\t")
files <- file.path(inputDir, paste0(metadata$LibraryName, ".genes.results"))
#get project name from input rsem data, or update this manually
project_name <- strsplit(basename(files[1]), split = "_")[[1]][1]
names(files) <- paste0(metadata$LibraryName)

# import data with metadata
txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
rownames(metadata) <- colnames(txi$counts)
metadata$Group <-factor(metadata$Group)
#metadata[order(metadata$Group, decreasing = FALSE), ]

#remove genes that have zero length
zero_length = (apply(txi$length, 1, min) == 0)
txi$length = txi$length[!zero_length,]
txi$abundance = txi$abundance[!zero_length,]
txi$counts = txi$counts[!zero_length,]

print("STEP 2: creating QC plots")

#_________________QC_________________##
dir.create(outputDir)
dir.create(paste0(outputDir, "QC/"))

### QC: PCA
qc.dds <- DESeqDataSetFromTximport(txi, metadata, ~Group)
qc.dds <- estimateSizeFactors(qc.dds)
# Transform counts for data visualization PCA
rld <- rlog(qc.dds, blind=TRUE)
#shorten name
colnames(rld) <- gsub("TEST_000._Ov_C_WT_", "", colnames(rld) )
#colnames(rld) <- gsub("_I5.*", "", colnames(rld) )
pdf(paste0(outputDir, "QC/", project_name, "_", "QC_PCA.pdf"))
plotPCA(rld, intgroup="Group")+
  geom_text(aes(label=name),vjust=2,check_overlap = TRUE,size = 2)
dev.off()

### QC: Correlation matrix
rld <- rlog(qc.dds, blind=TRUE)
#shorten name
colnames(rld) <- gsub("TEST_000._Ov_C_WT_", "", colnames(rld) )
#edit names for metadata too
metadata.qc <- metadata
rownames(metadata.qc) <- gsub("TEST_000._Ov_C_WT_", "", rownames(metadata.qc) )
# Extract the rlog matrix from the object
rld_mat <- assay(rld)
# Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function
# Plot heatmap using the correlation matrix and the metadata object
heat.colors <- RColorBrewer::brewer.pal(6, "Blues")
pdf(paste0(outputDir, "QC/", project_name, "_", "QC_correlation_matrix.pdf"))
pheatmap(rld_cor, annotation = metadata.qc, border_color=NA, fontsize = 10,
         fontsize_row = 10, height=30)

print("STEP 3: running DE analysis")

#_______________DESeq2_______________##
# run DESeq2
design <- ~Group
print(paste0("Selected design formula ", toString(design)))
dds <- DESeqDataSetFromTximport(txi, metadata, design) # tidy = TRUE)
dds$Group <- relevel(dds$Group, ref = "sensitive") #ref indicates the baseline group
dds <- DESeq(dds)

# get results from dds object
results <- results(dds, alpha=0.05)
results
summary(results)
# remove ensemble ID version number
rownames(results) <- sub('\\.[0-9]*$', '', rownames(results))

# export results
resOrdered <- results[order(results$padj),] # order them by pvalue
write.csv(as.data.frame(resOrdered), file=paste0(outputDir, project_name, "_differential_expression.csv"))

#create MA plots
# note from vignette
#It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with
#log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
pdf(paste0(outputDir, project_name, "_", "MA_plot_shrunken_LFC.pdf"))
res_shrunken <- lfcShrink(dds, coef="Group_resistant_vs_sensitive")
plotMA(res_shrunken, ylim = c(-5, 5))
dev.off()


#_______________OPTIONAL_______________##
# Optional steps
## visualize the results, the top line indicates what is used as baseline, the last item (sensitive) is the baseline
#head(results)
#summary(results)
## log2 foldchange vs normalized counts
#plotMA(results, ylim=c(-2,2))
## counts for a particular gene
#plotCounts(dds, gene=which.min(results$padj), intgroup="group")
## variance stabilization
#vst <- varianceStabilizingTransformation(dds)
#vsd <- assay(vst) # the vsd object contains the normalized data matrix which can be used in downstream analysis

print(paste0("Analysis complete, data is in ", outputDir))

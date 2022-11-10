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

cat("\n#### STEP 1: reading inputs ####\n\n")

# read inputs
inputDir <- opt$i
outputDir <- opt$o
metadata <- read.csv(opt$m, sep = "\t")
files <- file.path(inputDir, paste0(metadata$LibraryName, ".genes.results"))

#get project name from input rsem data, or update this manually
project_name <- strsplit(basename(files[1]), split = "_")[[1]][1]
names(files) <- paste0(metadata$LibraryName)

# import data with metadata
txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
rownames(metadata) <- colnames(txi$counts)
metadata$Group <-factor(metadata$Group)

#remove genes that have zero length
zero_length = (apply(txi$length, 1, min) == 0)
txi$length = txi$length[!zero_length,]
txi$abundance = txi$abundance[!zero_length,]
txi$counts = txi$counts[!zero_length,]

#prepare data for DE analysis
design <- ~Group
cat(paste0("\nSelected design formula: ", gsub(",", "", toString(design)),"\n"))
cat("Loading data to DESeq2\n")
dds <- DESeqDataSetFromTximport(txi, metadata, design) # tidy = TRUE) #this will be used in QC and analysis

cat("\n#### STEP 2: creating QC plots ####\n\n")

#_________________QC_________________##
dir.create(outputDir, showWarnings = FALSE)
dir.create(paste0(outputDir, "QC/"), showWarnings = FALSE)

### QC: PCA
qc.dds <- estimateSizeFactors(dds)
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

cat("\n#### STEP 3: running DE analysis ####\n\n")

#_______________DESeq2_______________##
# run DESeq2
dds$Group <- relevel(dds$Group, ref = "sensitive") #ref indicates the baseline group
dds <- DESeq(dds)

# get results from dds object
results <- results(dds, alpha=0.05)
cat("\n")
head(results, n = 0L)
cat("\n")
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
## counts for a particular gene
#plotCounts(dds, gene=which.min(results$padj), intgroup="group")
## variance stabilization
#vst <- varianceStabilizingTransformation(dds)
#vsd <- assay(vst) # the vsd object contains the normalized data matrix which can be used in downstream analysis

cat("\n#### Analysis complete ####\n\n")

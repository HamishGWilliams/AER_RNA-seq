## This Script will attempt to automate much of the work completed in a 
## RNA-seq analysis, including installing and loading all necessary 
## packages, importing, transforming, and analysing data; then producing
## figures and exporting the all to a single folder which can be accessed
#Uploaded to GitHub

## 0. Create a new Directory to put figures and plots into: ----
setwd("C:\\Users\\hamis\\OneDrive\\Documents\\PhD\\GitHub\\AER_RNA-seq")

dir.create(file.path(getwd(),"figures")) ## Change name of "figures" to 
                                         ## your given choice name.

## 1. Check if package is installed, and will install + Loads packages ----
# Create list of packages needed:
pkg_list<-c("tidyverse", "DESeq2", "tibble", "stats", "EDASeq",
            "ggplot2", "pheatmap", "dplyr", "tidyr", "ggfortify", "corrplot",
            "gprofiler2", "knitr", "gProfileR", "gage", "RUVSeq", "compGenomRData")

# Build function...
check_packages <- function(pkg_list) {
  for (pkg in pkg_list) {
    if (!require(pkg, character.only = TRUE)) {
      print(paste(pkg, "is not installed, installing package..."))
      install.packages(pkg)
    } else {
      suppressMessages(library(pkg, character.only = TRUE))
      print(paste(pkg, "is installed and loaded"))
    }
  }
}

# Having issues with Rtools installation
library("devtools")
find_rtools(T)
# Need to install Rtools42, compatible with R 4.2.1

# Execute Function
check_packages(pkg_list)
check_packages(pkg_list)

# Other packages which need to be forced
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("gage", force = TRUE)
library(gage)

## 2. Load in Count data ----
# Create path to count data & read into object:
                           #!! Creates a full path for your count data file
counts_file <- system.file("extdata/rna-seq/SRP029880.raw_counts.tsv", 
                           #!! Replace the file path to YOUR file Path 
                           package = "compGenomRData")
counts <- as.matrix(read.table(counts_file, header = T, sep = '\t'))
counts <- as.matrix(read.table("data/example/SRP029880.raw_counts.tsv", header = T, sep = '\t'))


## 3. Calculate CPM -> geneLengths -> rpkm -> rpk -> TPM -> V (variance) ----
{
cpm <- {apply(subset(counts, select = c(-width)), 2, 
              function(x) x/sum(as.numeric(x)) * 10^6)} #Compute CPM values for each sample
geneLengths <- {as.vector(subset(counts, select = c(width)))
} # create a vector of gene lengths 
rpkm <- {apply(X = subset(counts, select = c(-width)),
              MARGIN = 2, 
              FUN = function(x) {
                10^9 * x / geneLengths / sum(as.numeric(x))
              })
}
rpk <- {apply( subset(counts, select = c(-width)), 2, 
              function(x) x/(geneLengths/1000))} # find gene length normalized values 
tpm <- {apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)
} # normalize by the sample size using rpk values
V <- apply(tpm, 1, var) # Calculate Variance
selectedGenes <- {names(V[order(V, decreasing = T)][1:100])
} # Sort by variance, select top 100 genes
}
## Check each step for errors:
{
colSums(cpm) # = 10^6 for ALL
colSums(rpkm) # Different values
colSums(tpm) # = 10^6 for ALL
}

## 4. Figures & Plots ----

## Heatmap for clustered genes
{
pheatmap(tpm[selectedGenes,], scale = 'row', show_rownames = FALSE)

## We can also add annotation tracks to the clusters

coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv",
                            package = "compGenomRData")
colData <- read.table(coldata_file, header = T, sep = '\t', 
                      stringsAsFactors = TRUE)


# Generate and save Annotated HeatMap

pheatmap(tpm[selectedGenes,], scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = colData)

dev.copy(png, file = file.path(getwd(),"figures/example/clustered_genes_heatmap.png"))
dev.off()

dev.copy(svg, file = file.path(getwd(),"figures/example/clustered_genes_heatmap.svg"))
dev.off()
}
## PCA plot
{
# Transformations:
# transpose the matrix 
M <- t(tpm[selectedGenes,])
# transform the counts to log2 scale 
M <- log2(M + 1)
# compute PCA 
pcaResults <- prcomp(M)
# Plot PCA
autoplot(pcaResults, data = colData, colour = 'group')
summary(pcaResults)

}
## Correlation plots
{
correlationMatrix <- cor(tpm)
corrplot(correlationMatrix, order = 'hclust', 
         addrect = 2, addCoef.col = 'white', 
         number.cex = 0.7) 

dev.copy(png, file = file.path(getwd(),"figures/example/correlation_plot.png"))
dev.off()

dev.copy(svg, file = file.path(getwd(),"figures/example/correlation_plot.svg"))
dev.off()

# Transforming the correlation plot into a heatmap figure
# split the clusters into two based on the clustering similarity 
pheatmap(correlationMatrix,  
         annotation_col = colData, 
         cutree_cols = 2)

dev.copy(png, file = file.path(getwd(),"figures/example/corr_heatmap.png"))
dev.off()

dev.copy(svg, file = file.path(getwd(),"figures/example/corr_heatmap.svg"))
dev.off()
}

## 5. Differential Expression Analysis ----

#remove the 'width' column
countData <- as.matrix(subset(counts, select = c(-width)))
#define the experimental setup 
colData <- read.table(coldata_file, header = T, sep = '\t', 
                      stringsAsFactors = TRUE)
designFormula <- "~ group"

## Build initial DEseq matrix
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = as.formula(designFormula))

## Remove genes that have almost no information in any give samples
dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]

# Now perform Differential expression analysis
dds <- DESeq(dds)

#compute the contrast for the 'group' variable where 'CTRL' 
#samples are used as the control group. 
DEresults = results(dds, contrast = c("group", 'CASE', 'CTRL'))
#sort results by increasing p-value
DEresults <- DEresults[order(DEresults$pvalue),]

## 6. Diagnostics plots ----

## MA plot
{
# A scatter plot where the X-axis denotes the average normalized counts across
# samples and the y-axis denotes the log fold change in the given contrast.

# Most points will be on the horizontal 0 line as most genes are not
# differentially expressed
DESeq2::plotMA(object = dds, ylim = c(-5, 5))
  dev.copy(png, file = file.path(getwd(),"figures/example/MA_plot.png"))
  dev.off()
  dev.copy(svg, file = file.path(getwd(),"figures/example/MA_plot.svg"))
  dev.off()
}
## P-value distribution
{
# shows distribution of the raw p-values
# we expect to see a peak ~ low p-value and a uniform distribution above 0.1
#!! OTHERWISE, adjustment for multiple testing does not work, and results
#!! are NOT meaningful.
  {
  ggplot(data = as.data.frame(DEresults), aes(x = pvalue)) + 
    geom_histogram(bins = 100)
  dev.copy(png, file = file.path(getwd(),"figures/example/PValue_Distribution.png"))
  dev.off()
  dev.copy(svg, file = file.path(getwd(),"figures/example/PValue_Distribution.svg"))
  dev.off()
  }
}
## PCA plot
{
# checks the biological reproducibility of the sample replicates in a PCA plot
# or a heatmap.to plot PCA results, we need to extract the normalized counts 
# from the DEseqDataSet object. It is possible to colour the points in the 
# scatter plot by the variable of interest, which helps to see if the replicates
# cluster well.

# # extract normalized counts from the DESeqDataSet object
# countsNormalized <- DESeq2::counts(dds, normalized = TRUE)
# 
# # select top 500 most variable genes
# selectedGenes <- names(sort(apply(countsNormalized, 1, var), 
#                             decreasing = TRUE)[1:500])
# 
# plotPCA(countsNormalized[selectedGenes,], 
#         col = as.numeric(colData$group), adj = 0.5, 
#         xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.6))

##!! The above code doesnt seem to work, it comes up with the error:
# Error in (function (classes, fdef, mtable)  : 
# unable to find an inherited method for function ‘plotPCA’ for signature ‘"matrix"’

##!! It seems that the plotPCA function is unable to have 'matrix' type objects
##!! within it, so just need to use the plot below:
  {
  rld <- rlog(dds)
  DESeq2::plotPCA(rld, ntop = 500, intgroup = 'group') + 
    ylim(-50, 50) + theme_bw()
  dev.copy(png, file = file.path(getwd(),"figures/example/PCA_plot.png"))
  dev.off()
  dev.copy(svg, file = file.path(getwd(),"figures/example/PCA_plot.svg"))
  dev.off()
  }
}
## Relative Log Expression (RLE) plot:
{
# useful in finding out if the data at hand needs normalization.
# Runs a quick dignostic to be applied on the raw or normliazed count
# matrices to see if further processing is required.

#install.packages("EDASeq")
  {
    par(mfrow = c(1, 2))
    plotRLE(countData, outline=FALSE, ylim=c(-4, 4), 
            col=as.numeric(colData$group), 
            main = 'Raw Counts')
    plotRLE(DESeq2::counts(dds, normalized = TRUE), 
            outline=FALSE, ylim=c(-4, 4), 
            col = as.numeric(colData$group), 
            main = 'Normalized Counts')
    par(mfrow = c(1, 1))
  dev.copy(png, file = file.path(getwd(),"figures/example/RLE_plot.png"))
  dev.off()
  dev.copy(svg, file = file.path(getwd(),"figures/example/RLE_plot.svg"))
  dev.off()
  }


# Here the RLE plot is comprised of boxplots, where each box-plot = 
# distribution of the relative log expression of the genes expressed in the 
# corresponding sample
# Each genes expression is divided by the median expression value of that gene
# across all samples

# Ideally the boxplots are centered around the horizontal zero line and are
# as tightly distributed as possible (Risso, Ngai, Speed, et al. 2014)
# We can observe here how the noramlized dataset has mproved upon the raw 
# count data for all the samples
# However, in some cases, it is important to visualise the RLE plots in
# combination with other diagnostic plots such as PCA plots, heatmaps,
# and correlation plots to see if there is more unwanted variation in the data,
# which can be further accounted for using packages such as RUVseq and sva.
}

## 7. Functional Enrichment Analysis ----

## GO term analysis
{
# A commonly used tool to group genes by related functions is to do enrichment
# analysis. This groups genes into sets by shared functional terms.
# To make sure descriptions of functions are the same, initiatives such as 
# Gene Ontology Consortium have collate a list of Gene Ontology (GO) terms
# for each gene. 

# Probably the most common analysis appllied after a differential expression
# analysis. GO term analysis helps to quickly find out systematic changes that 
# can describe differences between groups of samples.

# Select genes that are significantly differentially expressed between case
# and control samples. We will select geens with adjusted p-values < 0.1
# and that show a -ve/+ve 2-fold chaange in case vs. control.

# extract differential expression results
DEresults <- results(dds, contrast = c('group', 'CASE', 'CTRL'))

#remove genes with NA values 
DE <- DEresults[!is.na(DEresults$padj),]
#select genes with adjusted p-values below 0.1
DE <- DE[DE$padj < 0.1,]
#select genes with absolute log2 fold change above 1 (two-fold change)
DE <- DE[abs(DE$log2FoldChange) > 1,]

#get the list of genes of interest
genesOfInterest <- rownames(DE)

#calculate enriched GO terms
goResults <- gprofiler(query = genesOfInterest, 
                       organism = 'hsapiens', 
                       src_filter = 'GO', 
                       hier_filtering = 'moderate')

## gprofiler2 Package - not included in guide:
# calculate enriched GO terms
goResults2 <- gost(query = genesOfInterest, 
                 organism = "hsapiens", ordered_query = FALSE, 
                 multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                 measure_underrepresentation = FALSE, evcodes = FALSE, 
                 user_threshold = 0.05, correction_method = "g_SCS", 
                 domain_scope = "annotated", custom_bg = NULL, 
                 numeric_ns = "", sources = NULL)

gostres_link <- gost(query = genesOfInterest, 
                     as_short_link = TRUE)
# https://biit.cs.ut.ee/gplink/l/E7naeUS5Su

gostplot(goResults2, capped = TRUE, interactive = TRUE)
# Use the interface in the Viewwer to save ->
}
## Gene set enrichment analysis
{
# valuable exploratory analysis tool that can associate systematic changes to a 
# high-level function rather than individual genes. 

#Let's define the first gene set as the list of genes from one of the
#significant GO terms found in the GO analysis. order go results by pvalue
goResults <- goResults[order(goResults$p.value),]
#restrict the terms that have at most 100 genes overlapping with the query
go <- goResults[goResults$overlap.size < 100,]
# use the top term from this table to create a gene set 
geneSet1 <- unlist(strsplit(go[1,]$intersection, ','))

#Define another gene set by just randomly selecting 25 genes from the counts
#table get normalized counts from DESeq2 results
normalizedCounts <- DESeq2::counts(dds, normalized = TRUE)
geneSet2 <- sample(rownames(normalizedCounts), 25)

geneSets <- list('top_GO_term' = geneSet1,
                 'random_set' = geneSet2)

# Using the defined gene sets, we'd like to do a group comparison between the 
# case samples with respect to the control samples

# use the normalized counts to carry out a GSEA. 
gseaResults <- gage(exprs = log2(normalizedCounts+1), 
                    ref = match(rownames(colData[colData$group == 'CTRL',]), 
                                colnames(normalizedCounts)), 
                    samp = match(rownames(colData[colData$group == 'CASE',]), 
                                 colnames(normalizedCounts)),
                    gsets = geneSets, compare = 'as.group')

# WE can see if there is significant up-regulation or down-regulation of the 
# gene set in the case group by accessing $greater or $less:
gseaResults$greater
gseaResults$less

# The random gene set shows no significant changes in regulation, but the 
# top GO term shows a significant up-regulation. Lets use a heatmap to see
# this:

# get the expression data for the gene set of interest
M <- normalizedCounts[rownames(normalizedCounts) %in% geneSet1, ]
# log transform the counts for visualization scaling by row helps visualizing
# relative change of expression of a gene in multiple conditions
{
pheatmap(log2(M+1), 
         annotation_col = colData, 
         show_rownames = TRUE, 
         fontsize_row = 8,
         scale = 'row', 
         cutree_cols = 2, 
         cutree_rows = 2)
dev.copy(png, file = file.path(getwd(),"figures/example/GSEA_Heatmap.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/example/GSEA_Heatmap.svg"))
dev.off()
}
}

## 8. Accounting for additional sources of variation ----

# Sometimes unforeseen variables may contribute to some of the variation in 
# our results (e.g. batch number, temperature of sample storage, etc...)

# First we use DESeq2 to account for possible sources of variation when they are
# known.

## Accounting for covariates using DESeq2
#!! USing a different dataset for this!!
{
counts_file <- system.file('extdata/rna-seq/SRP021193.raw_counts.tsv', 
                           package = 'compGenomRData')
colData_file <- system.file('extdata/rna-seq/SRP021193.colData.tsv', 
                            package = 'compGenomRData')

counts <- read.table(counts_file)
colData <- read.table(colData_file, header = T, sep = '\t', 
                      stringsAsFactors = TRUE)

# Look at how samples are clustered by calculating TPM as a heatmap:

#find gene length normalized values 
geneLengths <- counts$width
rpk <- apply( subset(counts, select = c(-width)), 2, 
              function(x) x/(geneLengths/1000))
#normalize by the sample size using rpk values
tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)

selectedGenes <- names(sort(apply(tpm, 1, var), 
                            decreasing = T)[1:100])
}
## Heatmap plot for covariates
{
pheatmap(tpm[selectedGenes,], 
         scale = 'row',
         annotation_col = colData, 
         show_rownames = FALSE)

dev.copy(png, file = file.path(getwd(),"figures/example/Covariate_heatmap.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/example/Covariate_heatmap.svg"))
dev.off()
# WE can see that library section variable is the dominating variable, rather
# than the 'diagnostic variable'
}
## Accounting for knwon Covariates
{
# Now we instruct DESeq2 to account for the library selection variable:
# remove the 'width' column from the counts matrix
countData <- as.matrix(subset(counts, select = c(-width)))
# set up a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = ~ LibrarySelection + group)

# RUn the DE
dds <- DESeq(dds)
# extract results
DEresults <- results(dds, contrast = c('group', 'CASE', 'CTRL'))
}

## 9. Accounting for Unknown Covariates with RUVseq ----

# Use tools such as RUVseq or sva to estimate potential sources of variation and
# clean up the counts table from those sources of variation. Later on, the 
# estimate covariates can be integrated into DESeq2's design formula.

# Use RUVseq first to diagnose the problem and solve it. 
    ## !! Using a New dataset for this part !!

counts_file <- system.file('extdata/rna-seq/SRP049988.raw_counts.tsv', 
                           package = 'compGenomRData')
colData_file <- system.file('extdata/rna-seq/SRP049988.colData.tsv', 
                            package = 'compGenomRData')

counts <- read.table(counts_file)
colData <- read.table(colData_file, header = T, 
                      sep = '\t', stringsAsFactors = TRUE)
# simplify condition descriptions
colData$source_name <- ifelse(colData$group == 'CASE', 
                              'EHF_overexpression', 'Empty_Vector')

# Making heatmaps of samples using TPM counts:
#find gene length normalized values 
geneLengths <- counts$width
rpk <- apply( subset(counts, select = c(-width)), 2, 
              function(x) x/(geneLengths/1000))
#normalize by the sample size using rpk values
tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)
selectedGenes <- names(sort(apply(tpm, 1, var), 
                            decreasing = T)[1:100])
{
pheatmap(tpm[selectedGenes,], 
         scale = 'row',
         annotation_col = colData, 
         cutree_cols = 2, 
         show_rownames = FALSE)
dev.copy(png, file = file.path(getwd(),"figures/example/Unknown_Covariate_heatmap.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/example/Unknown_Covariate_heatmap.svg"))
dev.off()
}
# CASE_5 clusters more closely to the control samples. Could be the result of
# some batch effect, or any other technical preparation steps. However, the 
# colData object doesnt contain any variables that we can use to pinpoint the
# exact cause of this. So we'll use RUVseq to estimate potential covariates
# to see if the clustering results can be improve

# Set up the experiment:

# remove 'width' column from counts
countData <- as.matrix(subset(counts, select = c(-width)))
# create a seqExpressionSet object using EDASeq package 
set <- newSeqExpressionSet(counts = countData,
                           phenoData = colData)

# Make a diagnostic RLE plot on the raw count data:

# make an RLE plot and a PCA plot on raw count data and color samples by group
{
par(mfrow = c(1,2))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=as.numeric(colData$group))
plotPCA(set, col = as.numeric(colData$group), adj = 0.5, 
        ylim = c(-0.7, 0.5), xlim = c(-0.5, 0.5))
dev.copy(png, file = file.path(getwd(),"figures/example/RLE_and_PCAPlot.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/example/RLE_and_PCAPlot.svg"))
dev.off()
par(mfrow = c(1,1))
}

## make RLE and PCA plots on TPM matrix 
{
par(mfrow = c(1,2))
plotRLE(tpm, outline=FALSE, ylim=c(-4, 4), col=as.numeric(colData$group))
plotPCA(tpm, col=as.numeric(colData$group), adj = 0.5, 
        ylim = c(-0.3, 1), xlim = c(-0.5, 0.5))
dev.copy(png, file = file.path(getwd(),"figures/example/RLE_PCA_TPMMatrix.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/example/RLE_PCA_TPMMatrix.svg"))
dev.off()
par(mfrow = c(1,1))
}

# Both RLE and PCA plots look better on normalized data, but still suggests the 
# necessity to further improvement, since case_5 sample still clusters with the
# conrol samples. We haven't yet accounted for the source of the variation.

## Removing unwanted variation from the data

## Using RUVg

# Using HKGs as negative controls and as a refernces to correct systematic
# biases in the data, we should be able to remove unwanted variation.
# Let’s use a list of ~500 house-keeping genes compiled here: 
# https://www.tau.ac.il/~elieis/HKG/HK_genes.txt.

#source for house-keeping genes collection:
#https://m.tau.ac.il/~elieis/HKG/HK_genes.txt
HK_genes <- read.table(file = system.file("extdata/rna-seq/HK_genes.txt", 
                                          package = 'compGenomRData'), 
                       header = FALSE)
# let's take an intersection of the house-keeping genes with the genes available
# in the count table
house_keeping_genes <- intersect(rownames(set), HK_genes$V1)

# Now will run RUVg() with different number of factors of unwanted variation.
# Plot the PCA after removing the unwante variation, we should see which k
# values, numbe of factors, produce better seperation between sample groups.

# now, we use these genes as the empirical set of genes as input to RUVg.
# we try different values of k and see how the PCA plots look 
{
par(mfrow = c(2, 2))
for(k in 1:4) {
  set_g <- RUVg(x = set, cIdx = house_keeping_genes, k = k)
  plotPCA(set_g, col=as.numeric(colData$group), cex = 0.9, adj = 0.5, 
          main = paste0('with RUVg, k = ',k), 
          ylim = c(-1, 1), xlim = c(-1, 1), )
  dev.copy(png, file = file.path(getwd(),"figures/example/RUVg_plot.png"))
  dev.off()
  dev.copy(svg, file = file.path(getwd(),"figures/example/RUVg_plot.svg"))
  dev.off()
  par(mfrow = c(1,1))
  
  }
}

# We can see that k = 1 is enough to separate the case vs. control. So now we
# re-run the RUVg() function with the HKGs to do more diagnostic plots:

# choose k = 1

set_g <- RUVg(x = set, cIdx = house_keeping_genes, k = 1)

## Diagnostic plots

# RLE plots
{
par(mfrow = c(1,2))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colData$group), main = 'without RUVg')
plotRLE(set_g, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colData$group), main = 'with RUVg')
dev.copy(png, file = file.path(getwd(),"figures/example/RLE_Adjusted.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/example/RLE_Adjusted.svg"))
dev.off()
par(mfrow = c(1,1))
}

# PCA plots
{
par(mfrow = c(1,2))
plotPCA(set, col=as.numeric(colData$group), adj = 0.5,
        main = 'without RUVg', 
        ylim = c(-1, 0.5), xlim = c(-0.5, 0.5))
plotPCA(set_g, col=as.numeric(colData$group), adj = 0.5, 
        main = 'with RUVg',
        ylim = c(-1, 0.5), xlim = c(-0.5, 0.5))
dev.copy(png, file = file.path(getwd(),"figures/example/PCA_adjusted.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/example/PCA_adjusted.svg"))
dev.off()
par(mfrow = c(1,1))
}

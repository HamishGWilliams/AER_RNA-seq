## This script will contain the code used to perform the DEA of the RNA-seq
## read data which has been aligned and output into the appropriate format
## for this analyses. 

# Set working directory:
setwd("C:\\Users\\hamis\\OneDrive\\Documents\\PhD\\GitHub\\AER_RNA-seq")

# Create folders for the output files: 
dir.create(file.path(getwd(),"figures")) ## Make figures folder
dir.create(file.path("figures","exp1data")) ## make exp1data folder

# Create separate folders for png and svg plots: 
dir.create(file.path("figures/exp1data","png_plots")) ## make folder for pngs
dir.create(file.path("figures/exp1data","svg_plots")) ## make folder for svgs

## Building, installing and loading required packages: ----

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

# Create list of packages needed:
pkg_list<-c("tidyverse", "DESeq2", "tibble", "stats", "EDASeq",
            "ggplot2", "pheatmap", "dplyr", "tidyr", "ggfortify", 
            "corrplot", "gprofiler2", "knitr", "gProfileR", 
            "gage", "RUVSeq", "compGenomRData")

# Execute function (x2):
check_packages(pkg_list)
check_packages(pkg_list)

av <- available.packages("DESeq2")
av[av[, "Package"] == pkg, ]

# Other packages which need to be forced
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("gage", force = TRUE)
BiocManager::install("DESeq2", force = TRUE)
BiocManager::install("RUVSeq", force = TRUE)
library(gage)

## Load in Count Data: ----
counts <- as.matrix(read.table("data/A_Equina/A_Equina_Counts_redo.tsv", header = T, sep = '\t'))
str(counts)
Counts_only <- subset(counts, select = c(-Chr,-Start,-End,-Strand,-Length))
str(Counts_only[2,])

countslong<- subset(counts[,c(1:5,13:26)])
countsshort <- subset(counts[,c(1:5,6:12,27:33)])

## Calculate Parameters: ----
{
# CPM
cpm <- {apply(subset(counts, select = c(-Chr,-Start,-End,-Strand,-Length)), 2, 
                function(x) as.numeric(x)/sum(as.numeric(x)) * 10^6)}
cpmlong <- {apply(subset(countslong, select = c(-Chr,-Start,-End,-Strand,-Length)), 2, 
              function(x) as.numeric(x)/sum(as.numeric(x)) * 10^6)}
cpmshort <- {apply(subset(countsshort, select = c(-Chr,-Start,-End,-Strand,-Length)), 2, 
              function(x) as.numeric(x)/sum(as.numeric(x)) * 10^6)}

counts_rownames <- as.vector(row.names(counts))
rownames(cpm) <- counts_rownames
rownames(cpmlong) <- counts_rownames
rownames(cpmshort) <- counts_rownames

# Genelengths
geneLengths <- as.vector(subset(counts, select = c(Length)))
geneLengths <- as.vector(as.numeric(geneLengths))
  
# RPKM
rpkm <- {apply(X = subset(counts, select = c(-Chr,-Start,-End,-Strand,-Length)),
                 MARGIN = 2, 
                 FUN = function(x) {
                   10^9 * as.numeric(x) / geneLengths / sum(as.numeric(x))
                 })
}
rownames(rpkm) <- counts_rownames

rpkmlong <- {apply(X = subset(countslong, select = c(-Chr,-Start,-End,-Strand,-Length)),
               MARGIN = 2, 
               FUN = function(x) {
                 10^9 * as.numeric(x) / geneLengths / sum(as.numeric(x))
               })
}
rownames(rpkmlong) <- counts_rownames

rpkmshort <- {apply(X = subset(countsshort, select = c(-Chr,-Start,-End,-Strand,-Length)),
               MARGIN = 2, 
               FUN = function(x) {
                 10^9 * as.numeric(x) / geneLengths / sum(as.numeric(x))
               })
}
rownames(rpkmshort) <- counts_rownames

# RPK
rpk <- {apply(subset(counts, select = c(-Chr,-Start,-End,-Strand,-Length)), 2, 
                 function(x) as.numeric(x)/(geneLengths/1000))}
rownames(rpk) <- counts_rownames

rpklong <- {apply(subset(countslong, select = c(-Chr,-Start,-End,-Strand,-Length)), 2, 
              function(x) as.numeric(x)/(geneLengths/1000))}
rownames(rpklong) <- counts_rownames

rpkshort <- {apply(subset(countsshort, select = c(-Chr,-Start,-End,-Strand,-Length)), 2, 
              function(x) as.numeric(x)/(geneLengths/1000))}
rownames(rpkshort) <- counts_rownames


# TPM
tpm <- apply(rpk, 2, function(x) as.numeric(x) / sum(as.numeric(x)) * 10^6)
rownames(tpm) <- counts_rownames

tpmlong <- apply(rpklong, 2, function(x) as.numeric(x) / sum(as.numeric(x)) * 10^6)
rownames(tpmlong) <- counts_rownames

tpmshort <- apply(rpkshort, 2, function(x) as.numeric(x) / sum(as.numeric(x)) * 10^6)
rownames(tpmshort) <- counts_rownames

# Variance
V <- apply(tpm, 1, var)
Vlong <- apply(tpmlong, 1, var)
Vshort <- apply(tpmshort, 1, var)
}

# Select top 100 genes: !! CHANGE THIS TO DESIRED NUMBER !!
selectedGenes <- names(V[order(V, decreasing = T)][1:100])
selectedGeneslong <- names(Vlong[order(Vlong, decreasing = T)][1:100])
selectedGenesshort <- names(Vshort[order(Vshort, decreasing = T)][1:100])


## Check each step for errors:
colSums(cpm) # = 10^6 for ALL
colSums(rpkm) # Different values
colSums(tpm) # = 10^6 for ALL

colSums(cpmlong) # = 10^6 for ALL
colSums(rpkmlong) # Different values
colSums(tpmlong) # = 10^6 for ALL

colSums(cpmshort) # = 10^6 for ALL
colSums(rpkmshort) # Different values
colSums(tpmshort) # = 10^6 for ALL

## 4. Figures & Plots ----

## Heatmap for clustered genes
pheatmap(tpm[selectedGenes,], scale = 'row', show_rownames = FALSE)
pheatmap(tpmlong[selectedGeneslong,], scale = 'row', show_rownames = FALSE)
pheatmap(tpmshort[selectedGenesshort,], scale = 'row', show_rownames = FALSE)


colData <- read.table("data/A_Equina/colData.tsv", header = T, sep = '\t', 
                      stringsAsFactors = TRUE)
colDatalong <- subset(colData[8:21,])
colDatashort <- subset(colData[c(1:7,22:28),])



pheatmap(tpm[selectedGenes,], scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = colData[c("Group")])
dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/clustered_genes_heatmap_all.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/exp1data/svg_plots/clustered_genes_heatmap_all.short"))
dev.off()

pheatmap(tpmlong[selectedGeneslong,], scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = colDatalong[c("Group")])
dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/clustered_genes_heatmap_long.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/exp1data/svg_plots/clustered_genes_heatmap_long.short"))
dev.off()

pheatmap(tpmshort[selectedGenesshort,], scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = colDatashort[c("Group")])
dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/clustered_genes_heatmap_short.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/exp1data/svg_plots/clustered_genes_heatmap_short.short"))
dev.off()

## PCA plot
{
  # Transformations:
  # transpose the matrix 
  M <- t(tpm[selectedGenes,])
  # transform the counts to log2 scale 
  M <- log2(M + 1)
  # compute PCA 
  pcaResults <- prcomp(M)
  str(M)
  # Plot PCA
  autoplot(pcaResults, data = colData, colour = 'Group')
  dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/PCA_plot.png"))
  dev.off()
  dev.copy(svg, file = file.path(getwd(),"figures/exp1data/svg_plots/PCA_plot.short"))
  dev.off()
  summary(pcaResults)
}


## Correlation plots
correlationMatrix <- cor(tpm)
corrplot(correlationMatrix, order = 'hclust', 
         addrect = 2, addCoef.col = 'white', 
         number.cex = 0.7) 

correlationMatrixlong <- cor(tpmlong)
corrplot(correlationMatrixlong, order = 'hclust', 
         addrect = 2, addCoef.col = 'white', 
         number.cex = 0.7) 

correlationMatrixshort <- cor(tpmshort)
corrplot(correlationMatrixshort, order = 'hclust', 
         addrect = 2, addCoef.col = 'white', 
         number.cex = 0.7) 

# Transforming the correlation plot into a heatmap figure
# split the clusters into two based on the clustering similarity 
pheatmap(correlationMatrix,  
         annotation_col = colData[c("Group")], 
         cutree_cols = 2)
dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/correlation_heatmap_all.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/exp1data/svg_plots/correlation_heatmap_all.short"))
dev.off()

pheatmap(correlationMatrixlong,  
         annotation_col = colDatalong[c("Group")], 
         cutree_cols = 2)
dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/correlation_heatmap_long.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/exp1data/svg_plots/correlation_heatmap_long.short"))
dev.off()

pheatmap(correlationMatrixshort,  
         annotation_col = colDatashort[c("Group")], 
         cutree_cols = 2)
dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/correlation_heatmap_short.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/exp1data/svg_plots/correlation_heatmap_short.short"))
dev.off()

## 5. Differential Expression Analysis ----

#remove the 'width' column
countData <- as.matrix(subset(counts, select = c(-Chr,-Start,-End,-Strand,-Length)))
countDatalong <- as.matrix(subset(countslong, select = c(-Chr,-Start,-End,-Strand,-Length)))
countDatashort <- as.matrix(subset(countsshort, select = c(-Chr,-Start,-End,-Strand,-Length)))

{apply(countData, 2, as.numeric)
sapply(countData, as.numeric)
class(countData) <- "numeric"
storage.mode(countData) <- "numeric"}

{apply(countDatalong, 2, as.numeric)
  sapply(countDatalong, as.numeric)
  class(countDatalong) <- "numeric"
  storage.mode(countDatalong) <- "numeric"}

{apply(countDatashort, 2, as.numeric)
  sapply(countDatashort, as.numeric)
  class(countDatashort) <- "numeric"
  storage.mode(countDatashort) <- "numeric"}

# define the experimental setup 
designFormula <- "~ Group"

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
DEresults = results(dds, contrast = c("Group", 'X', 'A', 'B', 'C'))
#sort results by increasing p-value
DEresults <- DEresults[order(DEresults$pvalue),]
## MA plot
# A scatter plot where the X-axis denotes the average normalized counts across
# samples and the y-axis denotes the log fold change in the given contrast.
# Most points will be on the horizontal 0 line as most genes are not
# differentially expressed
DESeq2::plotMA(object = dds, ylim = c(-5, 5))


## Build initial DEseq matrix
ddslong <- DESeqDataSetFromMatrix(countData = countDatalong, 
                              colData = colDatalong, 
                              design = as.formula(designFormula))
## Remove genes that have almost no information in any give samples
ddslong <- ddslong[ rowSums(DESeq2::counts(ddslong)) > 1, ]
# Now perform Differential expression analysis
ddslong <- DESeq(ddslong)
#compute the contrast for the 'group' variable where 'CTRL' 
#samples are used as the control group.
DEresultslong = results(ddslong, contrast = c("Group", 'C', 'B'))
#sort results by increasing p-value
DEresultslong <- DEresultslong[order(DEresultslong$pvalue),]
## MA plot
# A scatter plot where the X-axis denotes the average normalized counts across
# samples and the y-axis denotes the log fold change in the given contrast.
# Most points will be on the horizontal 0 line as most genes are not
# differentially expressed
DESeq2::plotMA(object = ddslong, ylim = c(-5, 5))
dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/MA_plot_long.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/exp1data/svg_plots/MA_plot_long.short"))
dev.off()

## Build initial DEseq matrix
ddsshort <- DESeqDataSetFromMatrix(countData = countDatashort, 
                              colData = colDatashort, 
                              design = as.formula(designFormula))
## Remove genes that have almost no information in any give samples
ddsshort <- ddsshort[ rowSums(DESeq2::counts(ddsshort)) > 1, ]
# Now perform Differential expression analysis
ddsshort <- DESeq(ddsshort)
#compute the contrast for the 'group' variable where 'CTRL' 
#samples are used as the control group.
DEresultsshort = results(dds, contrast = c("Group", 'X', 'A'))
#sort results by increasing p-value
DEresultsshort <- DEresultsshort[order(DEresultsshort$pvalue),]
## MA plot
# A scatter plot where the X-axis denotes the average normalized counts across
# samples and the y-axis denotes the log fold change in the given contrast.
# Most points will be on the horizontal 0 line as most genes are not
# differentially expressed
DESeq2::plotMA(object = ddsshort, ylim = c(-5, 5))
dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/MA_plot_short.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/exp1data/svg_plots/MA_plot_short.short"))
dev.off()

## P-value distribution

   # shows distribution of the raw p-values
   # we expect to see a peak ~ low p-value and a uniform distribution above 0.1
  #!! OTHERWISE, adjustment for multiple testing does not work, and results
   #!! are NOT meaningful.
ggplot(data = as.data.frame(DEresults), aes(x = pvalue)) + 
  geom_histogram(bins = 100)

ggplot(data = as.data.frame(DEresultslong), aes(x = pvalue)) + 
  geom_histogram(bins = 100)
dev.copy(png, file = file.path(getwd(),"figures/example/PValue_Distribution_long.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/example/PValue_Distribution_long.svg"))
dev.off()

ggplot(data = as.data.frame(DEresultsshort), aes(x = pvalue)) + 
  geom_histogram(bins = 100)
dev.copy(png, file = file.path(getwd(),"figures/example/PValue_Distribution_short.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/example/PValue_Distribution_short.svg"))
dev.off()

## PCA plot

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

rld <- rlog(dds)
DESeq2::plotPCA(rld, ntop = 500, intgroup = 'Group') + 
      ylim(-50, 50) + theme_bw()

## Relative Log Expression (RLE) plot:

  # useful in finding out if the data at hand needs normalization.
  # Runs a quick dignostic to be applied on the raw or normliazed count
  # matrices to see if further processing is required.
  
  #install.packages("EDASeq")
  
    par(mfrow = c(1, 2))
    plotRLE(countData, outline=FALSE, ylim=c(-4, 4), 
            col=as.numeric(colData$Group), 
            main = 'Raw Counts')
    plotRLE(DESeq2::counts(dds, normalized = TRUE), 
            outline=FALSE, ylim=c(-4, 4), 
            col = as.numeric(colData$Group), 
            main = 'Normalized Counts')
    par(mfrow = c(1, 1))

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

## 8. Accounting for additional sources of variation ----

# Sometimes unforeseen variables may contribute to some of the variation in 
# our results (e.g. batch number, temperature of sample storage, etc...)

# First we use DESeq2 to account for possible sources of variation when they are
# known.

## Accounting for covariates using DESeq2
#!! USing a different dataset for this!!
    
    # Look at how samples are clustered by calculating TPM as a heatmap:

    ## Heatmap plot for covariates
    colData <- read.table("data/A_Equina/colData.tsv", header = T, sep = '\t', 
                          stringsAsFactors = TRUE)
    
    
pheatmap(tpm[selectedGenes,], 
             scale = 'row',
             annotation_col = colData[c("Oil","Group","Length","Clone")], 
             show_rownames = FALSE)

pheatmap(tpmlong[selectedGeneslong,], 
         scale = 'row',
         annotation_col = colData[c("Group","Clone")], 
         show_rownames = FALSE)

pheatmap(tpmshort[selectedGenesshort,], 
         scale = 'row',
         annotation_col = colData[c("Oil","Clone")], 
         show_rownames = FALSE)


# define the experimental setup 
## !! The condition of interest should go at the end of the design formula, e.g. ~ subject + condition !! ##
designFormula <- "~ Clone + Group"

## Build initial DEseq matrix
vignette('DESeq2')
ddslong <- DESeqDataSetFromMatrix(countData = countDatalong, 
                                  colData = colDatalong, 
                                  design = as.formula(designFormula))
## Remove genes that have almost no information in any give samples
ddslong <- ddslong[ rowSums(DESeq2::counts(ddslong)) > 1, ]
# Now perform Differential expression analysis
ddslong <- DESeq(ddslong)
#compute the contrast for the 'group' variable where 'CTRL' 
#samples are used as the control group.
DEresultslong = results(ddslong, contrast = c("Group", 'C', 'B'))
#sort results by increasing p-value
DEresultslong <- DEresultslong[order(DEresultslong$pvalue),]
## MA plot
# A scatter plot where the X-axis denotes the average normalized counts across
# samples and the y-axis denotes the log fold change in the given contrast.
# Most points will be on the horizontal 0 line as most genes are not
# differentially expressed
DESeq2::plotMA(object = ddslong, ylim = c(-5, 5))
dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/MA_plot_Clone_long.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/exp1data/svg_plots/MA_plot_Clone_long.short"))
dev.off()

summary(DEresultslong)
sum(DEresultslong$padj < 0.1, na.rm=TRUE)
DEresultslong05 <- results(ddslong, alpha=0.05)
summary(DEresultslong05)
sum(DEresultslong05$padj < 0.05, na.rm=TRUE)
DElonglist <- as.vector(rownames(DEresultslong$padj < 0.05))

## Independent hypothesis weighting
#install.packages("IHW")
#BiocManager::install("IHW", force = TRUE)
library("IHW")
resIHW <- results(ddslong, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult

DESeq2::plotCounts(ddslong, gene=which.min(DEresultslong$padj), intgroup=c("Clone","Group"))

select <- order(rowMeans(counts(ddslong,normalized=TRUE)),
               decreasing=TRUE)[1:20]
df <- as.data.frame(colData(ddslong)[,c("Group","Clone")])
pheatmap(assay(ddslong)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/DElong_Heatmap.png"))
dev.copy(svg, file = file.path(getwd(),"figures/exp1data/svg_plots/DElong_Heatmap.svg"))
dev.off()

sampleDists <- dist(t(assay(ddslong)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(ddslong$Group, ddslong$Clone, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

## Bonferooni Multiple-testing Adjustemnt
# Create new pbonferonni variable - exactly the same as P-value
DEresultslong$pbonferroni <- DEresultslong$pvalue
# Adjust the values using the bonferonni method:
DEresultslong$pbonferroni<-p.adjust(DEresultslong$pbonferroni, method = "bonferroni", n = length(DEresultslong$pbonferroni))
# order the new pbonferonni values
DEresultslong <- DEresultslong[order(DEresultslong$pbonferroni),]
# Display P-vlue distribution:
hist(DEresultslong$pbonferroni)

## Banjamini and Hochberg Multiple-testing Adjustemnt
# Create new pbonferonni variable - exactly the same as P-value
DEresultslong$pfdr <- DEresultslong$pvalue
# Adjust the values using the bonferonni method:
DEresultslong$pfdr<-p.adjust(DEresultslong$pbonferroni, method = "fdr", n = length(DEresultslong$pbonferroni))
# order the new pbonferonni values
DEresultslong <- DEresultslong[order(DEresultslong$pfdr),]
# Display P-vlue distribution:
hist(DEresultslong$pfdr)

## Holm Multiple-testing Adjustemnt
# Create new pbonferonni variable - exactly the same as P-value
DEresultslong$pholm <- DEresultslong$pvalue
# Adjust the values using the bonferonni method:
DEresultslong$pholm<-p.adjust(DEresultslong$pbonferroni, method = "holm", n = length(DEresultslong$pbonferroni))
# order the new pbonferonni values
DEresultslong <- DEresultslong[order(DEresultslong$pholm),]
# Display P-vlue distribution:
hist(DEresultslong$pholm)

par(mfrow = c(1,2))
hist(DEresultslong$padj)
hist(DEresultslong$pvalue)
par(mfrow = c(1,1))

## Build initial DEseq matrix
ddsshort <- DESeqDataSetFromMatrix(countData = countDatashort, 
                                   colData = colDatashort, 
                                   design = as.formula(designFormula))
## Remove genes that have almost no information in any give samples
ddsshort <- ddsshort[ rowSums(DESeq2::counts(ddsshort)) > 1, ]
# Now perform Differential expression analysis
ddsshort <- DESeq(ddsshort)
#compute the contrast for the 'group' variable where 'CTRL' 
#samples are used as the control group.
DEresultsshort = results(ddsshort, contrast = c("Group", 'X', 'A'))
#sort results by increasing p-value
DEresultsshort <- DEresultsshort[order(DEresultsshort$pvalue),]
## MA plot
# A scatter plot where the X-axis denotes the average normalized counts across
# samples and the y-axis denotes the log fold change in the given contrast.
# Most points will be on the horizontal 0 line as most genes are not
# differentially expressed
DESeq2::plotMA(object = ddsshort, ylim = c(-5, 5))
dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/MA_plot_Clone_short.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/exp1data/svg_plots/MA_plot_Clone_short.short"))
dev.off()

summary(DEresultsshort)
sum(DEresultsshort$padj < 0.1, na.rm=TRUE)
DEresultsshort05 <- results(ddsshort, alpha=0.05)
summary(DEresultsshort05)
sum(DEresultsshort05$padj < 0.05, na.rm=TRUE)

## Independent hypothesis weighting
# BiocManager::install("IHW", force = TRUE)
library("IHW")
resIHW <- results(ddsshort, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult

DESeq2::plotCounts(ddsshort, gene=which.min(DEresultsshort$padj), intgroup=c("Clone","Group"))

select <- order(rowMeans(counts(ddsshort,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(ddsshort)[,c("Group","Clone")])
pheatmap(assay(ddsshort)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/DEshort_Heatmap.png"))
dev.copy(svg, file = file.path(getwd(),"figures/exp1data/svg_plots/DEshort_Heatmap.svg"))
dev.off()

sampleDists <- dist(t(assay(ddsshort)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(ddsshort$Group, ddslong$Clone, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

## P-value distribution

# shows distribution of the raw p-values
# we expect to see a peak ~ low p-value and a uniform distribution above 0.1
#!! OTHERWISE, adjustment for multiple testing does not work, and results
#!! are NOT meaningful.
ggplot(data = as.data.frame(DEresultslong), aes(x = pvalue)) + 
  geom_histogram(bins = 100)
dev.copy(png, file = file.path(getwd(),"figures/example/PValue_Distribution_long.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/example/PValue_Distribution_long.svg"))
dev.off()

ggplot(data = as.data.frame(DEresultsshort), aes(x = pvalue)) + 
  geom_histogram(bins = 100)
dev.copy(png, file = file.path(getwd(),"figures/example/PValue_Distribution_short.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/example/PValue_Distribution_short.svg"))
dev.off()

## Bonferooni Multiple-testing Adjustemnt
# Create new pbonferonni variable - exactly the same as P-value
DEresultsshort$pbonferroni <- DEresultsshort$pvalue
# Adjust the values using the bonferonni method:
DEresultsshort$pbonferroni<-p.adjust(DEresultsshort$pbonferroni, method = "bonferroni", n = length(DEresultsshort$pbonferroni))
# order the new pbonferonni values
DEresultsshort <- DEresultsshort[order(DEresultsshort$pbonferroni),]
# Display P-vlue distribution:
hist(DEresultsshort$pbonferroni)

## Banjamini and Hochberg Multiple-testing Adjustemnt
# Create new pbonferonni variable - exactly the same as P-value
DEresultsshort$pfdr <- DEresultsshort$pvalue
# Adjust the values using the bonferonni method:
DEresultsshort$pfdr<-p.adjust(DEresultsshort$pbonferroni, method = "fdr", n = length(DEresultsshort$pbonferroni))
# order the new pbonferonni values
DEresultsshort <- DEresultsshort[order(DEresultsshort$pfdr),]
# Display P-vlue distribution:
hist(DEresultsshort$pfdr)

## Holm Multiple-testing Adjustemnt
# Create new pbonferonni variable - exactly the same as P-value
DEresultsshort$pholm <- DEresultsshort$pvalue
# Adjust the values using the bonferonni method:
DEresultsshort$pholm<-p.adjust(DEresultsshort$pbonferroni, method = "holm", n = length(DEresultsshort$pbonferroni))
# order the new pbonferonni values
DEresultsshort <- DEresultsshort[order(DEresultsshort$pholm),]
# Display P-vlue distribution:
hist(DEresultsshort$pholm)

par(mfrow = c(1,2))
hist(DEresultsshort$padj)
hist(DEresultsshort$pvalue)
par(mfrow = c(1,1))


rldlong <- rlog(ddslong)
DESeq2::plotPCA(rldlong, ntop = 500, intgroup = "Group") + 
  ylim(-50, 50) + theme_bw()
DESeq2::plotPCA(rldlong, ntop = 500, intgroup = 'Clone') + 
  ylim(-50, 50) + theme_bw()


rldshort <- rlog(ddsshort)
DESeq2::plotPCA(rldshort, ntop = 500, intgroup = 'Group') + 
  ylim(-50, 50) + theme_bw()
DESeq2::plotPCA(rldshort, ntop = 500, intgroup = 'Clone') + 
  ylim(-50, 50) + theme_bw()


## Relative Log Expression (RLE) plot:

# useful in finding out if the data at hand needs normalization.
# Runs a quick dignostic to be applied on the raw or normliazed count
# matrices to see if further processing is required.

#install.packages("EDASeq")

par(mfrow = c(1, 2))
plotRLE(countDatalong, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colData$Group), 
        main = 'Raw Counts')
plotRLE(DESeq2::counts(ddslong, normalized = TRUE), 
        outline=FALSE, ylim=c(-4, 4), 
        col = as.numeric(colData$Group), 
        main = 'Normalized Counts')
par(mfrow = c(1, 1))
dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/RLE_plot_long.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/exp1data/svg_plots/RLE_plot_long.svg"))
dev.off()

par(mfrow = c(1, 2))
plotRLE(countDatashort, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colData$Group), 
        main = 'Raw Counts')
plotRLE(DESeq2::counts(ddsshort, normalized = TRUE), 
        outline=FALSE, ylim=c(-4, 4), 
        col = as.numeric(colData$Group), 
        main = 'Normalized Counts')
par(mfrow = c(1, 1))
dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/RLE_plot_short.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/exp1data/svg_plots/RLE_plot_short.svg"))
dev.off()
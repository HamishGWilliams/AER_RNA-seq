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

# Other packages which need to be forced
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("gage", force = TRUE)
library(gage)

## Load in Count Data: ----
counts <- as.matrix(read.table("data/A_Equina/Counts_Data.tsv", header = T, sep = '\t'))
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

pheatmap(tpmlong[selectedGeneslong,], scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = colDatalong[c("Group")])

pheatmap(tpmshort[selectedGenesshort,], scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = colDatashort[c("Group")])

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
  autoplot(pcaResults, data = colData, colour = 'group')
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

pheatmap(correlationMatrixlong,  
         annotation_col = colData[c("Group")], 
         cutree_cols = 2)

pheatmap(correlationMatrixshort,  
         annotation_col = colData[c("Group")], 
         cutree_cols = 2)

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

#define the experimental setup 
designFormula <- "~ Group"

## Build initial DEseq matrix
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = as.formula(designFormula))
ddslong <- DESeqDataSetFromMatrix(countData = countDatalong, 
                              colData = colDatalong, 
                              design = as.formula(designFormula))
ddsshort <- DESeqDataSetFromMatrix(countData = countDatashort, 
                              colData = colDatashort, 
                              design = as.formula(designFormula))

## Remove genes that have almost no information in any give samples
dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]
ddslong <- ddslong[ rowSums(DESeq2::counts(ddslong)) > 1, ]
ddsshort <- ddsshort[ rowSums(DESeq2::counts(ddsshort)) > 1, ]

# Now perform Differential expression analysis
dds <- DESeq(dds)
ddslong <- DESeq(ddslong)
ddsshort <- DESeq(ddsshort)

#compute the contrast for the 'group' variable where 'CTRL' 
#samples are used as the control group. 
DEresults = results(dds, contrast = c("group", 'X', 'A', 'B', 'C'))
DEresultslong = results(dds, contrast = c("Group", 'B', 'C'))
DEresultsshort = results(dds, contrast = c("Group", 'X', 'A'))
#sort results by increasing p-value
DEresults <- DEresults[order(DEresults$pvalue),]
DEresultslong <- DEresultslong[order(DEresultslong$pvalue),]
DEresultsshort <- DEresultsshort[order(DEresultsshort$pvalue),]

## MA plot

  # A scatter plot where the X-axis denotes the average normalized counts across
  # samples and the y-axis denotes the log fold change in the given contrast.
  
  # Most points will be on the horizontal 0 line as most genes are not
  # differentially expressed
  DESeq2::plotMA(object = dds, ylim = c(-5, 5))
  DESeq2::plotMA(object = ddslong, ylim = c(-5, 5))
  DESeq2::plotMA(object = ddsshort, ylim = c(-5, 5))
  
## P-value distribution

   # shows distribution of the raw p-values
   # we expect to see a peak ~ low p-value and a uniform distribution above 0.1
  #!! OTHERWISE, adjustment for multiple testing does not work, and results
   #!! are NOT meaningful.
ggplot(data = as.data.frame(DEresults), aes(x = pvalue)) + 
  geom_histogram(bins = 100)

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
DESeq2::plotPCA(rld, ntop = 500, intgroup = 'group') + 
      ylim(-50, 50) + theme_bw()

## Relative Log Expression (RLE) plot:

  # useful in finding out if the data at hand needs normalization.
  # Runs a quick dignostic to be applied on the raw or normliazed count
  # matrices to see if further processing is required.
  
  #install.packages("EDASeq")
  
    par(mfrow = c(1, 2))
    plotRLE(countData, outline=FALSE, ylim=c(-4, 4), 
            col=as.numeric(colData$group), 
            main = 'Raw Counts')
    plotRLE(DESeq2::counts(dds, normalized = TRUE), 
            outline=FALSE, ylim=c(-4, 4), 
            col = as.numeric(colData$group), 
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
    colData <- read.table("data/A_Equina/Sample_ID_query.tsv", header = T, sep = '\t', 
                          stringsAsFactors = TRUE)
    
    
    pheatmap(tpm[selectedGenes,], 
             scale = 'row',
             annotation_col = colData[c("group","Oil")], 
             show_rownames = FALSE)
    
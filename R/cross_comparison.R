# Load in Data
counts <- as.matrix(read.table("data/A_Equina/A_Equina_Counts_redo.tsv", header = T, sep = '\t'))

# Subset for A vs B
countsA_B <- subset(counts[,c(1:5,#metadata
                              6:12,#Group A
                              13:19)#GROUP B
                          ])

# Subset for A vs C
countsA_C<- subset(counts[,c(1:5,#metadata
                             6:12,#Group A
                             20:26)#Group C
                          ])

# Calculate parameters
cpmA_B <- {apply(subset(countsA_B, select = c(-Chr,-Start,-End,-Strand,-Length)), 2, 
                  function(x) as.numeric(x)/sum(as.numeric(x)) * 10^6)}
cpmA_C <- {apply(subset(countsA_C, select = c(-Chr,-Start,-End,-Strand,-Length)), 2, 
                   function(x) as.numeric(x)/sum(as.numeric(x)) * 10^6)}

counts_rownames <- as.vector(row.names(counts))
rownames(cpmA_B) <- counts_rownames
rownames(cpmA_C) <- counts_rownames

# Genelengths
geneLengths <- as.vector(subset(counts, select = c(Length)))
geneLengths <- as.vector(as.numeric(geneLengths))

# RPKM
rpkmA_B <- {apply(X = subset(countsA_B, select = c(-Chr,-Start,-End,-Strand,-Length)),
                   MARGIN = 2, 
                   FUN = function(x) {
                     10^9 * as.numeric(x) / geneLengths / sum(as.numeric(x))
                   })
}
rownames(rpkmA_B) <- counts_rownames

rpkmA_C <- {apply(X = subset(countsA_C, select = c(-Chr,-Start,-End,-Strand,-Length)),
                    MARGIN = 2, 
                    FUN = function(x) {
                      10^9 * as.numeric(x) / geneLengths / sum(as.numeric(x))
                    })
}
rownames(rpkmA_C) <- counts_rownames

# RPK
rpkA_B <- {apply(subset(countsA_B, select = c(-Chr,-Start,-End,-Strand,-Length)), 2, 
                  function(x) as.numeric(x)/(geneLengths/1000))}
rownames(rpkA_B) <- counts_rownames

rpkA_C <- {apply(subset(countsA_C, select = c(-Chr,-Start,-End,-Strand,-Length)), 2, 
                   function(x) as.numeric(x)/(geneLengths/1000))}
rownames(rpkA_C) <- counts_rownames

# TPM
tpmA_B <- apply(rpkA_B, 2, function(x) as.numeric(x) / sum(as.numeric(x)) * 10^6)
rownames(tpmA_B) <- counts_rownames

tpmA_C <- apply(rpkA_C, 2, function(x) as.numeric(x) / sum(as.numeric(x)) * 10^6)
rownames(tpmA_C) <- counts_rownames

# Variance
VA_B <- apply(tpmA_B, 1, var)
VA_C <- apply(tpmA_C, 1, var)


# Select top 100 genes: !! CHANGE THIS TO DESIRED NUMBER !!
selectedGenesA_B <- names(VA_B[order(VA_B, decreasing = T)][1:100])
selectedGenesA_C <- names(VA_C[order(VA_C, decreasing = T)][1:100])


## Check each step for errors:
colSums(cpmA_B) # = 10^6 for ALL
colSums(rpkmA_B) # Different values
colSums(tpmA_B) # = 10^6 for ALL

colSums(cpmA_C) # = 10^6 for ALL
colSums(rpkmA_C) # Different values
colSums(tpmA_C) # = 10^6 for ALL


## Figures & Plots 
# Heatmap for clustered genes
pheatmap(tpmA_B[selectedGenesA_B,], 
         scale = 'row', 
         show_rownames = FALSE,
         annotation_col = colData[c("Group", "Clone")])

pheatmap(tpmA_C[selectedGenesA_C,], 
         scale = 'row', 
         show_rownames = FALSE,
         annotation_col = colData[c("Group", "Clone")])

## PCA plots
# transpose the matrix 
MA_B <- t(tpmA_B[selectedGenesA_B,])
# transform the counts to log2 scale 
MA_B <- log2(MA_B + 1)
# compute PCA 
pcaResultsA_B <- prcomp(MA_B)
str(MA_B)
# Plot PCA [all]
autoplot(pcaResultsA_B, 
         data = colData, 
         colour = "Group",
         main = "PCA plot of A vs B")
summary(pcaResultsA_B)

# transpose the matrix 
MA_C <- t(tpmA_C[selectedGenesA_C,])
# transform the counts to log2 scale 
MA_C <- log2(MA_C + 1)
# compute PCA 
pcaResultsA_C <- prcomp(MA_C)
str(MA_C)
# Plot PCA [all]
autoplot(pcaResultsA_C, 
         data = colData, 
         colour = 'Group',
         main = "PCA plot of all groups")
summary(pcaResultsA_C)

countA_B <- as.matrix(subset(countsA_B, select = c(-Chr,-Start,-End,-Strand,-Length)))
countA_C <- as.matrix(subset(countsA_C, select = c(-Chr,-Start,-End,-Strand,-Length)))

{apply(countA_B, 2, as.numeric)
  sapply(countA_B, as.numeric)
  class(countA_B) <- "numeric"
  storage.mode(countA_B) <- "numeric"}

{apply(countA_C, 2, as.numeric)
  sapply(countA_C, as.numeric)
  class(countA_C) <- "numeric"
  storage.mode(countA_C) <- "numeric"}

colData <- read.table("data/A_Equina/colData.tsv", header = T, sep = '\t', 
                      stringsAsFactors = TRUE)
colDataA_B <- subset(colData[1:14,])
colDataA_C <- subset(colData[c(1:7,
                                 15:21),])

## Build initial DEseq matrix (A vs B)
designFormula <- "~ Group"

ddsA_B <- DESeqDataSetFromMatrix(countData = countA_B, 
                              colData = colDataA_B, 
                              design = as.formula(designFormula))
## Remove genes that have almost no information in any give samples
ddsA_B <- ddsA_B[ rowSums(DESeq2::counts(ddsA_B)) > 1, ]
# Now perform Differential expression analysis
ddsA_B <- DESeq(ddsA_B)
#compute the contrast for the 'group' variable where 'CTRL' 
#samples are used as the control group. 
DEresultsA_B = results(ddsA_B, contrast = c("Group", 'B', 'A'))
#sort results by increasing p-value
DEresultsA_B <- DEresultsA_B[order(DEresultsA_B$pvalue),]

DESeq2::plotMA(object = ddsA_B, 
               ylim = c(-5, 5),
               main = "MA plot [A vs B]")

# P-value Distribution
par(mfrow = c(1,2))
hist(DEresultsA_B$padj,
     main = "DEseq default adjustment [A vs B]",
     breaks = 150)
hist(DEresultsA_B$pvalue,
     main = "Raw P Values [A vs B]",
     breaks = 150)
par(mfrow = c(1,1))

## PCA plot
# extract normalized counts from the DESeqDataSet object
countsNormalizedA_B <- DESeq2::counts(ddsA_B, normalized = TRUE)

# select top 500 most variable genes
selectedGenes500A_B<- names(sort(apply(countsNormalizedA_B, 1, var), 
                                  decreasing = TRUE)[1:500])

DESeq2::plotPCA(countsNormalizedA_B[selectedGenes500A_B,], 
                col = as.numeric(colDataA_B$Group), 
                adj = 0.5, 
                xlim = c(-0.6, 0.6), 
                ylim = c(-0.6, 0.6),
                main = "PCA of Normalized counts [A vs B]")

## Relative Log Expression (RLE) plot:
par(mfrow = c(1, 2))
plotRLE(countA_B, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colDataA_B$Group), 
        main = 'Raw Counts [A vs B]')
plotRLE(DESeq2::counts(ddsA_B, normalized = TRUE), 
        outline=FALSE, ylim=c(-4, 4), 
        col = as.numeric(colDataA_B$Group), 
        main = 'Normalized Counts [A vs B]')
par(mfrow = c(1, 1))

# DEGs summary
summary(DEresultsA_B)
sum(DEresultsA_B$padj < 0.1, na.rm=TRUE)
DEresultsA_B05 <- results(ddsA_B, alpha=0.05)
summary(DEresultsA_B05)
sum(DEresultsA_B05$padj < 0.05, na.rm=TRUE)


## Build initial DEseq matrix (A vs C)
designFormula <- "~ Group"

ddsA_C <- DESeqDataSetFromMatrix(countData = countA_C, 
                                 colData = colDataA_C, 
                                 design = as.formula(designFormula))
## Remove genes that have almost no information in any give samples
ddsA_C <- ddsA_C[ rowSums(DESeq2::counts(ddsA_C)) > 1, ]
# Now perform Differential expression analysis
ddsA_C <- DESeq(ddsA_C)
#compute the contrast for the 'group' variable where 'CTRL' 
#samples are used as the control group. 
DEresultsA_C = results(ddsA_C, contrast = c("Group", 'C', 'A'))
#sort results by increasing p-value
DEresultsA_C <- DEresultsA_C[order(DEresultsA_C$pvalue),]

DESeq2::plotMA(object = ddsA_C, 
               ylim = c(-5, 5),
               main = "MA plot [A vs C]")

# P-value Distribution
par(mfrow = c(1,2))
hist(DEresultsA_C$padj,
     main = "DEseq default adjustment [A vs C]",
     breaks = 150)
hist(DEresultsA_C$pvalue,
     main = "Raw P Values [A vs C]",
     breaks = 150)
par(mfrow = c(1,1))

## PCA plot
# extract normalized counts from the DESeqDataSet object
countsNormalizedA_C <- DESeq2::counts(ddsA_C, normalized = TRUE)

# select top 500 most variable genes
selectedGenes500A_C<- names(sort(apply(countsNormalizedA_C, 1, var), 
                                 decreasing = TRUE)[1:500])

DESeq2::plotPCA(countsNormalizedA_C[selectedGenes500A_C,], 
                col = as.numeric(colDataA_C$Group), 
                adj = 0.5, 
                xlim = c(-0.6, 0.6), 
                ylim = c(-0.6, 0.6),
                main = "PCA of Normalized counts [A vs C]")

## Relative Log Expression (RLE) plot:
par(mfrow = c(1, 2))
plotRLE(countA_C, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colDataA_C$Group), 
        main = 'Raw Counts [A vs C]')
plotRLE(DESeq2::counts(ddsA_C, normalized = TRUE), 
        outline=FALSE, ylim=c(-4, 4), 
        col = as.numeric(colDataA_C$Group), 
        main = 'Normalized Counts [A vs C]')
par(mfrow = c(1, 1))

# DEGs summary
summary(DEresultsA_C)
sum(DEresultsA_C$padj < 0.1, na.rm=TRUE)
DEresultsA_C05 <- results(ddsA_C, alpha=0.05)
summary(DEresultsA_C05)
sum(DEresultsA_C05$padj < 0.05, na.rm=TRUE)

## Venn Diagram of A-X + A-B + A-C:
VD_Short <- na.omit(DEresultsshort[DEresultsshort$padj,])
VD_Short <- VD_Short[VD_Short$padj < 0.05,]

VD_A_B <- na.omit(DEresultsA_B)
VD_A_B <- VD_A_B[VD_A_B$padj < 0.05,]

VD_A_C <- na.omit(DEresultsA_C)
VD_A_C <- VD_A_C[VD_A_C$padj < 0.05,]

VD_Long <- na.omit(DEresultslong)
VD_Long <- VD_Long[VD_Long$padj < 0.05,]

VDList <- list('A-X'=rownames(VD_Short),
             'A-B'=rownames(VD_A_B),
             'A-C'=rownames(VD_A_C),
             'C-B'=rownames(VD_Long))

VDList <- list('A-X'=rownames(VD_Short),
               'C-B'=rownames(VD_Long))

VDList <- list('A-B'=rownames(VD_A_B),
               'A-C'=rownames(VD_A_C))

VDList <- list('A-X'=rownames(VD_Short),
               'A-B'=rownames(VD_A_B),
               'A-C'=rownames(VD_A_C))

ggVennDiagram(VDList) +
  ggplot2::scale_fill_gradient(low = "blue", high = "red")

ggVennDiagram(VDList,
              label_alpha = 0,
              label_color = ("black")) +
  ggplot2::scale_fill_gradient(low = "yellow", high = "red")

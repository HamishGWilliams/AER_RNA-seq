counts <- as.matrix(read.table("data/A_Equina/A_Equina_Counts_redo.tsv", header = T, sep = '\t'))
Counts_only <- subset(counts, select = c(-Chr,-Start,-End,-Strand,-Length))
countslong<- subset(counts[,c(1:5,13:26)])
countsshort <- subset(counts[,c(1:5,6:12,27:33)])
colData <- read.table("data/A_Equina/colData.tsv", header = T, sep = '\t', 
                      stringsAsFactors = TRUE)
colDatalong <- subset(colData[8:21,])
colDatashort <- subset(colData[c(1:7,22:28),])


designFormula <- "~ Group + Clone"
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
DEresultslong = results(ddslong, contrast = c("Group", 'B', 'C'))
#sort results by increasing p-value
DEresultslong <- DEresultslong[order(DEresultslong$pvalue),]

BiocManager::install("glmGamPoi", force = TRUE)
library("glmGamPoi")

ddslongtest <- DESeq(
  ddslong,
  test = "LRT",
  fitType = c("glmGamPoi"), # Fitting to a Gamma Poisson Distribution
  quiet = FALSE,
  full = design(object),
  reduced = TRUE,
  minReplicatesForReplace = 7,
  useT = FALSE,
  minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5,
  parallel = FALSE)

?DESeq

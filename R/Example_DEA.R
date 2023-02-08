## Load required packages:

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("alyssafrazee/RSkittleBrewer","ballgown",
                     "genefilter","dplyr","devtools"))
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

# Load phenotype data for samples
getwd()
pheno_data = read.csv("data/example/geuvadis_phenodata.csv")

# Read expression data calculated by StringTie
bg_chrX = ballgown(dataDir = "data/example/ballgown", samplePattern = "ERR", pData=pheno_data)

# Remove low abundance genes 
bg_chrX_filt = subset(bg_chrX,"rowVars(texpr(bg_chrX)) >1",genomesubset=TRUE)

# Identify transcripts that show significant DE
results_transcripts = stattest(bg_chrX_filt,
                               feature="transcript",covariate="sex",adjustvars =
                                 c("population"), getFC=TRUE, meas="FPKM")

# Identify Genes which are significant between groups
results_genes = stattest(bg_chrX_filt, feature="gene",
                         covariate="sex", adjustvars = c("population"), getFC=TRUE,
                         meas="FPKM")

# Add gene names and gene IDs to results_transcripts dataframe:
results_transcripts =
  data.frame(geneNames=ballgown::geneNames(bg_chrX_filt),
             geneIDs=ballgown::geneIDs(bg_chrX_filt), results_transcripts)

# Sort results by smallest P value to largest:
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)

# Write the results to a csv file that can be shared:
# Transcripts:
write.csv(results_transcripts, "data/example/chrX_transcript_results.csv",
          row.names=FALSE)
# Genes:
write.csv(results_genes, "data/example/chrX_gene_results.csv",
          row.names=FALSE) 

# Identify transcripts and genes with a q value <0.05:
subset(results_transcripts,results_transcripts$qval<0.05)
subset(results_genes,results_genes$qval<0.05)

tropical= c('darkorange', 'dodgerblue',
            'hotpink', 'limegreen', 'yellow')
palette(tropical)

# Show distrubution of gene abundances
fpkm = texpr(bg_chrX_filt, meas="FPKM")
fpkm = log2(fpkm+1)
boxplot(fpkm, col=as.numeric(pheno_data$sex), las=2, ylab='log2(FPKM+1)')

?texpr
str(pheno_data)

# Make plots of individual transcripts across samples:
ballgown::transcriptNames(bg_chrX)[12]
## 12
## "NM_012227"
ballgown::geneNames(bg_chrX)[12]
## 12
## "GTPBP6"
plot(fpkm[12,] ~ pheno_data$sex, border=c(1,2),
      main=paste(ballgown::geneNames(bg_chrX)[12],' : ',
                 ballgown::transcriptNames(bg_chrX)[12]),pch=19, xlab="Sex",
      ylab='log2(FPKM+1)')
points(fpkm[12,] ~ jitter(as.numeric(pheno_data$sex)),
        col=as.numeric(pheno_data$sex))

plotTranscripts(ballgown::geneIDs(bg_chrX)[1729], bg_chrX, main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))

plotMeans('MSTRG.56', bg_chrX_filt,groupvar="sex",legend=FALSE)




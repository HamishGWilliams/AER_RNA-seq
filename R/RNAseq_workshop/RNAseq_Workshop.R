## Metadata ----
# Date Created: 01/03/3
# Title: "RNAseq Workshop 'edgeR' coding"

## 0. Packages ----
# Build package checking function:
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

# BiocManager Installed packages:
BiocManager::install("edgeR")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Sc.sgd.db")

# Generate package list:
pkg_list<-c("BiocManager",
            "edgeR",
            "gplots",
            "RColorBrewer",
            "clusterProfiler",
            "org.Sc.sgd.db")

# Run Package Check function:
check_packages(pkg_list)
check_packages(pkg_list)


## 1. Reading in Data: ----
# Count Data:
countdata <- read.table("R/RNAseq_workshop/Europe_Malaysia_Counts.txt",
                        quote="",
                        sep="\t",
                        header=TRUE,
                        row.names=1)

# Column data:
sampleData <- read.table("R/RNAseq_workshop/samples.csv",
                         sep=",",
                         header = TRUE,
                         stringsAsFactors = TRUE)

# Combine 'counts' and 'columns' into DGEList Object:
data <- DGEList(counts = countdata,
                group = sampleData$Strain)
# Filter out low count genes:
keep <- rowSums(cpm(data)>1) >=2
data <- data[keep, ,keep.lib.sizes=FALSE]


## 2. Differential Expression: ----
# Normalize Data
data <- calcNormFactors(data)
# Estimate dispersion (uncertainty):
data <- estimateDisp(data)
# Differential Expression Execution:
et <- exactTest(data)
  # NOTE: For more complex study designs (e.g. paired), refer
  # to this command in the manual for 'edgeR': 
  # https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# Select genes (in this case, all -> specified by "Inf"):
sig <- topTags(et, n=Inf)
head(sig$table)

# Write table of DE results:
write.table(sig, "R/RNAseq_workshop/DEGs.txt",
            sep="\t")


## 3. Visualizing Raw Data and DE results: ----
# Plot MDS (i.e. PCA plot):
colours=c("red","dark green")
palette(colours)
plotMDS(data,
        col=as.numeric(sampleData$Strain),
        pch=as.numeric(sampleData$Strain),
        lwd = 3)
legend("topright",
       legend = c("Sample",
                  levels(sampleData$Strain)))

# Extract normalized counts:
norm <- cpm(data,
            normalized.lib.sizes = TRUE,
            log = TRUE)

# Plot Boxplot
colours=c("light blue","pink")
palette(colours)
boxplot(norm, ylab="Log Normalised Counts", las=2, 
        col=as.numeric(sampleData$Strain))

# Generate colour scale:
hmcol <- colorRampPalette(brewer.pal(9,"Blues"))(100)
sigcounts <- norm[rownames(sig$table)[sig$table$FDR<0.05 &
                                        abs(sig$table$logFC)>1.5],]

# Plot Heatmap:
heatmap.2(sigcounts,
          col=hmcol,
          density.info = "none",
          trace = "none",
          margins=c(7,4))

# Plot Volcano plot:
plot(sig$table$logFC, -log10(sig$table$FDR),
     ylab = "-log10 FDR",
     xlab = "LogFC",
     col=ifelse(sig$table$FDR < 0.05 &
                  sig$table$logFC>1.5 |
                  sig$table$FDR < 0.05 &
                  sig$table$logFC<(-1.5),
                "red","black"))

## 4. Functional Analysis: ----
# Extract results:
res <- as.data.frame(sig$table)
# Extract genes:
genes <- rownames(subset(res, FDR<=0.05))

# Over-representation of GO terms:
go <- enrichGO(genes, 
               org.Sc.sgd.db, 
               keyType = "GENENAME",
               ont = "ALL")
  # NOTE: the "org.Sc.sgd.db" is a YEAST specific genome wide
  # annotation package. Therefore it cannot be used with other
  # species of phylum of life. Refer to the "edgeR" manual
  # linked above for examples of how to perform GO.

# Write table of of significant GO terms:
write.table(go, "GO_enrich.txt",
            sep="\t",
            row.names=FALSE)

## Over-representation of the KEGG Pathways:
# Dictionary matching gene name to ENTREZ ID
dict <- select(org.Sc.sgd.db,
               keys=rownames(res),
               columns = c("ENTREZID",
                           "GENENAME"),
               keytype = "GENENAME")

# Make a list of significant ENTREZIDs for analysis:
genes_kegg <- dict$ENTREZID[match(genes,dict$GENENAME)]
kegg <- enrichKEGG(genes_kegg,
                   "sce",
                   keyType = "ncbi-geneid")

# Run Kegg Analysis:
kegg@result$geneID <- sapply(kegg@result$geneID,
                             function(x)
                               sapply (strsplit(x,"/"),
                                       function(y)
                                         paste(dict$GENENAME[match(y,
                                                                   dict$ENTREZID)],
                                               collapse = "/")))

# Write a table of significant KEGG pathways (q < 0.05):
write.table(kegg,
            "Pathway_enrich.txt",
            sep = "\t",
            row.names = FALSE)
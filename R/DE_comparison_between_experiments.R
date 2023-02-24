select <- order(rowMeans(counts(ddsshort,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(ddsshort)[,c("Group","Clone")])
pheatmap(assay(ddsshort)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/DEshort_Heatmap.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/exp1data/svg_plots/DEshort_Heatmap.svg"))
dev.off()

select <- order(rowMeans(counts(ddslong,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(ddslong)[,c("Group","Clone")])
pheatmap(assay(ddslong)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/DElong_Heatmap.png"))
dev.off()
dev.copy(svg, file = file.path(getwd(),"figures/exp1data/svg_plots/DElong_Heatmap.svg"))
dev.off()

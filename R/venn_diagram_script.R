if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")
install.packages("VennDiagram")
library("VennDiagram")

vendDEresultslong <- na.omit(DEresultslong)
vendDEresultslong <- vendDEresultslong[vendDEresultslong$padj < 0.1,]

vendDEresultsshort <- na.omit(DEresultsshort)
vendDEresultsshort <- vendDEresultsshort[vendDEresultsshort$padj < 0.1,]
?ggVennDiagram

ggVennDiagram(c(vendDEresultslong, vendDEresultsshort))
List <- list('Long'=rownames(vendDEresultslong),
             'Short'=rownames(vendDEresultsshort))
ggVennDiagram(List,
              label_alpha = 0,
              label_color = ("black")) +
  ggplot2::scale_fill_gradient(low = "white", high = "red")

dev.copy(png, file = file.path(getwd(),"figures/exp1data/png_plots/venndiagram_short_long.png"))
dev.off()

List <- list('Long'=rownames(DEresultslong),
             'Short'=rownames(DEresultsshort))
ggVennDiagram(List,
              label_alpha = 0,
              label_color = ("black")) +
  ggplot2::scale_fill_gradient(low = "white", high = "red")

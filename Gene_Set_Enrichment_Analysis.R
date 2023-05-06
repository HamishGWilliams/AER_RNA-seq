## Purpose of Script ----

## This script will contain all of the steps taken to perform Gene Set
## Enrichment Analysis on the Actinia eqiuna Differentially Expressed dataset

## From this analysis I am looking to with greater confidence assert terms
## of interest for further investigation, and as a way to collectively
## encapsulate the overall functional (phenotypic) changes that are occurring
## within Actinia equina in response to their exposure to Diesel in the 
## experiments.

## Downloading packages ----
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler",
                     force = TRUE)  # Have to use force to install

library(clusterProfiler)

browseVignettes("clusterProfiler")  # Not that useful, but provides links to 
                                    # 'some' resources of "use"


#!/usr/bin/env Rscript
#
# SampleAnalysis.R
#
# @author Rahul Dhodapkar <rahul.dhodapkar@yale.edu>
#

library(Seurat)
# import other required libraries here

################################################################################
## BUILD OUTPUT SCAFFOLDING
################################################################################
#
# Headers like these are not necessary, but I find they help to organize
# one's thoughts into coherent blocks, and help with reproducibility.
#
# This (optional) block creates some folders for any output files (calc/) 
# or figures (fig/) that you may generate during the course of the exploratory
# analysis.
#
# The results of these calls are assigned to a junk variable 's' simply to
# silence output for cleaner execusion.
#

s <- ifelse(!dir.exists("./fig"), dir.create("./fig"), FALSE)
s <- ifelse(!dir.exists("./calc"), dir.create("./calc"), FALSE)

################################################################################
## BEGIN ANALYSIS
################################################################################

retina.combined <- readRDS('./data/retina_combined.rds')

################################################################################
## BEGIN ANALYSIS
################################################################################
#
# Headers like these are not necessary, but I find they help to organize
# one's thoughts. Write a bit here about what you are trying/hoping to do.
#
# Find some samples of basic exploratory analysis below.
#

# Plot a UMAP labeled by cell type
DimPlot(retina.combined, label=T)

# Plot the expression of a specific gene on UMAP
FeaturePlot(retina.combined, 'rna_RLBP1')

# Plot the expression of a specific gene as a violin plot by cell type
VlnPlot(retina.combined, 'rna_RLBP1')

# Plot the expression of a specific gene split by disease and cell type
VlnPlot(retina.combined, 'rna_RLBP1', split.by='DiseaseGroup')


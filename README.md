# Wet AMD Prototyping Snapshot

This small repository is intended for use as a minimalist, easy-to-use
view of the wet AMD data.  These snapshots can be used to quickly check
in our snRNAseq data for different biological signals, and also to prototype
any downstream computational analysis compatible with Seurat objects.

A copy of this repository without data will also be made available on the
YCFNI github account for reference.

# Usage
To use this repository, you will first need to install the necessary command
line tools to allow Seurat to run on your computer.  

We have found that for macOS, you will need to install a few things
prior to attempting the installation. First you will need to install
the Xcode developer tools for mac. Please follow the instructions
[here](https://mac.install.guide/commandlinetools/index.html)

You can run:

    xcode-select -p

in a command prompt, to check if the command line utilities have been
properly installed.

You will also likely need to install Homebrew, a package manager for
mac. To check if homebrew is installed, type

    brew

on a command prompt (Terminal) window. If the command is not found,
you should install homebrew by running

    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

You will then need to install several other dependencies using homebrew
prior to running the installation from RStudio.  These commands should be
executed from a normal (Terminal) command prompt (bash) and not from the
R console.

    brew install openssl

Now in the R console, you may need to install.

    install.packages('httr')
    install.packages('igraph')
    install.packages('leiden')
    install.packages('Rtsne')
    install.packages('sctransform')
    install.packages('SeuratObject')
    install.packages('RcppEigen')

Please follow the
installation instructions for Seurat 
[here](https://satijalab.org/seurat/articles/install.html).

TL;DR
    
    install.packages('Seurat')

from the `R` console.

For a more unified graphical interface, with easy-to-save plots and
interactive visualizations for dataframes, you may wish to install and
use [RStudio](https://www.rstudio.com/)

# Getting Started

Before you run any analysis or generate any plots, you will need to set the
working directory for your RStudio (or R instance) to the directory on your
machine where this folder was unzipped (i.e. the parent folder that
contains *this* README document, on *your* computer).

To learn how to do this, please read the RStudio documentation
[here](https://support.rstudio.com/hc/en-us/articles/200711843-Working-Directories-and-Workspaces-in-the-RStudio-IDE). This can
also be done programmatically through the 
[`setwd`](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/getwd)
R command.

Once you have set up your environment, you are ready to load and analyze the
data.

    library(Seurat)

    retina.combined <- readRDS('./data/retina_combined.rds')

Your Seurat object, analogous to the `pbmc` object generated in the
Seurat [guided clustering tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
is now available to you as `retina.combined` within the R environment.

# Common Tasks

You can (and *should*) get much more detail on many of these tasks by reading
the excellent [documentation](https://rdrr.io/cran/Seurat/api/) from the Satija
lab on RDRR for the CRAN package.

For the below, we will assume that your Seurat object is called
`retina.combined`.  If you modified this, please update the function calls
accordingly.

## Plot a UMAP labeled by cell type

    DimPlot(retina.combined, label=T)
    
## Plot the expression of a specific gene on UMAP

Replace with the name of the gene you would like to plot, prefixed by
`rna_` to indicate that you would like to plot raw count data.
For example, to plot the raw recovered counts for the `RLBP1` transcript:

    FeaturePlot(retina.combined, 'rna_RLBP1')

## Plot the expression of a specific gene as a violin plot by cell type

Replace with the name of the gene you would like to plot, prefixed by
`rna_` to indicate that you would like to plot raw count data.
For example, to plot the raw recovered counts for the `RLBP1` transcript:

    VlnPlot(retina.combined, 'rna_RLBP1')

To generate side-by-side plots split by disease group, you may invoke:

    VlnPlot(retina.combined, 'rna_RLBP1', split.by='DiseaseGroup')

# Some useful notes and tips

- `Seurat` plotting utilities are built on top
of `ggplot2` and can sometimes be updated using the `ggplot2` framework if
required.

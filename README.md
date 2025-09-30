# Independent Hypothesis Weighting

<p align="center">
<img src="https://github.com/Bioconductor/BiocStickers/blob/devel/IHW/IHW.png?raw=true"  width="200">
</p>





Independent hypothesis weighting (IHW) is a multiple testing procedure that increases power compared to the method of Benjamini and Hochberg by assigning
data-driven weights to each hypothesis. The input to IHW is a two-column
table of p-values and covariates. The covariate can be any continuous-valued
or categorical variable that is thought to be informative on the statistical
properties of each hypothesis test, while it is independent of the p-value
under the null hypothesis. IHW is described in the following paper:

> Nikolaos Ignatiadis, Bernd Klaus, Judith Zaugg, Wolfgang Huber (2016)
*Data-driven hypothesis weighting increases detection power in genome-scale multiple testing.*
Nature methods, 13(7), 577-580.

The following paper describes the theoretical underpinning of the method:

> Nikolaos Ignatiadis, Wolfgang Huber (2021)
*Covariate powered cross-weighted multiple testing.*
Journal of the Royal Statistical Society: Series B (JRSS-B), 83(4), 720-751.


# Software availability

The package is available on  [Bioconductor](https://www.bioconductor.org/packages/devel/bioc/html/IHW.html), and may be installed as follows:

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("IHW")
```

The package can be installed as follows with `devtools` from the Github repository:

```R
devtools::install_github("nignatiadis/IHW")
```

The file contains most of the updated files for the IHW Binning Below are the descriptions for each folder:

R -method_new_binning_admm.R: the new method that we are using right now -ihw_convex.R: the first version IHW ADMM method that we desert -plots_new_naive.R: the naive version of the plot function that plots the weights of the IHW Binning method in method_new_binning_admm.R file -plots.R: the old plot function that plot the weights of IHW ADMM -weights.R: contains the weight function that transform the thresholds into weights

Real Dataset Result: this folder contains a rmd file Three_dataset_visualization.Rmd that apply IHW OLD and IHW Binning to three datasets described in the overleaf file

Simulation Results: this folder contains the simulation results described in the overleaf file. All the simulation files should be run on GPU and some of the loading packages and files in the R files might need to change the path to successfully run

-Simulation_ttest.R, Simulation_wass.R, Simulation_null.R: three simulation files that corresponds to the three simulation settings in the overleaf file (see the Simulation section) -Simulation_visualization.Rmd: a Rmarkdown file that plots the power and FDR control results of the two methods IHW OLD and IHW Binning for the three simulation settings (some of the load commands might need to change path in order to run)

-Simulation_nbins_runtime_comparison.R: compares the runtime of the two methods IHW OLD and IHW Binning with varying number of bins -Simulation_nrowx_runtime_comparison.R: compares the runtime of the two methods IHW OLD and IHW Binning with varying number of rows -Simulation_runtime_visualization.Rmd: a Rmarkdown file that plots the run time comparison

IHWpaper-master: contains R codes for IHWpaper

IHWold: contains codes for IHW OLD





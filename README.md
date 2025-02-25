# LipidTrend
"LipidTrend" is an R package that implements a permutation-based statistical test 
to identify significant differences in lipidomic features between groups. 
The test incorporates Gaussian kernel smoothing of region statistics to improve stability 
and accuracy, particularly when dealing with small sample sizes. 
This package also includes two plotting functions for visualizing significant tendencies in 
one-dimensional and two-dimensional feature data, respectively.

## Installation 
To install `LipidTrend`, ensure that you have R 4.5.0 or later installed 
(see the R Project at [http://www.r-project.org](http://www.r-project.org)) 
and are familiar with its usage.

`LipidTrend` package is available on Bioconductor repository
[http://www.bioconductor.org](http://www.bioconductor.org). 
Before installing `LipidTrend`, you must first install the core Bioconductor 
packages. If you have already installed them, you can skip this step.
```{r install_Bioconductor, eval=FALSE}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install()
```

Once the core Bioconductor packages are installed, you can proceed with 
installing `LipidTrend`.
```{r install_package, eval=FALSE}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("LipidTrend")
``` 

After the installation is complete, you’re ready to start using `LipidTrend`.

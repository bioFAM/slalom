[![Build Status](https://travis-ci.org/bioFAM/slalom.svg?branch=master)](https://travis-ci.org/bioFAM/slalom)   [![Documentation Status](https://readthedocs.org/projects/f-scLVM/badge/?version=latest)](http://f-scLVM.readthedocs.io/en/latest/?badge=latest)

# The factorial single-cell latent variable model (slalom)

## What is slalom
slalom is a scalable modelling framework for single-cell RNA-seq data that uses gene set annotations to dissect single-cell transcriptome heterogeneity, thereby allowing to identify biological drivers of cell-to-cell variability and model confounding factors.

## Philosophy

Observed heterogeneity in single-cell profiling data is multi-factorial. slalom provides an efficient framework for unravelling this heterogeneity by simultaneously inferring latent factors that reflect annotated factors from pathway databases, as well as unannotated factors that capture variation outside the annotation.
slalom builds on sparse factor analysis models, for which this implementation provides efficient approximate inference using Variational Bayes, allowing the application of slalom to very large datasets containing up to 100,000 cells.

## Implementation
We provide 2 implementation of the slalom model: an R/C++ implementation that is [available on Bioconductor](https://bioconductor.org/packages/devel/bioc/html/slalom.html) and a python implementation. Both implemetations implement the model described in the accompanying publication [1]. 

Software by Florian Buettner, Davis McCarthy and Oliver Stegle. 

## R Implmentation
The `slalom` R package is available Bioconductor, so the most reliable way
to install the package is to use the usual Bioconductor method:

```{R}
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("slalom")
```
The source code for the R package can be found in the [R_package folder](https://github.com/bioFAM/slalom/blob/master/R_package) of this repository.


## Python implementation
### Installation requirements python implementation:

slalom requires Python >= 2.7 or newer with
  - scipy, h5py, numpy, matplotlib, scikit-learn, re
  
slalom can be installed via pip with `pip install slalom`.
For best results, we recommend the [ANACONDA](https://anaconda.org) python distribution.


### How to use slalom?
The current software version should be considered as beta. More extensive documentation, tutorials and examples will be available soon. 

For an illustration of how slalom can be applied to mESC data considered in Buettner et al. [1], we have prepared a [notebook](http://nbviewer.ipython.org/github/bioFAM/slalom/blob/master/ipynb/f-scLVM.ipynb). Along with other notebooks, this illustrates example analyses/workflows with slalom that you can read, download and adapt for your own analyses. These notebooks can be viewed and downloaded from [here](http://nbviewer.ipython.org/github/bioFAM/slalom/blob/master/ipynb/) or [here](https://github.com/bioFAM/slalom/tree/master/ipynb).

Documentation of the code can be found [here](http://f-scLVM.readthedocs.io).
## References:

[1] Buettner, F.,Pratanwanich, N., Marioni, J., Stegle, O. Scalable latent-factor models applied to single-cell RNA-seq data separate biological drivers from confounding effects. [Submitted](http://biorxiv.org/content/early/2016/11/15/087775).




## License
See [Apache License (Version 2.0, January 2004)](https://github.com/bioFAM/slalom/blob/master/license.txt).

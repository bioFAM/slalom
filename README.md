[![Build Status](https://travis-ci.org/PMBio/f-scLVM.svg?branch=master)](https://travis-ci.org/PMBio/f-scLVM)   [![Documentation Status](https://readthedocs.org/projects/f-sclvm/badge/?version=latest)](http://f-sclvm.readthedocs.io/en/latest/?badge=latest)

#The factorial single-cell latent variable model (f-scLVM)


##What is f-scLVM?

f-scLVM is a scalable modelling framework for single-cell RNA-seq data that uses gene set annotations to dissect single-cell transcriptome heterogeneity, thereby allowing to identify biological drivers of cell-to-cell variability and model confounding factors.


Software by Florian Buettner and Oliver Stegle. For detail please see the accompanying publication [1]. 

##Philosophy

Observed heterogeneity in single-cell profiling data is multi-factorial. f-scLVM provides an efficient framework for unravelling this heterogeneity by simultaneously inferring latent factors that reflect annotated factors from pathway databases, as well as unannotated factors that capture variation outside the annotation.
f-scLVM builds on sparse factor analysis models, for which this implementation provides efficient approximate inference using Variational Bayes, allowing the application of f-scLVM to very large datasets containing up to 100,000 cells.

##Installation requirements:

f-scLVM requires Python 2.7 or newer with
  - scipy, h5py, numpy, matplotlib, scikit-learn, re
  
f-scLVM can be installed via pip with `pip install fscLVM`.
For best results, we recommend the [ANACONDA](https://anaconda.org) python distribution.


##How to use f-scLVM?
The current software version should be considered as beta. More extensive documentation, tutorials and examples will be available soon. 

For an illustration of how f-scLVM can be applied to mESC data considered in Buettner et al. [1], we have prepared a notebook that can be viewed [interactively](http://nbviewer.ipython.org/github/pmbio/f-scLVM/blob/master/ipynb/f-scLVM.ipynb).

Documentation of the code can be found [here](http://f-sclvm.readthedocs.io).
##References:

[1] Buettner, F.,Pratanwanich, N., Marioni, J., Stegle, O. Scalable latent-factor models applied to single-cell RNA-seq data separate biological drivers from confounding effects. [Submitted](http://biorxiv.org/content/early/2016/11/15/087775).




##License
See [Apache License (Version 2.0, January 2004)](https://github.com/PMBio/f-scLVM/blob/master/license.txt).

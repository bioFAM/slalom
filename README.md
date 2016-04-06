#f-scLVM


##What is f-scLVM?

scLVM is a scalable modelling framework for single-cell RNA-seq data that can be used to simultaneously dissect the observed heterogeneity into different sources, thereby allowing for the identification of specific drivers of cell-to-cell variability.


Software by Florian Buettner and Oliver Stegle. f-scLVM is explained in  detail in the accompanying publication [1].

##Philosophy

Observed heterogeneity in single-cell profiling data is multi-factorial. f-scLVM provides an efficient framework for unravelling this heterogeneity by simultaneously inferring latent factors reflecting a large number of potential sources of variability. scLVM builds on sparse factor analysis models. We implement an efficient approximate inference scheme which allows the application of f-scLVM to very large datasets containing up to 100,000 cells.

##Installation requirements:

f-scLVM requires Python 2.7 with
  - scipy, h5py, numpy, pylab

##How to use f-scLVM?
The current software version should be considered as beta. More extensive documentation, tutorials and examples will be available soon. 

For an illustration of how f-scLVM can be applied to mESC data considered in Buettner et al. [1], we have prepared a notebook that can be viewed [interactively](http://nbviewer.ipython.org/github/pmbio/scLVM2/blob/master/py/demo/ipynb/f-scLVM.ipynb).

#The factorial single-cell latent variable model (f-scLVM)


##What is f-scLVM?

f-scLVM is a scalable modelling framework for single-cell RNA-seq data that can be used to dissect and model single-cell transcriptome heterogeneity, thereby allowing to identify biological drivers of cell-to-cell variability and model confounding factors.


Software by Florian Buettner and Oliver Stegle. f-scLVM is explained in  detail in the accompanying publication [1]. 

##Philosophy

Observed heterogeneity in single-cell profiling data is multi-factorial. f-scLVM provides an efficient framework for unravelling this heterogeneity by simultaneously inferring latent factors reflecting a large number of potential sources of variability. scLVM builds on sparse factor analysis models. We implement an efficient approximate inference scheme which allows the application of f-scLVM to very large datasets containing up to 100,000 cells.

##Installation requirements:

f-scLVM requires Python 2.7 with
  - scipy, h5py, numpy, matplotlib, scikit-learn, re
  
f-scLVM can be installed via pip with `pip install fscLVM`.

##How to use f-scLVM?
The current software version should be considered as beta. More extensive documentation, tutorials and examples will be available soon. 

For an illustration of how f-scLVM can be applied to mESC data considered in Buettner et al. [1], we have prepared a notebook that can be viewed [interactively](http://nbviewer.ipython.org/github/pmbio/f-scLVM/blob/master/py/demo/ipynb/f-scLVM.ipynb).

Documentation of the code can be found [here](https://github.com/PMBio/f-scLVM/tree/master/doc/_build/html/index.html).
##References:

[1] Buettner, F.,Pratanwanich, N., Marioni, J., Stegle, O. Scalable latent-factor models applied to single-cell RNA-seq data separate biological drivers from confounding effects. Submitted 

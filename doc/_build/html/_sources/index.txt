fscLVM documentation
======================

f-scLVM is a scalable modelling framework for single-cell RNA-seq data that can be used to dissect and model single-cell transcriptome heterogeneity, thereby allowing to identify biological drivers of cell-to-cell variability and model confounding factors.

Software by Florian Buettner and Oliver Stegle. f-scLVM is explained in detail in the accompanying publication [1].


[1] Buettner, F.,Pratanwanich, N., Marioni, J., Stegle, O. Scalable latent-factor models applied to single-cell RNA-seq data separate biological drivers from confounding effects. Submitted

Tutorial
--------

All steps required to run f-scLVM are illustrated in a jupyter notebook that can be viewed `interactively <http://nbviewer.jupyter.org/github/pmbio/f-scLVM/blob/master/py/demo/ipynb/f-scLVM.ipynb>`_. A brief tutorial on how to load data, fit the model and plot the results is summarised in the following secitons.

******
Input  
******
f-scLVM expects a hdf file containing the normalised, log transformed gene expression data as well as a set of annotations. We provide an R script (LINK) that can be used to generate this input from a gene expression matrix, using annotation from REACTOME or MSigDB. 

Usage is also illustrated in the an ipython notebook. 

We provide a ``load_hdf5`` function for loading the input data form such hdf file.

********************************
Model initialisation and fitting
********************************

An f-scLVM model can be initialised using the ``intiFA`` convenience function. Arguments can be used to specify options, incuding number of unannotated factors (`nHidden`), minimum number of genes in a pathway (`minGenes`), whether to use the fast option by pruning genes (`pruneGenes`), noise model (`noise`) and the data directory (`data_dir`). Once a model is initialised, it can be fit using the ``iterate`` method.

********************************
Diagnostics, plotting and saving.
********************************

The ``printDiagnostics`` function can be used to print diagnositcs based ont eh number of genes included/excluded by the model.  If more than 100% of gene annotations are changed for at least one factor, the model should be re-fitted with sparse unannotated facotrs.



Loading data and model initialisation
-------------------------------------

.. autofunction:: fscLVM.load_txt
.. autofunction:: fscLVM.load_hdf5
.. autofunction:: fscLVM.initFA
.. autofunction:: fscLVM.preTrain



Model fitting
-------------------------------------
.. autoclass:: fscLVM.core.CSparseFA
    :members:


Plotting and saving results
----------------------------
.. autofunction:: fscLVM.plotTerms
.. autofunction:: fscLVM.plotFactors
.. autofunction:: fscLVM.plotFA
.. autofunction:: fscLVM.saveFA
.. autofunction:: fscLVM.dumpFA

Contents:

.. toctree::
   :maxdepth: 2

fscLVM documentation
======================

f-scLVM is a scalable modelling framework for single-cell RNA-seq data that can be used to dissect and model single-cell transcriptome heterogeneity, thereby allowing to identify biological drivers of cell-to-cell variability and model confounding factors.

Software by Florian Buettner and Oliver Stegle. f-scLVM is explained in detail in the accompanying publication [1].


[1] Buettner, F.,Pratanwanich, N., Marioni, J., Stegle, O. Scalable latent-factor models applied to single-cell RNA-seq data separate biological drivers from confounding effects. Submitted


Quickstart
----------

******
Input  
******
f-scLVM requires two input files, a gene expression file and an annotation file. The gene expression file is a text file containing the normalised, log-transformed gene expression matrix, with every row corresponding to a cell. Column names should be gene identifiers matching those in the annotation file (i.e. if gene symbols are used in the annotation file, column names should also be gene symbols). Row names are optional and can eg be a known covariate which is used for plotting.
The annotation file is a text file with every row containing the name of a gene set, followed by the gene identifiers annotated to that gene set. We recommend using annotations such as those published in the REACTOME database or the Molecular signature database (MSigDB) and  provide an annotation file containing annotations form the RECTOME database. The license of MSigDB does not permit redistribution of the raw annotation files. To use MSigDB annotations, please register at http://software.broadinstitute.org/gsea/msigdb, download the hallmark gene sets (gene symbols) and place the file in data folder.
Both text files can then be loaded using the ``load_text`` function.

NB: f-scLVM works best on a subset of highly avriable genes; these can be identified using a variance filter based on a mean-variance trend fitted using spike-in transcripts or endogenous genes. A step-by-step workflow on low-level processing of scRNA-seq data (including gene filtering) can be found here: https://f1000research.com/articles/5-2122/v2


********************************
Model initialisation and fitting
********************************

An f-scLVM model can be initialised using the ``initFA`` convenience function. Arguments can be used to specify options, incuding number of unannotated factors (`nHidden`), minimum number of genes in a pathway (`minGenes`), whether to use the fast option by pruning genes (`pruneGenes`), noise model (`noise`) and the data directory (`data_dir`). Once a model is initialised, it can be fit using the ``train`` method.

********************************
Diagnostics, plotting and saving.
********************************

The ``printDiagnostics`` function can be used to print diagnositcs based ont eh number of genes included/excluded by the model.  If more than 100% of gene annotations are changed for at least one factor, the model should be re-fitted with sparse unannotated facotrs.

Tutorial
--------

All steps required to run f-scLVM are illustrated in a jupyter notebook that can be viewed `interactively <http://nbviewer.jupyter.org/github/pmbio/f-scLVM/blob/master/ipynb/f-scLVM.ipynb>`_. 


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
.. autofunction:: fscLVM.plotRelevance
.. autofunction:: fscLVM.saveFA
.. autofunction:: fscLVM.dumpFA

Contents:

.. toctree::
   :maxdepth: 2

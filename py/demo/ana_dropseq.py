# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 13:37:12 2016

@author: flo
"""

import sys
import os
import ast
sys.path.append('./..')
import scipy as SP
import cPickle as pickle
import core.fscLVM as fscLVM
import core.utils as utils
import h5py
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector, StrVector
import rpy2.robjects as robjects
stats = importr('stats')
from rpy2.robjects.lib import ggplot2
grdevices = importr('grDevices')
cowplot = importr('cowplot')
gg2 = importr('ggplot2')
from sklearn.manifold import TSNE


data_dir = '../../data/'
out_base = './results'
dataset = 'retina'
it = 900

if __name__ == '__main__':
    YFile = h5py.File(os.path.join(data_dir, 'Y'+dataset+'.hdf5'), 'r')
    Y = YFile['Y'][:]
    dataFile = h5py.File(os.path.join(data_dir, dataset+'.hdf5'), 'r')
    
    
    clusterFile = h5py.File(os.path.join(data_dir, 'retina_clusters.hdf5'), 'r')   
    cluster = clusterFile['clusIDall'][:]
    out_name = os.path.join(out_base, dataset,'MSigDB_fast/resHidden399_Bias0_Known__nExpressed_Sort1e52k_it_'+str(it)+'.hdf5')
    res = utils.loadFA(out_name)    
    idx_genes  = SP.logical_and(SP.sum(dataFile['Pi'][:]>.5,0)>0, Y.mean(0)>0.)#SP.any(pi>.5,1)
    Yhet = Y[:,idx_genes]   
    
    
    utils.plotTerms( S=res['S'][:], alpha=res['alphaRaw'][:], terms=res['terms'][:])
    
    utils.plotFactors(0,1, X = res['S'][:],  lab=cluster, terms=res['terms'][:])
    
    idx_sample = SP.unique(SP.random.random_integers(0,Y.shape[0],20000))
    tsne = TSNE(2,angle=.8, n_iter=200)
    Xtsne = tsne.fit_transform(Yhet[idx_sample,:])
    
        
    
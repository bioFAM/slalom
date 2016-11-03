import sys
import os
import scipy as SP
sys.path.insert(0,'./../../')
import cPickle as pickle
from fscLVM.utils import *
import fscLVM.core as fscLVM
import pdb
import h5py

#specify where the hdf5 file is
data_dir = '../../../data/'
out_base = './../results/'

dFile = 'Buettneretal.hdf5'

data = load_hdf5(dFile, data_dir=data_dir)
#my hack to load old data file
I = 1.0*(data['Pi']>0.5).T
Y = data['Y']
terms = data['terms']

pdb.set_trace()
FA = initFA(Y, terms,I,noise='gauss', nHidden=3, minGenes=15)
#iterate
FA.iterate(nIterations=2000)
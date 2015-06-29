# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 12:59:02 2015
Run permutations to assess null distribtion of alphas
Then come up with heuristics to switch factors off (the non-relevant ones)
@author: flo
"""
import time
import io
import sys
import os
import pdb
import glob
sys.path.append('./..')
import scipy as SP
import pylab as PL
import pdb
import core.sparseFAard as sparseFA
from core.misc import *

python_cmd = 'python'
nthreads = 1
mem_thread = 10000
cluster_cmd = 'bsub -n %d -q research-rh6 -R "rusage[mem=%d,tmp=10000]" -M %d -o ./cluster_out' % (nthreads,mem_thread,mem_thread)

if __name__ == '__main__':
    if 'cluster' in sys.argv:
        for r in range(100):
            outfile = './out/res'+str(r)+'.hdf5'
            cmd = '%s %s %s %s' % (cluster_cmd, python_cmd, sys.argv[0],out_file)
            print cmd
            os.system(cmd)
            time.sleep(sleep)
    #1. load data
    elif 'collect' in sys.argv:
        dir_name_out = sys.argv[2]
        FL = glob.glob(os.path.join(dir_name_out,'res*.hdf5'))
        piList = list()
        alphaList = list()
        boundList = list()
        sample_res = {}
        for fn in FL:
            f = h5py.File(fn,'r')
            for key in f.keys():
                if key in sample_res.keys()==False:
                    sample_res[key]=list()
                sample_res[key].append(f[key][:])

        
    else: 
        out_file = sys.argv[1]
        terms = SP.loadtxt('../../R/terms_GESA.txt', dtype='str', skiprows=1)
        Y = SP.loadtxt('Ystd_hscGESA_sfMmus.csv')
        pi = SP.loadtxt('pistd_hscGESA_sfMmus.csv')
        idx_het = SP.loadtxt('../../R/idx_het_sfMmus.txt', skiprows=1, dtype='int') -1
        Y = Y[:,idx_het]
        pi = pi[idx_het,:]
        pi[pi>.5] =1.
        pi[pi<.5] =0.
        idx_aff = SP.logical_or(pi.sum(1)==0, Y.sum(0)==0.)
        Y = Y[:,~idx_aff]
        pi = pi[~idx_aff,:]
        terms = terms[pi.sum(0)>10]
        pi = pi[:,pi.sum(0)>10]
        
    
        #run "true terms"
        K = pi.shape[1]
    
        #data for sparseFA instance
        init={'init_data':sparseFA.CGauss(Y),'Pi':pi}
        sigmaOff = 1E-5
        sparsity = 'VB'
        nIterations = 2000
        #permutation move
        permutation_move = False
        #prior on noise level 
        priors = {'Eps': {'priors':[1,100]}}
        #how to initialize network?
        initType = 'pca'
        FAtrue = sparseFA.CSparseFA(components=K,sigmaOff=sigmaOff,sigmaOn=SP.ones(pi.shape[1])*1.0,sparsity=sparsity,nIterations=nIterations,permutation_move=permutation_move,priors=priors,initType=initType)
    
        FAtrue.init(**init)
        
        #do permutations
        Non = FAtrue.Non
        Nperm = 20

            
        
        piRand_ = SP.zeros(pi.shape)
        indRand_ = SP.array([SP.random.choice(range(pi.shape[1]),Non[i], replace=False) for i in range(len(Non))])
        for k in range(pi.shape[0]):
            piRand_[indRand_[k],k] = 1.
        init={'init_data':sparseFA.CGauss(Y),'Pi':piRand_}
        FArand_ = sparseFA.CSparseFA(components=K,sigmaOff=sigmaOff,sigmaOn=SP.ones(piRand_.shape[1])*1.0,sparsity=sparsity,nIterations=nIterations,permutation_move=permutation_move,priors=priors,initType=initType)
    
        FArand_.init(**init)
        FArand_.iterate()

        res = {}
        res['boundRand'] = FArand_.calcBound()        
        res['alphaRand'] = FArand_.Alpha.E1
        res['piRand'] = piRand
        
        fout = h5py.File(out_file,'w')
        smartDumpDictHdf5(res, fout)
        
            
    
    



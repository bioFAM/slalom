# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 12:59:02 2015
Run permutations to assess null distribtion of alphas
Then come up with heuristics to switch factors off (the non-relevant ones)
@author: flo
"""

import sys
sys.path.append('./..')
import scipy as SP
import pylab as PL
import pdb
import core.sparseFAard as sparseFA



if __name__ == '__main__':

    #1. load data
    terms = SP.loadtxt('/Users/flo/projects/Auto_Bionf/R/terms_GESA.txt', dtype='str', skiprows=1)
    Y = SP.loadtxt('Ystd_hscGESA_sfMmus.csv')
    pi = SP.loadtxt('pistd_hscGESA_sfMmus.csv')
    idx_het = SP.loadtxt('/Users/flo/projects/Auto_Bionf/R/idx_het_sfMmus.txt', skiprows=1, dtype='int') -1
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
    #FA.S.E1 = Sinit
    FAtrue.iterate()
    boundTrue = FAtrue.calcBound()
    
    alphaTrue = FAtrue.Alpha.E1
    print alphaTrue
    
    #do permutations
    Non = FAtrue.Non
    Nperm = 20
    piRand = list()
    alphaRand = list()
    boundRand = list()
        
    
    for r in range(20):
        piRand_ = SP.zeros(pi.shape)
        indRand_ = SP.array([SP.random.choice(range(pi.shape[1]),Non[i], replace=False) for i in range(len(Non))])
        for k in range(pi.shape[0]):
            piRand_[indRand_[k],k] = 1
        piRand.append(piRand_)
        init={'init_data':sparseFA.CGauss(Y),'Pi':piRand_}
        FArand_ = sparseFA.CSparseFA(components=K,sigmaOff=sigmaOff,sigmaOn=SP.ones(piRand_.shape[1])*1.0,sparsity=sparsity,nIterations=nIterations,permutation_move=permutation_move,priors=priors,initType=initType)
    
        FArand_.init(**init)
        #FA.S.E1 = Sinit
        FArand_.iterate()
        print r 'folds finished'
        boundRand.append(FArand_.calcBound())
        
        alphaRand.append(FArand_.Alpha.E1)
        
    
    



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
from matplotlib import pyplot as plt

def plot_heatMap(FA,thre=10, row_labels=None):        
    thre=10
    data=SP.corrcoef(FA.S.E1[:,FA.Alpha.E1<thre].transpose())
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(data,  vmin=-1., vmax=1.)
#    plt.colorbar(fig)
    # put the major ticks at the middle of each cell
    ax.set_xticks(SP.arange(data.shape[0])+0.5, minor=False)
    ax.set_yticks(SP.arange(data.shape[1])+0.5, minor=False)
    
    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    if row_labels!=None:
        ax.set_xticklabels(row_labels, minor=False, rotation=75)
        ax.set_yticklabels(row_labels, minor=False)
    plt.show()



if __name__ == '__main__':
    #1. load data
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
    nIterations = 301
    #permutation move
    permutation_move = False
    #prior on noise level 
    priors = {'Eps': {'priors':[1,100]}}
    #how to initialize network?
    initType = 'prior'
    FAtrue = sparseFA.CSparseFA(components=K,sigmaOff=sigmaOff,sigmaOn=SP.ones(pi.shape[1])*1.0,sparsity=sparsity,nIterations=nIterations,permutation_move=permutation_move,priors=priors,initType=initType)

    FAtrue.init(**init)
    #FA.S.E1 = Sinit
    FAtrue.iterate(tolerance=1e-3)
    boundTrue = FAtrue.calcBound()
    
    alphaTrue = FAtrue.Alpha.E1
    print alphaTrue

    #plot
    thre=10
    
    row_labels = terms[alphaTrue<thre]
    plot_heatMap(FAtrue,thre=10,row_labels=row_labels)
    
    
    #recursively eliminate factors
    #start with the one with hightst alpha
    FAlist = list()
    termsList = list()    
    boundList = list()
    tol=1e-3
    minFac=3    
    FAcurr = FAtrue
    alpha_ = FAcurr.Alpha.E1
    terms_ = terms
    inds = range(len(terms))
    deltaBound=tol
    numF = 0
    while(len(alpha_)>minFac and deltaBound>=tol):
        FAold = FAcurr
        bound_ = FAcurr._bound
        alpha_ = FAcurr.Alpha.E1
        terms_ = terms_[inds]
        FAlist.append(FAold)
        termsList.append(terms_)
        boundList.append(bound_)
        print terms_[SP.argsort(alpha_)]  
        
        inds = SP.setdiff1d(range(len(alpha_)),SP.argmax(alpha_))#SP.array([1,9])
        permutation_move = False
        init_factors = {}
        init_factors['S'] = FAold.S.E1[:,inds]
        init_factors['W'] = FAold.W.E1[:,inds]
        init2={'init_data':sparseFA.CGauss(Y),'Pi':pi[:,inds], 'init_factors':init_factors}
        init_type = 'pca'
        FAcurr = sparseFA.CSparseFA(components=len(inds),sigmaOff=sigmaOff,sigmaOn=SP.ones(len(inds))*1.0, sparsity=sparsity,nIterations=nIterations,permutation_move=permutation_move,priors=priors,initType=initType)    
        FAcurr.init(**init2)
        FAcurr.iterate(tolerance=1e-10)
        deltaBound = FAcurr._bound - bound_
        numF+=1
        print '%i factors excluded.' % (numF)

        
 
#    #do permutations
#    Non = FAtrue.Non
#    Nperm = 20
#    piRand = list()
#    alphaRand = list()
#    boundRand = list()
#        
#    
#    for r in range(20):
#        piRand_ = SP.zeros(pi.shape)
#        indRand_ = SP.array([SP.random.choice(range(pi.shape[1]),Non[i], replace=False) for i in range(len(Non))])
#        for k in range(pi.shape[0]):
#            piRand_[indRand_[k],k] = 1
#        piRand.append(piRand_)
#        init={'init_data':sparseFA.CGauss(Y),'Pi':piRand_}
#        FArand_ = sparseFA.CSparseFA(components=K,sigmaOff=sigmaOff,sigmaOn=SP.ones(piRand_.shape[1])*1.0,sparsity=sparsity,nIterations=nIterations,permutation_move=permutation_move,priors=priors,initType=initType)
#    
#        FArand_.init(**init)
#        #FA.S.E1 = Sinit
#        FArand_.iterate()
#        print r 'folds finished'
#        boundRand.append(FArand_.calcBound())
#        
#        alphaRand.append(FArand_.Alpha.E1)
#        
#    
#    
#
#

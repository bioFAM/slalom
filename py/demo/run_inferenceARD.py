"""
run inference on tody data
"""

import sys
sys.path.append('./..')
import scipy as SP
import pylab as PL
import pdb
import h5py
import core.sparseFAvem as sparseFA



if __name__ == '__main__':

#    #1. load data
#    terms = SP.loadtxt('/Users/flo/projects/Auto_Bionf/R/terms_GESA.txt', dtype='str', skiprows=1)
#    Y = SP.loadtxt('Ystd_hscGESA_sfMmus.csv')
#    pi = SP.loadtxt('pistd_hscGESA_sfMmus.csv')
#    idx_het = SP.loadtxt('/Users/flo/projects/Auto_Bionf/R/idx_het_sfMmus.txt', skiprows=1, dtype='int') -1
#    Y = Y[:,idx_het]
#    pi = pi[idx_het,:]
#    pi[pi>.5] =1.
#    pi[pi<.5] =0.
#    idx_aff = SP.logical_or(pi.sum(1)==0, Y.sum(0)==0.)
#    Y = Y[:,~idx_aff]
#    pi = pi[~idx_aff,:]
#    terms = terms[pi.sum(0)>10]
#    pi = pi[:,pi.sum(0)>10]
#    
#    
#    
#    Y = SP.loadtxt('Ystd_hsc9_sfMmus.csv')
#    pi = SP.loadtxt('pistd_hsc9_sfMmus.csv')
#    idx_het = SP.loadtxt('/Users/flo/projects/Auto_Bionf/R/idx_het_sfMmus.txt', skiprows=1, dtype='int') -1
    #Y = Y[:,idx_het]
    #pi = pi[idx_het,:]
    
    
    dataFile = h5py.File('../../../data/zeisel_microglia.hdf5', 'r')
    pi = dataFile['Pi'][:].T
    Y = SP.log2(dataFile['Yhet'][:].T+1)
    terms = dataFile['terms'][:]    
    
    pi[pi>.5] =0.98
    pi[pi<.5] =1e-10
    
    
    #idx_aff = SP.logical_or(pi.sum(1)==0, Y.sum(0)==0.)
 #   pi = pi[~idx_aff,:]
#    #Y = Y[:,~idx_aff]
 #   pi = pi[:,pi.sum(0)>10]
    #Y = SP.loadtxt('Ytoy2.csv')
    #pi = SP.loadtxt('pitoy2.csv')

    K = pi.shape[1]

    #data for sparseFA instance
    init={'init_data':sparseFA.CGauss(Y),'Pi':pi}
    sigmaOff = 1E-5
    sparsity = 'VB'
    nIterations = 3000
    #permutation move
    permutation_move = False
    #prior on noise level 
    priors = {'Eps': {'priors':[1,100]}}
    #how to initialize network?
    initType = 'pcaRand'
    FA = sparseFA.CSparseFA(components=K,sigmaOff=sigmaOff,sigmaOn=SP.ones(pi.shape[1])*1.0,sparsity=sparsity,nIterations=nIterations,permutation_move=permutation_move,priors=priors,initType=initType)

    FA.init(**init)
    #FA.S.E1 = Sinit
    FA.iterate(forceIterations=True)
    FA.calcBound()
    
    alpha = FA.Alpha.E1
    print alpha
    #look at variance explained and weigh by numnber of annotated genes
    # gfa_trimmed2$D/(alpha[,GF_i]*(gfa_trimmed2$datavar - gfa_trimmed2$D/gfa_trimmed2$tau))
#    Ion = FA.W.C[:,:,1]>.5
#    rel_contrib = SP.zeros(K)
#    for k in range(K):
#        rel_contrib[k] = FA.Non[k]/(alpha[k]*(Y[:,Ion[:,k]].var(0)-Ion[Ion[:,k],:].sum(1)/FA.Eps.E1[Ion[:,k]]).sum())    
##        residual[k] = (Y[:,Ion[:,k]].var(0)-1./FA.Eps.E1[Ion[:,k]]).sum()
#        
#        
#
#    #data for sparseFA instance
#    inds = SP.where(alpha<11)[0]#SP.array([1,9])
#    permutation_move = False
#    init_factors = {}
#    init_factors['S'] = FA.S.E1[:,inds]
#    init_factors['W'] = FA.W.E1[:,inds]
#    init2={'init_data':sparseFA.CGauss(Y),'Pi':pi[:,inds], 'init_factors':init_factors}
#    init_type = 'data'
#    FA2 = sparseFA.CSparseFA(components=len(inds),sigmaOff=sigmaOff,sigmaOn=SP.ones(len(inds))*1.0, sparsity=sparsity,nIterations=nIterations,permutation_move=permutation_move,priors=priors,initType=initType)
#
#    FA2.init(**init2)
##    FA2.S.E1 = FA.S.E1[:,inds]
##    FA2.W.E1 = FA.W.E1[:,inds]    
#    
#    FA2.iterate()
#    print FA2.calcBound()
#    print FA2.Alpha.E1
#    
#
#    #results:
#    #factor activities
#    X = FA.S.E1
#    #weights
#    W = FA.W.E1
#    #sparsity pattern
#    Z = FA.W.C[:,:,1]
#
#    #SP.savetxt('X_inferenceSTDGESA_sfMmus.csv',X)
#    #SP.savetxt('W_inferenceSTDGESA_sfMmus.csv',W)
#    #SP.savetxt('Z_inferenceSTDGESA_sfMmus.csv',Z)
##    
##    SP.savetxt('X_inferenceGFA.csv',X)
##    SP.savetxt('Z_inferenceGFA.csv',Z)
##    SP.savetxt('W_inferenceGFA.csv',W)
##    
#    import h5py
#    #fres['Z'] = Z
#    fres= h5py.File('./resFAtoy.h5py', 'w')
#    fres['X'] = X
#    fres['W'] = W
#    fres.close()
#    
#
#    
#    

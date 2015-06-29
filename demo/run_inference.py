"""
run inference on tody data
"""

import sys
sys.path.append('./..')
import scipy as SP
import pylab as PL
import pdb
import core.sparseFA as sparseFA



if __name__ == '__main__':

    #1. load data
    #Y = SP.loadtxt('Ystd_hscGESA_sfMmus.csv')
    #pi = SP.loadtxt('pistd_hscGESA_sfMmus.csv')

    Y = SP.loadtxt('Ytoy.csv')
    pi = SP.loadtxt('pitoy.csv')
    
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

    K = pi.shape[1]

    #data for sparseFA instance
    init={'init_data':sparseFA.CGauss(Y),'Pi':pi, }
    sigmaOff = 1E-5
    sparsity = 'VB'#'EPV'
    nIterations = 50
    #permutation move
    permutation_move = False
    #prior on noise level 
    priors = {'Eps': {'priors':[1,100]}}
    #how to initialize network?
    initType = 'prior'
    FA = sparseFA.CSparseFA(schedule=['S','Eps','W'], components=K,sigmaOff=sigmaOff,sparsity=sparsity,nIterations=nIterations,permutation_move=permutation_move,priors=priors,initType=initType)

    FA.init(**init)
    FA.iterate()



    #data for sparseFA instance
    inds = SP.array([0,3])
    permutation_move = False
    init={'init_data':sparseFA.CGauss(Y),'Pi':pi, 'S':Zinit}
    FA = sparseFA.CSparseFA(components=len(inds),sigmaOff=sigmaOff,sparsity=sparsity,nIterations=nIterations,permutation_move=permutation_move,priors=priors,initType=initType)

    FA.init(**init)
    FA.iterate()
    FA.calcBound()
    
    

    #results:
    #factor activities
    X = FA.S.E1
    #weights
    W = FA.W.E1
    #sparsity pattern
    Z = FA.W.C[:,:,1]

    #SP.savetxt('X_inferenceSTDGESA_sfMmus.csv',X)
    #SP.savetxt('W_inferenceSTDGESA_sfMmus.csv',W)
    #SP.savetxt('Z_inferenceSTDGESA_sfMmus.csv',Z)
    
    SP.savetxt('X_inferencetoy.csv',X)
    SP.savetxt('Z_inferencetoy.csv',Z)
    SP.savetxt('W_inferencetoy.csv',W)
    
    import h5py
    #fres['Z'] = Z
    fres= h5py.File('./resFAtoy.h5py', 'w')
    fres['X'] = X
    fres['W'] = W
    fres.close()
    

    
    

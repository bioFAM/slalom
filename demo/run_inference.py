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
    Y = SP.loadtxt('Y.csv')
    pi = SP.loadtxt('pi.csv')

    K = pi.shape[1]

    #data for sparseFA instance
    init={'init_data':sparseFA.CGauss(Y),'Pi':pi}
    sigmaOff = 1E-5
    sparsity = 'EPV'
    nIterations = 100
    #permutation move
    permutation_move = False
    #prior on noise level 
    priors = {'Eps': {'priors':[1,100]}}
    #how to initialize network?
    initType = 'prior'
    FA = sparseFA.CSparseFA(components=K,sigmaOff=sigmaOff,sparsity=sparsity,nIterations=nIterations,permutation_move=permutation_move,priors=priors,initType=initType)

    FA.init(**init)
    FA.iterate()



    #results:
    #factor activities
    X = FA.S.E1
    #weights
    W = FA.W.E1
    #sparsity pattern
    Z = FA.Z.E1

    SP.savetxt('X_inference.csv',X)
    SP.savetxt('Z_inference.csv',Z)
    SP.savetxt('W_inference.csv',W)
    
    

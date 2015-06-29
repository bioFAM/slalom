"""create_toy_data
- create toy data for factor analysis infernece and store as CSV files
"""

import sys
import scipy as SP
import numpy.random as random

import scipy.io
import os
import pylab as PL
import pdb


if __name__ =='__main__':
    #number of samples
    N = 100
    #number of genes
    G = 50
    #number of factors
    K = 10
    #sparcity
    sparse = 0.01

    #uncertainty on the prior network to be generated
    FPR = 0.01
    FNR = 0.01

    #noise level
    eps = 0.1

    #1. create random netowrk prior
    Ion = (random.rand(G,K)<sparse)
    Ioff = ~Ion

    #2. create pi matrix with prior probability of a link
    pi  = SP.zeros([G,K])
    pi[Ion] = 1-FNR
    pi[Ioff] = FPR

    #3. sample form pi true network
    Z = 1.0*(random.rand(G,K)<pi)

    #4. create weight matrix
    W = random.randn(G,K)
    W *= Z

    #5. create factor matrix
    X = random.randn(N,K)

    Y = SP.dot(X,W.T) + eps*random.randn(N,G)

    #6. store results
    SP.savetxt('Yp.csv',Y)
    SP.savetxt('Xp.csv',X)
    SP.savetxt('Wp.csv',W)
    SP.savetxt('pip.csv',pi)
    SP.savetxt('Zp.csv',Z)
    

    
    # try toy data for GFA: look at a few overlappig categories
    # decide a few being on (with few genes to make it hard) and one with all genes on as "baseline"
    
     #number of samples
    N = 80
    #number of genes
    G = 500
    #number of factors
    
    #sparcity
    sparse = SP.hstack([SP.arange(0.05,0.55,0.1),0.9])
    K = len(sparse)

    #noise level
    eps = 0.1

    #1. create random netowrk prior
    Ion = SP.zeros((G,K))
    #pi  = SP.zeros([G,K])
    for k in range(K):
        onk_ = random.rand(G,)<sparse[k]
        if k==3:
            onk_[SP.where(Ion[:,1]==True)[0]]=True
        Ion[:,k] = onk_
        #pi[:,k] = sparse[k]

    Ion = (Ion>0)
    Ioff = ~Ion

    #2. create pi matrix with prior probability of a link
    pi  = SP.zeros([G,K])
    pi[Ion] = 1.0
    pi[Ioff] = 0.0
    pi = SP.hstack([pi, (random.rand(G,1)<.3)*1.0,(random.rand(G,1)<.7)*1.0,(random.rand(G,1)<.1)*1.0, pi[:,2:3]])


    #3. sample form pi true network
    Z = 1.0*Ion

    #4. create weight matrix
    W = random.randn(G,K)
    W[:,0] = W[:,0]/5.0
    #W[:,sparse==.11] = W[:,sparse==.11]/2.0
    W *= Z

    #5. create factor matrix
    X = random.randn(N,K)

    Y = SP.dot(X,W.T) + eps*random.randn(N,G)
   
    #6. store results
    SP.savetxt('Ytoy2.csv',Y)
    SP.savetxt('Xtoy2.csv',X)
    SP.savetxt('Wtoy2.csv',W)
    SP.savetxt('pitoy2.csv',pi)
    SP.savetxt('Ztoy2.csv',Z) 
    SP.savetxt('sparsetoy2.csv',sparse) 

    #fres['Z'] = Z
    fres= h5py.File('./toymodel.h5py', 'w')
    fres['X'] = X
    fres['W'] = W
    fres.close()


    
    
       

    


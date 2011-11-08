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
    SP.savetxt('Y.csv',Y)
    SP.savetxt('X.csv',X)
    SP.savetxt('W.csv',W)
    SP.savetxt('pi.csv',pi)
    SP.savetxt('Z.csv',Z)

    
    
       

    


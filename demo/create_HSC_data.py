"""create_test_data from HSCs
- create input data for factor analysis infernece and store as CSV files
"""

import sys
import scipy as SP
import h5py



if __name__ =='__main__':
    #load data from R
    #fdata = h5py.File('/Users/flo/projects/Auto_Bionf/R/dataHSCstd_GESA_sfMmus.hdf5', 'r')
    fdata = h5py.File('/Users/flo/projects/Auto_Bionf/R/dataHSCstd_9_sfMmus.hdf5', 'r')
    Y = fdata['Yscale'][:].transpose()
    Y[SP.isnan(Y)] = 0
    idx = []
    for key in fdata['group_list'].keys():
        idx.append(fdata['group_list'][key][:]-1)

        
    initGFA = SP.loadtxt('/Users/flo/projects/Auto_Bionf/R/probZ_1015.txt', skiprows=1)
    idx_het = SP.loadtxt('/Users/flo/projects/Auto_Bionf/R/idx_het_sfMmus.txt', skiprows=1, dtype='int32')-1

    #number of samples
    N = Y.shape[0]
    #number of genes
    G = Y.shape[1]
    #number of factors
    K = 2#len(idx)

    #uncertainty on the prior network to be generated
    FPR = 0.01
    FNR = 0.01


    #1. create netowrk prior
    Ion = SP.zeros((G,K))
    
    for k in range(2):        
        for idx_ in SP.where(initGFA[:,k]>0.15)[0]:
            Ion[idx[idx_],k] = 1
    Ion = (Ion==1)
    
    Ioff = ~Ion

    #2. create pi matrix with prior probability of a link
    pi  = SP.zeros([G,K])
    pi[Ion] = 1-FNR
    pi[Ioff] = FPR
    Y = Y[:,idx_het]
    pi = pi[idx_het,:]    
    
    Sinit = SP.loadtxt('/Users/flo/projects/Auto_Bionf/R/gfaZ_1015.txt', skiprows=1)
    Sinit = Sinit/Sinit.std(0)
  


    #6. store results
    SP.savetxt('Ystd_hscGFA_sfMmus.csv',Y)
    SP.savetxt('pistd_hscGFA_sfMmus.csv',pi)

    


    
    
       

    


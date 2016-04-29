# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 14:46:40 2016

@author: flo
"""
#import core.fscLVM as sparseFA

from sklearn.decomposition import RandomizedPCA,PCA
import h5py
import scipy as SP
import matplotlib as mpl
import matplotlib.lines as mlines
mpl.use('Agg')
import pylab as plt
import os
import brewer2mpl

data_dir = '../../../data/'
out_base = './../results/'



def mad(X):
    median = SP.median(X, axis=0)
    return SP.median(abs(X-median), axis=0)
    
def secdev(x):
    return 1/(2*SP.pi)*(SP.exp(-x*x/2.))*(x*x-1)
    
    
def saveFA(FA, idx=None, Ycorr=None, Ycorr_lat = None, F=None, Fcorr=None):
    MAD = mad(FA.S.E1)
    alpha02 = (MAD>.5)*(1/(FA.Alpha.E1))
    out_file = h5py.File(FA.out_name+'_it_'+str(FA.iterationCount)+'.hdf5','w')    
    out_file['alphaRaw'] = FA.Alpha.E1
    out_file['alpha02'] = alpha02  
    out_file['W'] = FA.W.E1
    out_file['Eps'] = FA.Eps.E1
    out_file['S'] = FA.S.E1        
    out_file['Gamma'] = FA.W.C[:,:,0]
    out_file['pi'] = FA.Pi
    out_file['terms'] = FA.terms
    if idx != None:
        out_file['idx'] = idx
    if Ycorr != None:
        out_file['Ycorr'] = Ycorr
    if Ycorr_lat != None:
        out_file['Ycorr_lat'] = Ycorr_lat 
    if Fcorr != None:
        out_file['Fcorr'] = Fcorr
    if F != None:
        out_file['F'] = F               
    out_file.close()    
    
    
    
def loadFA(out_name):
    out_file = h5py.File(out_name,'r')    
    res = {}   
    for key in out_file.keys():    
        res[key] = out_file[key]    
    return res
    
    
def plotFactors(idx1, idx2,FA=None, X = None,  lab=None, terms=None, cols=None, isCont=True,madFilter=0.5):
    if FA!=None:
        print FA
        MAD = mad(FA.S.E1)
        alpha = (MAD>madFilter)*(1/(FA.Alpha.E1))
        idxF = SP.argsort(-alpha)    
        X1 = FA.S.E1[:,idxF[idx1]]
        X2 = FA.S.E1[:,idxF[idx2]]
    else:
        X1 = X[:,idx1]
        X2 = X[:,idx2]        
    
    if isCont==False:
        uLab = SP.unique(lab)  
        if cols==None:
            try:
                import brewer2mpl
            except ImportError:
                print 'Specify colors using the cols argument or install the brewer2mpl module'
            bmap=brewer2mpl.get_map('Paired', 'Qualitative', len(uLab))
            cols = bmap.hex_colors         
        pList=list()
        for i in range(len(X1)):
            pList.append(plt.plot(X1[i], X2[i], '.',color=cols[SP.where(lab[i]==uLab)[0]]))
        plt.xlabel(terms[idxF[idx1]])
        plt.ylabel(terms[idxF[idx2]])
        lList=list()
        for i in range(len(uLab)):
            lList.append( mlines.Line2D([], [], color=cols[i], marker='.',
                              markersize=7, label=uLab[i], linewidth=0))     
        plt.legend(handles=lList)
    else:
        plt.scatter(X1, X2, c=lab, s=20)
        plt.xlabel(terms[idxF[idx1]])
        plt.ylabel(terms[idxF[idx2]])
    plt.show()


    
def plotTerms(FA=None, S=None, alpha=None, terms=None, madFilter=.5):
    assert terms!=None
#        print 'terms need to be same length as relevance score'    
    if FA!=None:
        S = FA.S.E1
        alpha = FA.Alpha.E1
    MAD = mad(FA.S.E1)
    alpha = (MAD>madFilter)*(1/(alpha))

                 
    idx_sort = SP.argsort(terms)
    Y = alpha[idx_sort]
    X =SP.arange(len(alpha))#[idx_sort]
    plt.plot(X, Y, '.',markersize=10)
    plt.xticks(X, terms[idx_sort], size='small', rotation='vertical')    
    plt.ylabel("Relevance score")
    plt.show()
    


def getF(Y,FA=None, S=None, W=None, C=None, use_latent=False):
    #assert Y.shape[1] == FA.W.E1.shape[0] and Y.shape[0] == FA.W.E1.shape[0]
    if FA != None:
        S = FA.S.E1
        W = FA.W.E1
        C = FA.W.C[:,:,0]  

    isExpressed = (Y>0)*1.
    F = SP.dot(S,(C*W).transpose())
    F[isExpressed==0] = (F - 4.*(1./(1.+SP.exp(-F))))[isExpressed==0]  
    return F


def regressOut(Y,idx, FA=None, S=None, W=None, C=None, use_latent=False):
    #assert Y.shape[1] == FA.W.E1.shape[0] and Y.shape[0] == FA.W.E1.shape[0]
    if FA != None:
        S = FA.S.E1
        W = FA.W.E1
        C = FA.W.C[:,:,0]        
    idx = SP.array(idx)  
    isOn =  (C>.5)*1.0    
    if use_latent==False:
 
        Ycorr = Y-SP.dot(S[:,idx], (isOn[:,idx]*W[:,idx]).T)
    else:
        idx_use = SP.setxor1d(SP.arange(S.shape[1]),idx)
        Ycorr = SP.dot(S[:,idx_use], (isOn[:,idx_use]*W[:,idx_use]).T)    
    return Ycorr
    
    

def vcorrcoef(X,y):
    Xm = SP.reshape(SP.mean(X,axis=1),(X.shape[0],1))
    ym = SP.mean(y)
    r_num = SP.sum((X-Xm)*(y-ym),axis=1)
    r_den = SP.sqrt(SP.sum((X-Xm)**2,axis=1)*SP.sum((y-ym)**2))
    r = r_num/r_den
    return r

 
def getIlabel(order, Y, terms, pi,init_factors=None):
    assert (order in ['preTrain', 'PCA'])


    if order=='preTrain':
        assert init_factors!=None
        Ilabel = preTrain(Y, terms, pi,init_factors)
        return Ilabel
    else:
        PCs = SP.zeros((Y.shape[0], pi.shape[1]))
        for k in pi.shape[1]:
            pca = PCA(n_components=1)
            pca.fit_transform(Y[:,pi[:,k]>.5])
            PCs[:,k] = pca.score[:,0]

        X  = pca.fit_transform(Y)
        nFix = (SP.where(terms=='hidden')[0]).min()+len(SP.where(terms=='hidden')[0])
        MPC = abs(vcorrcoef(PCs.T,X.T))[nFix:]
        IpiRev = SP.argsort(MPC.ravel())      
        Ilabel = range(len(terms))
        Ilabel[nFix:] = IpiRev+nFix
        return Ilabel
        

def preTrain(Y, terms, pi00,init_factors, nFix=None, initType='pcaRand'):
    import core.fscLVM as fscLVM

    pi = pi00.copy()
    K = pi.shape[1]

    #data for sparseFA instance    
    pi[pi>.2] =0.99
    pi[pi<.8] =1e-8
    

    init={'init_data':fscLVM.CGauss(Y),'Pi':pi,'init_factors':init_factors}
    sigmaOff = 1E-3
    sparsity = 'VB'

    #prior on noise level 
    priors = {'Eps': {'priors':[1E-3,1E-3]}}
    #how to initialize network?
    #initType = 'pcaRand'
    terms0=terms
    pi0=pi.copy()
    FA0 = fscLVM.CSparseFA(components=K,sigmaOff=sigmaOff,sigmaOn=SP.ones(pi.shape[1])*1.0,sparsity=sparsity,nIterations=50,permutation_move=False,priors=priors,initType=initType)
    FA0.init(**init)
    if nFix==None:
        nFix = FA0.nKnown+FA0.nLatent
                        
#Fit PCA        
    pca = PCA(n_components=1)
    pca.fit(FA0.Z.E1)
    X = pca.transform(FA0.Z.E1)


#Sort by correlation to PC1    
    MPC = abs(vcorrcoef(FA0.initS[:,SP.argsort(FA0.W.Ilabel)].T,X.T))[nFix:]
    Ipi = SP.argsort(-MPC.ravel())
    IpiRev = SP.argsort(MPC.ravel())


    mRange = range(FA0.components)
    mRange[nFix:] = Ipi+nFix
    
    mRangeRev = range(FA0.components)
    mRangeRev[nFix:] = IpiRev+nFix

#Run model for 50 iterations         
    pi = pi0[:,mRange]
    terms = terms0[mRange]     
    init={'init_data':fscLVM.CGauss(Y),'Pi':pi,'init_factors':init_factors}
    FA = fscLVM.CSparseFA(components=K,sigmaOff=sigmaOff,sigmaOn=SP.ones(pi.shape[1])*1.0,sparsity=sparsity,nIterations=50,permutation_move=False,priors=priors,initType=initType)
    FA.shuffle=True
    FA.init(**init) 
    for j in range(50):
        FA.update()      
        

    #Run reverse model for 50 iterations         
    pi = pi0[:,mRangeRev]
    terms = terms0[mRangeRev]
    init={'init_data':fscLVM.CGauss(Y),'Pi':pi,'init_factors':init_factors}
    FArev = fscLVM.CSparseFA(components=K,sigmaOff=sigmaOff,sigmaOn=SP.ones(pi.shape[1])*1.0,sparsity=sparsity,nIterations=50,permutation_move=False,priors=priors,initType=initType)
    FArev.shuffle=True
    FArev.init(**init) 
    #FArev.iterate(forceIterations=True, nIterations=nIterations)
    for j in range(50):
        FArev.update() 
            
    
    
    IpiM = (0.5*FArev.Alpha.E1[SP.argsort(mRangeRev)][nFix:]+.5*FA.Alpha.E1[SP.argsort(mRange)][nFix:]).argsort()    
    Ilabel = SP.hstack([SP.arange(nFix),IpiM+nFix])

    return Ilabel
    
def load_data(dFile, annotation='MSigDB', minGenes=15, nHidden=3, doFast=True, FNR=0.001, noise='gauss'):
    data={}
    if doFast==False:
        out_dir = os.path.join(out_base,  dFile.split('.')[0],annotation)
    else:
        out_dir = os.path.join(out_base,  dFile.split('.')[0],annotation+'_fast')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    data['out_dir']=out_dir
    dataFile = h5py.File(os.path.join(data_dir, dFile), 'r')

    terms = dataFile['terms'][:]#[50:]
    pi = dataFile['Pi'][:].T#[:,50:]

    Y = dataFile['Yhet'][:].T


    #exclude small terms
    terms = terms[SP.sum(pi>.5,0)>minGenes]
    pi = pi[:,SP.sum(pi>.5,0)>minGenes]

    #Set FNR
    pi[pi<.1] =FNR

    #fast option?
    if doFast==True:
        idx_genes  = SP.logical_and(SP.sum(pi>.5,1)>0, Y.mean(0)>0.)#SP.any(pi>.5,1)
        Y = Y[:,idx_genes]
        pi = pi[idx_genes,:]


    #center data for Gaussian observation noise
    if noise=='gauss':
        Y-=SP.mean(Y,0)

    #include hidden variables
    terms = SP.hstack([SP.repeat('hidden',nHidden), terms])
    pi = SP.hstack([SP.ones((Y.shape[1],nHidden))*.99,pi])

    data['Y'] = Y
    data['pi'] = pi
    data['terms'] = terms

    return data

def addKnown(init_factors,dFile,data, idx_known=None):
    dataFile = h5py.File(os.path.join(data_dir, dFile), 'r')
    if len(idx_known)>0:
        known_names = dataFile['known_names'][:][idx_known]
        if len(dataFile['Known'][:].shape)>1:
            known = dataFile['Known'][:].T[:,idx_known]
        else:
            known = dataFile['Known'][:][:,SP.newaxis]
        #known -= known.mean(0)
        #known /= known.std(0)
        data['terms'] = SP.hstack([ known_names,data['terms']])
        pi = SP.hstack([SP.ones((data['Y'].shape[1],len(idx_known)))*.99,data['pi']])
        data['pi'] = pi
        init_factors['Known'] = known
        init_factors['iLatent'] = init_factors['iLatent'] + len(idx_known)

    return init_factors, data

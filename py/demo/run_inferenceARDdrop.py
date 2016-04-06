"""
run inference on tody data
"""
import sys
import os
import ast
sys.path.append('./..')
import scipy as SP
import pdb
import cPickle as pickle
#import core.sparseFAknown as sparseFA
from core.utils import *

import core.utils as utils
import core.sparseFAdropFast2 as sparseFA

import matplotlib as mpl
import matplotlib.lines as mlines
mpl.use('Agg')
import pylab as plt
import h5py
import brewer2mpl
from sklearn.decomposition import PCA
from sklearn.decomposition import RandomizedPCA
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector, StrVector
import rpy2.robjects as robjects
stats = importr('stats')
from rpy2.robjects.lib import ggplot2
grdevices = importr('grDevices')
cowplot = importr('cowplot')
gg2 = importr('ggplot2')


nthreads=15
os.environ['MKL_NUM_THREADS'] = str(nthreads)
os.environ['MKL_DOMAIN_NUM_THREADS'] = 'MKL_ALL='+str(nthreads)+', MKL_BLAS='+str(nthreads)
os.environ['MKL_DYNAMIC'] = 'true'


def mad(X):
    median = SP.median(X, axis=0)
    return SP.median(abs(X-median), axis=0)
    
def plotFac(FA, idx1, idx2, lab=None, terms=None, cols=None, isCont=False):
    X1 = FA.S.E1[:,idx1]
    X2 = FA.S.E1[:,idx2]
    if isCont==False:
        uLab = SP.unique(lab)  
        if cols==None:
             bmap=brewer2mpl.get_map('Paired', 'Qualitative', len(uLab))
             cols = bmap.hex_colors         
        pList=list()
        for i in range(len(X1)):
            pList.append(plt.plot(X1[i], X2[i], '.',color=cols[SP.where(lab[i]==uLab)[0]]))
        plt.xlabel(terms[idx1])
        plt.ylabel(terms[idx2])
        lList=list()
        for i in range(len(uLab)):
            lList.append( mlines.Line2D([], [], color=cols[i], marker='.',
                              markersize=7, label=uLab[i], linewidth=0))     
        plt.legend(handles=lList)
    else:
        plt.scatter(X1, X2, c=lab, s=20)
        plt.xlabel(terms[idx1])
        plt.ylabel(terms[idx2])        
    #plt.savefig('./Tcells_scatter.pdf')
    
def vcorrcoef(X,y):
    Xm = SP.reshape(SP.mean(X,axis=1),(X.shape[0],1))
    ym = SP.mean(y)
    r_num = SP.sum((X-Xm)*(y-ym),axis=1)
    r_den = SP.sqrt(SP.sum((X-Xm)**2,axis=1)*SP.sum((y-ym)**2))
    r = r_num/r_den
    return r
        
data_dir = '../../data/'
out_base = './results/'


if __name__ == '__main__':

    if 'cluster' in sys.argv:
        dFile = sys.argv[2]
        anno = sys.argv[3]
        nIterations = int(sys.argv[4])
        nHidden = int(sys.argv[5])
        idx_known = ast.literal_eval(sys.argv[6])
        doFast=bool(int(sys.argv[7]))
        if len(sys.argv)>7:
            idxCol = ast.literal_eval(sys.argv[8])
        else:
            idxCol=None
    else:
        
        dFile = 'retina.hdf5'
        anno = 'REACTOME'
        nHidden = 3
        idx_known = []#SP.arange(8,18)#SP.hstack([SP.arange(1,4),SP.arange(8,18)])
        nIterations = 6
        idxCol=[]
        doFast=True


 
            
    minGenes = 20
    dataFile = h5py.File(os.path.join(data_dir, dFile), 'r')
    if anno=='REACTOME':
        terms = dataFile['termsR'][:]  
        pi = dataFile['PiR20'][:].T
    else:
        terms = dataFile['terms'][:]#[50:]
        pi = dataFile['Pi'][:].T#[:,50:]
        
    pi[pi>.5] =0.99
    pi[pi<.5] =1e-08   
    
    if dFile=='Tcell.hdf5':
        Y = SP.log10(dataFile['Yhet'][:].T+1)
    elif dFile=='zeisel_microgliaR.hdf5':
        Y = SP.log2(dataFile['Yhet'][:].T+1)
    elif dFile=='retina.hdf5':
        YFile = h5py.File(os.path.join(data_dir, 'Yretina.hdf5'), 'r')
        Y = YFile['Y'][:]
    elif dFile=='retinaCounts.hdf5':
        YFile = h5py.File(os.path.join(data_dir, 'YretinaCounts.hdf5'), 'r')
        Y = YFile['Y'][:]
    else:
        Y = dataFile['Yhet'][:].T
    
    
    terms = terms[SP.sum(pi>.5,0)>minGenes]
    pi = pi[:,SP.sum(pi>.5,0)>minGenes]
    if 1:
        Yvar = Y.var(0)
        Yvar_idx = SP.argsort(-Yvar)
        Y = Y[:,Yvar_idx[:5000]]
        pi = pi[Yvar_idx[:5000],:]
     
    if doFast==True:
        idx_genes  = SP.logical_and(SP.sum(pi>.5,1)>0, Y.mean(0)>0.)#SP.any(pi>.5,1)
        Y = Y[:,idx_genes]
        pi = pi[idx_genes,:]

    terms = terms[SP.sum(pi>.5,0)>minGenes]
    pi = pi[:,SP.sum(pi>.5,0)>minGenes]        
    
    
     
    terms = SP.hstack([SP.repeat('hidden',nHidden), terms])
    pi = SP.hstack([SP.ones((Y.shape[1],nHidden))*.99,pi])

    init_factors = {}
    
    isBias=1
    if len(idx_known)>0:
        known_names = dataFile['known_names'][:][idx_known]
        if len(dataFile['Known'][:].shape)>1:
            known = dataFile['Known'][:].T[:,idx_known]
        else:
            known = dataFile['Known'][:][:,SP.newaxis]
        known -= known.mean(0)
        known = known+SP.random.random_sample(known.shape)/100
        known /= known.std(0)
        terms = SP.hstack([ known_names,terms])
        pi = SP.hstack([SP.ones((Y.shape[1],len(idx_known)))*.5,pi])
        init_factors['Known'] = known      
        if isBias==1:
            terms = SP.hstack([ 'bias',terms])
            pi = SP.hstack([SP.ones((Y.shape[1],1))*(1.-1e-10),pi])
            init_factors['Known'] =SP.hstack([SP.ones((Y.shape[0],1)), known])        
        
    else:
        known_names = '0'
        if isBias==1:
            terms = SP.hstack([ 'bias',terms])
            pi = SP.hstack([SP.ones((Y.shape[1],1))*(1.-1e-10),pi])
            init_factors['Known'] =SP.ones((Y.shape[0],1)) 
        #known_names = 'bias'
            
    if doFast==False:      
        out_dir = os.path.join(out_base,  dFile.split('.')[0],anno)
    else:
        out_dir = os.path.join(out_base,  dFile.split('.')[0],anno+'_FastScMax50EpsPC_pi15')

                
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    out_file = 'resHidden'+str(nHidden)+'99_Bias'+str(isBias)+'_Known__'+'__'.join(known_names)+'_Sort1e52k.hdf5'
    out_name = os.path.join(out_dir, out_file)
    print out_name  

    init_factors['iLatent'] = SP.where(terms=='hidden')[0]



    #Y-=SP.mean(Y,0)
    Y/=Y.max()
    nFix = len(SP.where(terms=='hidden')[0])+isBias
    K = pi.shape[1]
    IonPi = pi>.5
    varExpl = SP.zeros((K-nFix))
    for k in range(K-nFix):
        pca = RandomizedPCA()
        Xpc = pca.fit_transform(Y[:,IonPi[:,k+nFix]])[:,0:1]
        varExpl[k] = pca.explained_variance_ratio_[0]



#    Ipi = SP.argsort(varExpl.ravel())
    IpiRev = SP.argsort(-varExpl.ravel())
    mRangeRev = range(K)
    mRangeRev[nFix:] = (IpiRev+nFix)
    #mRangeRev = mRangeRev[:50]

    pi = pi[:,mRangeRev]
    terms = terms[mRangeRev]

    K = pi.shape[1]
    #data for sparseFA instance


    init={'init_data':sparseFA.CGauss(Y),'Pi':pi,'init_factors':init_factors}
    sigmaOff = 1E-3
    sparsity = 'VB'

    #permutation move
    permutation_move = False
    #prior on noise level
    priors = {'Eps': {'priors':[1E-3,1E-3]}}
    #how to initialize network?
    init={'init_data':sparseFA.CGauss(Y),'Pi':pi,'init_factors':init_factors}
    FA = sparseFA.CSparseFA(components=K,sigmaOff=sigmaOff,sigmaOn=SP.ones(pi.shape[1])*1.0,
                            sparsity=sparsity,nIterations=nIterations,
                            permutation_move=permutation_move,priors=priors,initType='pcaRand',minTerms=50, minIter=100)
    FA.shuffle=True
    FA.out_name = '.'.join(out_name.split('.')[0:2])
    FA.terms = terms

    FA.init(**init)
    FA.preTrain=1
    FA.doUpdate[50:] = 0
    for it in range(FA.nIterations):
        if it==100:
            FA.doUpdate[50:] = 1
        FA.update()
        FA.iterationCount+=1
             #calc reconstruction error
        if (SP.mod(FA.iterationCount,100)==0):
            alpha = FA.Alpha.E1/FA.S.E1.var(0)
            print SP.vstack([pi.sum(0)[FA.W.Ilabel][SP.argsort(alpha)],FA.W.C[:,:,0].sum(0)[SP.argsort(alpha)]])
            Zr = SP.dot(FA.S.E1[:,FA.doUpdate>=0],FA.W.E1.T[FA.doUpdate>=0,:])
            Zd = FA.Z.E1-Zr
            error = (Zd**2).mean()
            print "reconstruction error: %f" % (error)
            saveFA(FA)


=======
            
    
    Ion = FA.W.C[:,:,1]>.5
    NonInf = SP.sum(FA.W.C[:,:,0]>.5,0)
    IonInf = FA.W.C[:,:,0]>.5
    rel_contrib = SP.zeros(K)
    for k in range(K):
        rel_contrib[k] = NonInf[k]/(FA.Alpha.E1[k]*(Y[:,IonInf[:,k]].var(0)-1./FA.Eps.E1[IonInf[:,k]]).sum())
        #rel_contrib[k] = NonInf[k]/(FA.Alpha.E1[k]*(Y[:,IonInf[:,k]].var(0)).sum())

    alpha02 = (MAD>.05)*(1/(FA.Alpha.E1))

    print alpha02
    print SP.vstack([pi.sum(0)[FA.W.Ilabel][SP.argsort(alpha)],FA.W.C[:,:,0].sum(0)[SP.argsort(alpha)]])


    idxF = SP.argsort(-alpha)#[0:10]
    idxF0 = SP.argsort(-alpha0)#[0:10]
    idxF2 = SP.argsort(-alpha2)#[0:10]
    idxF02 = SP.argsort(-alpha02)#[0:10]

    plt_idx1 = 0
    plt_idx2 = 1

    #plot
    dataf = robjects.DataFrame({'alpha': FloatVector(alpha02),'terms':terms})
    if type(idxCol) == list and len(idxCol)==1: idxCol=idxCol[0]

    if SP.isscalar(idxCol):
        cc = dataFile['Known'][:].T[:,idxCol]
        datafs = robjects.DataFrame({'col': FloatVector(cc), 'X1': FloatVector(FA.S.E1[:,idxF02[plt_idx1]]),
                            'X2': FloatVector(FA.S.E1[:,idxF02[plt_idx2]])})
        colName = dataFile['known_names'][:][idxCol]
    elif type(idxCol) == list and len(idxCol)==2:
        cc = dataFile['Known'][:].T[:,idxCol[0]]+2*dataFile['Known'][:].T[:,idxCol[1]]
        S2 = FA.S.E1[:,idxF02[plt_idx2]]
        S1 = FA.S.E1[:,idxF02[plt_idx1]]
        datafs = robjects.DataFrame({'col': StrVector.factor(StrVector(cc)),
                                     'X1': FloatVector(S1), 'X2': FloatVector(S2)})
        colName = "Cell Cycle"
    elif type(idxCol) == list and len(idxCol)>2:
        cc = dataFile['Known'][:].T[:,idxCol[0]]
        for i in range(len(idxCol)):
            cc = cc+(i+1)*dataFile['Known'][:].T[:,idxCol[i]]
        S2 = FA.S.E1[:,idxF02[plt_idx2]]
        S1 = FA.S.E1[:,idxF02[plt_idx1]]
        datafs = robjects.DataFrame({'col': StrVector.factor(StrVector(cc)),
                                     'X1': FloatVector(S1), 'X2': FloatVector(S2)})
        colName = "Cell Type"
    else:
        S2 = FA.S.E1[:,idxF02[plt_idx2]]
        S1 = FA.S.E1[:,idxF02[plt_idx1]]
        datafs = robjects.DataFrame({'X1': FloatVector(S1), 'X2': FloatVector(S2)})
        colName=None



    p = ggplot2.ggplot(dataf)+ggplot2.aes_string(x='terms', y='alpha')+ggplot2.geom_point()+ \
        ggplot2.theme(**{'axis.title.y': ggplot2.element_text(angle=90,size = gg2.rel(1.2)),
        'axis.title.x': ggplot2.element_text(size = gg2.rel(1.2)),
        'axis.text.x': ggplot2.element_text(size=6, angle=90)})+\
        ggplot2.labs(x="Processes", y='Relevance')#+ theme(axis.title.y = 'element_text(size = rel(1.2))')+theme(axis.title.x = element_text(size = rel(1.2)))
    p.save(out_name+'_relevance.pdf', height=4, width=7)

       
    labels=str(tuple(['G1', 'S', 'G2M']))
    if colName=="Cell Cycle":
        pScatter = ggplot2.ggplot(datafs)+\
            ggplot2.aes_string(x='X1', y='X2', colour="factor(col, labels=c%s)" % labels) +ggplot2.geom_point()+\
            ggplot2.theme(**{'axis.title.y': ggplot2.element_text(angle=90,size = gg2.rel(.9)),
            'axis.title.x': ggplot2.element_text(size = gg2.rel(.9))})+\
            ggplot2.labs(x=terms[idxF02[plt_idx1]], y=terms[idxF02[plt_idx2]])+\
            ggplot2.scale_colour_manual(name=colName, values=StrVector(['#1b9e77', '#d95f02', '#7570b3']))
    elif colName=="Cell Type":
        labels=str(tuple(dataFile['known_names'][:][idxCol]))
        pScatter = ggplot2.ggplot(datafs)+ggplot2.aes_string(x='X1', y='X2', colour="factor(col, labels=c%s)" % labels)+ggplot2.geom_point()+\
            ggplot2.theme(**{'axis.title.y': ggplot2.element_text(angle=90,size = gg2.rel(.8)),
            'axis.title.x': ggplot2.element_text(size = gg2.rel(.8)),
            'axis.text.x': ggplot2.element_text(size = gg2.rel(.9))})+\
            ggplot2.labs(x=terms[idxF02[plt_idx1]], y=terms[idxF02[plt_idx2]])+\
            ggplot2.scale_colour_discrete(name =colName)

    elif type(colName)!=None and colName!=None:
        pScatter = ggplot2.ggplot(datafs)+ggplot2.aes_string(x='X1', y='X2', colour='col')+ggplot2.geom_point(alpha=0.4)+\
            ggplot2.theme(**{'axis.title.y': ggplot2.element_text(angle=90,size = gg2.rel(.8)),
            'axis.title.x': ggplot2.element_text(size = gg2.rel(.8)),
            'axis.text.x': ggplot2.element_text(size = gg2.rel(.9))})+\
            ggplot2.labs(x=terms[idxF02[plt_idx1]], y=terms[idxF02[plt_idx2]])+\
            ggplot2.scale_colour_continuous(name =colName)
    else:
         pScatter = ggplot2.ggplot(datafs)+ggplot2.aes_string(x='X1', y='X2')+ggplot2.geom_point()+\
            ggplot2.theme(**{'axis.title.y': ggplot2.element_text(angle=90,size = gg2.rel(.8)),
            'axis.title.x': ggplot2.element_text(size = gg2.rel(.8)),
            'axis.text.x': ggplot2.element_text(size = gg2.rel(.9))})+\
            ggplot2.labs(x=terms[idxF02[plt_idx1]], y=terms[idxF02[plt_idx2]])

    pScatter.save(out_name+'_'+terms[idxF02[plt_idx1]]+'_'+terms[idxF02[plt_idx2]]+'_scatter.pdf', height=4, width=7)

    out_file = h5py.File(out_name+'.hdf5','w')
    out_file['alphaRaw'] = FA.Alpha.E1
    out_file['alpha'] = alpha
    out_file['alpha2'] = alpha2
    out_file['alpha02'] = alpha02
    out_file['alpha0'] = alpha0
    out_file['W'] = FA.W.E1
    out_file['Eps'] = FA.Eps.E1
    out_file['S'] = FA.S.E1
    out_file['Gamma'] = FA.W.C[:,:,0]
    out_file['pi'] = pi
    out_file['terms'] = terms
    out_file.close()
    #pickle.dump(FA, open(out_name+'_FA.pickle', 'w'))

    


  
>>>>>>> 564e8fac635e2a87da2c01e8881e2a3f8a474fd8





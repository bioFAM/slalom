"""
run inference on tody data
"""
import sys
import os
import ast
sys.path.append('./..')
import scipy as SP
import cPickle as pickle
import core.fscLVM as fscLVM
from core.utils import *
import h5py
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
        
        dFile = 'Buettneretal_sfERCC.hdf5'
        anno = 'MSigDBtest'
        nHidden = 3
        idx_known = []
        nIterations = 100
        idxCol=[0,1]
        doFast=True


 
            
    minGenes = 15
    dataFile = h5py.File(os.path.join(data_dir, dFile), 'r')
    if anno=='REACTOME':
        terms = dataFile['termsR'][:]  
        pi = dataFile['PiR20'][:].T
    else:
        terms = dataFile['terms'][:]#[50:]
        pi = dataFile['Pi'][:].T#[:,50:]
        
   
    Y = dataFile['Yhet'][:].T
    if dFile=='Tcell.hdf5':
        Y = SP.log10(dataFile['Yhet'][:].T+1)
    elif dFile=='zeisel_microgliaR.hdf5':
        Y = SP.log2(dataFile['Yhet'][:].T+1)
    
    Y-=SP.mean(Y,0)
    
    terms = terms[SP.sum(pi>.5,0)>minGenes]
    pi = pi[:,SP.sum(pi>.5,0)>minGenes]
    
    if doFast==True:
        idx_genes  = SP.logical_and(SP.sum(pi>.5,1)>0, Y.mean(0)>0.)#SP.any(pi>.5,1)
        Y = Y[:,idx_genes]
        pi = pi[idx_genes,:]
    
     
    terms = SP.hstack([SP.repeat('hidden',nHidden), terms])
    pi = SP.hstack([SP.ones((Y.shape[1],nHidden))*.99,pi])

    init_factors = {}
    

    if len(idx_known)>0:
        known_names = dataFile['known_names'][:][idx_known]
        if len(dataFile['Known'][:].shape)>1:
            known = dataFile['Known'][:].T[:,idx_known]
        else:
            known = dataFile['Known'][:][:,SP.newaxis]
        known -= known.mean(0)
        known /= known.std(0)
        terms = SP.hstack([ known_names,terms])
        pi = SP.hstack([SP.ones((Y.shape[1],len(idx_known)))*.5,pi])
        init_factors['Known'] = known             
    else:
        known_names = '0'

            
    if doFast==False:      
        out_dir = os.path.join(out_base,  dFile.split('.')[0],anno)
    else:
        out_dir = os.path.join(out_base,  dFile.split('.')[0],anno+'_fast')
                
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    out_file = 'resHidden'+str(nHidden)+'_Known__'+'__'.join(known_names)+'.hdf5'
    out_name = os.path.join(out_dir, out_file)
    print out_name  

    terms0=terms.copy()
    pi0=pi.copy()
    init_factors['iLatent'] = SP.where(terms=='hidden')[0]

    
    Ilabel = preTrain(Y, terms, pi,init_factors) 
    print terms0[Ilabel]

    
    #Ilabel = mRangeRev
    pi = pi0[:,Ilabel]
    #Set FNR
    pi[pi<.1] =1e-3
    terms = terms0[Ilabel] 
    init={'init_data':fscLVM.CGauss(Y),'Pi':pi,'init_factors':init_factors}
    priors = {'Eps': {'priors':[1E-3,1E-3]}}
    FA = fscLVM.CSparseFA(components=pi.shape[1], nIterations=nIterations)            
    FA.shuffle=True
    FA.init(**init) 
    FA.iterate(forceIterations=True, nIterations=nIterations)
        
    MAD = mad(FA.S.E1)
    alpha02 = (MAD>.5)*(1/(FA.Alpha.E1))
    idxF02 = SP.argsort(-alpha02)
    
    
#plot 
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

#save    
    out_file = h5py.File(out_name+'.hdf5','w')    
    out_file['alphaRaw'] = FA.Alpha.E1
    out_file['alpha02'] = alpha02  
    out_file['MAD'] = MAD  
    out_file['W'] = FA.W.E1
    out_file['Eps'] = FA.Eps.E1
    out_file['S'] = FA.S.E1        
    out_file['Gamma'] = FA.W.C[:,:,0]
    out_file['pi'] = pi
    out_file['terms'] = terms
    out_file.close()
    pickle.dump(FA, open(out_name+'_FA.pickle', 'w'))
    


  





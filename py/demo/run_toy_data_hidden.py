"""create_toy_data
- create toy data for factor analysis infernece and store as CSV files
"""

import sys
sys.path.append('./..')
#scLVM_BASE = '../../../../scLVM/'
scLVM_BASE = '../../../scLVM/'
sys.path.insert(1,scLVM_BASE)
sys.path.insert(2, scLVM_BASE +'..')
import scipy as SP
import numpy.random as random
import core.sparseFAvem as sparseFA
import scipy.io
import os
import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
import pdb
import sys
import os
import h5py
import cPickle as pickle
import glob
from sklearn.metrics import roc_curve, auc
from sklearn.metrics.pairwise import pairwise_distances

python_cmd = 'python'
# settings
nthreads = 1
sleep = 0
mem_thread = 600
basename_out = './outSimRandHidden2FPTP'

cluster_cmd = 'bsub -n %d -q research-rh6 -R "rusage[mem=%d]" -M %d -o ./cluster_out' % (nthreads,mem_thread,mem_thread)



if __name__ =='__main__':
    if 'cluster' in sys.argv:
        nSim = 100        
        N = int(sys.argv[2])
        G = int(sys.argv[3])


        for iSim in range(nSim):
            out_dir = os.path.join(basename_out,'Sim_'+str(iSim))
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            in_file = os.path.join(out_dir,'data.h5py')
            #simulate data    
            #sparcity
            sparse = SP.hstack([SP.arange(0.1,0.55,0.1),0.6])
            K = len(sparse)
        
            #noise level
            eps = 0.1        
        #    1. create random netowrk prior      
            Ion = SP.zeros((G,3*K))>0
            #pi  = SP.zeros([G,K])
            for k in range(K):
                onk_ = random.rand(G,)<sparse[k]
                Ion[:,k] = onk_
        
            for k in range(K):
                onk_ = random.rand(G,)<sparse[k]
                Ion[:,k+K] = onk_

            for k in range(K):
                onk_ = random.rand(G,)<sparse[k]
                Ion[:,k+2*K] = onk_
        
            #onany_ = SP.where(SP.sum(Ion[:,:K],1)>0)[0]        
            #for k in SP.arange(K):
            #    idxon_ = random.choice(onany_, round(sparse[k]*G), replace=False)
            #    Ion[idxon_,k+2*K] = True
        
            Ion = (Ion>0)
            Ioff = ~Ion
        
            #2. create pi matrix with prior probability of a link
            Pi  = SP.zeros(Ion.shape)
            Pi[Ion] = 1.-1E-2
            Pi[Ioff] = 1E-5
            for k in range(Pi.shape[1]):
		TP = SP.where(Pi[:,k]>.5)[0]
		TN = SP.where(Pi[:,k]<.5)[0]
                FP = SP.random.choice(TN, size = min([round(len(TN)*.01), 20]))
                FN = SP.random.choice(TP, size = min([round(len(TP)*.01), 20]))
                Pik = Pi[:,k]
                Pik[FP] = 1.-1E-2
                Pik[FN] = 1E-5
                Pi[:,k] = Pik 
            #3. sample form pi true network    
                
            Z = 1.0*Ion[:,(2*K):]        
            #4. create weight matrix
            W = random.randn(G,K)
            W *= Z
        
            #5. create factor matrix
            X = random.randn(N,K)
            
            Khidden = 2
            IonHidden = SP.zeros((G,Khidden))>0
            sparseHidden = SP.arange(0.2,0.8,0.05)
            #pi  = SP.zeros([G,K])
            for k in range(Khidden):
                onk_ = random.rand(G,)<sparseHidden[SP.random.randint(0,len(sparseHidden))]
                IonHidden[:,k] = onk_
            Zhidden = 1.0*IonHidden        
            #4. create weight matrix
            Whidden = random.randn(G,Khidden)
            Whidden *= Zhidden
        
            #5. create factor matrix
            Xhidden = random.randn(N,Khidden)            
            
            Pi = SP.hstack([SP.zeros((G,Khidden+1))+.51,Pi])
            IonTerm = SP.arange(2*K,3*K)+Khidden+1
            
            
            Y = SP.dot(SP.hstack([X, Xhidden]),SP.hstack([W, Whidden]).T) + eps*random.randn(N,G)
            #data for sparseFA instance        
            Y-=SP.mean(Y,0)
            Y/=SP.std(Y,0)
            
            data_file = h5py.File(in_file, 'w')
            data_file['Y'] = Y
            data_file['Pi'] = Pi
            data_file['X'] = X
            data_file['W'] = W
            data_file['FN'] = FN
            data_file['FP'] = FP
            data_file['TN'] = TN
            data_file['TP'] = TP
            data_file['Ion'] = Ion
            data_file['Xhidden'] = Xhidden
            data_file['WHidden'] = Whidden
            data_file['IonHidden'] = IonHidden            
            data_file['Nhidden'] = Khidden+1
            data_file['IonTerm'] = IonTerm
            data_file.close()

            #cmd = '%s %s %s %s %s' % (python_cmd, sys.argv[0],'scLVM',in_file, out_dir)            
            cmd = '%s %s %s %s %s %s' % (cluster_cmd, python_cmd, sys.argv[0],'scLVM',in_file, out_dir)
            print cmd            
            #os.system(cmd)
            
            #cmd = '%s %s %s %s %s' % (python_cmd, sys.argv[0],'sscLVM',in_file,out_dir)
            cmd = '%s %s %s %s %s %s' % (cluster_cmd, python_cmd, sys.argv[0],'sscLVM',in_file,out_dir)
            print cmd
            #cmd = '%s %s %s %s' % (python_cmd, sys.argv[0],fn,out_file)
            os.system(cmd)

        
        
        
    elif 'collect' in sys.argv:
        dir_name_out = sys.argv[2]
        DL = glob.glob(dir_name_out)
        
        tprDict = {}
        fprDict = {}
        aucDict = {}	
        metrDict = {}	
        isOnList = []
        isWOnList = []
        isOnList = []
        corrList = []
        corrPCList = []
        metrics = ['varCompMean','var','cv2','alpha','rel_contrib']     
        for metric in metrics:
            fprDict[metric] = []
            tprDict[metric] = []
            aucDict[metric] = []
            metrDict[metric] = []

        for root, dirs, files in os.walk(dir_name_out):
            for dir_i in dirs:
                dfile = h5py.File(os.path.join(dir_name_out,dir_i, 'data.h5py'),'r')
                scLVMfile = h5py.File(os.path.join(dir_name_out,dir_i, 'out_scLVM.h5py'),'r')
                sscLVMfile = h5py.File(os.path.join(dir_name_out,dir_i, 'out_sscLVM.h5py'),'r')
                K = dfile['X'][:].shape[1]
                Nhidden = dfile['Nhidden'][()]
                idxOn = dfile['IonTerm'][:]-Nhidden
                    
                cnt = 3
                for metric in metrics:
                    if cnt<3:
                        metr = scLVMfile[metric][:]
                    else:
                        metr = sscLVMfile[metric][:][Nhidden:]
                    if metric=='alpha':
                        metr=1.0/metr
                             
                    isOn = SP.zeros((len(metr)))
                    isOn[idxOn] = 1.0
                    if cnt==3:
                        isOnList.append(isOn)                        
                        _corr = []
                        _corrPC = []
                        for i in range(len(idxOn)):
                            #_idx+=Nhidden
                            _corr.append(SP.corrcoef(sscLVMfile['S'][:][:,idxOn[i]+Nhidden], dfile['X'][:][:,i])[0,1])
                            _corrPC.append(SP.corrcoef(scLVMfile['Xlist'][:][idxOn[i]].ravel(), dfile['X'][:][:,i])[0,1])
                        _corrHidden = abs(1-pairwise_distances( sscLVMfile['S'][:][:,:Nhidden].T,dfile['Xhidden'][:].T, metric='correlation'))
                        corrList.append(SP.hstack([SP.diag(_corrHidden[SP.argmax(_corrHidden,0)]),_corr]))
                        corrPCList.append(SP.array(_corrPC))
                        
                    fpr, tpr, _ = roc_curve(isOn,metr)
                    fprDict[metric].append(fpr)
                    tprDict[metric].append(tpr)
                    metrDict[metric].append(metr)
                    aucDict[metric].append(auc(fpr, tpr))
                    cnt+=1

        aucDictAll = {} 
        fprDictAll = {} 
        tprDictAll = {} 
        for metric in metrics:
            fpr, tpr, _ = roc_curve(SP.hstack(isOnList),SP.hstack(metrDict[metric]))            
            aucDictAll[metric] = auc(fpr, tpr)
            tprDictAll[metric] = tpr 
            fprDictAll[metric] = fpr
        leg_str = list()
        for metric in metrics: 
            plt.plot(fprDictAll[metric], tprDictAll[metric], '-')
            leg_str.append(metric+', AUC = '+str(SP.round_(aucDictAll[metric],3)))
        plt.xlabel('FPR')
        plt.ylabel('TPR')
        plt.legend(leg_str, loc='lower-right')
        ax = plt.gca()
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        plt.savefig('./ROC_tes.pdf',bbox_inches='tight') 

        
        
    elif 'scLVM' in sys.argv:
        data_file = h5py.File(sys.argv[2], 'r')
        out_file = h5py.File(os.path.join(sys.argv[3],'out_scLVM.h5py'),'w')
        Y = data_file['Y'][:]
        Pi = data_file['Pi'][:]
        Nhidden = data_file['Nhidden'][()]
        from sklearn.decomposition import RandomizedPCA
        from scLVM import scLVM
     
        IonPi = (Pi>.5)[:,Nhidden:]
        K = IonPi.shape[1]
        Klist = list()
        Xlist = list()
        cv2 = SP.zeros((K))
        var = SP.zeros((K))
        for k in range(K):
            pca = RandomizedPCA(n_components=1)
            Xpc = pca.fit_transform(Y[:,IonPi[:,k]])                 
            Xlist.append(Xpc)
            Klist.append(SP.dot(Xpc, Xpc.T))
            cv2[k] = SP.mean((SP.std(Y[:,IonPi[:,k]],1)/SP.mean(Y[:,IonPi[:,k]],1))**2)
            var[k] = SP.mean(SP.var(Y[:,IonPi[:,k]],1))
            
        sclvm = scLVM(Y,tech_noise=SP.ones((Y.shape[1]))*1e-20)
        sclvm.varianceDecomposition(K=Klist,i0=0,i1=Y.shape[1], maxiter=2)
        varComp = sclvm.getVarianceComponents()[0][:,:K]
        varCompMean = SP.mean(varComp,0)
        
        out_file['varComp'] = varComp
        out_file['varCompMean'] = varCompMean
        out_file['var'] = var
        out_file['cv2'] = cv2
        out_file['Xlist'] = Xlist
        out_file.close()
        
        
        
        
    elif 'sscLVM' in sys.argv:
        out_dir = sys.argv[3]
        data_file = h5py.File(sys.argv[2], 'r')
        out_file = h5py.File(os.path.join(sys.argv[3],'out_sscLVM.h5py'),'w')
        Y = data_file['Y'][:]
        Pi = data_file['Pi'][:]
         #number of samples  
        init={'init_data':sparseFA.CGauss(Y),'Pi':Pi}
        sigmaOff = 1E-5
        sparsity = 'VB'
        nIterations = 5000
        #permutation move
        permutation_move = False
        #prior on noise level 
        priors = {'Eps': {'priors':[1,10]}}
        #how to initialize network?
        initType = 'pcaRand'# 'pca'
    
        K = Pi.shape[1]
        
        FAtrue = sparseFA.CSparseFA(components=K,sigmaOff=sigmaOff,\
        sigmaOn=SP.ones(Pi.shape[1])*1.0,sparsity=sparsity,nIterations=nIterations,permutation_move=permutation_move,priors=priors,initType=initType)
        FAtrue.init(**init)
    
        FAtrue.iterate(forceIterations=True)
        rel_contrib = SP.zeros(K)
        NonInf = SP.sum(FAtrue.W.C[:,:,0]>.5,0)
        IonInf = FAtrue.W.C[:,:,0]>.5
        for k in range(K):
            rel_contrib[k] = NonInf[k]/(FAtrue.Alpha.E1[k]*(Y[:,IonInf[:,k]].var(0)-1./FAtrue.Eps.E1[IonInf[:,k]]).sum())    
            
            
        pickle.dump(FAtrue, open(os.path.join(out_dir, 'FA.pickle'),'w'))
        out_file['rel_contrib'] = rel_contrib     
        out_file['alpha'] = FAtrue.Alpha.E1
        out_file['W'] = FAtrue.W.E1
        out_file['Eps'] = FAtrue.Eps.E1
        out_file['S'] = FAtrue.S.E1        
        out_file['Pi'] = FAtrue.W.C[:,:,0]
        out_file.close()
        
        
        
        plt.plot(SP.arange(len(FAtrue.Alpha.E1)), 1/FAtrue.Alpha.E1,'*')
        plt.savefig(os.path.join(out_dir, 'rel_plot.pdf'))

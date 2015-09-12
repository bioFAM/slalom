"""create_toy_data
- create toy data for factor analysis infernece and store as CSV files
"""

import sys
sys.path.append('./..')
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

python_cmd = 'python'
# settings
nthreads = 1
sleep = 0
mem_thread = 1000
basename_out = './outSim'

cluster_cmd = 'bsub -n %d -q research-rh6 -R "rusage[mem=%d]" -M %d -o ./cluster_out' % (nthreads,mem_thread,mem_thread)



if __name__ =='__main__':
    if 'cluster' in sys.argv:
        nSim = 200        
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
        
            onany_ = SP.where(SP.sum(Ion[:,:K],1)>0)[0]        
            for k in SP.arange(K):
                idxon_ = random.choice(onany_, round(sparse[k]*G), replace=False)
                Ion[idxon_,k+2*K] = True
        
            Ion = (Ion>0)
            Ioff = ~Ion
        
            #2. create pi matrix with prior probability of a link
            Pi  = SP.zeros(Ion.shape)
            Pi[Ion] = 1.-1E-2
            Pi[Ioff] = 1E-5

            #3. sample form pi true network    
            IonTerm = SP.arange(2*K,3*K)    
            Z = 1.0*Ion[:,(2*K):]        
            #4. create weight matrix
            W = random.randn(G,K)
            W *= Z
        
            #5. create factor matrix
            X = random.randn(N,K)
            Y = SP.dot(X,W.T) + eps*random.randn(N,G)
            #data for sparseFA instance        
            Y-=SP.mean(Y,0)
            Y/=SP.std(Y,0)
            
            data_file = h5py.File(in_file, 'w')
            data_file['Y'] = Y
            data_file['Pi'] = Pi
            data_file['X'] = X
            data_file['W'] = W
            data_file['Ion'] = Ion
            data_file['IonTerm'] = IonTerm
            data_file.close()

            #cmd = '%s %s %s %s %s' % (python_cmd, sys.argv[0],'scLVM',in_file, out_dir)            
            cmd = '%s %s %s %s %s %s' % (cluster_cmd, python_cmd, sys.argv[0],'scLVM',in_file, out_dir)
            print cmd            
            os.system(cmd)
            
            #cmd = '%s %s %s %s %s' % (python_cmd, sys.argv[0],'sscLVM',in_file,out_dir)
            cmd = '%s %s %s %s %s %s' % (cluster_cmd, python_cmd, sys.argv[0],'sscLVM',in_file,out_dir)
            print cmd
            #cmd = '%s %s %s %s' % (python_cmd, sys.argv[0],fn,out_file)
            os.system(cmd)

        
        
        
    elif 'collect' in sys.argv:
        pass
        
        
    elif 'scLVM' in sys.argv:
        data_file = h5py.File(sys.argv[2], 'r')
        out_file = h5py.File(os.path.join(sys.argv[3],'out_scLVM.h5py'),'w')
        Y = data_file['Y'][:]
        Pi = data_file['Pi'][:]
        from sklearn.decomposition import RandomizedPCA
        from scLVM import scLVM
     
        IonPi = Pi>.5
        K = Pi.shape[1]
        Klist = list()
        cv2 = SP.zeros((K))
        var = SP.zeros((K))
        for k in range(K):
            pca = RandomizedPCA(n_components=1)
            Xpc = pca.fit_transform(Y[:,IonPi[:,k]])                 
            Klist.append(SP.dot(Xpc, Xpc.T))
            cv2[k] = SP.mean((SP.std(Y[:,IonPi[:,k]],1)/SP.mean(Y[:,IonPi[:,k]],1))**2)
            var[k] = SP.mean(SP.var(Y[:,IonPi[:,k]],1))
            
        sclvm = scLVM(Y,tech_noise=SP.ones((Y.shape[1]))*1e-20)
        sclvm.varianceDecomposition(K=Klist,i0=0,i1=Y.shape[1])
        varComp = sclvm.getVarianceComponents()[0][:,:K]
        varCompMean = SP.mean(varComp,0)
        
        out_file['varComp'] = varComp
        out_file['varCompMean'] = varCompMean
        out_file['var'] = var
        out_file['cv2'] = cv2
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

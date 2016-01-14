 # f-scLVM
# variational sparse factor analysis model with spike and slab prior


from vbfa import *
import scipy as SP
#from mlib.utils import *
#from mlib.plot import *
import copy
import pdb
from sklearn.decomposition import PCA

def mad(X):
    median = SP.median(X, axis=0)
    return SP.median(abs(X-median), axis=0)
    
def vcorrcoef(X,y):
    Xm = SP.reshape(SP.mean(X,axis=1),(X.shape[0],1))
    ym = SP.mean(y)
    r_num = SP.sum((X-Xm)*(y-ym),axis=1)
    r_den = SP.sqrt(SP.sum((X-Xm)**2,axis=1)*SP.sum((y-ym)**2))
    r = r_num/r_den
    return r


class CNodeAlphasparse(AGammaNode):
    def __init__(self,net,prior=[1E-3,1E-3]):
        AGammaNode.__init__(self,dim=[net.components],prior=prior)
        


class CNodeEpsSparse(CNodeEps):
    """Extensions of CNodeEps that can handle fixed prior expectaitons"""
    def __init__(self,net,prior=S.array([100,1])):
            CNodeEps.__init__(self,net,prior)
                
        
        

class CNodeSsparse(AVGaussNode):
    def __init__(self,net,prec=1):
#        CNodeS.__init__(self,net,prec=1)
        AVGaussNode.__init__(self,dim=[net._N,net.components],cdim=1)
        self.diagSigmaS = SP.ones((net._N,net.components))

        

class CNodeWsparse(CNodeW):
    """Abstract CnodeWsparse basclass
    - this ensures consitent handling of alternative CNodeW nodes; in particular the sparsity prior
    - main application : permutation move of factors"""
    def __init__(self,net,**kwargin):
        #call base class initialisation
        #CNodeW.__init__(self,net,**kwargin)
        #create indicator arrays for moment calculation
      
        self.Pi = zeros([net.Pi.shape[0],net.Pi.shape[1],2])
        self.Pi[:,:,0] =net.Pi
        self.Pi[:,:,1] = 1.-net.Pi


        self.C    = self.Pi.copy()
        #self.C[self.Pi>.8] = .9
        #self.C[self.Pi<.1]= .1  
        #labeling of factors in case we use permutation move
        self.Ilabel = SP.arange(net.components)
        


class CNodeWsparseVB(CNodeWsparse):
    def __init__(self,net,prec=1):
        CNodeWsparse.__init__(self,net,prec=prec)
        #variable initialisation in CNodeWsparse
        self.sigma2 = (1.0/prec)*SP.ones((net._D, net.components))
        self.E1 = SP.randn(net._D, net.components)
        self.E2diag = SP.zeros((net._D, net.components))#is calculated in update for W
        pass


class CSparseFA(AExpressionModule):
    '''CVBFA(AExpressionModule)
    - Variational Bayesian Factor analysis module'''


    def getDefaultParameters(self):
        """return a hash with default parameter value for this BayesNet"""
        #
        dp = AExpressionModule.getDefaultParameters(self)
        dp['initType'] = 'pcaRand'
        dp['nIterations'] = 2000
        #VB best
        dp['schedule'] = ['W','S','Eps']
#        #experiments(EP)
        dp['schedule'] = ['Eps','S','W']
        dp['schedule'] = ['W','S','Eps','Alpha']
        dp['permutation_move'] = False
        dp['shuffle']=True
        dp['components'] = 5
        dp['priors']     = {}
        dp['sigmaOff']   = 1E-3
        dp['sigmaOn']    = S.ones(dp['components'])*1.0
        dp['sparsity']   = 'VB'
        dp['iterative']  = True
        dp['randperm']   = False
        return dp

    def getName(self,base_name='sparseFA'):
        """return a name summarising the  main parameters"""
        print "Not implemented yet"


    def iterate(self, nIterations=None, forceIterations=None, tolerance=None, minIterations=10):
        '''iterate(nIteations=None,forceIterations=None)
        - perform nIterations; per default(None) parameters are tken from local intsance settings
        '''
        #forceIterations=True
        L.debug('fscLVM iterate')                
        
        iter=0
        if tolerance is None: tolerance = self.tolerance
        if nIterations is None: nIterations = self.nIterations
        if forceIterations is None: forceIterations = self.forceIterations
        Ion = (self.W.C[:,:,0]>.5)*1.
        Zr = S.dot(self.S.E1,(self.W.E1.T*Ion.T))
        Zd = self.Z.E1-Zr
        error = (Zd**2).mean()
        self.calcBound()
        if SP.mod(iter,1)==0:
            print "reconstruction error: %f lower bound: %f" % (error,  self._bound)
        errorOld = error

        for iter in range(nIterations):
            #self.iterationCount+=1
            t = time.time();
            #for node in self.schedule:
            #    self.updateNode(node)
                #pdb.set_trace()
            self.update()
            if self.permutation_move==True:
                #greedily select factors
                print "Not implemented for this model"

              
            self.iterationCount+=1
            #calc reconstruction error
            Zr = S.dot(self.S.E1,self.W.E1.T)
            Zd = self.Z.E1-Zr
            error = (Zd**2).mean()

            if (SP.mod(iter,100)==0):
                print "reconstruction error: %f lower bound: %f" % (error,  self._bound) 

                if (abs(errorOld - error) < tolerance) and not forceIterations and iter>minIterations:
                    print 'Converged after %i iterations' % (iter)
                    break

            L.info("Iteration %d: time=%.2f bound=%f" % (iter,time.time() - t, self._bound))
            errorOld = error

        return error

        
        #1. itearte
        AExpressionModule.iterate(self,*argin,**kwargin)
        #2. debug
        pass
    
    def logProbkk(self,k,l):
        """evaluate the probability of C_k and Pi_l"""
        pp = self.W.C[:,k,:]*self.W.Pi[:,l,:]
        lpp = SP.log(pp.sum(axis=1))
        return lpp.sum()
        
    def updateS(self,m):
        M = self.components
        if m>=self.nKnown:
            setMinus = SP.int_(SP.hstack([range(M)[0:m],range(M)[m+1::]]))

            #update S
            SW_sigma = (self.W.C[:, m,0]*self.W.E1[:, m])*self.Eps.E1
            SW2_sigma = (self.W.C[:, m,0]*(self.W.E2diag[:, m]))*self.Eps.E1              
            setMinus = SP.int_(SP.hstack([range(M)[0:m],range(M)[m+1::]]))
             
            b0 = SP.dot(self.S.E1[:,setMinus],(self.W.C[:, setMinus,0]*self.W.E1[:, setMinus]).transpose())
            b=SP.dot(b0,SW_sigma)
            #alphaSm = SP.dot(SP.ones(self.Z.E1.shape),SW2_sigma)
            alphaSm = SP.sum(SW2_sigma, 0);
            barmuS = SP.dot(self.Z.E1,SW_sigma) - b
            self.S.diagSigmaS[:,m] = 1./(1 + alphaSm)     
            self.S.E1[:,m] = barmuS/(1. + alphaSm)      
            
            #keep diagSigmaS
            self.Eps.diagSigmaS[m] = SP.sum(self.S.diagSigmaS[:,m])
        else:
            SW2_sigma = (self.W.C[:, m,0]*(self.W.E2diag[:, m]))*self.Eps.E1 
            alphaSm = SP.sum(SW2_sigma, 0)
            self.S.diagSigmaS[:,m] = 1./(1 + alphaSm)     
            
    def updateW(self,m):
        M = self.components
        logPi = SP.log(self.Pi[:,m]/(1-self.Pi[:,m]))            
        sigma2Sigmaw = (1.0/self.Eps.E1)*self.Alpha.E1[m]

                   
        setMinus = SP.int_(SP.hstack([range(M)[0:m],range(M)[m+1::]]))
        SmTSk = SP.sum( SP.tile(self.S.E1[:,m:m+1],(1, M-1))*self.S.E1[:,setMinus], 0)
        SmTSm = SP.dot(self.S.E1[:,m].transpose(),self.S.E1[:,m]) + self.S.diagSigmaS[:,m].sum()
#            SmTSm = SP.dot(SmTSm.transpose(),OnesND)
        
                   
#            b = SP.sum( (self.W.C[:, setMinus,0]*self.W.E1[:, setMinus])*(SmTSk.transpose()), 1)  
        b = SP.dot( (self.W.C[:, setMinus,0]*self.W.E1[:, setMinus]),(SmTSk.transpose()))                         
        diff = SP.dot(self.S.E1[:,m].transpose(),self.Z.E1) - b
        
        SmTSmSig = SmTSm + sigma2Sigmaw
        
        #update C and W 
        
        u_qm = logPi + 0.5*SP.log(sigma2Sigmaw) - 0.5*SP.log(SmTSmSig) + (0.5*self.Eps.E1)*((diff**2)/SmTSmSig)
        self.W.C[:, m,0] = 1./(1+SP.exp(-u_qm))

        self.W.C[:,m,1] = 1-self.W.C[:,m,0]
        self.W.E1[:, m] = (diff/SmTSmSig)                                #q(w_qm | s_qm=1), q=1,...,Q
        self.W.sigma2[:, m] = (1./self.Eps.E1)/SmTSmSig
        self.W.E2diag[:,m] = self.W.E1[:,m]**2 + self.W.sigma2[:,m]  

            
    def updateAlpha(self,m):
            #update Alpha
        Ewdwd = SP.sum(self.W.C[:,m,0]*self.W.E2diag[:,m])
        self.Alpha.a[m] = self.Alpha.pa + 0.5*Ewdwd        
        self.Alpha.b[m] = self.Alpha.pb + SP.sum(self.W.C[:,m,0])/2.0             
        self.Alpha.E1[m] = self.Alpha.b[m]/self.Alpha.a[m]
        
        
    def updateEps(self):
                #update Eps

        SW_sigma  = self.W.C[:,:,0]*self.W.E1
        SW2_sigma  = self.W.C[:,:,0]*self.W.E2diag
                    
        muSTmuS = SP.dot(self.S.E1.transpose(),self.S.E1)#  + self.S.diagSigmaS
        muSTmuS0 = muSTmuS - SP.diag(SP.diag(muSTmuS))# SP.dot(muSTmuS.transpose(),OnesND)

        t1 = SP.sum(SW_sigma*SP.dot(self.Z.E1.transpose(),self.S.E1), 1)
        t2 = SP.sum(SW2_sigma*SP.tile(SP.diag(muSTmuS).T + self.Eps.diagSigmaS,(self._D,1)), 1) 
        t3 = SP.sum( SP.dot(SW_sigma,muSTmuS0)*SW_sigma, 1)

        self.Eps.E1 = 1./((.5*(self.ZZ  + (-2*t1 + t2 + t3))+self.Eps.pa)/(0.5*self._N+ self.Eps.pb))
        self.Eps.E1[self.Eps.E1>1E5]=1E5

    def update(self):
        M = self.components
        self.Eps.diagSigmaS = SP.zeros((M,))
        mRange = range(M)
        nFix = self.nKnown+self.nLatent
        if self.shuffle==True and self.iterationCount>0:
            mRange[nFix:] = SP.random.permutation(mRange[nFix:])

        for m in mRange:                     
            self.updateW(m)
            self.updateAlpha(m)
            self.updateS(m) 
        self.updateEps()       
        pass
            


    def __init__(self,init_data=None,E1=None,E2=None,**parameters):
        """create the object"""
        #handle setting of parameters via Bayesnet constructor
        ABayesNet.__init__(self,parameters=parameters)
        #priors for the various components:
        if(not self.priors.has_key('Alpha')): self.priors['Alpha']={'priors': [1E-3,1E-3]}
        if(not self.priors.has_key('Eps')):   self.priors['Eps']={'priors':  [1E-3,1E-3]}
        
        self.dataNode=None
        if init_data is None and E1 is not None:
            init_data = CGauss(E1=E1,E2=E2)
        if init_data is not None:
            self.init(init_data)

    def init(self,init_data,Pi=None, init_factors=None):
        """initialize the sparse factor analysis instance"""

        #1. check input data
        #AGAussNode is defined in ExpresisonNet
        #expr Y ~ N(\mu= expr, \sigma = 0)
        if not isinstance(init_data,AGaussNode):
            raise Exception("initialization is only possible from a GaussNode")
        self.Z = CNodeZ(node=init_data)
        #datanode hold the data
        self.dataNode = self.Z
        
        #known factors
        if init_factors!=None and init_factors.has_key('Known'):
            self.nKnown = init_factors['Known'].shape[1]
            self.Known = init_factors['Known']
            assert self.Known.shape[0] == self.Z.E1.shape[0]
            self.nHidden = self.components-self.nKnown
            self.iHidden = SP.arange(self.nHidden,self.nHidden+self.nKnown)
            if init_factors.has_key('Intr'):
                self.nKnown = init_factors['Known'].shape[1]
                self.Known = init_factors['Known']
                assert self.Known.shape[0] == self._N
                self.nHidden = self.components-self.nKnown
                self.iHidden = SP.arange(self.nHidden,self.nHidden+self.nKnown)
        else:
            self.nHidden = self.components
            self.nKnown = 0
            self.iHidden = list()
            
        if init_factors!=None and init_factors.has_key('iLatent'):
            self.iLatent = init_factors['iLatent']
            self.nLatent = len(init_factors['iLatent'])
        else:
            self.nLatent = 0
            
            
        #interaction
        
        #Pi is prior probability of link for genes x factors
        self.Pi = Pi
        self.Non = (self.Pi>.5).sum(0)
        # set dimensionality of the data
        [self._N, self._D] = self.Z.E1.shape
        self.ZZ = SP.zeros((self._D,))
        for d in range(self._D):
            self.ZZ[d] = SP.sum(self.Z.E1[:,d]*self.Z.E1[:,d], 0)
        if self.Pi is not None:
            assert self.Pi.shape == (self._D,self.components)

        #which approximation to use: EP/VB?
        if self.sparsity=='EPV':
            W_node = CNodeWsparseEPV(self)
        elif self.sparsity=='EPVp':
            W_node = CNodeWsparseEPV(self)
            self.permutation_move = True
        elif self.sparsity=='EP':
            W_node = CNodeWsparseEP(self)
        elif self.sparsity=='VB':
            W_node = CNodeWsparseVB(self)
        else:
            W_node = CNodeW(self)
            
#        self.nodes = {'S':CNodeS(self),'W':W_node,'Eps':CNodeEpsFix(self,self.priors['Eps']['priors'])}
        self.nodes = {'S':CNodeSsparse(self),'W':CNodeWsparseVB(self), 'Alpha':CNodeAlphasparse(self,self.priors['Alpha']['priors']),'Eps':CNodeEpsSparse(self,self.priors['Eps']['priors'])}
        for n in self.nodes.keys(): setattr(self,n,self.nodes[n])

        #pca initialisation
        Ion = None
        if self.initType == 'pca':
            Ion = random.rand(self.Pi.shape[0],self.Pi.shape[1])<self.Pi
            self.W.C[:,:,0] = self.Pi
            self.W.C[:,:,0][self.W.C[:,:,0]<=.2] = .1
            self.W.C[:,:,0][self.W.C[:,:,0]>=.8] = .9
            for k in range(self.components):
                sv = linalg.svd(self.Z.E1[:,Ion[:,k]], full_matrices = 0);
                [s0,w0] = [sv[0][:,0:1], S.dot(S.diag(sv[1]),sv[2]).T[:,0:1]]
                v = s0.std(axis=0)
                s0 /= v;
                w0 *= v;
                self.S.E1[:,k] = s0.ravel()
                self.W.E1[Ion[:,k],k] = w0.ravel()
                self.W.E1[~Ion[:,k],k]*=self.sigmaOff
                self.S.diagSigmaS[:,k] = 1./2
        if self.initType == 'pcaRand':
            Ion = self.Pi>.5#random.rand(self.Pi.shape[0],self.Pi.shape[1])<self.Pi
            self.W.C[:,:,0] = self.Pi
            self.W.C[:,:,0][self.W.C[:,:,0]<=.1] = .1
            self.W.C[:,:,0][self.W.C[:,:,0]>=.9] = .9
            for k in range(self.nHidden):
                k+=self.nKnown
                if Ion[:,k].sum()>0:
                    #pdb.set_trace()
                    sv = linalg.svd(self.Z.E1[:,Ion[:,k]], full_matrices = 0);
                    [s0,w0] = [sv[0][:,0:2], S.dot(S.diag(sv[1]),sv[2]).T[:,0:2]]#PC 1 or 2???
                    v = s0.std(axis=0)
                    s0 /= v;
                    w0 *= v;
                    #pdb.set_trace()
                    self.W.E1[:,k] = SP.sqrt(1./self.components)*SP.randn(self._D)        
                    #self.W.E1[Ion[:,k],k] = w0[:,0].ravel()
                    #self.W.E1[~Ion[:,k],k]*=self.sigmaOff

                    self.S.E1[:,k] =0.0* SP.randn(self._N)+(s0[:,0])#+(s0[:,1])
                    
                else:
                    self.S.E1[:,k] = random.randn(self._N,)
                    self.W.E1[:,k] = SP.sqrt(1./self.components)*SP.randn(self._D)
                    

                self.S.diagSigmaS[:,k] = 1./2
                
            
            if self.nKnown>0:
                for k in SP.arange(self.nKnown):
                    self.W.E1[:,k] = SP.sqrt(1./self.components)*SP.randn(self._D)
                    self.S.diagSigmaS[:,k] = 1./2
                self.S.E1[:,SP.arange(self.nKnown)] =  self.Known
            if self.nLatent>0:
                #pca = PCA(n_components=1)
                #pca.fit(self.Z.E1)
                #X = pca.transform(self.Z.E1)
#                covM=SP.cov(self.Z.E1)
#                EV = linalg.eigh(covM)[1][:,(self._N-self.nLatent):(self._N)]
#                EV-=EV.mean(0)
#                EV/=EV.std(0)
                #self.S.E1[:,self.iLatent]=X+.1*SP.randn(self._N,self.nLatent)
                for iL in self.iLatent:
                    self.S.E1[:,iL] = random.randn(self._N,)#+X[:,0]#EV[:,self.nLatent-1]
            if 0:
#                covM=SP.cov(self.Z.E1)
#                EV = linalg.eigh(covM)[1][:,(self._N-1):(self._N)]
#                EV-=EV.mean(0)
#                EV/=EV.std(0)
                pca = PCA(n_components=1)
                pca.fit(self.Z.E1)
                X = pca.transform(self.Z.E1)

                nFix = self.nKnown+self.nLatent
                MPC = abs(vcorrcoef(self.S.E1[:,nFix:].T,X.T))
                Ipi = SP.argsort(-MPC.ravel())
                #pdb.set_trace()
                self.W.Ilabel = SP.hstack([SP.arange(nFix),Ipi+nFix])
                #update the prior Pi
                self.Pi[:,nFix:] = self.Pi[:,Ipi+nFix]
                self.S.E1[:,nFix:] = self.S.E1[:,Ipi+nFix]
                self.W.Pi[:,nFix:] = self.W.Pi[:,Ipi+nFix]                
                self.W.C[:,nFix:,:] = self.W.C[:,Ipi+nFix,:]   
            if 0:

                nFix = self.nKnown+self.nLatent
                FA0 = self
                alphaList=list()
                MADList=list()
                I0 = SP.arange(self.components-nFix)
                for i in range(5):
                    FA=FA0
#                    Ipi = SP.random.permutation(I0)
#                    FA.W.Ilabel = SP.hstack([SP.arange(nFix),Ipi+nFix])
#                #update the prior Pi
#                    FA.Pi[:,nFix:] = FA.Pi[:,Ipi+nFix]
#                    FA.S.E1[:,nFix:] = FA.S.E1[:,Ipi+nFix]
#                    FA.W.Pi[:,nFix:] = FA.W.Pi[:,Ipi+nFix]                
#                    FA.W.C[:,nFix:,:] = FA.W.C[:,Ipi+nFix,:]  
                    for j in range(5):
                        FA.update()
                    MAD = mad(FA.S.E1)
                    alphaList.append(FA.Alpha.E1)
                    MADList.append(MAD)
                Ipi = SP.argsort(SP.mean(SP.array(alphaList),0)[nFix:])
                #pdb.set_trace()
                self.W.Ilabel = SP.hstack([SP.arange(nFix),Ipi+nFix])
                #update the prior Pi
                self.Pi[:,nFix:] = self.Pi[:,Ipi+nFix]
                self.S.E1[:,nFix:] = self.S.E1[:,Ipi+nFix]
                self.W.Pi[:,nFix:] = self.W.Pi[:,Ipi+nFix]                
                self.W.C[:,nFix:,:] = self.W.C[:,Ipi+nFix,:] 
                print self.terms[self.W.Ilabel]
                pdb.set_trace() 
                        
                    
                
            self.initS = self.S.E1.copy()
                
        elif self.initType == 'greedy':
            self.S.E1 = random.randn(self._N,self.components)
            self.W.E1 = random.randn(self._D,self.components)
            Ion = (self.Pi>0.5)
            self.W.E1[~Ion]*= self.sigmaOff
            for k in range(Ion.shape[1]):
                self.W.E1[Ion[:,k]]*=self.sigmaOn[k]

        elif self.initType == 'prior':
            Ion = random.rand(self.Pi.shape[0],self.Pi.shape[1])<self.Pi
            self.W.E1[~Ion]*=self.sigmaOff
            for k in range(Ion.shape[1]):
                self.W.E1[Ion[:,k],k]*=self.sigmaOn[k]
        elif self.initType == 'on':
            for k in range(Ion.shape[1]):
                self.W.E1[:,k]*=self.sigmaOn[k]
        elif self.initType == 'random':
            for k in range(self.Pi.shape[1]):
                self.S.diagSigmaS[:,k] = 1./2
                self.S.E1[:,k] = 2*SP.randn(self._N)
            self.W.E1 = SP.randn(self._D, self.Pi.shape[1])
            self.W.C[:,:,0] = self.Pi
            self.W.C[:,:,0][self.W.C[:,:,0]<=.2] = .1
            self.W.C[:,:,0][self.W.C[:,:,0]>=.8] = .9
            self.initS = self.S.E1.copy()
        elif self.initType == 'data':
            assert ('S' in init_factors.keys())
            assert ('W' in init_factors.keys())
#            Ion = init_factors['Ion']
            Sinit = init_factors['S']
            Winit = init_factors['W']
            self.W.C[:,:,0] = self.Pi
            self.W.C[:,:,0][self.W.C[:,:,0]<=.2] = .1
            self.W.C[:,:,0][self.W.C[:,:,0]>=.8] = .9
            for k in range(self.components):
                self.S.E1[:,k] = Sinit[:,k]
                self.W.E1[:,k] = Winit[:,k]
                self.S.diagSigmaS[:,k] = 1./2




    #calculate the variational bound; currently not used
    def calcBound(self):
        L.debug('CVBFA calcBound')
        #print "Not used"
        #self._bound = ABayesNet.calcBound(self)
#        self.bound=0.
#        
#        
#        #p(data|..)
#        #OLI: is this right? the last term should be <-1/2*tau(D-W*x)^{2}>
#        #try: here we recyle the calculation made in the update Eps:
#        #Bx = -self._N*self._D/2.0*S.log(2*pi) + self._N/2.0*self.Eps.E2.sum() - sum(self.Eps.E1*(self.Eps.a-self.Eps.pa))
#        #Bx = -self._N*self._D/2.0*S.log(2*pi) + self._N/2.0*self.Eps.E2.sum() + sum(self.Eps.E1*self.Eps.pa-self.Eps.b)
#
#
#        #note : trace (S.cov) comes from the fact that the 2nd moment of S is not just S.E1**2 but + cov!
#        #KL q(S)/P(S)
#        
#        #orig
#        #Bss= -self._N/2.0*logdet(self.S.cov) - self._N/2.0*trace(eye(self.components)-self.S.cov) + 0.5*(self.S.E1**2).sum()
#
#        #KL q(W)/p(W|alpha)
#        Bww = 0
#        for k in range(self.components):
#            Bww-=0.5*self.Non[k]*(special.digamma(self.Alpha.b[k])-S.log(self.Alpha.a[k]))
#        #Bww= -self._D/2.0*sum(special.digamma(self.Alpha.b)-S.log(self.Alpha.a))
#
# #       for d in range(self._D):
# #           Bww = Bww - 1/2.0*( logdet(self.W.cov[d,:,:]) + trace(eye(self.components)-dot(self.W.E2[d,:,:],diag(self.Alpha.E1))))
#
#        #self._bound = self._bound + Bx - Bss - Bww
#        self._boundLOG.append(self._bound)
#        L.debug('CVBFA bound = %.2f'%self._bound)

        return self._bound






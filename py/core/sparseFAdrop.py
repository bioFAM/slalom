# SparseFA
# variational VB/EP sparse factor analysis model
# this class implements either a standard variational inference procedure for a sparse model or a hybrid approach.
import sys
sys.path.append('/Users/flo/projects/Auto_Bionf/scLVM2/py/core')
from vbfa import *
import scipy as SP
from scipy.stats import truncnorm
#from mlib.utils import *
#from mlib.plot import *
import copy
import pdb


class CNodeAlphasparse(AGammaNode):
    def __init__(self,net,prior=[1E-3,1E-3]):
        AGammaNode.__init__(self,dim=[net.components],prior=prior)
        
    def update(self,net):
        W = net.W
        #
        #WmuWmuT = SP.zeros((Q,K,K))
        #for q in range(Q):
        #    WmuWmuT[q,:,:] = SP.dot(W.E1[q:q+1,:].transpose(), W.E1[q:q+1,:])
        #pdb.set_trace()
        #WmuWmuT = WmuWmuT.sum(axis=0)
        Ewdwd = SP.sum(W.C[:,:,0]*W.E2diag,axis=0)
        #self.a[:] = self.pa + 0.5*diag(Ewdwd+WmuWmuT)
#        (sum( vardist.Gamma(:,m).*(vardist.muW(:,m).^2+vardist.sigma2W(:,m)) ))/2;
        self.a[:] = self.pa + 0.5*Ewdwd
        
        self.b[:] = self.pb + SP.sum(W.C[:,:,0],axis=0)/2.0
        #pdb.set_trace()
        #update expectation values
        AGammaNode.update(self)

class CNodeEpsSparse(CNodeEps):
    """Extensions of CNodeEps that can handle fixed prior expectaitons"""
    def __init__(self,net,prior=S.array([100,1])):
            CNodeEps.__init__(self,net,prior)
            self.E1 = SP.minimum(self.E1, 1/4.)

    def update(self,net):
    #add terms 	tiled_A = np.tile(np.resize(A, [1, D, K]), [N, 1, 1])
        self.E1 = SP.repeat(net.kappa, net._D)
      
                
        
        

class CNodeSsparse(AVGaussNode):
    def __init__(self,net,prec=1):
#        CNodeS.__init__(self,net,prec=1)
        AVGaussNode.__init__(self,dim=[net._N,net.components],cdim=1)
        self.diagSigmaS = SP.ones((net._N,net.components))

    def update(self,net=None):
        """update(net):
        where net is the "base networ, i.e. sparseFA or similar"""
        W     = net.W
        M = net.components      
        for m in xrange(M):       
            SW_sigma = (W.C[:, m,0]*W.E1[:, m])*net.Eps.E1
            SW2_sigma = (W.C[:, m,0]*(W.E2diag[:, m]))*net.Eps.E1   
            alphaSm = SP.dot(SP.ones(net.Z.E1.shape),SW2_sigma)
            setMinus = SP.int_(SP.hstack([range(M)[0:m],range(M)[m+1::]]))
            #diagSigmaSm = 1./(1 + SW2_sigma)           
            
            #SmTSk = SP.tile(self.E1[:,m:m+1],(1,M-1))*self.E1[:,setMinus]           
            #SmTSm = self.E1[:,m]*self.E1[:,m] + diagSigmaSm           
            #sigma2Sigmaw = (1.0/net.Eps.E1)/W.sigma2[m]
             
            b0 = SP.dot(self.E1[:,setMinus],(W.C[:, setMinus,0]*W.E1[:, setMinus]).transpose())
            b=SP.dot(b0,SW_sigma)
            barmuS = SP.dot(net.Z.E1,SW_sigma) - b
            self.diagSigmaS[:,m] = 1./(1 + alphaSm)     
            self.E1[:,m] = barmuS/(1. + alphaSm)
#           pdb.set_trace()
            #Eps=net.Eps
            #W = net.W
            #self.cov = linalg.inv(p)
            #p=eye(net.components)+tensordot(Eps.E1,W.E2,[0,0])
            #matrix   = dot( dot(self.cov,W.E1.T), diag(Eps.E1) )
            #self.E1  = dot(net.Z.E1,matrix.T)
            #AVGaussNode.update(self)
        pass

class CNodeWsparse(CNodeW):
    """Abstract CnodeWsparse basclass
    - this ensures consitent handling of alternative CNodeW nodes; in particular the sparsity prior
    - main application : permutation move of factors"""
    def __init__(self,net,**kwargin):
        #call base class initialisation
        CNodeW.__init__(self,net,**kwargin)
        #create indicator arrays for moment calculation
      
        self.Pi = zeros([net.Pi.shape[0],net.Pi.shape[1],2])
        self.Pi[:,:,0] =net.Pi
        self.Pi[:,:,1] = 1.-net.Pi
 #       self.prec = S.hstack([net.sigmaOff**(-2),net.sigmaOn**(-2)])
        #log prior

        self.C    = self.Pi.copy()
        #self.C[self.Pi>.8] = .9
        #self.C[self.Pi<.1]= .1  
        #labeling of factors in case we use permutation move
        self.Ilabel = SP.arange(net.components)

    def update(self,net=None):
        def logProbkk(k,l):
            """evaluate the probability of C_k and Pi_l"""
            pp = C[:,k,:]*Pi[:,l,:]
            lpp = SP.log(pp.sum(axis=1))
            return lpp.sum()

        if (net is None) or (net.permutation_move==False): 
            return
        print "pong"
        #do factor permutation if active
        #use the marignal indicators to calculate this; I thik they contain all we need; however we need to divide out the prior

        pass



        
class CNodeWsparseVEM(CNodeWsparse):
    def __init__(self,net,prec=1.):
        CNodeWsparse.__init__(self,net,prec=prec)
        self.sigma2 = (1.0/prec)*SP.ones((net._D, net.components))
        self.E1 = SP.randn(net._D, net.components)
        self.E2diag = SP.zeros((net._D, net.components))#is calculated in update for W
        #variable initialisation in CNodeWsparse
        #self.sigma2 = (1.0/prec)*SP.ones((net._D, net.components))
        #self.E2diag = SP.zeros((net._D, net.components))
        #    self.E2diag[d,:] = SP.diag(self.E2[d,:,:])
        #for d in range(net._D):

    def update(self,net=None):
        if(net==None):
            AVGaussNode.update(self)
            print "pong VEM"
            return
        #need to replace Y*q(S) with q(XS) in diff
        M = net.components    
        for m in xrange(M):
            logPi = SP.log(net.Pi[:,m]/(1-net.Pi[:,m]))
            setMinus = SP.int_(SP.hstack([range(M)[0:m],range(M)[m+1::]]))
            SmTSk = SP.tile(net.S.E1[:,m:m+1],(1,M-1))*net.S.E1[:,setMinus];
            SmTSk = SP.dot(SmTSk.transpose(),SP.ones((net._N,net._D)))
            SmTSm = net.S.E1[:,m]*net.S.E1[:,m] + net.S.diagSigmaS[:,m]
            SmTSm = SP.dot(SmTSm.transpose(),SP.ones((net._N,net._D)))
            
            sigma2Sigmaw = (1.0/self.Eps.E1)*self.Alpha.E1[m]
            
            b = SP.sum( (self.C[:, setMinus,0]*self.E1[:, setMinus])*(SmTSk.transpose()), 1)                         
            diff = SP.dot(net.S.E1[:,m].transpose(),net.Z.E1) - b.transpose()
            SmTSmSig = SmTSm + sigma2Sigmaw
            u_qm = logPi + 0.5*SP.log(sigma2Sigmaw) - 0.5*SP.log(SmTSmSig) + (0.5*net.Eps.E1)*((diff**2)/SmTSmSig)
            self.C[:, m,0] = 1./(1+SP.exp(-u_qm))

            self.C[:,m,1] = 1-self.C[:,m,0]
            self.E1[:, m] = (diff/SmTSmSig)                                #q(w_qm | s_qm=1), q=1,...,Q
            self.sigma2[:, m] = (1./net.Eps.E1)/SmTSmSig
            self.E2diag[:,m] = self.E1[:,m]**2 + self.sigma2[:,m]




class CNodeWsparseVB(CNodeWsparse):
    def __init__(self,net,prec=1):
        CNodeWsparse.__init__(self,net,prec=prec)
        #variable initialisation in CNodeWsparse
        pass

    def update(self,net=None):
        if(net==None):
            AVGaussNode.update(self)
            return
        S     = net.S
        Eps   = net.Eps
        M = S.E2.sum(axis=0)
        tauOff = net.sigmaOff**(-2)
        tauOn = net.Alpha.E1

        #integrated procedure: save but does not add anything:
        if 1:
            #1. calcualte the sufficient statistics of the incoming messages
            #we only store the diagonal stuff here
            E1_in  = zeros_like(self.E1)
            prec_in= zeros_like(self.E1)
            for d in xrange(net._D):
                pin   =  Eps.E1[d]*M
                pprior = diag(tauOn*self.C[d,:,1] + self.C[d,:,0]*(tauOff))
                ptotal = pin + pprior
                self.cov[d,:,:] = linalg.inv(ptotal)
                self.E1[d,:] = dot(self.cov[d,:,:],Eps.E1[d]*dot(S.E1.T,net.Z.E1[:,d]))
            #update moments
            CNodeW.updateE2(self,net)

#        CNodeWsparse.update(self,net)
        pass
    




class CSparseFA(AExpressionModule):
    '''CVBFA(AExpressionModule)
    - Variational Bayesian Factor analysis module'''


    def getDefaultParameters(self):
        """return a hash with default parameter value for this BayesNet"""
        #
        dp = AExpressionModule.getDefaultParameters(self)
        dp['initType'] = 'greedy'
        dp['nIterations'] = 20
        #VB best
        dp['schedule'] = ['S','Eps','W']
#        #experiments(EP)
        dp['schedule'] = ['Eps','S','W']
        dp['schedule'] = ['W','S','Eps','Alpha']
        dp['permutation_move'] = False
        dp['shuffle']=True
        dp['components'] = 5
        dp['priors']     = {}
        dp['sigmaOff']   = 0.1
        dp['sigmaOn']    = S.ones(dp['components'])*1.0
        dp['sparsity']   = 'VB'
        dp['iterative']  = True
        dp['randperm']   = False
        return dp

    def getName(self,base_name='sparseFA'):
        """return a name summarising the  main parameters"""
        sparsity = copy.copy(self.sparsity)
        if sparsity=='EPV':
            sparsity = 'VB/EP'
        if self.permutation_move:
            base_name = base_name + 'p'
        sigma_str = '%.1e' % (self.sigmaOff**2)
        exponent = sigma_str[4::]
        sigma_str = "$\sigma^{2}_{0}=10^{%s}$"% (exponent)
        #translate the sigma
        name = "%s %s, %s" % (base_name,sparsity,sigma_str)
        return name


    def iterate(self, nIterations=None, forceIterations=None, tolerance=None, minIterations=10):
        '''iterate(nIteations=None,forceIterations=None)
        - perform nIterations; per default(None) parameters are tken from local intsance settings
        '''
        #forceIterations=True
        L.debug('SparseFA iterate')
        iter=0
        if tolerance is None: tolerance = self.tolerance
        if nIterations is None: nIterations = self.nIterations
        if forceIterations is None: forceIterations = self.forceIterations
        Zr = S.dot(self.S.E1,self.W.E1.T)
        Zd = self.Z.E1-Zr
        error = (Zd**2).mean()
        self.calcBound()
        if SP.mod(iter,1)==0:
            print "reconstruction error: %f lower bound: %f" % (error,  self._bound)
        LB = 0

        for iter in range(nIterations):
            #self.iterationCount+=1
            t = time.time();
            #for node in self.schedule:
            #    self.updateNode(node)
                #pdb.set_trace()
            self.update()
            self.iterationCount+=1
            #calc reconstruction error
            xi = SP.dot(self.S.E1,self.W.E1.T)
            #Zr = (SP.random.poisson(rate_fun(xi)))
            Zd = self.Z.E1-xi
            error = (Zd**2).mean()
            #self.calcBound()
            if SP.mod(iter,100)==0:
                print "reconstruction error: %f lower bound: %f" % (error,  self._bound) 

            if (abs(LB - self._bound) < tolerance) and not forceIterations and iter>minIterations:
                print 'Converged after %i iterations' % (iter)
                break

            L.info("Iteration %d: time=%.2f bound=%f" % (iter,time.time() - t, self._bound))
            LB = self._bound

        return self._bound

        
        #1. itearte
        AExpressionModule.iterate(self,*argin,**kwargin)
        #2. debug
        pass

    def updateS(self,m):
        M = self.components
        if m>=self.nKnown:        
            YmeanX = self.meanX
            setMinus = SP.int_(SP.hstack([range(M)[0:m],range(M)[m+1::]]))
            #diagSigmaSm = 1./(1 + SW2_sigma)           
            
            #SmTSk = SP.tile(self.E1[:,m:m+1],(1,M-1))*self.E1[:,setMinus]           
            #SmTSm = self.E1[:,m]*self.E1[:,m] + diagSigmaSm           
            #sigma2Sigmaw = (1.0/net.Eps.E1)/W.sigma2[m]
            #update S
            SW_sigma = (self.W.C[:, m,0]*self.W.E1[:, m])*self.Eps.E1
            SW2_sigma = (self.W.C[:, m,0]*(self.W.E2diag[:, m]))*self.Eps.E1             
             
            b0 = SP.dot(self.S.E1[:,setMinus],(self.W.C[:, setMinus,0]*self.W.E1[:, setMinus]).transpose())
            b=SP.dot(b0,SW_sigma)
            alphaSm = SP.dot(SP.ones(self.Z.E1.shape),SW2_sigma)
            barmuS = SP.dot(YmeanX,SW_sigma) - b# - SP.dot(self.bias,SW_sigma)
            self.S.diagSigmaS[:,m] = 1./(1 + alphaSm)     
            self.S.E1[:,m] = barmuS/(1. + alphaSm)    
            
        else:
            SW2_sigma = (self.W.C[:, m,0]*(self.W.E2diag[:, m]))*self.Eps.E1 
            alphaSm = SP.sum(SW2_sigma, 0)
            self.S.diagSigmaS[:,m] = 1./(1 + alphaSm) 
            
    def updateW(self,m):
        M = self.components
        logPi = SP.log(self.Pi[:,m]/(1-self.Pi[:,m]))            
        sigma2Sigmaw = (1.0/self.Eps.E1)*self.Alpha.E1[m]
        YmeanX = self.meanX
                   
        setMinus = SP.int_(SP.hstack([range(M)[0:m],range(M)[m+1::]]))
        SmTSk = SP.tile(self.S.E1[:,m:m+1],(1,M-1))*self.S.E1[:,setMinus];
        SmTSk = SP.dot(SmTSk.transpose(),SP.ones((self._N,self._D)))
        SmTSm = self.S.E1[:,m]*self.S.E1[:,m] + self.S.diagSigmaS[:,m]
        SmTSm = SP.dot(SmTSm.transpose(),SP.ones((self._N,self._D)))
        
                   
        b = SP.sum( (self.W.C[:, setMinus,0]*self.W.E1[:, setMinus])*(SmTSk.transpose()), 1)                         
        diff = SP.dot(self.S.E1[:,m].transpose(),YmeanX) - b# - SP.dot(self.S.E1[:,m].T,self.bias)
        SmTSmSig = SmTSm + sigma2Sigmaw
        
        #update C and W 
        u_qm = logPi + 0.5*SP.log(sigma2Sigmaw) - 0.5*SP.log(SmTSmSig) + (0.5*self.Eps.E1)*((diff**2)/SmTSmSig)

        self.W.C[:, m,0] = 1./(1+SP.exp(-u_qm))
        self.W.C[:, m,0]
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
        M = self.components
        SW_sigma  = self.W.C[:,:,0]*self.W.E1
        SW2_sigma  = self.W.C[:,:,0]*self.W.E2diag
        muSTmuS = self.S.E1*self.S.E1  + self.S.diagSigmaS
        muSTmuS = SP.dot(muSTmuS.transpose(),self.isExpressed)
        t1 = SP.sum(SW_sigma*SP.dot(self.Z.E1.transpose(),self.S.E1), 1)
        t2 = SP.sum(SW2_sigma.transpose()* muSTmuS,0)
        t3 = SP.zeros((self._D,))
        for m in range(M):
            for m1 in SP.arange(m+1,M):
                tt = ( (self.W.C[:, m1,0]*self.W.E1[:, m1])*SW_sigma[:, m] )
                t3 = t3 + tt*SP.dot((self.S.E1[:,m1]*self.S.E1[:,m]).transpose(),self.isExpressed)
        #pdb.set_trace()
        self.Eps.E1 = 1./((self.ZZ  + (-2*t1 + t2 + 2*t3))/self.numExpressed)        
                
        self.Eps.E1[self.Eps.E1<1/4.]=1/4.    
        
        #EITHER: Standardize (expressed) genes to unit variance  and set kappa to 0.25
        #OR: use min(0.25, Eps.E1) (!) - maybe with standardization...
        #OR use Poisson noise model and then max(kappa poisson, 0.25) todo on the train..then check with toy data...
        
            
    def update(self):
        M = self.components
        mRange = range(M)
        nFix = self.nKnown+self.nLatent
        if self.shuffle==True and self.iterationCount>0:
            mRange[nFix:] = SP.random.permutation(mRange[nFix:])
            #mRange[nFix:] = SP.random.permutation(mRange[nFix:])
        #print mRange[0:4]
        for m in mRange:                      
            self.updateW(m)                                              
            self.updateS(m)                                                             
       #  for m in mRange:
            self.updateAlpha(m) 
            
        self.updateEps()                
            
                
        #update xi folowing Seeger et al 2012, AISTATS
        Xi = SP.dot(self.S.E1,(self.W.C[:, :,0]*self.W.E1).transpose())
#        if(SP.any(SP.isnan(Xi))):
#            pdb.set_trace()
        
        self.meanX[self.isExpressed==0] = (Xi - 4.*(1./(1.+SP.exp(-Xi))))[self.isExpressed==0]
       
        pass
            


    def __init__(self,init_data=None,E1=None,E2=None,**parameters):
        """create the object"""
        #handle setting of parameters via Bayesnet constructor
        ABayesNet.__init__(self,parameters=parameters)
        #priors for the various components:
        if(not self.priors.has_key('Alpha')): self.priors['Alpha']={'priors': [1E-3,1E-3]}
        if(not self.priors.has_key('Eps')):   self.priors['Eps']={'priors': [1,100]}
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
        self.kappa = 1./4.0# + 0.17*self.Z.E1.max()
        
        self.isExpressed = (self.Z.E1>0)*1.

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
        
        
        
        
        #Pi is prior probability of link for genes x factors
        self.Pi = Pi
        self.Non = (self.Pi>.5).sum(0)
        # set dimensionality of the data
        [self._N, self._D] = self.Z.E1.shape
        self.ZZ = SP.zeros((self._D,))
        
        for d in range(self._D):
            self.ZZ[d] = SP.sum(self.Z.E1[self.isExpressed[:,d]==1.,d]*self.Z.E1[self.isExpressed[:,d]==1.,d], 0)
        self.numExpressed = SP.sum(self.Z.E1,0)
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
            
        #self.nodes = {'S':CNodeS(self),'W':W_node,'Eps':CNodeEpsFix(self,self.priors['Eps']['priors'])}
        self.nodes = {'S':CNodeSsparse(self),'W':CNodeWsparseVEM(self), 'Alpha':CNodeAlphasparse(self,self.priors['Alpha']['priors']),'Eps':CNodeEpsSparse(self,self.priors['Eps']['priors'])}
        for n in self.nodes.keys(): setattr(self,n,self.nodes[n])

        #pca initialisation
        Ion = None
        Zstd = self.Z.E1
        Zstd-=Zstd.mean(0)
        if self.initType == 'pca':
            Ion = random.rand(self.Pi.shape[0], self.Pi.shape[1]) < self.Pi
            for k in range(self.components):
                sv = linalg.svd(self.Z.E1[:, Ion[:, k]], full_matrices=0)
                s0, w0 = sv[0][:, 0:1], S.dot(S.diag(sv[1]), sv[2]).T[:, 0:1]
                v = s0.std(axis=0)
                s0 /= v
                w0 *= v
                self.S.E1[:, k] = s0.ravel()
                self.W.E1[Ion[:, k], k] = w0.ravel()
                self.W.E1[(~Ion[:, k], k)] *= self.sigmaOff
                
        if self.initType == 'pcaRand':
           
            Ion = self.Pi>.5
            self.W.C[:,:,0] = self.Pi
            self.W.C[:,:,0][self.W.C[:,:,0]<=.2] = .1
            self.W.C[:,:,0][self.W.C[:,:,0]>=.8] = .9
            self.W.C[:,:,1] = 1.-self.W.C[:,:,0]
            for k in range(self.nHidden):
                k+=self.nKnown
                if Ion[:,k].sum()>0:
                    sv = linalg.svd(Zstd[:, Ion[:, k]], full_matrices = 0);
                    [s0,w0] = [sv[0][:,0:2], S.dot(S.diag(sv[1]),sv[2]).T[:,0:2]]#PC 1 or 2???
                    v = s0.std(axis=0)
                    s0 /= v;
                    w0 *= v;
                    self.W.E1[:,k] = SP.sqrt(1./self.components)*SP.randn(self._D)        
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
            #pdb.set_trace()        
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
        elif self.initType == 'data':
            assert ('S' in init_factors.keys())
            assert ('W' in init_factors.keys())
#            Ion = init_factors['Ion']
            Sinit = init_factors['S']
            Winit = init_factors['W']
            for k in range(self.components):
                self.S.E1[:,k] = Sinit[:,k]
                self.W.E1[:,k] = Winit[:,k]
#                self.W.E1[~Ion[:,k],k]*=self.sigmaOff
        #Xi = SP.dot(self.S.E1,(self.W.C[:, :,0]*self.W.E1).transpose() )
        self.meanX = self.Z.E1.copy()#Xi - self.fprime(Xi, self.Z.E1)/self.kappa
        self.initS = self.S.E1.copy()
        self.initW = self.W.E1.copy()

        

        #update moments
        if 0:
            self.S.update(self)
            self.W.update(self)



    #calculate the variational bound:
    def calcBound(self):
        L.debug('CVBFA calcBound')
        #self._bound = ABayesNet.calcBound(self)
        self.bound=0.
        
        
        #p(data|..)
        #OLI: is this right? the last term should be <-1/2*tau(D-W*x)^{2}>
        #try: here we recyle the calculation made in the update Eps:
        #Bx = -self._N*self._D/2.0*S.log(2*pi) + self._N/2.0*self.Eps.E2.sum() - sum(self.Eps.E1*(self.Eps.a-self.Eps.pa))
        #Bx = -self._N*self._D/2.0*S.log(2*pi) + self._N/2.0*self.Eps.E2.sum() + sum(self.Eps.E1*self.Eps.pa-self.Eps.b)


        #note : trace (S.cov) comes from the fact that the 2nd moment of S is not just S.E1**2 but + cov!
        #KL q(S)/P(S)
        
        #orig
        #Bss= -self._N/2.0*logdet(self.S.cov) - self._N/2.0*trace(eye(self.components)-self.S.cov) + 0.5*(self.S.E1**2).sum()

        #KL q(W)/p(W|alpha)
#        Bww = 0
#        for k in range(self.components):
#            Bww-=0.5*self.Non[k]*(special.digamma(self.Alpha.b[k])-S.log(self.Alpha.a[k]))
#        #Bww= -self._D/2.0*sum(special.digamma(self.Alpha.b)-S.log(self.Alpha.a))
#
#        #for d in range(self._D):
#        #    Bww = Bww - 1/2.0*( logdet(self.W.cov[d,:,:]) + trace(eye(self.components)-dot(self.W.E2[d,:,:],diag(self.Alpha.E1))))
#
#        #self._bound = self._bound + Bx - Bss - Bww
#        self._boundLOG.append(self._bound)
#        L.debug('CVBFA bound = %.2f'%self._bound)

        return self._bound

    def fprime(self, xi, y):
        return pi_fun(xi)*(1-y/rate_fun(xi))

    def getPrediction(self):
        L.info('CVBFA getPrediction')

        #make sure we always produce a prediction even if not intialized
        if self.iterationCount==0:
            return CGauss(E1=S.array([0]),prec=S.array([0]))
        p = dot(self.S.E1,self.W.E1.T)
        E1 = real(p)
        prec = ones(shape = self.Z.E1.shape)*self.Eps.E1

        return CGauss(E1=E1,prec=prec)


    def residuals(self):
        L.info('CVBFA residuals')
        return self.Z.E1 - self.getPrediction().E1


def rate_fun(x):
    rate =  SP.log(1+SP.exp(x))
    rate[rate<1e-10]=1e-10
    return rate
    
def pi_fun(x):
    return 1.0/(1.+SP.exp(-x))



def calc_accuracy(Z,X):
    """calc predictive accuracy
    Z: responsibilties
    X: truth"""
    Nok = (Z==X).sum()
    return double(Nok)/X.size
    

def accuracy_experiment(vbfa,Xt,nIblock =5):
    """helper function to run an accuracy experiment
    vbfa: vbfa object to use (needs to be initialised)
    Xt: true network for evaluation
    nIblock: number of itertions per block (2)
    """
    nIterations = vbfa.nIterations
    #number of blocks?
    nblocks     = nIterations/nIblock

    PL.ion()
    R = S.zeros([nblocks,10])
    for i in xrange(nblocks):
        vbfa.iterate(nIterations=nIblock)
        #get the prdictions and calculate precisions
        Z = vbfa.W.C[:,:,0]
        # calc accuracy
        auc = ROC.auroc(labels=Xt,predictions=Z)
        acc = calc_accuracy(Z>0.5,Xt)
        fon = double(Z.sum())/Z.size
        if 0:
            #DEBUGGING
            E1_in = vbfa.W.E1_in
            prec_in = vbfa.W.prec_in
            PL.figure(1)
            PL.clf()
            PL.hist(E1_in,100)
            PL.figure(2)
            PL.clf()
            PL.hist(prec_in,100)
            PL.draw()
            raw_input()        
        R[i,0] = acc
        R[i,1] = auc
        R[i,2] = fon
        print "accuracy: %f (fraction on on %f)" % (acc,fon)
    return R
    


if __name__ =='__main__':
    #load data and perform sparse FA:
    data_base = './../Data/DataUAI/Toy_486_20_20'
#    data_base = './../Data/DataUAI/Toy_1000_60_20'
    import scipy
    import scipy.io
    import os
    import logging as L
    import pylab as PL
    L.getLogger().setLevel(L.INFO)

    random.seed(5)


    #data:
    data_train = os.path.join(data_base,'X_PI_Y_AlphaY_BetaY.mat')
    data_true  = os.path.join(data_base,'SubnetData.mat')
    data      = scipy.io.loadmat(data_train)
    data_true = scipy.io.loadmat(data_true)
    #prior
    Pi        = data['Pi']
    Y         = data['Y'].T
    #true data
    truth     = data_true['SubnetData']
    Xt        = (truth.X>0)
    Pit       = truth.Pi

    if 0:
        X = random.randn(20,5)
        W = random.randn(5,50)
        Y = S.dot(X,W)
        components = 5

    components = Pi.shape[1]
    

    if 0:
        Y = Y[0:15,:]
        vbfa2 = CVBFA(components=components,initType='random')
        vbfa2.init(CGauss(Y))
        vbfa2.iterate(nIterations=20)

        
    nIterations=100
    sigmaOff   =0.1

    #make a differnet problem: subset well determined, rest unknown:
    
#    X      = Pi>0.5
#    Pi[X] = 0.5
#    Pi[~X] = 0.5
    #1. VB
#    Pi[:,:] = 1.0
    if 1:
        vbfaVB = CSparseFA(components=components,sigmaOff=sigmaOff,sparsity='VB',nIterations=nIterations)
        vbfaVB.init(CGauss(Y),Pi=Pi)
        v2 = copy.copy(vbfaVB)
        vbfaVB.init(CGauss(Y),Pi=Pi)
        RVB    = accuracy_experiment(vbfaVB,Xt)

    #2. EP
    if 0:
        vbfaEP = CSparseFA(components=components,sigmaOff=sigmaOff,sparsity='EP',nIterations=nIterations)
        vbfaEP.init(CGauss(Y),Pi=Pi)
        REP    = accuracy_experiment(vbfaEP,Xt)

    #3. EPV
    if 0:
        vbfaEPV = CSparseFA(components=components,sigmaOff=sigmaOff,sparsity='EPV',nIterations=nIterations)
        vbfaEPV.init(CGauss(Y),Pi=Pi)
        REPV    = accuracy_experiment(vbfaEPV,Xt)
        
    
    if 0:

        #how many blocks and iterations per block
        nIterations = 4
        nIterations_block = 2
        R = S.zeros([nIterations,3])
        for i in xrange(nIterations):
            vbfa.iterate(nIterations=nIterations_block)
            #get the prdictions and calculate precisions
            Z = vbfa.W.C[:,:,1]
            # calc accuracy
            auc = ROC.auroc(labels=Xt,predictions=Z)
            acc = calc_accuracy(Z>0.5,Xt)
            fon = double(Z.sum())/Z.size
            R[i,0] = acc
            R[i,1] = auc
            R[i,2] = fon
            print "accuracy: %f (fraction on on %f)" % (acc,fon)
        #calculat area under the curve

    

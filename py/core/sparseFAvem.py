# SparseFA
# variational VB/EP sparse factor analysis model
# this class implements either a standard variational inference procedure for a sparse model or a hybrid approach.

from vbfa import *
import scipy as SP
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


    def update(self,net):
#           SW_sigma  = vardist.Gamma.*vardist.muW;
#           SW2_sigma = vardist.Gamma.*(vardist.muW.^2 + vardist.sigma2W); 
#           %muSTmuS = vardist.muS'*vardist.muS;
#           muSTmuS = vardist.muS.*vardist.muS  + vardist.diagSigmaS;
#           muSTmuS = muSTmuS'*model.indexPresent;
#           
#           %muSTmuS0 = muSTmuS - diag(diag(muSTmuS));
#           t1 = sum(SW_sigma.*(model.Y'*vardist.muS), 2);
#           %t2 = sum(SW2_sigma.*repmat(diag(muSTmuS)' + diagSigmaS,Q,1), 2); 
#           t2 = sum(SW2_sigma'.* muSTmuS,1); 
#           %t3 = sum( (SW_sigma*muSTmuS0).*SW_sigma, 2);
#           
#           t3 = zeros(1, model.Q);
#           t3 = t3 + tt.*((vardist.muS(:,m1).*vardist.muS(:,m))'*model.indexPresent);
        M = net.components
        SW_sigma  = net.W.C[:,:,0]*net.W.E1
        SW2_sigma  = net.W.C[:,:,0]*net.W.E2diag
        
        muSTmuS = net.S.E1*net.S.E1  + net.S.diagSigmaS
        muSTmuS = SP.dot(muSTmuS.transpose(),SP.ones((net._N,net._D)))
        t1 = SP.sum(SW_sigma*SP.dot(net.Z.E1.transpose(),net.S.E1), 1)
        t2 = SP.sum(SW2_sigma.transpose()* muSTmuS,0)
        t3 = SP.zeros((net._D,))
        for m in range(M):
            for m1 in SP.arange(m+1,M):
                tt = ( (net.W.C[:, m1,0]*net.W.E1[:, m1])*SW_sigma[:, m] )
                t3 = t3 + tt*SP.dot((net.S.E1[:,m1]*net.S.E1[:,m]).transpose(),SP.ones((net._N,net._D)))

        self.E1 = 1./((net.ZZ  + (-2*t1 + t2 + 2*t3))/net._N)
      
                
        
        

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
        #variable initialisation in CNodeWsparse
        self.sigma2 = (1.0/prec)*SP.ones((net._D, net.components))
        self.E2diag = SP.zeros((net._D, net.components))
        for d in range(net._D):
            self.E2diag[d,:] = SP.diag(self.E2[d,:,:])

    def update(self,net=None):
        if(net==None):
            AVGaussNode.update(self)
            print "pong VEM"
            return

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
            self.iterationCount+=1
            t = time.time();
            #for node in self.schedule:
            #    self.updateNode(node)
                #pdb.set_trace()
            self.update()
            #calc reconstruction error
            Zr = S.dot(self.S.E1,self.W.E1.T)
            Zd = self.Z.E1-Zr
            error = (Zd**2).mean()
            self.calcBound()
            if SP.mod(iter,500)==0:
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

    def update(self):
        M = self.components
        #OnesND = SP.ones((self._N,self._D))
        diagSigmaS = SP.zeros((M,))
        for m in range(M):
            logPi = SP.log(self.Pi[:,m]/(1-self.Pi[:,m]))
            sigma2Sigmaw = (1.0/self.Eps.E1)*self.Alpha.E1[m]
            
            setMinus = SP.int_(SP.hstack([range(M)[0:m],range(M)[m+1::]]))
#            SmTSk = SP.dot(SP.tile(self.S.E1[:,m:m+1],(1,M-1)).transpose(),self.S.E1[:,setMinus])
            #SmTSk = SP.dot(SmTSk.transpose(),OnesND)
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
            
            #update Alpha
            Ewdwd = SP.sum(self.W.C[:,m,0]*self.W.E2diag[:,m])
            self.Alpha.a[m] = self.Alpha.pa + 0.5*Ewdwd        
            self.Alpha.b[m] = self.Alpha.pb + SP.sum(self.W.C[:,m,0])/2.0     
            
            self.Alpha.E1[m] = self.Alpha.b[m]/self.Alpha.a[m]
            
            
            #update S
            SW_sigma = (self.W.C[:, m,0]*self.W.E1[:, m])*self.Eps.E1
            SW2_sigma = (self.W.C[:, m,0]*(self.W.E2diag[:, m]))*self.Eps.E1   
            
            setMinus = SP.int_(SP.hstack([range(M)[0:m],range(M)[m+1::]]))
            #diagSigmaSm = 1./(1 + SW2_sigma)           
            
            #SmTSk = SP.tile(self.E1[:,m:m+1],(1,M-1))*self.E1[:,setMinus]           
            #SmTSm = self.E1[:,m]*self.E1[:,m] + diagSigmaSm           
            #sigma2Sigmaw = (1.0/net.Eps.E1)/W.sigma2[m]
             
            b0 = SP.dot(self.S.E1[:,setMinus],(self.W.C[:, setMinus,0]*self.W.E1[:, setMinus]).transpose())
            b=SP.dot(b0,SW_sigma)
            #alphaSm = SP.dot(SP.ones(self.Z.E1.shape),SW2_sigma)
            alphaSm = SP.sum(SW2_sigma, 0);
            barmuS = SP.dot(self.Z.E1,SW_sigma) - b
            self.S.diagSigmaS[:,m] = 1./(1 + alphaSm)     
            self.S.E1[:,m] = barmuS/(1. + alphaSm)      
            
            #keep diagSigmaS
            diagSigmaS[m] = SP.sum(self.S.diagSigmaS[:,m])
        
        #update Eps
        M = self.components
        SW_sigma  = self.W.C[:,:,0]*self.W.E1
        SW2_sigma  = self.W.C[:,:,0]*self.W.E2diag
             
        
        muSTmuS = SP.dot(self.S.E1.transpose(),self.S.E1)#  + self.S.diagSigmaS
        muSTmuS0 = muSTmuS - SP.diag(SP.diag(muSTmuS))# SP.dot(muSTmuS.transpose(),OnesND)

        t1 = SP.sum(SW_sigma*SP.dot(self.Z.E1.transpose(),self.S.E1), 1)
        t2 = SP.sum(SW2_sigma*SP.tile(SP.diag(muSTmuS).T + diagSigmaS,(self._D,1)), 1) 
        t3 = SP.sum( SP.dot(SW_sigma,muSTmuS0)*SW_sigma, 1)

        self.Eps.E1 = 1./((self.ZZ  + (-2*t1 + t2 + t3))/self._N) 
        self.Eps.E1[self.Eps.E1>1E6]=1E6
        #pdb.set_trace()
        
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
        self.nodes = {'S':CNodeSsparse(self),'W':CNodeWsparseVEM(self), 'Alpha':CNodeAlphasparse(self,self.priors['Alpha']['priors']),'Eps':CNodeEpsSparse(self,self.priors['Eps']['priors'])}
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
                self.W.E1[:,k] = SP.sqrt(1./self.components)*SP.randn(self._D)                
                self.S.diagSigmaS[:,k] = 1./2
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
        Bww = 0
        for k in range(self.components):
            Bww-=0.5*self.Non[k]*(special.digamma(self.Alpha.b[k])-S.log(self.Alpha.a[k]))
        #Bww= -self._D/2.0*sum(special.digamma(self.Alpha.b)-S.log(self.Alpha.a))

 #       for d in range(self._D):
 #           Bww = Bww - 1/2.0*( logdet(self.W.cov[d,:,:]) + trace(eye(self.components)-dot(self.W.E2[d,:,:],diag(self.Alpha.E1))))

        #self._bound = self._bound + Bx - Bss - Bww
        self._boundLOG.append(self._bound)
        L.debug('CVBFA bound = %.2f'%self._bound)

        return self._bound



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

    

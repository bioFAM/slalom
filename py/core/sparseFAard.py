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
        Ewdwd = 0
        Ewdwd = W.E2.sum(axis=0)
        self.a[:] = self.pa + 0.5*diag(Ewdwd)
        self.b[:] = self.pb + net.Non/2.0
        #pdb.set_trace()
        #update expectation values
        AGammaNode.update(self)

class CNodeEpsFix(CNodeEps):
    """Extensions of CNodeEps that can handle fixed prior expectaitons"""
    def __init__(self,net,prior=S.array([100,1])):
        if SP.isscalar(prior):
            self.fixE1=prior
            CNodeEps.__init__(self,net)
        else:
            self.fixE1 = None
            CNodeEps.__init__(self,net,prior)


    def update(self,net):
        if self.fixE1 is not None:
            self.E1[:] = self.fixE1
        else:
            CNodeEps.update(self,net)


class CNodeWsparse(CNodeW):
    """Abstract CnodeWsparse basclass
    - this ensures consitent handling of alternative CNodeW nodes; in particular the sparsity prior
    - main application : permutation move of factors"""
    def __init__(self,net,**kwargin):
        #call base class initialisation
        CNodeW.__init__(self,net,**kwargin)
        #create indicator arrays for moment calculation
        self.Pi = zeros([net.Pi.shape[0],net.Pi.shape[1],2])
        self.Pi[:,:,1] =net.Pi
        self.Pi[:,:,0] = 1.-net.Pi
        self.prec = S.hstack([net.sigmaOff**(-2),net.sigmaOn**(-2)])
        #log prior
        self.lpC1 = log(net.Pi+1e-10)
        self.lpC0 = log(1-net.Pi+1e-10)
        self.C    = self.Pi.copy()
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
        #do factor permutation if active
        #use the marignal indicators to calculate this; I thik they contain all we need; however we need to divide out the prior
        C = self.C/self.Pi
        Pi = self.Pi
        #normalise
        Cs = (C+1E-6).sum(axis=2)
        C[:,:,0]/=Cs
        C[:,:,1]/=Cs
        #todo: make this faster
        #now evaluate the probability of C under the (network) prior
        M = SP.zeros([net.components,net.components])
        for k in xrange(net.components):
            for l in xrange(net.components):
                M[k,l] = logProbkk(k,l)
        print "pong"
        
        #greedily select factors
        K = random.permutation(net.components)
        K = SP.arange(net.components)
        F = SP.arange(net.components)
        Ipi = SP.zeros(net.components,dtype='int')
        for k in K:
            #get beset one
            Ibest  = F[M[k,F].argmax()]
            Ipi[k] = Ibest
            #remove from list
            F  = SP.setdiff1d(F,[Ibest])
        #keep track of the changes also
        self.Ilabel = self.Ilabel[Ipi]
        #update the prior Pi
        self.Pi = self.Pi[:,Ipi,:]
        #and the precalculated log versions:
        self.lpC1 = self.lpC1[:,Ipi]
        self.lpC0 = self.lpC0[:,Ipi]
        pass



        


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

        ##analogous to the EP approximation, slower
        if 0:
            #1. calcualte the sufficient statistics of the incoming messages
            #we only store the diagonal stuff here
            E1_in  = zeros_like(self.E1)
            prec_in= zeros_like(self.cov)
            for d in range(net._D):
                pin   =  Eps.E1[d]*M
                ptotal= pin + 0
                cov   = linalg.inv(ptotal)
                prec_in[d,:,:] = pin
                E1_in[d,:]   = dot(cov, Eps.E1[d]*dot(S.E1.T,net.Z.E1[:,d]))
            self.E1_in = E1_in
            self.prec_in = prec_in
            for d in range(net._D):
                pprior = diag(1.0*self.C[d,:,1] + self.C[d,:,0]*(tauOff))
                ptotal = prec_in[d,:,:] + pprior
                ctotal = linalg.inv(ptotal)
                C      = dot(ctotal,prec_in[d,:,:])
                self.E1[d,:]    = dot(C,E1_in[d,:])
                self.cov[d,:,:] = ctotal
            CNodeW.updateE2(self,net)

        if 0:
            #update the state of the indicators
            E2 = self.E2.diagonal(axis1=1,axis2=2)
            #copy log priors
            L0 = self.lpC0.copy()
            L1 = self.lpC1.copy()
            #add log likelihood terms
            L0+= 0.5*log(tauOff/(2*pi)) - 0.5*tauOff*E2
            L1+= 0.5*log(tauOn/(2*pi))    - 0.5*tauOn*E2
            self.C[:,:,0] = exp(L0) + 1E-6
            self.C[:,:,1] = exp(L1) + 1E-6
            self.C[:,:,1]/= self.C.sum(axis=2)
            self.C[:,:,0] = 1-self.C[:,:,1]
        #call base calss for permutation moves etc.
        CNodeWsparse.update(self,net)
        pass
        






class CNodeWsparseEP(CNodeW):
    """EP mixture prior sparseness node"""

    def GaussPDF(self,x,mu,prec):
        """Gauss PDF with precisions"""
        return exp(self.LGaussPDF(x,mu,prec))

    def LGaussPDF(self,x,mu,prec):
        """Gauss PDF with precisions"""
        logPDF = 0.5*log(prec/(2*pi)) - 0.5 * prec*(x-mu)**2
        return logPDF
    
    
    def __init__(self,net,prec=1):
        CNodeW.__init__(self,net,prec)
        #create spase for responsibilities and set them to prior
        #ZE is expectation value: 
        #Z are both options for normalisation etc.
        self.C  = zeros([net.Pi.shape[0],net.Pi.shape[1],2])
        #create indicator arrays for moment calculation
        self.Pi = zeros([net.Pi.shape[0],net.Pi.shape[1],2])
        self.Pi[:,:,1] =net.Pi
        self.Pi[:,:,0] = 1-net.Pi
        self.prec = S.array([net.sigmaOff**(-2),net.sigmaOn**(-2)])
        #log prior
        self.lpC1 = log(net.Pi)
        self.lpC0 = log(1-net.Pi)        
        self.C    = self.Pi.copy()
        pass

    def update(self,net=None):
        if(net==None):
            AVGaussNode.update(self)
            return
        
        S     = net.S
        Eps   = net.Eps
        M = S.E2.sum(axis=0)
        prec_off = net.sigmaOff**(-2)

        #1. calcualte the sufficient statistics of the incoming messages
        #we only store the diagonal stuff here
        E1_in  = zeros_like(self.E1)
        prec_in= zeros_like(self.E1)
        for d in range(net._D):
            pin   =  Eps.E1[d]*M
            ptotal= pin + 0
            cov   = linalg.inv(ptotal)
            prec_in[d,:] = pin.diagonal()
            E1_in[d,:]   = dot(cov, Eps.E1[d]*dot(S.E1.T,net.Z.E1[:,d]))

        ###DEBUG
        self.prec_in = prec_in
        self.E1_in   = E1_in
        #fix retrospectively
        if 0:
            for d in range(net._D):
                prec = prec_in[d,:] + 1
                cov  = diag(prec**(-1))
                self.cov[d,:,:] = cov
                #update E
                self.E1[d,:] = 1/prec*prec_in[d,:]*E1_in[d,:]
            CNodeW.updateE2(self,net)
            return
        #calc EP sparseness prior
        #calculate EP update for W_{g,k}
        
        #temporary array
        M = zeros([E1_in.shape[0],E1_in.shape[1],3])
        #Z0: normalisation
        #Z1: mean
        #Z2: variance
        for i in range(2):
            #total precision (for product of gaussian == sum)
            prec_product = prec_in + self.prec[i]
            mean_product = 1.0/prec_product*prec_in*E1_in
            #0. calculate scale for normalisation:
            #(convolution of Gaussians)
            prec_conv     = (1.0/prec_in + 1.0/self.prec[i])**(-1)
            M0       = self.Pi[:,:,i]*(self.GaussPDF(0,E1_in,prec_conv)+1E-6)
            #M0 is also what is needed to calculate the responsibilities
            #cache: 
            self.C[:,:,i] = M0
            M1       = M0*mean_product
            #variance
            M2       = M0*(1/prec_product + mean_product**2)
            M[:,:,0]+= M0
            M[:,:,1]+= M1
            M[:,:,2]+= M2
        #now calculate final mean and variance
        #devide by overall scale
        M[:,:,1:3] = M[:,:,1:3]/M[:,:,0:1]
        self.E1 = M[:,:,1]
        #variance
        V2 = M[:,:,2]-self.E1**2
        #update moments
        for d in range(net._D):
            self.cov[d,:,:] = diag(V2[d,:])
        CNodeW.updateE2(self,net)

        #indicators:
        # pi  = 1/s * pi^{0}* \int_{w} \normal{w|w_in,\tau_in}\normal{w|0,\tau_{c}
        # so we can get those from M0 already
        #the quick way
        self.C[:,:,:] = self.C[:,:,:]/M[:,:,0:1]
        return

        ##SLOW but save:
        for i in range(2):
            prec_conv     = (1.0/prec_in + 1.0/self.prec[i])**(-1)
            #prior
            self.C[:,:,i] = self.Pi[:,:,i]
            #message from incoming Gaussian
            self.C[:,:,i]*=self.GaussPDF(E1_in,0,prec_conv) + 1E-6
        #normalize
        self.C[:,:,1]/= self.C.sum(axis=2)
        #self.C.sum(axis=2) must be the same then M[:,:,0]!
        #we can move this calculation up into the first for loop for efficiency (once it works)
        self.C[:,:,0] = 1-self.C[:,:,1]
            
        return
        
        
        #update the state of the indicators
        E2 = self.E2.diagonal(axis1=1,axis2=2)
        #copy log priors
        L0 = self.lpC0.copy()
        L1 = self.lpC1.copy()
        #add log likelihood terms
        L0+= 0.5*log(tauOff/(2*pi)) - 0.5*tauOff*E2
        L1+= 0.5*log(1.0/(2*pi))    - 0.5*1.0*E2
        self.C[:,:,0] = exp(L0) + 1E-6
        self.C[:,:,1] = exp(L1) + 1E-6
        self.C[:,:,1]/= self.C.sum(axis=2)
        if 0:
            self.C[:,:,1] = self.pC
            self.C[:,:,0] = 1-self.pC
            #the likelihood ratio from obsevations
            self.C[:,:,0]*= exp(0.5*log(tauOff/2*pi))*exp(-0.5*tauOff*E2)
            self.C[:,:,1]*= exp(0.5*log(1.0/2*pi))*exp(-0.5*1.0*E2)
            self.C[:,:,1]/= self.C.sum(axis=2)
            self.C[:,:,0] = 1-self.C[:,:,1]
        #print isnan(self.E1).any()
        #print isnan(self.E2).any()
        #print isnan(self.C[:,:,1]).any()
        pass





class CNodeWsparseEPV(CNodeWsparse):
    """EP mixture prior sparseness node in vectorised form"""

    def GaussPDF(self,x,mu,prec):
        """Gauss PDF with precisions"""
        return exp(self.LGaussPDF(x,mu,prec))

    def LGaussPDF(self,x,mu,prec):
        """Gauss PDF with precisions"""
        logPDF = 0.5*log(prec/(2*pi)) - 0.5 * prec*(x-mu)**2
        return logPDF
    
    
    def __init__(self,net,prec=1):
        CNodeWsparse.__init__(self,net,prec=prec)
        #all variables needed are initialised in CNodeWsparse.__init__
        self.g    = zeros([2,net._D,net.components])
        self.damp = ones( [net._D,net.components])

        self.iterative=net.iterative
        pass


       

    def update(self,net=None):
        if(net==None):
            AVGaussNode.update(self)
            return

        def n2mode(x):
            """conversion nat. parameter and back"""
            return array([x[0]/x[1],1/x[1]])

        def n2modeV(x):
            """conversion nat. parameter and back"""
            rv = zeros_like(x)
            rv[0] = x[0]/x[1]
            rv[1] = 1.0/x[1]
            return rv

        #define convenience variables
        S     = net.S
        Eps   = net.Eps
        M = S.E2.sum(axis=0)
 
        #1. calcualte the sufficient statistics of the incoming messages
        #get full covariances from the VB approximation
        E1_in  = zeros_like(self.E1)
        prec_in= zeros_like(self.cov)
        cov_in = zeros_like(self.cov)
        PE    =  zeros_like(self.E1)
        for d in range(net._D):
            pin   =  Eps.E1[d]*M
            ptotal= pin + 0
            cov   = linalg.inv(ptotal)
            cov_in[d,:,:] = cov
            prec_in[d,:,:] = pin
            E1_in[d,:]     = dot(cov, Eps.E1[d]*dot(S.E1.T,net.Z.E1[:,d]))
            #also calculate mean prec*e1_in
            PE[d,:] = dot(prec_in[d,:,:],E1_in[d,:])

        #1. update W
        #updates of W factorise in Wd slices, but do this updates in parallel
        t  = E1_in
        K  = cov_in
        #1. set up approximation
        Sigma = K.copy()
        KI    = prec_in

        #2. the local site functions
        #note: these are cached between iterations
        g     = self.g
        if 1:
            g[:,:,:] = 0
        mu = zeros([net._D,net.components])
        for d in range(net._D):
            Sigma[d,:,:] = linalg.inv(KI[d,:,:] + diag(g[1,d,:]))
            mu[d,:]      = dot(Sigma[d],g[0,d,:]) + dot(Sigma[d],PE[d,:])
        #current approx. mean, dot as this can be cached..
        self.damp[:,:] = 1.0
        g2 = g.copy()
        
        #update iterations
        ##now iterate through all dimensions, only perform one single update
        for niter in range(1):
            if net.randperm:
                perm = random.permutation(net.components)
            else:
                perm = arange(net.components)
            for k in perm:
                #cavity as natural parameter representation
                cav_np = n2modeV(array([mu[:,k],Sigma[:,k,k]]))-g[:,:,k]
                #ensure we don't have negative variances. good idea?
                cav_np[1] = abs(cav_np[1])
                cav_mp = n2modeV(cav_np)

                #1. update W
                #calc expectation values cavity distribution*sparseness prior
                M = zeros([3,net._D])
                for i in range(2):
                    #prec for product (cav_np[1] is prec)
                    prec_product = cav_np[1]+self.prec[i]     # prec_cav + prec_prior(c) 
                    mean_product  = 1.0/prec_product*cav_np[0] # 1/prec_product*(prec_cav*mean_cav)
                    #convoultion parameters
                    prec_conv    = (1.0/cav_np[1] + 1.0/self.prec[i])**(-1)  # (1/prec + 1/prec)**-1
                    #scale
                    M0 = self.Pi[:,k,i]*(self.GaussPDF(0,cav_mp[0],prec_conv)+1E-6)
                    #keep scale for calculation of indicator
                    self.C[:,k,i] = M0
                    #1. moment
                    M1 = M0*mean_product
                    M2 = M0*(1/prec_product + mean_product**2)
                    M[0]+=M0
                    M[1]+=M1
                    M[2]+=M2
                    if isnan(M).any():
                        pdb.set_trace()
                        pass
                #calc resulting approx. for site function
                mu_g = M[1]/M[0]
                v2_g = M[2]/M[0]-mu_g**2
                #update the site function
                gn = n2modeV(array([mu_g,v2_g]))-cav_np
                #TODO: include damping etc.
                #variance difference for rank-1 updates
                ds2= gn[1]-g[1,:,k]
                dg = gn-g[:,:,k]
                g[:,:,k] = g[:,:,k] + self.damp[:,k]*dg
#                g[:,:,k] = gn
                #2. update Z
                #this is simply normalisation now
                self.C[:,k,:] = self.C[:,k,:]/M[0:1].T

                if not self.iterative:
                    #if not iterative, don't update Sigma
                    continue
                #TODO: rank-one updates
                Sigma2= Sigma
                #updates of sigma and mu again iterated over all elements
                for d in range(net._D):
                    if 1:
                    #rank one update
                        Sigma[d] = Sigma[d] - ds2[d]/(1+ds2[d]*Sigma[d,k,k])*outer(Sigma[d,:,k],Sigma[d,k,:])
                        if 0:
                            try:
                                Csigma = linalg.cholesky(Sigma[d])
                            except linalg.LinAlgError:
                                print "outch"
                                Sigma[d]=Sigma2[d]
                                g[:,d,k] = g2[:,d,k]
                    else:
                        Sigma[d] = linalg.inv(KI[d] + diag(g[1,d,:]))
                    #update mu:
                    mu[d]   = dot(Sigma[d],g[0,d,:]) + dot(Sigma[d],PE[d])
            #::end for k
            
            #do full updates at the end of every sweep
            #if no iterative learning, update covariances here
            for d in range(net._D):
                Sigma[d] = linalg.inv(KI[d] + diag(g[1,d,:]))
                #update mu:
                mu[d]   = dot(Sigma[d],g[0,d,:]) + dot(Sigma[d],PE[d])
            #update copy of sites
            g2 = g.copy()
        #update approximation in data structure
        self.E1 = mu
        self.cov = Sigma
        CNodeW.updateE2(self,net)
        #call the update of the CsparseFA
        CNodeWsparse.update(self,net)
        return





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
        dp['schedule'] = ['S','Eps','W','Alpha']
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
                
        if tolerance is None: tolerance = self.tolerance
        if nIterations is None: nIterations = self.nIterations
        if forceIterations is None: forceIterations = self.forceIterations
        LB = 0

        for iter in range(nIterations):
            self.iterationCount+=1
            t = time.time();
            for node in self.schedule:
                self.updateNode(node)
            #self.calcBound()
            #calc reconstruction error
            Zr = S.dot(self.S.E1,self.W.E1.T)
            Zd = self.Z.E1-Zr
            error = (Zd**2).mean()
            self.calcBound()
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
        self.nodes = {'S':CNodeS(self),'W':W_node,'Eps':CNodeEps(self,self.priors['Eps']['priors']), 'Alpha':CNodeAlphasparse(self,self.priors['Alpha']['priors'])}
        for n in self.nodes.keys(): setattr(self,n,self.nodes[n])

        #pca initialisation
        Ion = None
        if self.initType == 'pca':
            Ion = random.rand(self.Pi.shape[0],self.Pi.shape[1])<self.Pi
            for k in range(self.components):
                sv = linalg.svd(self.Z.E1[:,Ion[:,k]], full_matrices = 0);
                [s0,w0] = [sv[0][:,0:1], S.dot(S.diag(sv[1]),sv[2]).T[:,0:1]]
                v = s0.std(axis=0)
                s0 /= v;
                w0 *= v;
            self.S.E1[:,k] = s0.ravel()
            self.W.E1[Ion[:,k],k] = w0.ravel()
            self.W.E1[~Ion[:,k],k]*=self.sigmaOff
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
                self.W.E1[Ion[:,k]]*=self.sigmaOn[k]
        elif self.initType == 'on':
            for k in range(Ion.shape[1]):
                self.W.E1[:,k]*=self.sigmaOn[k]
        elif self.initType == 'data':
            assert ('S' in keys(init_factors))
            assert ('W' in keys(init_factors))
#            Ion = init_factors['Ion']
            Sinit = init_factors['S']
            Winit = init_factors['W']
            for k in range(self.components):
                self.S.E1[:,k] = Sinit[:,k]
                self.W.E1[:,k] = Winit[:,k]
#                self.W.E1[~Ion[:,k],k]*=self.sigmaOff

        #update moments
        if 0:
            self.S.update()
            self.W.update()
            #self.W.updateE2(self)
        if 1:
            #do we have data to initialise the cov. of W?
            if Ion is not None:
                for d in xrange(self._D):
                    for k in range(Ion.shape[1]):
                        self.W.E1[Ion[:,k]]*=self.sigmaOn[k]
                    pprior = diag(self.sigmaOn**(-2)*Ion[d,:]+(~Ion[d,:])*self.sigmaOff**(-2))
                    self.W.cov[d,:,:] = linalg.inv(pprior)
                pass
            self.S.update()
            self.W.update()
            #self.Eps.update(self)


    #calculate the variational bound:
    def calcBound(self):
        L.debug('CVBFA calcBound')
        self._bound = ABayesNet.calcBound(self)

        
        
        #p(data|..)
        #OLI: is this right? the last term should be <-1/2*tau(D-W*x)^{2}>
        #try: here we recyle the calculation made in the update Eps:
        Bx = -self._N*self._D/2.0*S.log(2*pi) + self._N/2.0*self.Eps.E2.sum() - sum(self.Eps.E1*(self.Eps.a-self.Eps.pa))
        #Bx = -self._N*self._D/2.0*S.log(2*pi) + self._N/2.0*self.Eps.E2.sum() + sum(self.Eps.E1*self.Eps.pa-self.Eps.b)


        #note : trace (S.cov) comes from the fact that the 2nd moment of S is not just S.E1**2 but + cov!
        #KL q(S)/P(S)
        
        #orig
        Bss= -self._N/2.0*logdet(self.S.cov) - self._N/2.0*trace(eye(self.components)-self.S.cov) + 0.5*(self.S.E1**2).sum()

        #KL q(W)/p(W|alpha)
        Bww = 0
        for k in range(self.components):
            Bww-=0.5*self.Non[k]*(special.digamma(self.Alpha.b[k])-S.log(self.Alpha.a[k]))
        #Bww= -self._D/2.0*sum(special.digamma(self.Alpha.b)-S.log(self.Alpha.a))

        for d in range(self._D):
            Bww = Bww - 1/2.0*( logdet(self.W.cov[d,:,:]) + trace(eye(self.components)-dot(self.W.E2[d,:,:],diag(self.Alpha.E1))))

        self._bound = self._bound + Bx - Bss - Bww
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
        Z = vbfa.W.C[:,:,1]
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

    

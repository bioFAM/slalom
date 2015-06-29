##BayesNet
# abstract Basis classes for Variational Bayesian models:

import scipy as S
import logging as L
import numpy.random as random
import numpy.linalg as linalg
import scipy.special as special
import time
import pdb



##TODO
# 1) we should agree whether node updates take in E1/E2 or not
#    currently we have a wild mix (see ABernoulli, AGauss)
# 2) message are current Gauss; separate these out form nodes?


def savelog(x):
    return S.log(x+1E-10)

class ANode(object):
    ''' Abstract base class for all components of a probabilistic approximation model '''


    def __init__(self):
        self._bound = 0;


    # entropy of the approximating distribution
    def entropy(self):
        L.warn("Entropy not implemented")
        return None

    # calculate variational lower bound
    def calcBound(self, net):
        L.warn("Bound not implemented");
        return self._bound;


    # update the node
    def update(self):
        L.warn("Update not implemented!")


class ADirichletNode(ANode):
    ''' Class for Dirichlet node for mixture models
    The prior is the concentration parameter '''


    def __init__(self, dim=[1], prior=0.1):
        ''' Initialise dimension and prior of the node ''' 
        L.debug('ADirichletNode __init__')
        ANode.__init__(self)
        self._prior = prior
        self._a0= S.array(prior)
        self.a = S.zeros(dim)
        self.E1= S.zeros(dim)
        self.lnE=S.zeros(dim)

        self.a[:] = self.a0
        ADirichletNode.update(self)

    def update(self):
        ''' update expectation value '''
        L.debug('ADirichletNode update')
        S = self.a.sum(axis=-1) # add up responsibilities
        shape = list(S.shape)
        shape.append(1)
        self.E1 = self.a*(1/S).reshape(shape) # normalise over last dimension to calculate E1 and lnE
        self.lnE= special.psi(self.a)-special.psi(self.a.sum(axis=-1)).reshape(shape)

        
    def calcBound(self, net):
        L.debug('ADirichletNode calcBound')
        return ((self.a0 - 1)*self.lnE).sum() + special.gammaln(self.a0.sum()) - special.gammaln(self.a0).sum() + self.entropy()

        
    def entropy(self):
        L.debug('ADirichletNode entropy')
        return -special.gammaln(self.a.sum()) + special.gammaln(self.a).sum() + (self.a.sum() - self.a.shape[0])*special.digamma(self.a.sum()) - ((self.a - 1)*special.digamma(self.a)).sum()



class ABernoulliNode(ANode):
    ''' Abstract base class for Bernoulli node for mixture models '''
    def __init__(self, dim=[2],prior=S.array([0.5,0.5]), E1=None, init='prior'):
        L.debug('ABernoulliNode __init__')
        ANode.__init__(self)        
        self._prior = prior
        
        if E1 is None:
            if init == 'prior': E1 = prior[1]*S.ones(dim)
            elif init == 'rnd': E1 = random.random_sample(dim)
            elif init == 'rndass': # random assignment
                E1 = S.zeros(dim);
                for i in range(dim[0]): E1[i,random.randint(0,dim[1])] = 1            
            elif init == 'rndsoftass': # random soft assignment
                E1 = S.ones(dim)/(10.*(dim[1] - 1));
                for i in range(dim[0]): E1[i,random.randint(0,dim[1])] = 0.9
            self.E1 = E1
            self.E2 = E1

    #update simply normalises the distribution.
    #Assume normalise around dimension 0
    def update(self,P=None):
        """update(P=None)
        update ABernoulliNode
        - if P is assumed to contain unnormalized ratios"""
        L.debug('ABernoulliNode update')
        if P is not None:
            S = P.sum(axis=-1)
        shape = list(S.shape)
        shape.append(1)
        P *= (1/(S+1E-10)).reshape(shape)
        #expectatoin value is the p(1) which is assumed to tbe in the last column; hence this construction...
        self.E1 = P.swapaxes(0,-1)[1]
        self.E2 = self.E1
        
    
    def entropy(self):
        L.debug('ABernoulliNode entropy')
        #return -(self.E1*S.log(self.E1 + 1E-10)).sum()
        return -(self.E1*savelog(self.E1) + (1-self.E1)*savelog(1-self.E1))


    def calcBound(self,net):
        L.debug('ABernoulliNode calcBound')
        return (self._prior[1]*self.E1 + self._prior[0]*(1-self.E1) + self.entropy()).sum()
        



class AGammaNode(ANode):


    def __init__(self, dim=[1], prior=[1E-3,1E-3], E1=None):
        L.debug('AGammaNode __init__')
        ANode.__init__(self)

        #save prior:
        self._prior = prior
        self.pa = prior[0]
        self.pb = prior[1]
        # self._dim = dim

        self.a  = S.ones(dim)*self.pa
        self.b  = S.ones(dim)*self.pb

        AGammaNode.update(self) # updates E1, E2

        # Manually set E1 if needed
        if(E1 is not None): self.E1 = E1


    def update(self):
        L.debug('AGammaNode update')

        self.E1 = self.b/self.a
        self.E2 = -S.log(self.a) + special.digamma(self.b)
        self.lnE = special.digamma(self.a) - S.log(self.b)


    def entropy(self):
        L.debug('AGammaNode entropy')
        #I think this is wrong
        #return (self.a - S.log(self.b) + special.gammaln(self.a) + (1. - self.a)*special.digamma(self.a)).sum()
        return (self.b - S.log(self.a) + special.gammaln(self.b) + (1-self.b)*special.digamma(self.b)).sum()
    


    def calcBound(self, net):
        L.debug('AGammaNode calcBound')

        #Ithink this is wrong
        #b1=self.b.shape[0]*self.pa*S.log(self.pb) + ((self.pa - 1)*(special.digamma(self.a) - S.log(self.b)) - self.pb*(self.a/self.b) - special.gammaln(self.pa)).sum() + self.entropy()

        b1= (-special.gammaln(self.pb) + self.pb*S.log(self.pa) + (self.pb-1)*(S.log(self.b)-S.log(self.a)) - self.pa*(self.b/self.a)).sum() + self.entropy()
        #b2=(-self.b*S.log(self.a) + self.pb*S.log(self.pa) + special.gammaln(self.b) - special.gamma(self.pb) - self.pa*(self.E1) + self.b -(self.b-self.pb)*(special.digamma(self.b)-S.log(self.a))).sum()
        return b1
        pass





class AVGaussNode(ANode):
    ''' Vector-Gauss-Node which sits in two plates
    dim is supposed to be (N,d) where d is the vectorized gauss-node dimension and N is the data-dimension'''



    def __init__(self,dim=[1,1],cdim=S.nan,prior=[0,1]):
        L.debug('AVGaussNode __init__')
        ANode.__init__(self)        
#        self._dim = dim
#        self._cdim = cdim
        self._prior = prior
        if(S.isnan(cdim)):  cdim = dim[0]

        self.E1 = prior[0] + random.standard_normal(size=dim)
        self.E2 = S.zeros([dim[0],dim[1],dim[1]])
        self.cov= S.zeros([cdim,dim[1],dim[1]])

        for c in range(cdim): self.cov[c,:,:] = prior[1]*S.eye(dim[1])
        AVGaussNode.update(self)


    def update(self):
        L.debug('AVGaussNode update')
        self.E2[:,:,:] = self.cov
        for n in range(self.E2.shape[0]):
            self.E2[n,:,:] = self.E2[n,:,:] + S.outer(self.E1[n,:],self.E1[n,:])

    def calcBound(self, net):
        return 0



class AGaussNode(ANode):
    ''' Completely separated Gauss node
        We usually use this representation to model the dataset itself '''



    def __init__(self,node=None,dim=[1,1],prior=[0,1],E1=None,E2=None,prec=None):
        ANode.__init__(self)
        #copy constructor
        if node is not None:
            self.E1 = node.E1.copy()
            self.E2 = node.E2.copy()
            self.cov= node.cov.copy()
            return
        self._prior=prior
        if E1 is None:  E1 = prior[0]*random.standard_normal(size=dim)
        AGaussNode.update(self, E1=E1, E2=E2, prec=prec)

        
    def update(self,E1,E2=None,cov=None,prec=None):
        """update(E2,E2,cov,prec)
        -updaes the state of the Gaussnode from E1 and a representation of the second moment"""
        self.E1 = E1
        #make sure this works also without any second moment
        if (E2 is None) and (cov is None) and (prec is None):
            cov = S.array(1E-6)
        if E2 is not None:
            self.E2 = E2
            self.cov = self.E2 - self.E1**2
        elif cov is not None:
            self.cov = cov
            self.E2 = E1**2 + cov
        elif prec is not None:
            self.cov = 1.0/prec
            self.E2 = E1**2 + self.cov


    def getPrec(self):
        '''return the precision of this node from current moments'''
        return 1.0/self.cov

    def getMean(self):
        '''return the mean'''
        return self.E1

    def getVariance(self):
        '''return the variance'''
        return self.cov
        

class CGauss(AGaussNode):
    '''class for general Gaussian messages
    commonly this is used to represent predictions and messages between expressionmodule'''
    def __init__(self,E1=None,E2=None,cov=None,prec=None):
        if E1 is not None:
            #in caes the CGauss is only constructed from mean
            self.update(E1=E1,E2=E2,cov=cov,prec=prec)

    
class ABayesNet(ANode):
    ''' Abstract Bayes net class. Consists of nodes, schedule of updates, tolerance for convergence, number of iterations'''
    


    def getDefaultParameters(self):
        '''getDefaultParameters()
        - return a dictionary with defaultparameters of the BayesNet class'''
        dp = {}
        dp['nIterations'] = 1
        dp['tolerance']   = 1E-4
        dp['forceIterations'] = False
        dp['schedule'] = []
        return dp

    def __init__(self,parameters=None,nodes={}, schedule=[]):
        ''' Initialise the Bayes net. Requires a map of names to updatable objects (of type AUnit, called nodes for conciseness), and an update schedule of the node names
        - this also takes care of default parameter handling etc.
        '''
        L.debug('ABayesNet init')
        #0. set default parameters
        dp = self.getDefaultParameters()
        #update with specified parameters
        dp.update(parameters)
        #set all parameters as member variables
        for param in dp.keys():
            self.__setattr__(param,dp[param])
        
        #set a default tollerance for convergence
        if len(self.schedule)==0:
            self.schedule = nodes.keys()
        self.nodes = nodes
        self.iterationCount = 0
        self._bound = 0
        self._boundLOG = []


    def update(self,**kwargs):
        '''update of this node per'''
        L.debug('ABayesNet update')
        self.iterate(**kwargs)
        

    def updateNode(self,node):
        """update node (node)
        - this is mainly to allow flexibility by overriding this method"""
        if node in self.nodes.keys():
            #the update has the net as an argument (self!)
            self.nodes[node].update(self)
        else:
            raise Exception("node %s in schedule but not in nodes" % node)

        
    
    def iterate(self, nIterations=None, forceIterations=None):
        '''iterate(nIteations=None,forceIterations=None)
        - perform nIterations; per default(None) parameters are tken from local intsance settings
        '''
        
        L.debug('ABayesNet iterate')
                
        if nIterations is None: nIterations = self.nIterations
        if forceIterations is None: forceIterations = self.forceIterations
        LB = 0

        for iter in range(nIterations):
            self.iterationCount+=1
            t = time.time();
            for node in self.schedule:
                self.updateNode(node)
            self.calcBound()

            if (abs(LB - self._bound) < self.tolerance) and not forceIterations:
                L.info('Converged')
                break

            L.info("Iteration %d: time=%.2f bound=%f" % (iter,time.time() - t, self._bound))
            LB = self._bound

        return self._bound


    def meanLogProb(self):
        L.debug('ABayesNet meanLogProb')
        return 0;



    def calcBound(self):
        L.debug('ABayesNet calcBound')
        self._bound = self.meanLogProb()
        for node in self.nodes:
            self._bound += getattr(self,node).calcBound(self)

        return self._bound

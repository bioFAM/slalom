# Copyright(c) 2016, The f-scLVM developers (Florian Buettner, Oliver Stegle)
#
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

# VBFA
# variational BAyesian Factor Analysis

from numpy import *
import numpy.random as random
import numpy.linalg as linalg
import scipy.special as special
import sys
import pickle
import logging as L
import scipy as S
import pdb
#use everything from BayesNet
from .bayesnet import *
from .expressionnet import AExpressionModule


#node which captures the data which might have a 1./2. moment as well:
class CNodeZ(AGaussNode):
    pass
#    def __init__(self,E1, E2 = None, cov = None, prec = None):
#        AGaussNode.update(self, E1, E2, cov, prec)



class CNodeS(AVGaussNode):
    def __init__(self,net,prec=1):
        AVGaussNode.__init__(self,dim=[net._N,net.components],cdim=1)

    def update(self,net=None):
        """update(net):
        where net is the "base networ, i.e. sparseFA or similar"""
        if(net==None):
            AVGaussNode.update(self)
            return

        W = net.W
        Eps=net.Eps
        p=eye(net.components)+tensordot(Eps.E1,W.E2,[0,0])
        self.cov = linalg.inv(p)
        matrix   = dot( dot(self.cov,W.E1.T), diag(Eps.E1) )
        self.E1  = dot(net.Z.E1,matrix.T)
        AVGaussNode.update(self)
        pass



class CNodeW(AVGaussNode):
    def __init__(self,net,prec=1):
        AVGaussNode.__init__(self,dim=[net._D,net.components])
        #obsolete due to update in Eps
        #self.E2W = zeros([net._D,net._N])

    def updateE2(self,net=None):
        if net is None:
            AVGaussNode.update(self)
            return
        self.E2[:,:,:] = self.cov

        Ss = net.S.E2.sum(axis=0).T
        
        for d in range(self.E2.shape[0]):
            self.E2[d,:,:] = self.E2[d,:,:] + S.outer(self.E1[d,:],self.E1[d,:])
        #this is obsolete due to update in Eps
        #self.E2W_[d] = (self.E2[d,:,:]*Ss).sum()
        pass



    def update(self,net=None):
        if(net==None):
            
            self.updateE2()
            return
        
        S     = net.S
        Alpha = net.Alpha
        Eps   = net.Eps
        M = S.E2.sum(axis=0)
            
        for d in range(net._D):
            p = diag(Alpha.E1) + Eps.E1[d]*M
            self.cov[d,:,:] = linalg.inv(p)
            self.E1[d,:] = dot(self.cov[d,:,:],Eps.E1[d]*dot(S.E1.T,net.Z.E1[:,d]))
            
        #for 2. moment-calculation
        #AVGaussNode.update(self)
        self.updateE2(net)



class CNodeEps(AGammaNode):
    def __init__(self,net,prior=[100,1]):
        AGammaNode.__init__(self,dim=[net._D],prior=prior)

    def update(self,net):
        S   = net.S
        W   = net.W
        Z   = net.Z.E1
        # take second moment into account.            
        self.b[:] = self.pb + net._N/2.0
        #set a:
        TD = Z*tensordot(S.E1,W.E1,[1,1])
        #linearise in dimension space to handle high dimensions better:
        if 0:
            TD2= zeros([net._D,net._N])
            for d in range(net._D):
                TD2[d,:] = tensordot(W.E2[d,:,:],S.E2,([1],[1])).trace(axis1=0,axis2=2).sum()
        else:
            Ss = net.S.E2.sum(axis=0).T
            TD2 = (W.E2*Ss).sum(axis=1).sum(axis=1)


        t  = net.Z.E2.sum(axis=0) - 2*TD.sum(axis=0) + TD2     
        
        self.a[:] = self.pa + 0.5*t
            #the non-compressed version of this is below:
            #for d in range(net._D):
            #    t=0
            #    for n in range(net._N):
            #        t = t +Z[n,d]**2 -2*Z[n,d]*dot(W.E1[d,:],S.E1[n,:]) + trace(dot(W.E2[d,:,:],S.E2[n,:,:]))
            #    self.a[d] = self.a[d] + 0.5*t
            #

            #update expectation values:
        AGammaNode.update(self)
            


class CNodeAlpha(AGammaNode):
    def __init__(self,net,prior=[1E-3,1E-3]):
        AGammaNode.__init__(self,dim=[net.components],prior=prior)
        
    def update(self,net):
        W = net.W
        Ewdwd = 0
        Ewdwd = W.E2.sum(axis=0)
        self.a[:] = self.pa + 0.5*diag(Ewdwd)
        self.b[:] = self.pb + net._D/2.0
        #update expectation values
        AGammaNode.update(self)

          

class CVBFA(AExpressionModule):
    '''CVBFA(AExpressionModule)
    - Variational Bayesian Factor analysis module'''


    def getDefaultParameters(self):
        """return a hash with default parameter value for this BayesNet"""
        #
        dp = AExpressionModule.getDefaultParameters(self)
        dp['initType'] = 'pca'
        dp['nIterations'] = 20
        dp['schedule'] = ['S','W','Alpha','Eps']
        dp['components'] = 5
        dp['priors']  = {}
        dp['name_str'] = {}
        return dp



    def __init__(self,init_data=None,E1=None,E2=None,**parameters):
        """create the object"""
        #handle setting of parameters via Bayesnet constructor
        ABayesNet.__init__(self,parameters=parameters)
        #priors for the various components:
        if('Alpha' not in self.priors): self.priors['Alpha']={'priors': [1E-3,1E-3]}
        if('Eps' not in self.priors):   self.priors['Eps']={'priors': [1,100]}
        self.dataNode=None
        if init_data is None and E1 is not None:
            init_data = CGauss(E1=E1,E2=E2)
        if init_data is not None:
            self.init(init_data)


    def init(self,init_data,Pi=None):
        if not isinstance(init_data,AGaussNode):
            raise Exception("initialization is only possible from a GaussNode")
        self.Z = CNodeZ(node=init_data)
        self.dataNode = self.Z

        # set dimensionality of the data
        [self._N, self._D] = self.Z.E1.shape
        

        #add the new nodes - to be replaced by XML init:
        self.nodes = {'S':CNodeS(self),'W':CNodeW(self),'Eps':CNodeEps(self,self.priors['Eps']['priors']),'Alpha':CNodeAlpha(self,self.priors['Alpha']['priors'])}
        for n in list(self.nodes.keys()): setattr(self,n,self.nodes[n])

        #pca initialisation
        if self.initType == 'pca':
            sv = linalg.svd(self.Z.E1, full_matrices = 0);
            [s0,w0] = [sv[0][:,0:self.components], S.dot(S.diag(sv[1]),sv[2]).T[:,0:self.components]]
            v = s0.std(axis=0)
            s0 /= v;
            w0 *= v;
            self.S.E1 = s0
            self.W.E1 = w0
            self.W.update()
            self.S.update()
        else:
            print("random init")
            self.S.E1 = random.randn(self._N,self.components)
            self.W.E1 = random.randn(self._D,self.components)
            self.S.update()
            self.W.update()
            self.W.updateE2(self)
        self._K = net.W.E1.shape[1]


    def getName(self):
        """return a name summarising the  main parameters"""
        name = "VBFA: %s C=%d" % (self.name_str,self.components)
        return name

    def iterate(self, nIterations=None, forceIterations=None):
        '''iterate(nIteations=None,forceIterations=None)
        - perform nIterations; per default(None) parameters are tken from local intsance settings
        '''
        forceIterations=True
        L.debug('SparseFA iterate')
                
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
            error = ((Zd)**2).mean()
            print("reconstruction error: %f" % (error))
            

            if (abs(LB - self._bound) < self.tolerance) and not forceIterations:
                L.info('Converged')
                break

            L.info("Iteration %d: time=%.2f bound=%f" % (iter,time.time() - t, self._bound))
            LB = self._bound

        return self._bound




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

        #KL q(W)/p(W)
        Bww= -self._D/2.0*sum(special.digamma(self.Alpha.b)-S.log(self.Alpha.a))

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

        

def logdet(M):
    UC = linalg.cholesky(M)
    return 2*sum(S.log(diag(UC)))

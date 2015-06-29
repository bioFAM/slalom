"""expressionNet
- Aditive models for Gene Expression

AExpressionModule(ABayesNet)
- prototype class for any exprssoin Module used by CSumExpressionNet

CSumExpressionNet(AExpressionModule)
- Implementation of ExpressionModule which models gene expression data by a sum of influences
"""

from bayesnet import *
import pdb
import logging as L
import time
import mxml


class AExpressionModule(ABayesNet):
    ''' A determinant of gene expression; has a Gaussian node that can make predictions and a Gamma precision node for the Gaussian
    The update schedule first updates the mean of the Gaussian node based on the incoming message, then everything else, and then Gaussian node again to prepare the outgoing message. The data node and variance node are updated regardless of the rest of the specified schedule
    -dataNode holds the local datanode which is overwritten by updateDataNode
    '''
    

    
    def __init__(self, E1, E2):
        L.debug('AExpressionModule __init__')
        # read initial data, initialise net - default behaviour is to create and populate data node
        self.dataNode = AGaussNode(E1=E1, E2=E2)


    def updateDataNode(self, message):
        if not isinstance(message,AGaussNode):
            raise Exception('updateDataNode requires a message of type AGaussNode')
        L.debug('updateDataNode')

        #0. calculate the current pseudo dataset from up and down messages
        prediction = self.getPrediction()
        #store E1 and prec
        pE1  = prediction.getMean()
        pPrec= prediction.getPrec() 
        mE1  = message.getMean()
        mPrec= message.getPrec()
        precNew = pPrec + mPrec
        E1New   = (pE1*pPrec + mE1*mPrec)/precNew
        self.dataNode.update(E1=E1New,prec=precNew)


    def getPrediction(self):
        L.warn('AExpressionModule: getPrediction needs to be implemented')
        return None



###############CSumExpressionNet#################

class CNodeExpr(AGaussNode):
    def __init__(self, **kwargin):
        AGaussNode.__init__(self,**kwargin)

class CNodeEps(AGammaNode):
    def __init__(self, **kwargin):
        AGammaNode.__init__(self,**kwargin)


    def update(self,net):
        '''update function for noise level of expressionnet'''
        #data
        data = net.dataNode
        #the sum from all modules
        msum  = net.sumDataNode
        #calculate residuals using uncertainties:
        t = data.E2 - 2*data.E1*msum.E1 + msum.E2
        t = t.sum(axis=0)
        self.a[:] = self.pa + 0.5*t
        self.b[:] = self.pb + net.N/2.0
        #update expectation values etc.
        AGammaNode.update(self)
        
        
        


class CNodeSum(AGaussNode):
    '''SumNode handling the incoming messages from all the ExpressionModule'''
    def __init__(self,**kwargin):
        AGaussNode.__init__(self,**kwargin)

    def update(self,net=None):
        '''update()
        -update the sum node based on the current state of all expression modules'''
        mean = 0
        variance = 0
        for module in net.expressionModules:
            prediction = net.expressionModules[module].getPrediction()
            mean = mean + prediction.E1
            variance=variance + prediction.getVariance()
        #update the GaussNode
        AGaussNode.update(self,E1=mean,cov=variance)


    
class CSumExpressionNet(AExpressionModule):
    ''' Class for a Bayesian linear model of gene expression'''


    def __init__(self, initData=None,xml=None,**parameters):
        L.debug('CExpressionNet __init__')

        if xml is None:
            L.error('CSumExpressionNet can only be constructed from an CXML object')
        self.xml = xml

        self.expressionModules = {}
        #get xml parameters
        xml_parameters = self.xml.getParameters()
        #setting of parameters and handling of defaultparameters is done via BayesNet constructor
        ABayesNet.__init__(self,parameters=xml_parameters)

        #dimension of expression levels
        self.N    = self.expr.shape[0]
        self.D    = self.expr.shape[1]
        L.info('Loaded expression data set %dx%d' % (self.N,self.D))

        #create nodes
        #get node definition from xml
        nodes = self.xml.getNodes()
        #we could cycle through the nodes but this might be overkill.
        #let's just get the node parameters we care about
        self.Eps       = CNodeEps(dim=self.D,prior=nodes['Eps']['priors'])
        self.dataNode  = AGaussNode(dim=[self.N,self.D])
        #Data for SumExpressionNet is observed and hence arb. precise
        self.dataNode.update(E1=self.expr,prec=1E6)

        #EPS node is updated...
        self.nodes['Eps'] = self.Eps
        #set expressoin variables as very certain dataNode
        ##create expressionModules
        modules = self.xml.getModels()
        #update CNodeSum
        self.sumDataNode = CNodeSum()
        for module in modules:
            name           = module.getAttribute('name')
            class_name     = module.getAttribute('class')
            class_instance = mxml.CXml(xml=module).createModelClass()
            self.expressionModules[name] = class_instance
        #initialize modules in the order specfied in schedule
        for node in self.schedule:
            if node not in self.expressionModules.keys():
                continue
            module = self.expressionModules[node]
            init_data      = self.messageModule(node)
            module.init(init_data)
            #perform one iteration to update bound and be able to produce valid predictions
            module.update(nIterations=1)
            #update data Node now including the initialization of the newly added module
            self.sumDataNode.update(self)
        pass

    def messageModule(self,module):
        '''update the dataNode of a particular module'''
        # this function assumes that the dataNode of the expressionNet is up2date and we hence can
        # calculate the residual message going to a particular node
        # m_E2: from noise node
        E1 = self.dataNode.E1.copy()
        E1-= self.sumDataNode.E1
        # add the mdoules contribution again as we sent the messag to this module
        E1+=self.expressionModules[module].getPrediction().E1
        #take expectatin value from Eps and repeat to match shape
        prec = S.repeat(self.Eps.E1.reshape([1,-1]),E1.shape[0],axis=0)
#        print "warning: no 2nd moement in messageModule"
        message = CGauss(E1=E1,prec=prec)
#        message = CGauss(E1=E1)
        
        return message

    def meanLogProb(self):
        """calculate man log prob for the variational bound"""
        L.debug('ABayesNet meanLogProb')
        #the constant factors from the N x D gaussians
        mlm = 0.5*self.N*S.log(self.Eps.E2) - 0.5*self.N*S.log(2*S.pi)
        #now the expectation value of the squred diffrence: this has alrady been calculated in the Eps node!
        mlm = mlm + (self.Eps.a-self.Eps.pa)
        return mlm.sum()

    def calcBound(self):
        self._bound = ABayesNet.calcBound(self)
        for module in self.expressionModules:
            self._bound += self.expressionModules[module]._bound
        L.debug('ExpressionNet bound = %.2f'%self._bound)
        self._boundLOG.append(self._bound)
        return self._bound


    def updateNode(self,node,calc_bound=False,**kwargs):
        """override updateNode method as we have to treat expressionModules and normal nodes
        separately"""
        if node in self.expressionModules:
            #expressionModules:
            data_message = self.messageModule(node)
            #send upward message
            self.expressionModules[node].updateDataNode(data_message)
            #updat the module
            self.expressionModules[node].update(**kwargs)
            #always update sumNode as it is a helper only
            self.sumDataNode.update(self)
        elif node in self.nodes.keys():
            #if "normal node (Eps)"
            self.nodes[node].update(self)
        if calc_bound:
            return self.calc_bound()


if __name__ == '__main__':
    L.debug("temporary Main for expressionnet.py")
    

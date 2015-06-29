"""xml
- class to standardize how different models access the xml interfacing (parameter handling, calss creation etc.)
"""

import pdb
from xml.dom import minidom
import logging as L
import os.path
import string
import cPickle
import scipy as S
#models which need to be made public to the xml parsing function as we create classes in here


class CXml(object):
    """Helper class to handle xml data structures as needed for experiments"""

    def __init__(self,xml=None,xml_file=None):
        #open xml file
        if xml is None:
            self.xml = self.openXml(xml_file)
        else:
            self.xml = xml


    ##node creation functions etc
    def createModelClass(self):
        """create a class from the next model definition in self.xml"""
        if not self.xml.nodeName=="model":
            xml = self.xml.getElementsByTagName("model")[0]
        else:
            xml = self.xml
        _xml = CXml(xml)
        #get All parameters
        parameters = _xml.getParameters()
        #get node information
        nodes      = _xml.getNodes()
        #do we pass o the xml object?
        fwdxml     = xml.getAttribute("xml")
        class_name = xml.getAttribute("class")
        if fwdxml:
            eval_str = "%s(priors=nodes,xml=_xml,**parameters)" % class_name
        else:
            eval_str = "%s(priors=nodes,**parameters)" % class_name
        return eval(eval_str)
        pass

    ##functions to open/close search&replace in XML
    def openXml(self,xml_file):
        """open xml_file and replace parameters with IDs as in kwargin"""
        #1. open XML file
        xml = minidom.parse(xml_file)
        #find experiment node (take first one)
        xml = xml.getElementsByTagName('experiment')[0]
        #2. replace Ids in the parametersection of the experiment
        return xml

    def getElementsByTagName(self,tagName,level=None):
        """like standard getElementsByTagName but restrict to a certain level of depth"""
        def parentNode(node,n):
            if n>0:
                return parentNode(node.parentNode,n-1)
            else:
                return node.parentNode

        elements = self.xml.getElementsByTagName(tagName)
        if level is not None:
            elements = filter(lambda x: parentNode(x,level)==self.xml, elements)
        return elements

    def getModels(self):
        """return all models underneath the current xml node"""
        return self.getElementsByTagName('model',level=1)

    def getParameters(self,key="name",parse=True):
        """return the parameters of an xml model structure(key: key of the attributes, parse: True/False if true attributes are parsed, i.e. eval evaluated etc."""
        params = self.getElementsByTagName('param',1)
        rv = {}
        for param in params:
            value = param.getAttribute('value')
            if parse:
                ptype = param.getAttribute('type')
                if(param.getAttribute('eval')):
                    value = eval(value)
                elif(ptype=='matrix'):
                    value = self.parseMatrixParameter(value)
                elif(ptype=='double'):
                    value = S.double(value)
                elif(ptype=='int'):
                    value = S.int32(value)
                elif(ptype=='str'):
                    #no action for string
                    pass
                else:
                    raise Exception("Invalid Attribute exception attribute %s has no type or eval!" % param)
            rv[str(param.getAttribute(key))]=value
        return rv

    def parseMatrixParameter(self,value):
        """get a file parameter and automate the handling of pickled files, references etc"""
        if os.path.exists(value):
            if (string.lower(os.path.splitext(value)[1])=='.pickle'):
                #load file and return the value
                return cPickle.load(open(value,'rb'))
            else:
                L.error('gut a file with unsupported extension')
        else:
            #assume it its a global varible
            return globals()[value]

    def getNodes(self):
        """get all Nodes in the current xml structure"""
        nodes = self.getElementsByTagName('node',1)
        rv = {}
        for node in nodes:
            node_name = node.getAttribute('class')
            #get parameters
            node_parameters = CXml(node).getParameters()
            rv[node_name] = node_parameters
        return rv
        

    def getParameter(self,name,key="name"):
        parameters = self.getParameters(key=key)
        return parameters[name]


        

    def replaceXML(self,parameters):
        """replace entries with specific IDs from kwargin"""
        #2. replace all ID nodes with the kwargin parameters if they match
        params = self.xml.getElementsByTagName('param')
        for param in params:
            ID = param.getAttribute('id')
            if ID=='':
                continue
            elif ID  in parameters.keys():
                value = parameters[ID]
                param.setAttribute('value',value)
            pass



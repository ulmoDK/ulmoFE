import numpy as np

from Elements import BeamElem


class BeamMesh:
    
    def __init__(self,elemType="beam",):
        self.meshElems = None
        self.elemType = elemType
        self.mesh = None
        if self.elemType == "beam":
            self.ndof=2
        
    
    
    

    def nodeDOFs(self,nodeNum):
        
        return [2*nodeNum, 2*nodeNum+1]
    
    def beamNodes(self,elemNum):
        
        return [elemNum,elemNum+1]
        
    
    def createBeamElem(self,elemNum,length,EI,rho,A):
        
        nodes  = self.beamNodes(elemNum)
        
        dof = []
        for node in nodes:
            for d in self.nodeDOFs(node):
                dof.append(d)
           
        
        #BE=BeamElem(length,EI,rho,A)
        
        #elem = {"nodes":  nodes, "dofs" : dof, "K" : BE.K, "M" : BE.M}
        elem = {"nodes":  nodes, "dofs" : dof, "length" : length, "EI" : EI, "rho" : rho, "A" : A}
        
        
        return elem
        
        
        
    def createBeamMesh(self,numElems,length,EI,rho,A):
        
        #elemLength = beamLength / numElems
        numNodes = numElems+1
        numDOFs =  numNodes*2
        
        meshElems = []#
        for elem in range(numElems):
            meshElems.append(self.createBeamElem(elem,length,EI,rho,A))
            
                
        self.mesh = {"mesh" : meshElems, "numNodes" : numNodes, "numElems" : numElems, "numDOFs" : numDOFs }
        
        return self.mesh
    
  
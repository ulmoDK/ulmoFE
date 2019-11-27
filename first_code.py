# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 23:16:28 2019

@author: ulmo
"""

import numpy as np
from scipy.linalg import eigh



class BeamElem:
    
    def __init__(self,length=1.0,EI=1.0,rho=1.0,A=1.0):
        self.length = length
        self.EI = EI
        self.rho=rho
        self.A=A
        
        self.K = self.Kelem()
        self.M = self.Melem()
        
        
    def Kelem(self):
        
        
        L = self.length
        l=L
        
        """
        K = self.EI/L**3 * np.array([[12, 6*l, -12, 6*l],
				  [6*l, 4*l*l, -6*l, 2*l*l],
				  [-12, -6*l, 12, -6*l],
				  [6*l, 2*l*l, -6*l, 4*l*l]])
        """
        K = self.EI/L**3 * np.matrix([[12, 6*L, -12, 6*L],
                                    [6*L, 4*L**2, -6*L, 2*L**2],
                                    [-12, -6*L, 12, -6*L],
                                    [6*L, 2*L**2, -6*L, 4*L**2]])
        
        return K
    
    
    def Melem(self):
        
        L = self.length
        
        M = self.rho*self.A*L/420 * np.matrix([ [156,22*L, 54, -13*L],
                                                [22*L, 4*L**2, 13*L, -3*L**2],
                                                [54, 13*L, 156, -22*L],
                                                [-13*L, -3*L**2, -22*L, 4*L**2]])
        
        return M



        
    





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
           
        
        BE=BeamElem(length,EI,rho,A)
        
        elem = {"nodes":  nodes, "dofs" : dof, "K" : BE.K, "M" : BE.M}
        
        
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
    
    
    
class Assembly:
    
    def __init__(self,mesh):

        self.mesh = mesh
        self.K = self.createGlobal("K")
        self.M = self.createGlobal("M")
        #print(self.K)
        
    def createGlobal(self,GM="K"):
        G = np.zeros((self.mesh["numNodes"]*2,self.mesh["numNodes"]*2))
        
        
        localDofs = list(range(4))
        for elem in self.mesh["mesh"]:
            globalDofs = elem["dofs"]
            for globalRow,localRow in zip(globalDofs,localDofs):
                for globalCol,localCol in zip(globalDofs,localDofs):
                    G[globalRow,globalCol] += elem[GM][localRow,localCol]
                    
        return G
        

    
    
if __name__ == "__main__":
    #BE = BeamElem()
    #BE.Kelem()
    #print(BE.K)
    #print(BE.M)
    
    EI = 1.0
    A = 1.0
    rho = 1.0
    numElems =6000
    length = 1.0/numElems
    
    BM = BeamMesh()
    BM.createBeamMesh(numElems,length,EI,rho,A)
    
    
    
    #print(BM.meshElems)
    
    A = Assembly(BM.mesh)
    
    M = A.M
    K = A.K
    #print(K)
    #print(M)
    print(BM.mesh["mesh"][0]["K"])
    
    restrained_dofs = [1, 0, -2, -1]
    # remove the fixed degrees of freedom
    for dof in restrained_dofs:
        for i in [0,1]:
            M = np.delete(M, dof, axis=i)
            K = np.delete(K, dof, axis=i)
    
    
    evals, evecs = eigh(K,M)
    frequencies = np.sqrt(evals)
    print(frequencies[0])
    print(4.730041**2)
	









































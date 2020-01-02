# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 23:16:28 2019

@author: ulmo
"""

import numpy as np
#from scipy.linalg import eigh
#from scipy.linalg import eig

from Elements import BeamElem
from Mesh import BeamMesh




class Assembly:
    
    def __init__(self,mesh):

        self.mesh = mesh
        self.K = self.createGlobal("K")
        self.M = self.createGlobal("M")
        
        
        
        
        
    def createGlobal(self,GM="K"):
        
        G = np.zeros((self.mesh["numNodes"]*2,self.mesh["numNodes"]*2))
        
        
        localDofs = list(range(4))
        for elem in self.mesh["mesh"]:
            BE=BeamElem(elem["length"],elem["EI"],elem["rho"],elem["A"])
            elem["K"]=BE.K
            elem["M"]=BE.M
            globalDofs = elem["dofs"]
            for globalRow,localRow in zip(globalDofs,localDofs):
                for globalCol,localCol in zip(globalDofs,localDofs):
                    G[globalRow,globalCol] += elem[GM][localRow,localCol]
            elem["K"]=None # The matrix is set to zero to save memory
            elem["M"]=None # The matrix is set to zero to save memory                   
        return G
        


    def imposeBC(self,bcNodes):

        
        for node in bcNodes:
            self.K[node,:]=0.0
            self.K[:,node]=0.0
            self.K[node,node]=1.0
            
            self.M[node,:]=0.0
            self.M[:,node]=0.0
            
        
        
        return self





    def createLoadVector(self,loadNodes,loadMagnitude):
        
        
        if len(loadNodes) != len(loadMagnitude):
            raise ValueError("loadNodes and loadMagnitude must both be lists, and must be the same length")
        
        self.loadVector = np.zeros([self.K.shape[0],1])
        
        for node,magnitude in zip(loadNodes,loadMagnitude):
            self.loadVector[node]=magnitude

        return self.loadVector



    
if __name__ == "__main__":
 
    # Staticproperties    
    EI = 1.0
    A = 1.0
    rho = 1.0
    numElems = 50
    length = 1.0/numElems
    
    # Create mesh and elements
    BM = BeamMesh()
    BM.createBeamMesh(numElems,length,EI,rho,A)
    
    # Assemble and impose B.C.'s
    A = Assembly(BM.mesh)    
    restrained_dofs = [-1, 0, -2, 1]
    A.imposeBC(restrained_dofs)

    # solve
    #Solve(A,"eigen_frequency")
    
    """
    evals, evecs = eig(A.K,A.M)
    evals.sort()
    frequencies = np.sqrt(evals)
    print(frequencies[0].real)
    print(4.730041**2)
	"""
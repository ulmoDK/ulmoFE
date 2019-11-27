# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 23:16:28 2019

@author: ulmo
"""

import numpy as np
from scipy.linalg import eigh

from Elements import BeamElem
from Mesh import BeamMesh

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
	






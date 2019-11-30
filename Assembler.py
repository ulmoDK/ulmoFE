# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 23:16:28 2019

@author: ulmo
"""

import numpy as np
from scipy.linalg import eigh
from scipy.linalg import eig

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
        


    def imposeBC(self,bcNodes):
        
        
        for node in bcNodes:
            self.K[node,:]=0.0
            self.K[:,node]=0.0
            self.K[node,node]=1.0
            
            self.M[node,:]=0.0
            self.M[:,node]=0.0
            
        
        
        return self
    
    
    
    
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
    restrained_dofs = [1, 0, -2, -1]
    A.imposeBC(restrained_dofs)

    # solve
    evals, evecs = eig(A.K,A.M)
    evals.sort()
    frequencies = np.sqrt(evals)
    print(frequencies[0].real)
    print(4.730041**2)
	






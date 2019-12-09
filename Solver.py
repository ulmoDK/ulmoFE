# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 22:34:37 2019

@author: ulmo
"""

import numpy as np
from scipy.linalg import eigh
from scipy.linalg import eig

from Elements import BeamElem
from Mesh import BeamMesh
from Assembler import Assembly




class Solve:
    
    def __init__(self,Assembly,method):
        
        self.K = Assembly.K
        self.M = Assembly.M
        
        if method == "eigen_frequency":
            self.eigenFreq()
            self.printLowFreq()
    
    

    def eigenFreq(self):
        
        evals, self.evecs = eig(A.K,self.M)
        evals.sort()
        self.freq = np.sqrt(evals)
        
        
        return self
    
    def printLowFreq(self):
        
        print("The lowest frequency is {}".format(frequencies[0].real))
    
    
    
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
    Solve(A,"eigen_frequency")
    
    """
    evals, evecs = eig(A.K,A.M)
    evals.sort()
    frequencies = np.sqrt(evals)
    print(frequencies[0].real)
    print(4.730041**2)
	"""
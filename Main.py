import numpy as np
from scipy.linalg import eigh

from Elements import BeamElem
from Mesh import BeamMesh
from Assembler import Assembly
from Solver import Solve

# to-do: create doc-strings

if __name__ == "__main__":
    
    # Properties for beam elements (in this case all beam elements are created equal)    
    EI = 1.0
    A = 1.0
    rho = 1.0
    numElems = 50
    length = 1.0/numElems
    
    
    BM = BeamMesh()
    BM.createBeamMesh(numElems,length,EI,rho,A)
    
    
    
    
    
    # Assemble and impose B.C.'s
    A = Assembly(BM.mesh)
    A = Assembly(BM.mesh)    
    restrained_dofs = [-1, 0, -2, 1]
    A.imposeBC(restrained_dofs)

    # solve
    Solution=Solve(A,"eigen_frequency")
    print(4.730041**2)

	

    print("\nThe relative error is {} %\b".format(100*abs(Solution.freq[0].real-4.730041**2)/(4.730041**2)))
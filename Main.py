import numpy as np
from scipy.linalg import eigh

from Elements import BeamElem
from Mesh import BeamMesh
from Assembler import Assembly

# to-do: create doc-strings

if __name__ == "__main__":
    #BE = BeamElem()
    #BE.Kelem()
    #print(BE.K)
    #print(BE.M)
    
    EI = 1.0
    A = 1.0
    rho = 1.0
    numElems =50
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
	


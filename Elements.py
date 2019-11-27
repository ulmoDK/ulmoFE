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



 
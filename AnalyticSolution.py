# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 23:15:33 2020

@author: ulmo
"""


import numpy as np

def CantileverBeam(numPoints,length,EI,P):
    """Solution to cnatilever beam with point load at the free end"""
    
    Displacement = np.zeros([numPoints,1])
    
    x = np.linspace(0.0, length, num=numPoints)
    print(x)
    for i in range(x.shape[0]):
        Displacement[i] = P*length*x[i]**2 /(6*EI) *(3 - x[i]/length)
    return Displacement
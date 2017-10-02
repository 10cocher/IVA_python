# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:11:12 2017

@author: Emmanuel
"""
import numpy as np
from prec import myfloat

def oneRefl(zrefl,Rcoeff,dz,dh,nz,nx):
    nh = np.int(1)
    izrefl = np.int(np.floor(zrefl/dz))
    #
    xi = np.zeros((nz,nx,nh), dtype=myfloat, order='F')
    xi[izrefl,:,:] = Rcoeff/(dz*dh)
    return xi

# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 12:57:55 2017

@author: Emmanuel
"""

"""
prec.py
-------
single precision: np.float32 np.complex64
double precision: np.float64 np.complex128
"""

import numpy as np

myfloat = np.dtype(np.float32) ; mycomplex = np.dtype(np.complex64)
#myfloat = np.dtype(np.float64) ; mycomplex = np.dtype(np.complex128)


def quad(num):
    return np.float32(num)
    #return np.float64(num)

def complexquad(num):
    return np.complex64(num)
    #return np.complex128(num)
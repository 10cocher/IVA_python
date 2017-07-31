# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 17:25:56 2017

@author: Emmanuel
"""

import numpy as np

def SmoothSin(a):
        """ Creates a smoothing function monotically mapping the interval
        [0,1] to  [0,1] with a zero-derivative at x=0 and x=1
        ----------
        Parameters
        a, float
        ----------
        """
        return np.square( np.sin( np.pi*np.float(a)/2. ) ) 
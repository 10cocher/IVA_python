# -*- coding: utf-8 -*-
"""
Created on Sun Jul  2 16:26:49 2017

@author: Emmanuel
"""

import numpy as np
import matplotlib.pyplot as plt
import src_wavelet as src_wvlt

if __name__=="__main__":
    #
    plt.close("all")
    #
    zmin=0
    zmax=600
    dz=5
    
    tmin=0
    tmax=0.5
    dt=0.001
    
    tsrc=0.4
    fmax=40
    
    nz=np.int(np.floor((zmax-zmin)/dz))
    nt=np.int(np.floor((tmax-tmin)/dt))
    nsrc2=np.int(np.floor(tsrc/dt))
    nsrc=np.int(2)*nsrc2 + np.int(1)
    ntt=nt+nsrc2
    
    sax=np.linspace(-tsrc,tsrc,num=nsrc,dtype=float)
    
    print("nz = %s nt = %s" % (nz, ntt))
    print('nz={} nt={}'.format(nz,ntt))
    
    
    
    src=src_wvlt.defsrc(nsrc,fmax,dt)
    
    src_wvlt.defsrcdcnv(src,dt)
    
    '''
    #plt.figure(1)
    plt.plot(sax,src)
    plt.xlabel('time (s)')
    plt.title('source wavelet')
    plt.show()
    '''
    

# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 11:44:48 2017

@author: Emmanuel
"""

import numpy as np
from prec import quad, myfloat
#import matplotlib.pyplot as plt
#import mymaths


def prop_1D(S,v,dz,dt):
    ntt = S.shape[0]
    nz  = S.shape[1]
    nze = 5*nz
    #
    vext = np.ones((nze,), dtype=myfloat, order='F')
    Sext = np.zeros((ntt,nze), dtype=myfloat, order = 'F')
    Rext = np.zeros_like(Sext)
    #
    Sext[:,2*nz:3*nz] = S
    vext[2*nz:3*nz] = np.squeeze(v)
    vext[0   :2*nz] = v[0,0]
    vext[3*nz:nze ] = v[nz-1,0]

    #
    for it in np.arange(1,ntt-1,dtype=int):
        if (it < 5) or (it > ntt-5):
            print(it,ntt-2)
        for iz in np.arange(0,nze, dtype=int):
            izm = np.amax([0,iz-1])
            izp = np.amin([nze-1,iz+1])
            #
            Rext[it+1,iz] = 2.*Rext[it,iz] - Rext[it-1,iz] \
                + np.square( vext[iz]*dt/dz ) * \
                (Rext[it,izm] - 2.*Rext[it,iz] + Rext[it,izp]) \
                + Sext[it,iz] * np.square( dt*vext[iz] )
    #
    R = Rext[:,2*nz:3*nz]
    return R

def calcS0(src,zsrc,isrc,dz,dx,noff,ntt,nz,nxx):
    nsrc   = src.shape[0]
    izsrc  = np.int(np.floor(zsrc/dz))
    maxoff = (noff - np.int(1))/np.int(2)
    ispos  = np.int( np.minimum(isrc, maxoff + 1) - 1)
    #
    S0 = np.zeros((ntt,nz,nxx), dtype=myfloat, order='F')
    S0[0:nsrc,izsrc,ispos] = src/(dx*dz)
    return S0

def verifCond(dz, dx, dt, fmax, vmin, vmax, mode1D2D):
    """ checks if current parameters satisfy stability and dispertion conditions
    Parameters
    ----------
    dz  : spatial sampling in depth  (meter)
    dx  : spatial sampling laterraly (meter)
    dt  : time sampling (meter)
    fmax: maximum frequency in source wavelet
    vmin: min value of velocity model
    vmax: max value of velocity model
    mode1D2D : int, 1D or 2D case
    """
    # Initialisation
    okDispersionZ = False
    okDispersionX = False
    okStability   = False
    # Compute maximum admissible parameters
    dzmax = vmin / (quad(10)*fmax)
    dxmax = vmin / (quad(10)*fmax)
    if mode1D2D == np.int(1):
        dtmax = dz / vmax # in 1D
    if mode1D2D == np.int(2):
        dtmax = dz / (vmax * np.sqrt(2)) # in 2d
    print('----------------------------------------------------------')
    # Dispersion
    print('Dispersion condition...')
    if (dz <= dzmax):
        okDispersionZ = True
        print ('   dz = %5.3f m  ; dzmax = %5.3f m  --- ok!' %(dz,dzmax))
    else:
        print ('   dz = %5.3f m  ; dzmax = %5.3f m  --- !!! problem !!!' %(dz,dzmax))
    if (mode1D2D == np.int(2)):
        if (dx <= dxmax):
            okDispersionX = True
            print ('   dx = %5.3f m  ; dxmax = %5.3f m  --- ok!' %(dx,dxmax))
        else:
            print ('   dx = %5.3f m  ; dxmax = %5.3f m  --- !!! problem !!!' %(dx,dxmax))
    # Stability
    print('Stability condition...')
    if (dt <= dtmax):
        okStability = True
        print ('   dt = %5.3f ms ; dtmax = %5.3f ms --- ok!' %(dt*1000,dtmax*1000))
    else:
        print ('   dt = %5.3f ms ; dtmax = %5.3f ms --- !!! problem !!!' %(dt*1000,dtmax*1000))
    #
    print('----------------------------------------------------------')
    ok = okStability and okDispersionZ and (okDispersionX or (mode1D2D == np.int(1)))
    return ok

# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 23:03:09 2017

@author: Emmanuel
"""

import numpy as np
from prec import quad
import matplotlib.pyplot as plt

def comparNuage1D(ax,v1,v2):
    min1 = np.amin(v1)
    min2 = np.amin(v2)
    max1 = np.amax(v1)
    max2 = np.amax(v2)
    minmin = np.amin([min1, min2])
    maxmax = np.amax([max1, max2])
    #
    ax.plot([minmin, maxmax],[0, 0],'k')
    ax.plot([0, 0],[minmin, maxmax],'k')
    ax.plot([minmin, maxmax],[minmin, maxmax],'k')
    #
    ax.grid(color='gray')
    ax.plot(v1,v2,'b+')
    ax.axis('equal')
    ax.axis([minmin, maxmax, minmin, maxmax])

def comparNuage2D(ax,v1,v2):
    n1 = v1.shape[0]
    n2 = v1.shape[1]
    w1 = v1.reshape((n1*n2,))
    w2 = v2.reshape((n1*n2,))
    comparNuage1D(ax,w1,w2)

def comparNuage3D(ax,v1,v2):
    n1 = v1.shape[0]
    n2 = v1.shape[1]
    n3 = v1.shape[2]
    w1 = v1.reshape((n1*n2*n3,))
    w2 = v2.reshape((n1*n2*n3,))
    comparNuage1D(ax,w1,w2)

def plot_wvfld_1D(ax,P,tmin,tmax,zmin,zmax):
    limP  = np.amax(np.abs(P))
    #
    im = ax.imshow(np.transpose(P), clim=((-limP,limP)), cmap='seismic', \
                     aspect='auto', extent=[tmin,tmax,zmax,zmin])
    plt.colorbar(im)
    ax.set_xlabel('time (s)'); ax.set_ylabel('depth (m)')
    ax.grid(color='gray')

def plot_Pobs_2D(ax, P, tmin, tmax, offmin, offmax, nsrc, cut=True):
    nsrc2 = np.int((nsrc-1)/2)
    #
    if cut:
        itmin = nsrc2
        tmin2 = quad(0)
    else:
        itmin = 0
        tmin2 = tmin
    #
    limP = np.amax(np.abs(P[itmin:,:]))
    #
    im = ax.imshow(P[itmin:,:], clim=((-limP,limP)), cmap='seismic', \
                   aspect='auto', extent=[offmin,offmax,tmax,tmin2])
    plt.colorbar(im)
    ax.set_xlabel('offset (m)'); ax.set_ylabel('time (s)')
    ax.grid(color='gray')
    #

def plot_xi_zx(ax,xi,xmin,xmax,zmin,zmax):
    xi_zx = np.squeeze(xi)
    limxi = np.amax(np.abs(xi_zx))
    #
    im = ax.imshow(xi_zx, clim=((-limxi,limxi)), cmap='seismic', \
                   aspect='equal', extent=[xmin,xmax,zmax,zmin])
    plt.colorbar(im)
    ax.set_xlabel('x-position (m)'); ax.set_ylabel('depth (m)')
    ax.grid(color='gray')

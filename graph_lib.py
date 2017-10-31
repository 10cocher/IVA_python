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

def plot_Pobs(ax, P, tax, tmin, tmax, offmin, offmax, nsrc, mode1D2D, cut=True):
    nsrc2 = np.int((nsrc-1)/2)
    # removes negative times if cut==True
    if cut:
        itmin = nsrc2
        tmin2 = quad(0)
    else:
        itmin = 0
        tmin2 = tmin
    # 1D case
    if (mode1D2D == np.int(1)):
        ax.plot(tax[itmin:],np.squeeze(P[itmin:,:,:]))
        ax.set_xlabel('time (s)')
        ax.set_title('Pobs')
    # 2D case
    if (mode1D2D == np.int(2)):
        limP = np.amax(np.abs(P[itmin:,:]))
        im = ax.imshow(P[itmin:,:], clim=((-limP,limP)), cmap='seismic', \
                       aspect='auto', extent=[offmin,offmax,tmax,tmin2])
        plt.colorbar(im)
        ax.set_xlabel('offset (m)'); ax.set_ylabel('time (s)')
    #
    ax.grid(color='gray')

def plot_trace(ax, P, tax, isrc, ioff, nsrc, cut=True, mylabel=''):
    nsrc2 = np.int((nsrc-1)/2)
    if cut:
        itmin = nsrc2
    else:
        itmin = 0
    #
    trace = np.squeeze(P[:,ioff])
    ax.plot(trace[itmin:],tax[itmin:],label=mylabel)
    ax.set_ylabel('time (s)')
    ax.grid(color='gray')

def plot_xi_1D(ax, xi, zax, title='xi'):
    ax.plot(zax,np.squeeze(xi))
    ax.set_xlabel('depth (m)')
    ax.set_title(title)
    ax.grid(color='gray')

def plot_xi_zx(ax, xi, xmin, xmax, zmin, zmax, title='xi'):
    #
    nh    = xi.shape[2]
    nh2   = np.int((nh-1)/2)
    xi_zx = np.squeeze(xi[:,:,nh2])
    #
    limxi = np.amax(np.abs(xi))
    #
    im = ax.imshow(xi_zx, clim=((-limxi,limxi)), cmap='seismic', \
                   aspect='equal', extent=[xmin,xmax,zmax,zmin])
    plt.colorbar(im)
    ax.set_xlabel('x-position (m)'); ax.set_ylabel('depth (m)')
    ax.grid(color='gray')
    ax.set_title(title + ' at h=0 m')

def plot_xi_CIG(ax, xi, xax, hmin, hmax, zmin, zmax, title='xi'):
    #
    nx     = xi.shape[1]
    nx2    = np.int( np.floor(nx/2) )
    xi_CIG = np.squeeze(xi[:,nx2,:])
    #
    limxi = np.amax(np.abs(xi))
    #
    im_CIG = ax.imshow(xi_CIG, clim=((-limxi,limxi)), cmap='seismic', \
                   aspect='equal', extent=[hmin,hmax,zmax,zmin])
    plt.colorbar(im_CIG)
    ax.set_xlabel('h (m)'); ax.set_ylabel('depth (m)')
    ax.grid(color='gray')
    ax.set_title(title + ': CIG at x=%5.1f m'%(xax[nx2]))

def plot_vel(ax, v, vmin, vmax, zax, zmin, zmax, xmin, xmax, mode1D2D,\
             title='velocity model'):
    if mode1D2D == np.int(1):
        ax.plot(zax,np.squeeze(v))#, label='vmod')
        #plt.legend()#handles=[line_vmod,line_vini])
        ax.set_xlabel('depth (m)')
    if mode1D2D == np.int(2):
        im = ax.imshow(v, clim=((vmin,vmax)), cmap='seismic', \
                     aspect='equal', extent=[xmin,xmax,zmax,zmin])
        plt.colorbar(im)
        ax.set_xlabel('x-position (m)')
        ax.set_ylabel('depth (m)')
    ax.set_title(title)
    ax.grid(color='gray')

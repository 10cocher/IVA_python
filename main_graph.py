# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 19:17:12 2017

@author: Emmanuel
"""

import numpy as np
from prec import myfloat
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import graph_lib as glib
import pickle
import sys

if __name__=="__main__":
    plt.close("all")
    #%% main parameters
    #%% stores param in a python-readable file
    folder = 'out/'
    f = open(folder + 'param.pckl', 'rb')
    nz, nx, nh, nt, ntt, nsrc, nstotal, ns, nrcv, noff,\
    dz, dx, dh, dt, dsint,\
    zmin, zmax, xmin, xmax, hmin, hmax, tmin, tmax, tsrc,\
    vmin, vmax, mode1D2D,\
    = pickle.load(f)
    f.close()
    #%%
    maxoff = (noff-1)/2
    offmax = dx * maxoff
    offmin = -offmax
     #%% print all dimensions
    print('----------------------------------------------------------')
    print('     nz = %4i ; dz = %5.3f m' %(nz,dz))
    print('     nx = %4i ; dx = %5.3f m' %(nx,dx))
    print('     nh = %4i ; dh = %5.3f m' %(nh,dh))
    print('     nt = %4i ; dt = %5.3f ms (positive time samples)' %(nt,dt*1000))
    print('    ntt = %4i (neg. and pos. time samples)' %(ntt))
    print('   nsrc = %4i (number of time samples for src wavelet)' %(nsrc))
    print('nstotal = %4i ; ds = %2i (total number of sources in the acquisition)' %(nstotal,dsint))
    print('     ns = %4i (number of sources treated by this processor)' %(ns))
    print('   nrcv = %4i' %(nrcv))
    print('   noff = %4i' %(noff))
    print('----------------------------------------------------------')
    #%% axies for graphics
    zax = np.linspace(zmin,  zmax-dz, num=nz,   dtype=myfloat)
    xax = np.linspace(xmin,  xmax-dx, num=nx,   dtype=myfloat)
    hax = np.linspace(hmin,  hmax, num=nh,   dtype=myfloat)
    tax = np.linspace(tmin,  tmax, num=ntt,  dtype=myfloat)
    sax = np.linspace(-tsrc, tsrc, num=nsrc, dtype=myfloat)
    offax = np.linspace(offmin, offmax, num=noff, dtype=myfloat)

    #%%
    def createStack(P,limP,ax,offmin,offmax,zmin,zmax,tax):
        ntt  = P.shape[0]
        nz   = P.shape[1]; nz2   = np.int(np.floor(nz/2))
        noff = P.shape[2]; noff2 = np.int(np.floor(noff/2))
        ims = [] #
        for it in np.arange(0, ntt, dtype=int):
            im = ax.imshow(np.squeeze(P[it,:,:]), clim=((-limP,limP)), cmap='seismic', \
                     aspect='auto', extent=[offmin,offmax,zmax,zmin], animated=True)
            title = ax.annotate('it=%3i t=%6.4f' %(it,tax[it]),(nz2,noff2)) # add text
            #title = ax.annotate(it,(5,5))
            #im.set_title('t = %5.3f' %(tax[it]))
            ims.append([im,title])
        return ims
    #%% subrtouine to read files
    def readWavefield(folder,name):
        P = np.load(folder + name + '.npy')
        limP = np.amax(np.abs(P))
        return P, limP
    def readXi(folder,name):
        return readWavefield(folder,name)
    #%% Observed data and velocity models
    isrc = 128
    if mode1D2D == np.int(1):
        Pobs, limPobs = readWavefield(folder,'Pobs')
    else:
        Pobs, limPobs = readWavefield(folder,'Pobs_%06i'%(isrc))
    #
    vini = np.load(folder + 'vini.npy')
    vmod = np.load(folder + 'vmod.npy')
    #
    fig = plt.figure('vmod and ximod')
    #
    ax1 = fig.add_subplot(221)
    glib.plot_Pobs(ax1, Pobs, tax, tmin, tmax, offmin, offmax, nsrc, mode1D2D, True)
    #
    ax2 = fig.add_subplot(222)
    glib.plot_vel(ax2, vmod, vmin, vmax, zax, zmin, zmax, xmin, xmax, mode1D2D, 'vmod')
    #
    ax3 = fig.add_subplot(223)
    glib.plot_vel(ax3, vini, vmin, vmax, zax, zmin, zmax, xmin, xmax, mode1D2D, 'vini')
    #glib.plot_vmod(ax2, vmod, vini, zax)
    #ax2.plot(zax,np.squeeze(vmod),zax,np.squeeze(vini))
    #
    #%% compares adjoint and inverse reflectivity
    #
    xiadj, limxiadj = readWavefield(folder,'xiadj')
    xiinv, limxiinv = readWavefield(folder,'xiinv')
    #
    if mode1D2D == np.int(1):
        Fxiadj, limP = readWavefield(folder,'Fxiadj')
        Fxiinv, limP = readWavefield(folder,'Fxiinv')
    else:
        Fxiadj, limP = readWavefield(folder,'Fxiadj_%06i'%(isrc))
        Fxiinv, limP = readWavefield(folder,'Fxiinv_%06i'%(isrc))
    #
    if mode1D2D==np.int(1):
        fig2 = plt.figure('test adjoint/inverse')
        ax2 = fig2.add_subplot(222)
        ax4 = fig2.add_subplot(224)
        ax1 = fig2.add_subplot(221)
        glib.plot_xi_1D(ax2, xiadj, zax, 'xi adjoint')
        glib.plot_xi_1D(ax4, xiinv, zax, 'xi inverse')
        #
        ax1.plot(tax,np.squeeze(Pobs)  ,label='Pobs')
        ax1.plot(tax,np.squeeze(Fxiadj),label='F xiadj')
        ax1.plot(tax,np.squeeze(Fxiinv),label='F xiinv')
        ax1.legend()
        ax1.grid(color='gray')
        #
    else:
        fig2 = plt.figure('test adjoint/inverse')
        ax1 = fig2.add_subplot(221)
        glib.plot_xi_zx( ax1, xiadj, xmin, xmax, zmin, zmax, title='xiadj')
        ax2 = fig2.add_subplot(222)
        glib.plot_xi_CIG(ax2, xiadj, xax, hmin, hmax, zmin, zmax, title='xiadj')
        ax3 = fig2.add_subplot(223)
        glib.plot_xi_zx( ax3, xiinv, xmin, xmax, zmin, zmax, title='xiinv')
        ax4 = fig2.add_subplot(224)
        glib.plot_xi_CIG(ax4, xiinv, xax, hmin, hmax, zmin, zmax, title='xiinv')
        #
        #
        #ioff = np.int((noff-1)/2)
        ioff = np.int(np.floor(3*noff/4))
        #
        fig3 = plt.figure('test demigrate')
        ax1 = fig3.add_subplot(221)
        #ax2 = fig3.add_subplot(222)
        #ax3 = fig3.add_subplot(223)
        glib.plot_trace(ax1, Pobs, tax, isrc, ioff, nsrc, cut=True, mylabel='Pobs')
        glib.plot_trace(ax1, Fxiadj, tax, isrc, ioff, nsrc, cut=True, mylabel='Fxiadj')
        glib.plot_trace(ax1, Fxiinv, tax, isrc, ioff, nsrc, cut=True, mylabel='Fxiinv')
        ax1.legend()
        #ax3 = fig2.add_subplot(2,2,3)
        #glib.plot_wvfld_1D(ax3,np.squeeze(sw),tmin,tmax,zmin,zmax)
        #
        #ax4 = fig2.add_subplot(2,2,4)
        #glib.plot_wvfld_1D(ax4,np.squeeze(rw),tmin,tmax,zmin,zmax)

    #%%
    isrc = 128
    if False:
        Pobs, limPobs = readWavefield(folder,'Pobs_%06i'%(isrc))
        #
        fig = plt.figure('Pobs at ix=%i m' %(isrc))
        ax1 = fig.add_subplot(111)
        glib.plot_Pobs(ax1, Pobs, tax, tmin, tmax, offmin, offmax, nsrc, mode1D2D, True)
    #
    #%%
    if False:
        ximod1 = np.squeeze(np.load(folder + 'ximod_1.npy'))
        ximod2 = np.squeeze(np.load(folder + 'ximod_2.npy'))
        fig1 = plt.figure(num='ximod')
        #
        ax1 = fig1.add_subplot(2,1,1)
        glib.plot_xi_zx(ax1,ximod1,xmin,xmax,zmin,zmax)
        #
        ax2 = fig1.add_subplot(2,1,2)
        glib.plot_xi_zx(ax2,ximod2,xmin,xmax,zmin,zmax)
    #%%
    if False:
        S0, limS0 = readWavefield(folder,'S0')
        P0, limP0 = readWavefield(folder,'P0')
        S1, limS1 = readWavefield(folder,'S1')
        P1, limP1 = readWavefield(folder,'P1')
        #
        print('shape(S0)=',S0.shape)
        print('shape(P0)=',P0.shape)
        print('shape(S1)=',S1.shape)
        print('shape(P1)=',P1.shape)
    #%% 2D anim
    if False:
        fig56 = plt.figure(num='compar S0 3D')
        ax1    = fig56.add_subplot(111)
        imsP = createStack(P1,limP1,ax1,offmin,offmax,zmin,zmax,tax)
        ani1 = animation.ArtistAnimation(fig56, imsP, interval=5, blit=True,
                                    repeat_delay=1000)
        plt.show
    #%% 1D wvfld
    if False:
         #%%
        fig2 = plt.figure(num='S0 P0')
        #
        ax1 = fig2.add_subplot(1,4,1)
        glib.plot_wvfld_1D(ax1,np.squeeze(S0),tmin,tmax,zmin,zmax)
        #
        ax2 = fig2.add_subplot(1,4,2)
        glib.plot_wvfld_1D(ax2,np.squeeze(P0),tmin,tmax,zmin,zmax)
        #
        ax3 = fig2.add_subplot(1,4,3)
        glib.plot_wvfld_1D(ax3,np.squeeze(S1),tmin,tmax,zmin,zmax)
        #
        ax4 = fig2.add_subplot(1,4,4)
        glib.plot_wvfld_1D(ax4,np.squeeze(P1),tmin,tmax,zmin,zmax)
        #
        #ax3 = fig2.add_subplot(1,3,3)
        #glib.comparNuage2D(ax3,np.squeeze(S0_P),np.squeeze(S0_F))
    #%%
    sys.exit('debug')
    #%%
    fig3 = plt.figure(num='compar S0 3D')
    #
    ax1 = fig3.add_subplot(1,1,1)
    glib.comparNuage3D(ax1,S0_P,S0_F)


    #%%
    sys.exit('debug')
    #%% graphics
    if False:
        #xax = np.arange(0, nx, dtype=int)
        print('shape(xax) =',xax.shape)
        print('shape(sOK) =',sOK.shape)
        #print('shape(staptotal) =',staptotal.shape)
        #
        #
        fig1 = plt.figure(num='src wavelet')
        #
        ax1 = fig1.add_subplot(1, 1, 1)
        ax1.grid(color='gray')
        ax1.plot(xax,sOK)
        ax1.set_xlabel('x index')
        ax1.set_title('sOK')
        #
        #ax2 = fig1.add_subplot(2, 2, 2)
        #ax2.grid(color='gray')
        #ax2.plot(xax,staptotal)
        #ax2.set_xlabel('x index')
    #%%
    fig1 = plt.figure(num='source acquisition')
    #
    ax1 = fig1.add_subplot(2,1,1)
    ax1.plot(xax,sOK)
    ax1.grid(color='gray')
    #%%
    S = np.load('S.npy')
    P = np.load('P.npy')
    S1 = np.load('S1.npy')
    P1 = np.load('P1.npy')
    #%%
    fig2 = plt.figure(num='vext')
    #
    ax1 = fig2.add_subplot(2,3,1)
    glib.plot_wvfld_1D(ax1,S,tmin,tmax,zmin,zmax)
    #
    ax2 = fig2.add_subplot(2,3,2)
    glib.plot_wvfld_1D(ax2,S1,tmin,tmax,zmin,zmax)
    #
    ax4 = fig2.add_subplot(2,3,4)
    glib.plot_wvfld_1D(ax4,P,tmin,tmax,zmin,zmax)
    #
    ax5 = fig2.add_subplot(2,3,5)
    glib.plot_wvfld_1D(ax5,P1,tmin,tmax,zmin,zmax)
    #
    ax3 = fig2.add_subplot(2,3,3)
    glib.comparNuage2D(ax3,S,S1)
    #
    ax6 = fig2.add_subplot(2,3,6)
    glib.comparNuage2D(ax6,P,P1)

    #%%
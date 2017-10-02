# -*- coding: utf-8 -*-
"""
Created on Sun Jul  2 16:26:49 2017

@author: Emmanuel Cocher
"""


do_mpi = True

if do_mpi:
    from mpi4py import MPI
from prec import myfloat, quad
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import src_wavelet as src_wvlt
import acquisition
import propag_python as prop_py
import refmodels
import sys
#import monmodule
import propag as prop
#import f2py2e
import time
import pickle

if __name__=="__main__":
    #
    #plt.close("all")
    #%% MPI
    #
    if do_mpi:
        comm = MPI.COMM_WORLD
        fcomm = comm.py2f()
        rank = comm.Get_rank()
        nb_procs = comm.Get_size()
        #
        #print('type(comm)=',type(comm))
        #print('Hello, I am processor ranked %01i/%01i with fcomm=%010i' %(rank,nb_procs,fcomm))
        #time.sleep(5)
        #comm.Barrier()
        #
    else:
        print('working without MPI...')
        comm = 0
        fcomm = 1
        rank = np.int(0)
        nb_procs = np.int(1)
    #print( 'comm=',comm)
    #print('fcomm=',fcomm)
    #%%
    def mystop(do_mpi, comm, rank):
        if do_mpi:
            print('I''m number %i and I''m disconnecting' %(rank))
            #comm.Abort()
        sys.exit('debug')
    #%% main parameters
    folder = 'out/'
    mode1D2D = np.int(2) # (1) = 1D ; (2) = 2D
    MethAcq  = np.int(0) # see acquisition.okReceiver for details
    MethAcq2 = np.int(0) # (0)=surface acqusition ; (1)=VSP
    #%% Propagation code
    MethProp = np.int(1) # (0)=rays ; (1)=finite differences
    npml     = np.int(30)
    #%%
    MethXi   = np.int(0) # (0) xi = deltac/2c0 ; (1) xi = 2deltac/c0^3
    #%%
    MethJxi  = np.int(0) # (0) ||P+M-Pobs-Mobs|| ; (1) ||P-Pobs|| + am*||M-Mobs||
    am = 20.
    #%% tapers
    # taper on the extended reflectivity model during wavefield propagation
    MethTap  = np.int(1) # (0) = no taper ; (1) = taper according to MethTap2
    MethTap2 = np.int(7) # conversion from the three binary digits ZXH
    ptap = np.array([50., 130., 50., 200., 60.], dtype=myfloat)
    """ details :
        ptap[0] (meters) (z) zeros  from z=0.      to z=ptap[0]
        ptap[1] (meters) (z) smooth from z=ptap[0] to z=ptap[1]
        ptap[2] (meters) (z) smooth from z=zmax    to z=ptap[2]
        ptap[3] (meters) (x) width of apodisation in x
        ptap[4] (meters) (h) width of apodisation in h
    """
    # taper on the extended reflectivity for the regularisation function
    MethTpRg = np.int(0) # (0) = no taper ; (1) taper according to MethTpR2
    MethTpR2 = np.int(0) # same as MethTap2
    # taper on receivers
    MethTapD = np.int(3) # (0) = no taper ; (1) see inside prop.f90 for details
    # taper on sources
    MethTapS = np.int(0) # (0) = no taper ; (1) see acquisition.defTapSources for details
    # taper on data
    MethMR   = np.int(0) # (0) = no taper ; (1) see inside prop.f90 for details
    Mtap = np.array([1., 0.1, 120., 0.05, 3.0e-4], dtype=myfloat)
    """ details :
        Mtap[0] time exposant
        Mtap[1] width of apodisation along time axis (in seconds)
        Mtap[2] width of apodisation along receiver axis (in meters)
        Mtap[3] t0 for MethMR = 2
        Mtap[4] at for MethMR = 2
    """
    #%%
    MethQ    = np.int(1)
    MethDcnv = np.int(1)
    #%%
    PMeths = np.array([MethJxi , MethXi, MethTap ,MethTap2, MethMR,
                       MethTapD, MethQ , MethDcnv,MethAcq2], dtype=int)
    #%% spatial dimensions
    zmin = quad(   0) ; zmax=quad( 600) ;  dz=quad(5) ; nz=np.int(np.floor((zmax-zmin)/dz))
    xmin = quad(   0) ; xmax=quad(1000) ;  dx=quad(5) ; nx=np.int(np.floor((xmax-xmin)/dx))
    hmin = quad(-180) ; hmax=-hmin      ;  dh=dx      ; nh=np.int(np.floor((hmax-hmin)/dh))
    if mode1D2D == np.int(1):
        xmin = quad(0) ; xmax = quad(0) ; dx = quad(1) ; nx = np.int(1)
        hmin = quad(0) ; hmax = quad(0) ; dh = quad(1) ; nh = np.int(1)
    #%% temporal dimension
    tmin = quad(0) ; tmax=quad(0.5) ; dt=quad(0.001) ; nt=np.int(np.floor((tmax-tmin)/dt))
    #%% simulation of wave propagation around source bewteen -maxoffR and +maxoffR
    maxoffR = quad(400) # in meters
    if mode1D2D == np.int(1):
        maxoffR = quad(0)
    maxoff = np.int(np.floor(maxoffR/dx))
    noff   = 2*maxoff + 1
    #%% acquisition : position of sources and receivers
    # vertical position of sources
    izsrc = np.int(1) ;  zsrc = izsrc * dz
    # lateral position of sources
    # typically every 'ds' between xSlim[0] and xSlim[1] and between xSlim[2] and xSlim[3]
    xSlim = np.array([450., 450., 450., 465.], dtype=myfloat)
    ds = np.int(1)
    sOKs = False # forces symmetric sources repartition
    sOKf = False # forces sources to be at xSlim[1] and xSlim[2]
    # lateral positions of the sources on the grid
    sOK, nstotal = acquisition.okSource(xSlim,dx,nx,ds,sOKs,sOKf,mode1D2D)
    # repartition of sources between processors
    slist, stap = acquisition.repartSource(sOK,nb_procs,rank)
    ns = slist.shape[0]
    ## careful : slist starts at 0, maybe a problem when passing to Fortran...
    # receivers position
    xRlim = np.array([-400., 400.], dtype=myfloat)
    zVSP  = np.array([ 100., 200.], dtype=myfloat)
    if MethAcq2==0 :
        nrcv = noff
    if MethAcq2==1 :
        nrcv = nz
    rOK = acquisition.okReceiver(xRlim,xSlim,zVSP,xmin,xmax,dz,dx,nrcv,slist,\
                                 MethAcq,MethAcq2,mode1D2D)
    if rank == np.int(0):
        np.save(folder + 'sOK.npy',sOK)
        np.save(folder + 'rOK.npy',rOK)
    #%% source wavelet
    tsrc = quad(0.4) # positive time in seconds
    fmax = quad(40)  # maximum frequency in Herz
    #
    nsrc2 = np.int( np.floor(tsrc/dt) )
    nsrc  = np.int(2)*nsrc2 + np.int(1)
    ntt   = nt + nsrc2
    #
    src     = src_wvlt.defsrc(nsrc, fmax, dt)
    srcdcnv = src_wvlt.defsrcdcnv(src, dt)
    #
    tmin = -nsrc2 * dt
    tmax = nt     * dt
    #%% vmin and vmax (in meter/second)
    vmin = quad(2450)
    vmax = quad(3050)
    #%% print all dimensions
    if rank == np.int(0):
        print('----------------------------------------------------------')
        print('     nz = %4i ; dz = %5.3f m' %(nz,dz))
        print('     nx = %4i ; dx = %5.3f m' %(nx,dx))
        print('     nh = %4i ; dh = %5.3f m' %(nh,dh))
        print('     nt = %4i ; dt = %5.3f ms (positive time samples)' %(nt,dt*1000))
        print('    ntt = %4i (neg. and pos. time samples)' %(ntt))
        print('   nsrc = %4i (number of time samples for src wavelet)' %(nsrc))
        print('nstotal = %4i ; ds = %2i (total number of sources in the acquisition)' %(nstotal,ds))
        print('     ns = %4i (number of sources treated by this processor)' %(ns))
        print('   nrcv = %4i' %(nrcv))
        print('   noff = %4i' %(noff))
        print('----------------------------------------------------------')
    #%% stores param in a python-readable file
    if rank == np.int(0):
        f = open(folder + 'param.pckl', 'wb')
        pickle.dump([nz, nx, nh, nt, ntt, nsrc, nstotal, ns, nrcv, noff,\
                     dz, dx, dh, dt, ds,\
                     zmin, zmax, xmin, xmax, hmin, hmax, tmin, tmax, tsrc,\
                     vmin, vmax, mode1D2D,\
                     ], f)
        f.close()
    #%% axies for graphics
    zax = np.linspace(zmin,  zmax, num=nz,   dtype=myfloat)
    xax = np.linspace(xmin,  xmax, num=nx,   dtype=myfloat)
    hax = np.linspace(hmin,  hmax, num=nh,   dtype=myfloat)
    tax = np.linspace(tmin,  tmax, num=ntt,  dtype=myfloat)
    sax = np.linspace(-tsrc, tsrc, num=nsrc, dtype=myfloat)
    #%% checks stability and dispersion conditions
    #if rank == np.int(0):
    ok = prop_py.verifCond(dz, dx, dt, fmax, vmin, vmax, mode1D2D)
    #if do_mpi:
    #comm.Bcast([ok, MPI.LOGICAL], root=0)
    #print('rank=%2i, ok=%l' %(rank,ok))
    if (not ok):
        sys.exit('stability or dispersion conditions is not satisfied by the current values of dz, dx and dt')
    #%% velocity models
    vmod = quad(3000) * np.ones((nz,nx), dtype=myfloat, order='F')
    vini = quad(2500) * np.ones((nz,nx), dtype=myfloat, order='F')
    #%% reflectivity models
    zrefl = quad(400)
    Rcoeff= quad(0.5)
    ximod = refmodels.oneRefl(zrefl,Rcoeff,dz,dh,nz,nx)
    ximodzero = False
    #%% outputs graphics
    if rank == np.int(0):
        np.save(folder + 'vmod.npy',vmod)
        np.save(folder + 'vini.npy',vini)
    #%% for the example
    isrc = np.int(maxoff+1)
    if mode1D2D == np.int(1):
        isrc = np.int(0)
    nxx = noff
    cfvmod = quad(4)/(vmod**2)
    #%%

    #%%
    t1 = time.time()
    Pobs, Pobszero = prop.linop.flin(ximod,vmod,src,rOK,Mtap,ptap,zsrc,\
                                     dx,dh,dz,dt,vmin,vmax,PMeths,ximodzero,\
                                    stap,slist,mode1D2D,fcomm, npml,noff,ntt)
    t2 = time.time()
    print('calcPobs took %8.3f seconds' %(t2-t1))
    for isrc in np.arange(0, ns, dtype=int):
        name = folder + 'Pobs_%06i' %(slist[isrc]) + '.npy'
        np.save(name,Pobs[:,:,isrc])
    mystop(do_mpi,comm,rank)
    #%%
    mystop(do_mpi,comm,rank)
    #%%
    #S0_P = prop_py.calcS0(src,zsrc,isrc,dz,dx,noff,ntt,nz,nxx);
    # with fortran
    #print('calcS0 using fortran')
    #t1 = time.time()
    #S0 = prop.propag.calcs0(src,zsrc,isrc+1,dz,dx,noff,ntt,nz,nxx);
    #t2 = time.time()
    #print('calcP0 using fortran')
    #t3 = time.time()
    #P0 = prop.propag.calcpropag(S0,vmod,vmin,vmax,dx,dz,dt,True,mode1D2D,\
    #   npml,noff,nz,ntt)
    #t4 = time.time()
    #print('calcS1 using fortran')
    #t5 = time.time()
    #S1 = prop.propag.calckp1(P0,ximod,cfvmod,dh,dt,False,2,True,mode1D2D)
    #t6 = time.time()
    #print('calcP1 using fortran')
    #t7 = time.time()
    #P1 = prop.propag.calcpropag(S1,vmod,vmin,vmax,dx,dz,dt,True,mode1D2D,\
    #   npml,noff,nz,ntt)
    #t8 = time.time()
    #print('S0=%6.4f\nP0=%6.4f\nS1=%6.4f\nP1=%6.4f' %(t2-t1,t4-t3,t6-t5,t8-t7))
    #np.save(folder + 'S0.npy',S0)
    #np.save(folder + 'P0.npy',P0)
    #np.save(folder + 'S1.npy',S1)
    #np.save(folder + 'P1.npy',P1)
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
    #%% Safeguards for the 1D case
    if mode1D2D == np.int(1):
        # removes tapers along x and h directions for extended model vectors
        ptap[3] = quad(0)
        ptap[4] = quad(0)
        if (MethTap2 >= 4):
            MethTap2 = 4
        else:
            MethTap2 = 0
            MethTap  = 0
        if (MethTpR2 >=4):
            MethTpR2 = 4
        else:
            MethTpR2 = 0
            MethTpRg = 0
        # no tapers on the acquisition
        MethTapS = 0
        if (MethTapD >= 2):
            MethTapD = MethTapD - 2
    #%%
    PMeths = np.array([MethJxi , MethXi, MethTap ,MethTap2, MethMR,
                       MethTapD, MethQ , MethDcnv,MethAcq2], dtype=int)
    #%% spatial dimensions
    zmin = quad(   0) ; zmax=quad( 500) ;  dz=quad(6) ; nz=np.int(np.floor((zmax-zmin)/dz))
    xmin = quad(   0) ; xmax=quad(1600) ;  dx=quad(6) ; nx=np.int(np.floor((xmax-xmin)/dx))
    hmin = quad(-180) ; hmax=-hmin      ;  dh=dx      ; nh=np.int(np.floor((hmax-hmin)/dh))
    if mode1D2D == np.int(1):
        xmin = quad(0) ; xmax = quad(0) ; dx = quad(1) ; nx = np.int(1)
        hmin = quad(0) ; hmax = quad(0) ; dh = quad(1) ; nh = np.int(1)
    #%% temporal dimension
    tmin = quad(0) ; tmax=quad(0.5) ; dt=quad(0.0013) ; nt=np.int(np.floor((tmax-tmin)/dt))
    #%% simulation of wave propagation around source bewteen -maxoffR and +maxoffR
    maxoffR = quad(600) # in meters
    if mode1D2D == np.int(1):
        maxoffR = quad(0)
    maxoff = np.int(np.floor(maxoffR/dx))
    noff   = 2*maxoff + 1
    #%% acquisition : position of sources and receivers
    # vertical position of sources
    izsrc = np.int(1) ;  zsrc = izsrc * dz
    # lateral position of sources
    # typically every 'ds' between xSlim[0] and xSlim[1] and between xSlim[2] and xSlim[3]
    xSlim = np.array([xmin, 450., 450., xmax], dtype=myfloat)
    dsint = np.int(8)
    if mode1D2D == np.int(1):
        dsint = np.int(1)
    ds = quad(dsint)
    sOKs = False # forces symmetric sources repartition
    sOKf = False # forces sources to be at xSlim[1] and xSlim[2]
    # lateral positions of the sources on the grid
    sOK, nstotal = acquisition.okSource(xSlim,dx,nx,dsint,sOKs,sOKf,mode1D2D)
    # repartition of sources between processors
    slist, stap = acquisition.repartSource(sOK,nb_procs,rank)
    ns = slist.shape[0]
    ## careful : slist starts at 0, maybe a problem when passing to Fortran...
    # receivers position
    xRlim = np.array([-600., 600.], dtype=myfloat)
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
        print('----------------------------------------------------------', flush=True)
        print('     nz = %4i ; dz = %5.3f m' %(nz,dz), flush=True)
        print('     nx = %4i ; dx = %5.3f m' %(nx,dx), flush=True)
        print('     nh = %4i ; dh = %5.3f m' %(nh,dh), flush=True)
        print('     nt = %4i ; dt = %5.3f ms (positive time samples)' %(nt,dt*1000), flush=True)
        print('    ntt = %4i (neg. and pos. time samples)' %(ntt), flush=True)
        print('   nsrc = %4i (number of time samples for src wavelet)' %(nsrc), flush=True)
        print('nstotal = %4i ; ds = %2i (total number of sources in the acquisition)' %(nstotal,dsint), flush=True)
        print('     ns = %4i (number of sources treated by this processor)' %(ns), flush=True)
        print('   nrcv = %4i' %(nrcv), flush=True)
        print('   noff = %4i' %(noff), flush=True)
        print('----------------------------------------------------------', flush=True)
    #%% stores param in a python-readable file
    if rank == np.int(0):
        f = open(folder + 'param.pckl', 'wb')
        pickle.dump([nz, nx, nh, nt, ntt, nsrc, nstotal, ns, nrcv, noff,\
                     dz, dx, dh, dt, dsint,\
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
    ok = prop_py.verifCond(dz, dx, dt, fmax, vmin, vmax, mode1D2D, rank)
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
    #isrc = np.int(maxoff+1)
    #if mode1D2D == np.int(1):
        #isrc = np.int(0)
    #nxx = noff
    #cfvmod = quad(4)/(vmod**2)
    #%%

    #%% Compute observed data
    t1 = time.time()
    Pobs, Pobszero = prop.linop.flin(ximod,ximodzero,vmod,src,rOK,Mtap,ptap,\
                                     zsrc,dx,dh,dz,dt,vmin,vmax,PMeths,\
                                     stap,slist,mode1D2D,fcomm, npml,noff,ntt)
    t2 = time.time()
    if rank == np.int(0):
        print('calcPobs took %8.3f seconds' %(t2-t1), flush=True)
    #%% Migration (using adjoint operator Fadj)
    t1 = time.time()
    xiadj, xizero = prop.linop.fadj(Pobs,Pobszero,vini,src,rOK,Mtap,ptap,zsrc,\
                                 dx,dh,dz,dt,vmin,vmax,PMeths,stap,slist,\
                                 mode1D2D,fcomm,npml,noff,nh)
    t2 = time.time()
    if rank == np.int(0):
        print('Fadj took %8.3f seconds' %(t2-t1), flush=True)
    #%% Migration (using adjoint operator Fdag)
    t1 = time.time()
    xiinv, xizero = prop.linop.fdag(Pobs,Pobszero,vini,srcdcnv,rOK,Mtap,ptap,zsrc,\
                                 dx,dh,dz,dt,ds,vmin,vmax,PMeths,stap,slist,\
                                 mode1D2D,fcomm,npml,noff,nh)
    t2 = time.time()
    if rank == np.int(0):
        print('Fdag took %8.3f seconds' %(t2-t1), flush=True)

    #%% Save results
    #--- Pobs ---
    for isrc in np.arange(0, ns, dtype=int):
        if mode1D2D == np.int(1):
            name = folder + 'Pobs' + '.npy'
            np.save(name,Pobs)
        if mode1D2D == np.int(2):
            name = folder + 'Pobs_%06i' %(slist[isrc]) + '.npy'
            np.save(name,Pobs[:,:,isrc])
    #--- xi adjoint ---
    if (rank == np.int(0)):
        np.save(folder + 'xiadj.npy',xiadj)
        np.save(folder + 'xiinv.npy',xiinv)
    #%% demigration
    print('Demigration...',flush=True)
    Fxiadj, Pzero = prop.linop.flin(xiadj,xizero,vini,src,rOK,Mtap,ptap,\
                                     zsrc,dx,dh,dz,dt,vmin,vmax,PMeths,\
                                     stap,slist,mode1D2D,fcomm, npml,noff,ntt)
    Fxiinv, Pzero = prop.linop.flin(xiinv,xizero,vini,src,rOK,Mtap,ptap,\
                                     zsrc,dx,dh,dz,dt,vmin,vmax,PMeths,\
                                     stap,slist,mode1D2D,fcomm, npml,noff,ntt)
    for isrc in np.arange(0, ns, dtype=int):
        if mode1D2D == np.int(1):
            np.save(folder + 'Fxiadj.npy',Fxiadj)
            np.save(folder + 'Fxiinv.npy',Fxiinv)
        if mode1D2D == np.int(2):
            name = folder + 'Fxiadj_%06i' %(slist[isrc]) + '.npy'
            np.save(name,Fxiadj[:,:,isrc])
            name = folder + 'Fxiinv_%06i' %(slist[isrc]) + '.npy'
            np.save(name,Fxiinv[:,:,isrc])
    #%%
    #mystop(do_mpi,comm,rank)
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

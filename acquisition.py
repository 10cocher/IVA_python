# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 15:20:43 2017

@author: Emmanuel
"""

import numpy as np
from prec import quad, myfloat
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt

def okSource(xSlim,dx,nx,ds,sOKs,sOKf,mode1D2D):
    """ Defines the lateral position of sources
    Parameters:
    ----------
    xSlim : float(4), sources are located between xSlim[0] and xSlim[1]
    AND between xSlim[2] and xSlim[3]
    dx : float, grid spacing (in direction x)
    nx : int, number of x-sample
    sOKs : boolean, if True, ensures a symmetric repartition of sources
    sOKf : boolean, if True, ensures a source at xSlim[1] and xSlim[2]
    mode1D2D : int, (1) 1D acquisition (2) 2D acquisition

    Returns:
        sOK (boolean, array(nx)) : if sOK(ix) is True, a Source is located at ix
        ns : number of sources in the acquisition
    """
    # in 1D, this is easy
    if mode1D2D == np.int(1):
        sOK = np.array([1], dtype=bool, order='F')
    # in 2D more work is needed
    if mode1D2D == np.int(2):
        eps = dx/quad(2)
        ixx = np.arange(0, nx, dtype=int)
        xx = dx * ixx
        # sources between xSlim[0] and xSlim[1] or xSlim[2] and xSlim[3]
        sOK = np.logical_and(\
               np.mod( ixx + 1 , ds) == 0 ,\
               np.logical_or(\
                np.logical_and( xx > xSlim[0] - eps , xx < xSlim[1] + eps ) ,\
                np.logical_and( xx > xSlim[2] - eps , xx < xSlim[3] + eps ) \
               )
              ) # ixx+1 for compliance with Fortran code
        # if sOKf, forces a source at xSlim[1] and xSlim[2]
        if sOKf:
            is1 = np.int( np.floor(xSlim[1]/dx))
            is2 = np.int( np.ceil (xSlim[2]/dx))
            sOK[is1] = True
            sOK[is2] = True
        # if sOKs, forces a symmetric (in x) repartition of sources
        nx2 = np.int( np.floor_divide(nx,quad(2)) )
        if sOKs:
            sOK[ nx-nx2 : nx ] = np.flip( sOK[0:nx2], 0 )
    # number of sources in the acquisition geometry
    nstotal = np.count_nonzero(sOK)
    return sOK, nstotal

def repartSource(sOK,nb_procs,rang):
    """ Allocates sources to processors
    Parameters:
        sOK : boolean array(nx), lateral position of sources
        nb_procs : number of processors
        rang : rank of the processor
    """
    nx = sOK.shape[0]
    nstotal = np.count_nonzero(sOK)
    dstotal = quad(1)/quad(nstotal)
    #
    staptotal = dstotal * np.ones((nx,), dtype=myfloat, order='F')
    #
    quotient, reste = np.divmod(nstotal, nb_procs)
    if (rang+1 > reste):
        ns = quotient
        ismin = reste*(ns+1) + (rang - reste)*ns #+1
        ismax = reste*(ns+1) + (rang - reste +1)*ns
    else:
        ns = quotient + 1
        ismin = rang    *ns #+ 1
        ismax = (rang+1)*ns
    #
    src_list    = np.arange(nx, dtype=int) + np.int(1) # +1 because Fortran indexing starts at 1
    src_vraies  = src_list[sOK]
    stap_vraies = staptotal[sOK]
    #
    slist = src_vraies[ismin:ismax]
    stap = stap_vraies[ismin:ismax]
    return slist, stap

def okReceiver(xRlim,xSlim,zVSP,xmin,xmax,dz,dx,nrcv,slist,MethAcq,MethAcq2,mode1D2D):
    """ Defines the position of receivers
    Parameters:
    ----------
    xRlim : float(2)
    xSlim : float(4)
    zVSP  : float(2) receivers located between zVSP[0] and zVSP[1]
    dz    : float, grid spacing (in direction z)
    dx    : float, grid spacing (in direction x)
    nrcv  : int, number of receivers
    xmin  : float, limit min of the domain along x-axis
    xmax  : float, limit max of the domain along x-axis
    MethAcq2 : int, (0) = surface acquisition
                    (1) = VSP acquisition
    MethAcq  : int, in the case of surface acquisition (MethAcq2=0)
              (0) receivers are located between xRlim[0] and xRlim[1]
              (1) receivers everywhere except between xSlim[1] and xSlim[2]
    mode1D2D : int, (1) 1D acquisition (2) 2D acquisition

    Returns:
        rOK (boolean, array(noff,ns)) : if rOK(ioff,isrc) is True, a receiver is
                                        located at position ioff for source isrc
    """
    ns = slist.shape[0]
    # in 1D, this is easy
    if mode1D2D == np.int(1):
        rOK = np.ones((nrcv,ns), dtype=bool, order='F')
        print('in 1D, shape(rOK)=',rOK.shape)
    # in 2D more work is needed
    if mode1D2D == np.int(2):
        rOK = np.zeros((nrcv,ns), dtype=bool, order='F')
        if MethAcq2 == np.int(0): # surface acquisition
            eps = dx/2.
            noff = nrcv
            maxoff = (noff-1)/2
            # defines pos_off = position in meters relative to the source
            pos_off = dx * np.arange(-maxoff, maxoff+1, dtype=int)
            for iis in np.arange(0, ns, dtype=int):
                # Initialization
                isrc    = slist[iis]
                # defines pos_x = absolute lateral position in meters
                pos_x   = dx * np.arange(isrc - maxoff, isrc + maxoff + 1, dtype=int)
                # rOK is false outside [xmin,xmax]
                rOKline = np.logical_and(pos_x > xmin - eps, pos_x < xmax + eps)
                # Case MethAcq=(0) receivers only between xRlim[0] and xRlim[1]
                if MethAcq == np.int(0):
                    rOKline = np.logical_and(rOKline,\
                                             np.logical_and(pos_off > xRlim[0]-eps,\
                                                            pos_off < xRlim[1]+eps)
                                             )
                # Case MethAcq=(1) receivers everywhere except between xSlim[1] and xSlim[2]
                if MethAcq == np.int(1):
                    rOKline = np.logical_and(rOKline,\
                                             np.logical_and(pos_x < xSlim[1]-eps,\
                                                            pos_x > xSlim[2]+eps)
                                             )
                # fill big array
                rOK[:,iis] = rOKline
        if MethAcq2 == np.int(1): # VSP acquisition
            nz = nrcv
            zax = dz * np.arange(0, nz, dtype=int)
            for iis in np.arange(0, ns, dtype=int):
                rOKline = np.logical_and(zax > zVSP[0]-eps,\
                                         zax < zVSP[1]+eps)
                rOK[:,iis] = rOKline
    return rOK
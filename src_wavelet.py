# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 13:58:55 2017

@author: Emmanuel Cocher
"""

import numpy as np
import matplotlib.pyplot as plt
import mymaths

def defsrc(ns,fmax,dt):
    """ Creates a Ricker wavelet
    Parameters
    ----------
    ns : int, number of (time) samples for the source
    fmax : maximum frequency of the Ricker
    dt : value of the time sample
    """
    ns2=(ns-np.int(1))/np.int(2)
    src=np.zeros((np.int(ns),), dtype=float, order='F')
    fc=fmax/2.5
    for i in np.arange(-ns2,ns2+1,dtype=int): 
        a1=i*dt*fc*np.pi
        a2=a1*a1
        idx=i+np.int(ns2)
        src[idx]=(1-2*a2)*np.exp(-a2)
    return src

def defsrcdcnv(src,dt):
    """ Outputs the deconvolved version of the input source wavelet
    Parameters
    ----------
    ns : number of (time) samples for the source
    fmax : maximum frequency of the Ricker
    dt : value of the time sample
    """
    ns = src.shape[0]
    ns2 = (ns-1)/2
    #
    fft_src = np.fft.fftshift(np.fft.fft(src)) # fft of src
    sqr_fft_src = np.multiply( fft_src , np.conjugate(fft_src) )
    eps = np.complex(np.amax(np.abs(sqr_fft_src))/100.)
    #
    fft_srcdcnv = np.divide( fft_src , sqr_fft_src + eps )
    srcdcnv = np.real( np.fft.ifft(np.fft.ifftshift(fft_srcdcnv)) )
    #
    smth = np.ones(ns, dtype=float, order='F')
    i1 = 4
    i2 = 8
    #smth[0:i1-1]=0.
    for i in np.arange(1, i1, dtype=int):
        smth[i-1]=0.
    for i in np.arange(i1, i2, dtype=int):
        smth[i-1] = mymaths.SmoothSin( np.float(i-i1)/np.float(i2-i1) )
    for i in np.arange(0, ns2-1, dtype=int):
        smth[ns-i-1] = smth[i]
    #
    srcdcnv = np.multiply(srcdcnv, smth) / (np.square(dt))
    if True:
        spectrum_src = np.absolute(fft_src) # spectrum of src
        spectrum_srcdcnv = np.absolute(fft_srcdcnv) # spectrum of srcdcnv
        sax = np.linspace(-ns2*dt,ns2*dt,num=ns,dtype=float)
        fax = np.linspace(1,ns,num=ns,dtype=int)
        #
        fig1 = plt.figure(num='src wavelet')
        #
        ax1 = fig1.add_subplot(2, 3, 1)
        ax1.grid(color='gray')
        ax1.plot(sax,src)
        ax1.set_xlabel('time (s)')
        #
        ax2 = fig1.add_subplot(2, 3, 2)
        ax2.grid(color='gray')
        ax2.plot(sax,srcdcnv)
        ax2.set_xlabel('time (s)')
        #
        ax2 = fig1.add_subplot(2, 3, 3)
        ax2.grid(color='gray')
        ax2.plot(fax,smth)
        ax2.set_xlabel('time (s)')
        #
        ax3 = fig1.add_subplot(2, 3, 4)
        ax3.grid(color='gray')
        ax3.plot(fax,spectrum_src)
        ax3.set_xlabel('frequency (Hz)')
        #
        ax4 = fig1.add_subplot(2, 3, 5)
        ax4.grid(color='gray')
        ax4.plot(fax,spectrum_srcdcnv)
        ax4.set_xlabel('frequency (Hz)')
        #
    return srcdcnv
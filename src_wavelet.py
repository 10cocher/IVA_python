# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 13:58:55 2017

@author: Emmanuel Cocher
"""

import numpy as np
from prec import quad, complexquad, myfloat
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
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
    src=np.zeros((ns,), dtype=myfloat, order='F')
    fc=fmax/quad(2.5)
    for i in np.arange(-ns2,ns2+1,dtype=int):
        a1=quad(i)*dt*fc*np.pi
        a2=a1*a1
        idx=i+np.int(ns2)
        src[idx]=(quad(1)-quad(2)*a2)*np.exp(-a2)
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
    eps = complexquad(np.amax(np.abs(sqr_fft_src))/quad(100))
    #
    fft_srcdcnv = np.divide( fft_src , sqr_fft_src + eps )
    srcdcnv = np.real( np.fft.ifft(np.fft.ifftshift(fft_srcdcnv)) )
    #
    smth = np.ones(ns, dtype=myfloat, order='F')
    i1 = 4
    i2 = 8
    #smth[0:i1-1]=0.
    for i in np.arange(1, i1, dtype=int):
        smth[i-1]=0.
    for i in np.arange(i1, i2, dtype=int):
        smth[i-1] = mymaths.SmoothSin( quad(i-i1)/quad(i2-i1) )
    for i in np.arange(0, ns2-1, dtype=int):
        smth[ns-i-1] = smth[i]
    #
    srcdcnv = np.multiply(srcdcnv, smth) / (np.square(dt))
    """
    if False:
        spectrum_src     = np.absolute(fft_src) # spectrum of src
        spectrum_srcdcnv = np.absolute(fft_srcdcnv) # spectrum of srcdcnv
        #
        src_conv  = np.convolve(src, srcdcnv, mode='same') * dt
        spec_prod = np.multiply(spectrum_src, spectrum_srcdcnv)
        #
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
        ax2.plot(sax,src_conv)
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
        ax5 = fig1.add_subplot(2, 3, 6)
        ax5.grid(color='gray')
        ax5.plot(fax,spec_prod)
        ax5.set_xlabel('frequency (Hz)')
        #
    """
    return srcdcnv
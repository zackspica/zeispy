#!/usr/bin/env python
import numpy as np
import math
import scipy.fftpack
from pyrocko.trace import moving_sum

def moving_avg(x,n):
    """Efficient  moving avg
    
    :type x: numpy.array
    :param x: the trace
    :type n: int
    :param n: window lenght #points
    :return: the new trace of same lenght
    """
    n = int(n)
    cx = x.cumsum()
    nn = len(x)
    y = np.zeros(nn, dtype=cx.dtype)
    y[n/2:n/2+(nn-n)] = (cx[n:]-cx[:-n])/n
    y[:n/2] = y[n/2]
    y[n/2+(nn-n):] = y[n/2+(nn-n)-1]
    return y


def spectral_normalization(data):
    """Compute simple Whitening of a trace.

    :type data: numpy.array
    :param data: the data
    :return: specrtral normalized data
    """
    spec = np.fft.rfft(data)
    spec /= np.abs(spec)
    data = np.fft.irfft(spec)
    del spec
    return data


def temp_norm_fast(data, N, safe=True):
    """ Computes the running-absolute-mean normalization of a the data, the trace.
    This function computes the running average of the absolute value of the 
    data in a normalization time window of fixed length (2N+1) and 
    weights the data at the center of the window by the inverse of this 
    average. As in Bensen et al., 2007. 
    A one-sample window (N = 0) is equivalent to one-bit normalization.

    :type data: numpy.array
    :param data: the data
    :type N: int  
    :param N: Number of point. Folowing Bensen2008: works well if N = half of the maximum 
        period of the bandpass filter but can vary considerabely and still produce similar 
        results. 2N+1 = width of the normalization window.
    :type safe: bool 
    :param safe: If False the function doesn't protect against division by zero but will 
        be about 4% run faster.
    :return: temporal normalized data
    """
    ydata = data.astype(np.float)
    mavg = moving_sum(np.abs(ydata), 2*N+1, mode='same')/(2*N+1)
    if safe:
        eps = mavg.max()*1e-9
        if eps == 0.0:
            eps = 1e-9
        ydata /= (mavg+eps)
    else:
        ydata /= mavg 
    return ydata


def corrconv(vb, va, eps=0.01):
    """ This function is used to get the cross correlation with deconvolution
    vb is the receiver and va is the virtual source 
    This function follow the processing of Nori.
    Note: prob of lenght... not full! Instable and must be improved 
    Last update 03/11/16 has to be rechecked. 

    :type vb: numpy array
    :param vb: data to correlate | receiver
    :type va: numpy array
    :param vb: data to correlate | virtual source
    :type eps: float 
    :param eps: epsilon used in the division
    :return: correlation of len(va)*2-1
    """
    l = len(va)*2 -1
    fva = np.fft.fft(va, l)#.astype(np.float32)
    fvb = np.fft.fft(vb, l)#.astype(np.float32)
 
    out = fvb * np.conj(fva)  
    #deno =  fva * np.conj(fva) + eps * np.mean(fva * np.conj(fva)) 
    bautocorr = fva * np.conj(fva)
    deno = np.maximum( bautocorr, eps * bautocorr.max() )
    
    out /= deno
    ccov = np.real( scipy.fftpack.fftshift( scipy.fftpack.ifft( out ) ) )
    
    del fva , fvb

    return ccov 


def correlation(vb, va, normalized=False):
    """Simple cross-correlation in the frequency domain 
    
    :type vb: numpy array
    :param vb: data to correlate | receiver
    :type va: numpy array
    :param vb: data to correlate | virtual source
    :return: correlation of len(va)*2-1
    """
    l = len(va)*2 -1
    fva = scipy.fftpack.fft(va, l)
    fvb = scipy.fftpack.fft(vb, l)

    cc =  fvb * np.conj(fva) 
    cc = np.real( scipy.fftpack.ifftshift( scipy.fftpack.ifft( cc ) ) )

    if normalized:
        # i have to dubble check that
        E = np.real(np.sqrt(np.mean(scipy.fftpack.ifft(fva))**2 * np.mean(scipy.fftpack.ifft(fvb))**2))
        cc /= E
        pass

    return cc
    
    
def crosscoherence(vb, va, eps=0.01, test=False):
    """Simple cross-coherence in the frequency domain 

    :type vb: numpy array
    :param vb: data to correlate | receiver
    :type va: numpy array
    :param vb: data to correlate | virtual source
    :type eps: float 
    :param eps: epsilon used in the division
    :return: correlation of len(va)*2-1
    """
    l = len(va)*2 - 1
    fva = scipy.fftpack.fft(va, l)
    fvb = scipy.fftpack.fft(vb, l)
    
    div = np.abs(fva) * np.abs(fvb) + eps * np.average(np.abs(fva) * np.abs(fvb))
    ccohe =  fvb * np.conj(fva) / div
    if test:
        if np.isnan(ccohe).any():
            raise FloatingPointError("/!\ Cross-Coherence: zeros encountered in division. Epsilon must be changed ")
    ccohe = np.real(scipy.fftpack.fftshift(scipy.fftpack.ifft(ccohe)))
    del div, fva, fvb

    return ccohe


def convolution(vb, va, npt=30., eps=1e-6):
    """ This function is used to get the cross correlation with deconvolution
    vb is the receiver and va is the virtual source
    npt is the number used for running window average
    Marine Denolle way (Yixiao)

    :type vb: numpy array
    :param vb: data to correlate | receiver
    :type va: numpy array
    :param vb: data to correlate | virtual source
    :type eps: float 
    :param eps: epsilon used in the division
    :return: correlation of len(va)*2-1
    """
    Len = len(va)*2. -1.
    fva = scipy.fftpack.fft(va, Len)
    fvb = scipy.fftpack.fft(vb, Len)

    fva2 = RWA(fva,npt); fvb2 = RWA(fvb, npt)
    p1 = np.ones(fva.size)*eps
    p2 = np.ones(fvb.size)*eps

    for i in xrange(int(Len)):
        p1[i] = np.sum(fva2[i:i+npt])/npt
        p2[i] = np.sum(fvb2[i:i+npt])/npt
    p1 += eps
    p2 += eps
    ANIRWF = fvb * np.conj(fva) / (p1*p2)    
    ANIRW = np.real( np.fft.fftshift(np.fft.ifft(ANIRWF)) )

    return ANIRW

    
def nextpow2(x):
    """
    Return the next power of 2 of `x`.
    :type x: int
    :param x: any value
    :rtype: int
    :returns: the next power of 2 of `x`
    """
    return np.ceil(np.log2(np.abs(x)))


def RWA(v,npt=10):
    """ This function is used to get the running window average
        should be improved
    """
    Len = int(len(v))
    npt = int(npt)
    v2 = np.zeros(Len + 2 * npt)# + np.zeros(Len + 2*npt)*0j
    
    v2[0:npt] = (v[0:npt])[::-1]
    v2[-npt:] = (v[Len-npt:])[::-1]
    v2[npt:npt+Len] = v
    v2 = np.abs(v2)

    return v2


if __name__ == "__main__":
	
	data1 = np.random.random((1000,))
	data2 = np.random.random((1000,))
	corr = corrconv(data2,data1,10)
	print np.mean(corr)
	plt.figure()
	plt.plot(corr)
	plt.show()


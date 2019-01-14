# -*- coding: utf-8 -*-
"""
Created on Weds 04 06 21:03:06 2016
@author: Zack Spica zspica@stanford.edu
Last Update: April 07
"""
import matplotlib.pyplot as plt
import numpy as np
from obspy.core import read, Trace, Stats
try:from obspy.signal.invsim import cosine_taper
except:from obspy.signal.invsim import cosTaper as cosine_taper
from zeispy.mycorr import *
from zeispy.stack import stack
#from obspy.signal import detrend






def getCodas(datapair, distpair, df, appvel, Xcoda, coda_len,\
                 detrend=True, onebit=False, plot=False):
    """
    Get the codas of a cross-correlation functions (C1)

    :param datapair: pair of traces A and B
    :type datapair: list of two `numpy.array`
    :param distpair: contain the two distance between S-A and S-B.
    :type distpair: list of two float
    :param df: sampling rate in Hz
    :type df: float
    :param appvel: Apparent velocity of the surface wave (typically 3 but can be much slower).
    :type appvel: float
    :param Xcoda: How many time the theoretical arriving of the coda (typically 2 as in Stehly et al., 2008)
    :type Xcoda: float
    :param coda_len: length of the coda window (typically = 800s).
                Can raise an error if too long or C1 too short.
    :type coda_len: int, float
    :param detrend: whether or not to detrend the C1. Use a polynomial 2nd order.
    :type detrend: bool
    :param onebit: whether or not to apply one bit normalization
    :type onebit: bool
    :param plot: whether or not to plot:
    :type plot: bool

    :rtype: list of :class:`numpy.array`
    :returns: list of the 4 codas in the following order: left A, right A, left B, right B.
    """

    mid = len(datapair[0])/2
    c1l, c1r = datapair[0][:mid][::-1], datapair[0][mid+1:]
    c2l, c2r = datapair[1][:mid][::-1], datapair[1][mid+1:]

    coda_len = int(coda_len * df)
    codaStartA = int((distpair[0] / appvel * Xcoda) * df )
    codaStartB = int((distpair[1] / appvel * Xcoda) * df )
    csmx = np.max([codaStartA, codaStartB])

    l1 = c1l[csmx:coda_len]
    r1 = c1r[csmx:coda_len]
    l2 = c2l[csmx:coda_len]
    r2 = c2r[csmx:coda_len]

    if onebit:
        l1 = np.sign(l1)
        r1 = np.sign(r1)
        l2 = np.sign(l2)
        r2 = np.sign(r2)

    l1 -= np.mean(l1)
    r1 -= np.mean(r1)
    l2 -= np.mean(l2)
    r2 -= np.mean(r2)

    #if plot:
    #    plt.plot(l1)
    #    plt.plot(r1)
    #    plt.plot(l2)
    #    plt.plot(r2)
    #    plt.show()

    return l1, r1, l2, r2


def corrCodas(a, b, shortWin, overlap, df=20, corr_method='crosscoherence', stack_method='linear'):

    """Moving window correlation function
    :rtype: numpy.array
    :param a, b: the codas to be correlated
    :type a, b: numpy.array
    :param shortWin: short time window for moving wind
    :type shortWin: int, float
    :param step: step to iterate (i.e., for overlap)
    :type step: int
    :rparam: the correlated function of lenght equal to shortWin*2-1.
    """
    N = len(a)
    c = np.zeros(shortWin*2-1)
    cp = cosine_taper(shortWin, 0.1)
    data = np.zeros(len(c))

    for i in xrange(0, N-shortWin, overlap):
        winA = a[i: i+shortWin]
        winB = b[i: i+shortWin]
        cp *= winA; cp*=winB

        if corr_method == 'crosscoherence':
            cc = crosscoherence(winB, winA)
        elif corr_method == 'correlation':
            cc = correlation(winB, winA)
        elif corr_method == 'convolution':
            cc = convolution(winB, winA)

        data = np.vstack((data,cc))

    corr = stack(data, stack_method=stack_method, df=df)

    return corr#/np.max(np.abs(c))



def getSlices(coda_len, shortWin, overlap, df):
    """
    get slices for moving window.

    :param coda_len: lenght of the coda window (typically = 800). Can raise an error if too long or C1 too short.
    :type coda_len: int, float
    :param shortWin:
    :type shortWin: int, float
    :param overlap: percentage of overlap during the moving window sliding (0-1).
    :type overlap: float
    :param df: sampling rate in Hz
    :type df: float.
    """
    begins, ends = [], []
    i = 0
    while i <=  (coda_len - shortWin):
        begins.append( int(i * df) )
        ends.append( int(i * df + shortWin * df) )
        i += int( shortWin * (1.0 - overlap) )
    print 'Slices: ', begins, ends
    return begins, ends


def polynomial(data, order):
    """
    Removes a polynomial trend from the data.

    :param data: The data to detrend. Will be modified in-place.
    :type data: :class:`numpy.ndarray`
    :param order: The order of the polynomial to fit.
    :type order: int
    !Modified from Obspy.
    """
    # Convert data if it's not a floating point type.
    if not np.issubdtype(data.dtype, float):
        data = np.require(data, dtype=np.float64)
    x = np.arange(len(data))
    fit = np.polyval(np.polyfit(x, data, deg=order), x)
    data -= fit
    return data


def simple(data):
    """
    Detrend signal simply by subtracting a line through the first and last
    point of the trace.

    :param data: Data to detrend, type numpy.ndarray.
    :return: Detrended data. Returns the original array which has been
        modified in-place if possible but it might have to return a copy in
        case the dtype has to be changed.
    !Function Obspy.
    """
    # Convert data if it's not a floating point type.
    if not np.issubdtype(data.dtype, float):
        data = np.require(data, dtype=np.float64)
    ndat = len(data)
    x1, x2 = data[0], data[-1]
    data -= x1 + np.arange(ndat) * (x2 - x1) / float(ndat - 1)
    return data

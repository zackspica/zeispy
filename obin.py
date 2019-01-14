import numpy as np
import os, time, glob
#from obspy.core import UTCDateTime
from pyrocko import model
from zeispy.binplot import *

def raw_open(file, dt=0.05, shape=None, test=False, plot=False, saveplot=None):
    """Open a binary or npy file. If it has no header it will 
    reshape it as a `numpy.ndarray`, generaly 2D matrix. If shape is not 
    specified, the function try to reshape the input array
    based on the lenght of the stations.txt in the folder.
    
    :param file: path and file name
    :type file: str
    :param dt: deltat
    :type dt: float
    :param shape: shape of the output `numpy.ndarray`
    :type shape: tulpe 
    :params test: wether or not to print the size of the array in order to difine shape
    :type test: bool
    :param plot: whether or not to plot the output
    :type plot: bool
    
    :rtype: :class:`numpy.ndarray`
    :returns: reshaped bin
    """
    if test == True:
        print np.shape(np.fromfile(file, dtype=np.float32))
        return
    if file.endswith('.npy'):
        bin = np.load(file)
        x = np.shape(bin)[0]
        y = np.shape(bin)[1]
        shape = (x, y)
    else:
        bin = np.fromfile(file, dtype=np.float32)
        if type(shape) is tuple:
            pass
        elif shape == None and os.path.isfile('stations.txt'):
            x = len(model.load_stations('stations.txt'))
            y = len(bin)/x
            shape = (x, y)
        else:
            raise ShapeError('No stations.txt found and no shape specified!')
        bin = np.reshape(bin, shape)

    if plot:
        raw_plotshot(shape[0], shape[1], bin, dt, show=plot, save=saveplot )
 
    return bin


def ropen(file, dt=0.05, shape=None):#, plot=False, saveplot=None):
    """Same as raw_open but return all arguments (shape, dt) for plot facilities.
    Open a binary or npy file. If it has no header it will 
    reshape it as a `numpy.ndarray`, generaly 2D matrix. If shape is not 
    specified, the function try to reshape the input array
    based on the lenght of the stations.txt in the folder.

    :param file: path and file name
    :type file: str
    :param dt: deltat
    :type dt: float
    :param shape: shape of the output `numpy.ndarray`
    :type shape: tulpe 
    #:param plot: whether or not to plot the output
    #:type plot: bool

    :rtype: :class:`numpy.ndarray`
    :returns: reshaped bin, shape1, shape2, dt
    """
    if file.endswith('.npy'):
        bin = np.load(file)
        #x = np.shape(bin)[0]
        #y = np.shape(bin)[1]
        #shape = (x, y)
    else:
        bin = np.fromfile(file, dtype=np.float32)
        if type(shape) is tuple:
            pass
        elif shape == None and os.path.isfile('stations.txt'):
            x = len(model.load_stations('stations.txt'))
            y = len(bin)/x
            shape = (x, y)
        else:
            raise Exception('No stations.txt found and no shape specified!')
        
        bin = np.reshape(bin, shape)
    #if plot:
    #    raw_plotshot(shape[0], shape[1], bin, dt, show=plot, save=saveplot )

    return bin#, shape[0], shape[1], dt   


def mergbins(bin1, bin2, dt=0.05, shape1=None, shape2=None):
    """ merge to sets of bins file together
        
    :param bin1: npy or bonary file with all the traces
    :type: str
    :param bin2: npy or bonary file with all the traces
    :type: str
    """
    fn1 = os.path.basename(bin1)[:-4]
    fn = fn1+'.merged'
    a = raw_open(bin1, dt=dt, shape=shape1)
    print fn1, 'has a shape of:', np.shape(a)
    b = raw_open(bin2, dt=dt, shape=shape2)
    fn2 = os.path.basename(bin2) 
    print fn2, 'has a shape of:', np.shape(b)
    #c = np.concatenate((a.T,b))
    c = np.concatenate((a, b))
    print '>>> The new file has a shape of:', np.shape(c)
    print '>>> Saving ', fn
    
    try: os.makedirs('bins/full')
    except: pass
    np.save('bins/full/%s'%fn,c.astype(np.float32))
    #os.system('mv bins/full/%s.npy bins/full/%s.bin'%(fn,fn))


def stack(nobin1, nobin2):
    pass
     

if __name__=='__main__':
    g = sorted(glob.glob('/data/beroza/zack/Groningen/si_average/*bin'))
    for nf in g: 
        print nf
        #fi = nf.split('/')[-1][:-4]
        #f2 = 'bins/Z/%s.bin'%fi
        #mergbins(nf, f2)
        raw_open(nf, dt=0.05, shape=(407,1953600/407), plot=True)
        break

import numpy as np
from numpy.ma import masked_where, masked_values
from matplotlib import cm, colors
import matplotlib.pyplot as plt
from obspy.signal.filter import *

def raw_plotshot(x, y, z, dt, cmap=cm.gray, show=True, save=None):
    """Plot pcolormesh of a binshot i.e. one station correlated with
    all the other int the list, ordered as in stations.txt.

    :param x: lenght of the x axis, typically bin.shape[0]
    :type x: int, float
    :param y: lenght of the y axis, typically bin.shape[1]
    :type y: int, float
    :param z: array to be plotted, the bin
    :type z: :class:`numpy.ndarray`
    :param cmap: colormap
    :type cmap: :class:`matplotlib.cm.cmap`
    :param show: whether or no to show the plot
    :type show: bool
    :params save: file basename for savinf. If exist the plot is saved
    :type save: str
    """
    X = np.arange(-y/2., (y/2.)) * dt
    Y = np.arange(x)
    np.meshgrid(X,Y)
    nz = masked_where(z==0., z)
    #norm = colors.Normalize(vmin=np.min(z), vmax=1, clip=False)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    #fmin =0.1 
    #fmax = 0.25
    #for i, u in enumerate(z):
    #    u = bandpass(u,fmin,fmax,1./dt)
    #    u /= np.max(u)
    #    z[i] = u
    #z /= np.max(z)
    ax.pcolormesh(X, Y, z, cmap=cmap, vmin=-0.0001, vmax=0.0001 )
    ax.set_xlim(X.min(),X.max()); ax.set_ylim(Y.min(),Y.max())
    #ax.set_xlim(5,60); ax.set_ylim(Y.min(),Y.max())
    ax.set_xlabel('lag time [$s$]')
    ax.set_ylabel('Station #')

    if save != None:
        filename = save+'.png'
        plt.savefig(filename, dpi=300, format='png')
    if show:
        plt.show()



def sorted_plotshot(x,y,z, dt, sourcesta, cmap=cm.gray, show=True, save=None):
    """to be continued""" 
    from pyrocko import model, orthodrome
    from obspy.geodetics.base import gps2dist_azimuth as gps2dist
    dists = []
    dists2 = []
    stations=model.load_stations('stations.txt')
    for a in stations:
        if a.station == sourcesta:
            break
  #  for s in stations[2602:]:
    for b in stations[:2602]:
        d = gps2dist(a.lat, a.lon, b.lat, b.lon)[0]/1000.
    #    di = orthodrome.distance_accurate50m(a, b)/1000.
#        print d[0],di
        dists.append(d)
    #    dists2.append(di)
#    dists = np.array(dists)[0] 
    
    np.set_printoptions(threshold=np.nan)
    #print  np.max(dists), len(dists)
    a = np.argsort(dists)
    #b = np.argsort(dists2)
    #print a[:100],b[:100]


    z = z[a]
    fmin = 0.1
    fmax = 0.25
    #for u in z:
    #    u = bandpass(u,fmin,fmax,1./dt)
    #    plt.plot(u)
    ##plt.plot(z[111])
    #    plt.show()
    np.save('tmp.npy', z) 
    raw_plotshot(x,y,z, dt, save='tmp')
    #3zip(*sorted(zip())

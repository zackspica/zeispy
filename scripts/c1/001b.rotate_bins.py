import numpy as np
import os, glob, time
from pyrocko import model
from obspy.signal.rotate import * 

path = 'bins/NL/'

days = np.arange(1,24,1)
stations = model.load_stations('stationsNL.txt')
a = np.loadtxt('angCorrectionG-array.txt', dtype={'names': ('names', 'ags'),'formats': ('S4','f')})
stag = a['names']
angles = a['ags']

#print angles

for d in days:
    
    E = np.load(path + 'daily.NL.2017.%03d.20sps.E.npy'%d)
    N = np.load(path + 'daily.NL.2017.%03d.20sps.N.npy'%d)
    for i, (e, n) in enumerate(zip(E,N)):
        starget = stations[i]
        for iss, s in enumerate(stag):
            if s == starget.station:
                ba = angles[iss]
                r, t = rotate_ne_rt(n, e, ba)
                E[i] = t
                N[i] = r
    print d
    np.save(path + 'daily.NL.2017.%03d.20sps.T.npy'%d, E)
    np.save(path + 'daily.NL.2017.%03d.20sps.R.npy'%d, N)

        
        
        


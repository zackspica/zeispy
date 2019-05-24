import os
from stack import stack
import glob
import numpy as np
from obspy import read

if __name__=='__main__':
    path = 'pc3/EE/'
   # ss = ['235716','235896', '236367']
    ss = [os.path.basename(s) for s in glob.glob('pc3/EE/*')]
    print ss
    for s1 in ss:
        g = glob.glob('%s%s/*'%(path, s1))
        print g
        stations = []
        for s in g:
            stations.append(os.path.basename(s))
        print stations
        mat = np.array([]) 
        for sta in stations:
            try:
                folder = '%s%s/%s/'%(path,s1,sta)
                g = glob.glob(folder+'/*')
                s2 =sta 
                print s1, s2
                for i, f in enumerate(g):
                    print f
                    tr = read(f, format='MSEED')[0]
                    data = tr.data
                    if i ==0:
                        mat = data
                    else:
                        mat = np.vstack((mat, data))
                cc = stack(mat)
                tr.data = cc
                try:
                    os.makedirs('cc_average_pc3/EE/%s'%(s1))
                except:pass
                tr.write('cc_average_pc3/EE/%s/%s.%s.mseed'%(s1,s1,s2), format='MSEED')
            except: pass

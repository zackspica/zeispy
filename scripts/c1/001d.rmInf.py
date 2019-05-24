import numpy as np
import glob, os
comps = ['Z','E','N']
for comp in comps:
    g = glob.glob('/data/beroza/zspica/BO/bins/full/*.%s.npy'%comp)
    for iff, f in enumerate(g):
        i=0
        #if iff <760: continue
        try:
            b = np.load(f)
        except:
            continue
        for il,  l in enumerate(b):
            try:
                if np.mean(l) > 10**6: 
                    b[il] = np.zeros(len(l))
                    i+=1
                if np.isinf(l[0]):
                    b[il] = np.zeros(len(l))
                    i+=1
                if np.isnan(l[0]):
                    b[il] = np.zeros(len(l))
                    i+=1
            except:continue
    #            i+=1
        print iff,i, f
        np.save(f, b)    

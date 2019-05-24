import numpy as np 
import os, glob

compos = ['Z', 'E', 'N']
for c in compos:
  g = glob.glob('bins/full/*%s.npy'%c)
  for f in g:
    print
    step = 1728000/4
    k=1
    bn = os.path.basename(f)
    
#    aname = 'bins/split/%s'%name
    
#    if os.path.isfile(aname):
#        print '>> removing %s'%f
#        os.system('rm -rf %s'%f)
#        continue
    #else:
    print '>> dealing with  %s'%f
    mat = np.load(f)# dtype=np.float32)
   
    print np.shape(mat)

    for i in np.arange(0,4,1):
            name = 'split'+bn[5:-12]+'.%02d.%s'%(i+1,bn[-11:])
            slice = mat[:,i*step:i*step+step]#   np.shape(mat)
            print '      >',f, '>>', name
            print np.shape(slice)
        #name = 'split_%03d'%k
            k+=1
            np.save('bins/split/%s'%name, slice.astype(np.float32)) 
            del slice
    del mat

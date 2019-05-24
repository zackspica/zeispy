from zeispy.onoribin import *
import os, glob

#path1 = sorted(glob.glob('/data/lawrence2/nnakata/Groningen/data/bin/2016/*Z*'))

comps = ['Z', 'E','N']
years = ['2017']
days = ['%03d'% d for d in np.arange(1,24,1)]

for compo in comps:
  if compo == 'E': compo2 = 'T'
  if compo == 'N': compo2 = 'R'
  if compo == 'Z': compo2 = 'Z'
#  for n in np.arange(2,9,1):
  for d in days :
        #the DH 
        #print '/data/lawrence2/zspica/Groningen/NL-TA/bins/NL/split/split.NL.2016.%s.%s.20sps.%s.npy'%(d, n, compo2)
        bin1 = glob.glob('/data/beroza/zspica/NLTA2/bins/NL/daily.NL.2017.%s.20sps.%s.npy'%(d, compo2))
        if bin1 == []:
            continue
        f1 = bin1[0]
        #print compo, yr, d
        #f1 = raw_open(bin1[0], test=True)
        try:
            f2 = glob.glob('/data/beroza/zspica/NLTA2/bins/TA/daily.TA2.2017.%s.20sps.%s.npy'%(d, compo))[0]
        except Exception as e:
            print e;continue
#        raw_open(f1)
#        raw_open(f2)
        outname = 'daily.TA.NL.2017.%s.20sps.%s.npy'%(d, compo)
        print '>> saving ' , outname
        #mergbins(f1,f2, shape1=(336,1728000), shape2=(10,1728000))
        mergbins(f2,f1, shape1=(415,216000), shape2=(336,216000), outname=outname)




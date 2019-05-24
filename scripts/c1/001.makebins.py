#!/usr/local/miniconda/bin/python

# copy of a020

import sys, glob
import numpy as np
from obspy import Trace, read, read_inventory
# from pylab import *
#import matplotlib.pyplot as plt
import math, os
import nnobspy
from pyrocko import model
from obspy.signal.invsim import corn_freq_2_paz

oridir = '/data/beroza/liuxin/NLphase2/'
bindir  = 'bins/full/'
f_prefilt=(0.02, 0.025, 9.5, 10.0)
resamp=20
nt = int(86400*resamp)
paz = corn_freq_2_paz(5.0, damp=0.67)
paz['sensitivity'] = 89.4  
#respdir = '../responses'
#inv = read_inventory('inv.xml')
#kjd
# ------------------------------------------------------------------------------
# Main
def Process(arg_list):

  # ------------------------------------------------------------
  # parameters
    stations = model.load_stations('stationsBO.txt')
    sday = int(arg_list[0])
    eday = int(arg_list[1])
    comp = arg_list[2]
    iyear = int(arg_list[3])
  # ------------------------------------------------------------
  
    #for iyear in range(2016,2017):
    for iday in range(sday,eday+1):
            bf = bindir+'daily.BO.'+str(iyear)+'.'+'%03d'%iday+'.'+str(resamp)+'sps.'+comp+'.npy'
      #if os.path.isfile(bf):
        #print '>> Exit', bf, 'exists'
        #exit()
      #binzfile = bindir+'/daily.TA.'+str(iyear)+'.'+'%03d'%iday+'.'+str(resamp)+'sps.Z.npy'
      #binnfile = bindir+'/daily.TA.'+str(iyear)+'.'+'%03d'%iday+'.'+str(resamp)+'sps.N.npy'
      #binefile = bindir+'/daily.TA.'+str(iyear)+'.'+'%03d'%iday+'.'+str(resamp)+'sps.E.npy'
      #binfiles = [binzfile, binnfile, binefile]
#      comp = ['HHZ','HHN','HHE']

      #for icomp, bf in enumerate(binfiles):   
        #if icomp <2:continue
            binarray = np.zeros((len(stations), nt))
            count = 0
            for ir, sta in enumerate(stations):
                if comp == 'Z': icomp = 0  
                if comp == 'N': icomp = 1  
                if comp == 'E': icomp = 2  
                #comp = sta.name.split('-')[icomp]
            #respfile = respdir+'/RESP.'+row[0]+'.'+row[1]+'.'+compnametmp[icomp]
            #mseedfile = oridir+'/%s/%03d/%s.%s..%s.D.%s.%03d.mseed'%(iyear,iday, sta.network, sta.station,comp[icomp], iyear,iday)
            #mseedfile = oridir+'/%03d/%s.%s.%s.mseed'%(iday, sta.network, sta.station,iday)
                mseedfile = glob.glob(oridir+'%03d/%s.%s.%03d.mseed'%(iday, sta.network, sta.station,iday))
                if mseedfile !=[]:
                    tr = read(mseedfile[0],format='MSEED')[icomp]
                    #print tr
              #st.resample(resamp)
                    #nnobspy.MergeTraces(st)
                    #for tr in st:
                    nnobspy.Detrend(tr)
                    nnobspy.CorrectStartSample(tr)
                    nnobspy.StartMidnight(tr)
                    tr.resample(resamp)
                   # print comp
                   # if ['E', '1'] in comp:
               # tr.resample(resamp)
                # la fonction resample cause de probleme et 
                # cree des gros spikes!!!
                #tr.decimate(5)#,strict_length=True)
                   # nnobspy.ResampDeconResp_RESPfile(st, inv, f_prefilt, resamp)
                    binarray[ir] = np.float32(tr.data)

                    count += nt

            
                print '>> ', ir, 'year: ',iyear,"; day: ",iday, "; sta:", sta.station, "; comp:", comp   
      
            np.save(bf, binarray)

            del binarray, bf


# ------------------------------------------------------------
if __name__ == "__main__":
  Process(sys.argv[1:])
  

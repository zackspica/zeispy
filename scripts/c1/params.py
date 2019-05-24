#!/usr/bin/env python
from pyrocko import model
from obspy.signal.invsim import corn_freq_2_paz
import os


def get_params():
   # Get Configuration
    params = {}
    ######################
    # params for BINFILE #
    ######################
    params['stations']  = model.load_stations('stationsBO.txt')
#    params['receivers']  = model.load_stations('receivers.txt')
    params['WORKDIR']               = os.getcwd()
    params['downsampling']          = False#True
    params['df']                    = 20. 
    params['dt']                    = 1./params['df']
    params['bin_duration']          = (3600.*24)#/4
    params['percentfill']           = 10
    params['freqmax']               = 9.8
    params['freqmin']               = 0.01
    params['stop']                  = 0 # important param but default = 0
    params['sizeout']               = 415 #size of the output matrix (sizein - stop)
    ###################
    # params for CORR #
    ###################
    params['rotation']              = False
    params['corrType']              = 'crosscoherence'
    params['analysis_duration']     = (3600.*24)#/24 # how much time to process as bunch. must stay like that not implemented to change now
    params['temp_norm']             = 0 # 0: removing eq with param 'clipping'; -1: 1-bit normalization; 1: Windsorinzing with param 'clipping'   
    params['clipping']              = 3. # clipping eq or windsorizing at 3 * std(trace)
    params['overlap']               = 0.4 
    params['corr_duration']         = 250. # slicing the 86400 in small windows
    params['npts']                  = params['corr_duration'] * params['df']
    params['nthreads']              = 12 
    params['outname']               = 'BO.BO' #net1.net2
    ##################
    # params for PAZ #
    ##################
    #"""/!\ hard coded in CorrelBinSingle.py. Sta a or sta b"""
    paz =corn_freq_2_paz(15.0, damp=0.57)
    #paz['sensitivity']=(76*16/16)/5.9605e-7
    params['paz'] = paz
    params['pre_filt'] = (0.05, 0.055, 9.0, 9.5)
    
    return params


        



if __name__=='__main__':
    params= get_params()
    print params
    #print filters

#EDF

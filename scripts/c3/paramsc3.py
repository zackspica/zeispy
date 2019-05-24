#!/usr/bin/env python
from pyrocko import model

def get_params():
    
    params = {}
    params['stations']      = model.load_stations('stationsBO.txt')
    params['sources']      = model.load_stations('stationsBO.txt')
    params['receivers']    = model.load_stations('stationsBO.txt')
    
    params['bindirA']                = '../NL-TA/bins/TA-NL/cc_average/'
    params['bindirB']                = '../NL-TA/bins/TA-NL/cc_average/'
    #params['bindirA']                = '../NLTA2/bins/cc_average/'
    #params['bindirB']                = '../NLTA2/bins/cc_average/'
    params['df']                    = 20.
    params['appvel']                = 1.2
    params['Xcoda']                 = 2. 
    params['coda_len']              = 410
    params['distThreshold']         = 100
    params['shortWin']              = 40
    params['overlap']               = 0.6 
    params['onebit']                = False # even smaller window if needed
    params['stack_method']          = 'linear'
    
    params['freqmax']               = 9.5#float(get_config(db, "preprocess_lowpass"))
    params['freqmin']               = 0.01#float(get_config(db, "preprocess_highpass"))
#    params['format_out']            = 'MSEED'
#    params['nthreads']              = 8
    
    return params


if __name__=='__main__':
    params= get_params()
    print params
    #print filters

#EDF

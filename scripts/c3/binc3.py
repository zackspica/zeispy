# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Created on Weds 04 06 21:03:06 2016
@author: Zack Spica zspica@stanford.edu
Last Update: Mar 18
"""
import os, glob, time
import sys, getopt
import numpy as np
from pyrocko import model, orthodrome
from obspy.core import read, Trace, Stats
from obspy.signal.filter import *
from params import *
from zeispy.obin import *
from zeispy.c3 import *
from zeispy.mycorr import nextpow2
from zeispy.stack import stack
import multiprocessing as mp


def get_bins_pairs(p, staTarget1, staTarget2, comp1, comp2):
    "This function is just to get the proper imputs. take 0.17 sec to compute"   
    for ia, asta in enumerate(p['sources']):
        if asta.station == staTarget1: 
            break
    for ib, bsta in enumerate(p['receivers']):#[ia+1:]):
        if bsta.station == staTarget2: 
            break
    A = asta.station
    B = bsta.station
    namepairA_B = A+'.'+B
    # get dist between pair
    distpairA_B = orthodrome.distance_accurate50m(asta, bsta)/1000.
    print  '>> Dist = %s km'%distpairA_B, 'for pair %s'%namepairA_B
    # get respective binfile
    bin1 = glob.glob('%s/*.%s.%s.*'%(p['bindirA'], A, comp1))
    bin2 = glob.glob('%s/*.%s.%s.*'%(p['bindirB'], B, comp2))
    # if both files exists open
    if bin1 !=[] and bin2 != []:
        bin1 = raw_open(bin1[0], dt=1./p['df'])#, shape=shape,)
        bin2 = raw_open(bin2[0], dt=1./p['df'])#, shape=shape,)
    else: 
        print '>> No files in common...'
        print '>> Exit!'
        sys.exit()

    return bin1, bin2, namepairA_B, asta, bsta, distpairA_B


def corr_c3_pair(p, datapair, distpair, combosta, ic):
    
    print '>> We look for ', combosta[0].station, combosta[1].station,\
            ' with ', combosta[2].station, ' | ', ic 

    l1, r1, l2, r2 = getCodas(datapair, distpair, p['df'], p['appvel'],\
        p['Xcoda'], p['coda_len'], detrend=True, onebit=p['onebit'] )
    #l1 = bandpass(l1, p['freqmin'], p['freqmax'], p['df'], zerophase=True)    
    #l2 = bandpass(l2, p['freqmin'], p['freqmax'], p['df'], zerophase=True)    
    #r1 = bandpass(r1, p['freqmin'], p['freqmax'], p['df'], zerophase=True)    
    #r2 = bandpass(r2, p['freqmin'], p['freqmax'], p['df'], zerophase=True)    
    c3l = corrCodas(l1, l2, int(p['shortWin']*p['df']), \
            int(p['overlap']*p['df']), stack_method=p['stack_method'])
    c3r = corrCodas(r1, r2, int(p['shortWin']*p['df']), \
            int(p['overlap']*p['df']), stack_method=p['stack_method'])

    return c3l, c3r, np.max(distpair)


def gdist(a, b):
    return int(orthodrome.distance_accurate50m(a, b) / 1000.)


def c3_from_bin(argv):
    """
    Compute C3 function between two stations using binfile.
    A binfile is contains all the correlations (C1) between
    a target station and all the other of the network. 
    All parameters are taken from params.py

    3 arguments: StaTarget1, staTarget2, depth
    param staTarget1-2: name of the stations to target (1 and 2)
    rtype: str
    param: depth: depth level of the starget stations - very spacific for Groningen. 
    rtype: float, int. 
    """
    p = get_params()
    try:
        opts, args = getopt.getopt(argv,"ha:b:c:d:v",["help", "staTarget1=", "staTarget2=", "comp1=", "comp2", "verbose"])
    except getopt.GetoptError,  e:
        print e
        print "binc3.py -a <staTarget1> -b <staTarget2> -c <comp1> -d <comp2> -v <verbose> -h <help>"
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print "binc3.py -a <staTarget1> -b <staTarget2> -c <comp1> -d <comp2> -v <verbose> -h <help>"
            sys.exit(2)
        elif opt in ('-a','--staTarget1'):
            staTarget1= arg
        elif opt in ('-b','--staTarget2'):
            staTarget2 = arg
        elif opt in ('-c','--comp1'):
            comp1 = arg
        elif opt in ('-d','--comp2'):
            comp2 = arg
    
    c3 = np.zeros(int(p['shortWin']*p['df']) * 2 - 1 ) # len of output
    datac3 = np.zeros( len(c3) )
    bin1, bin2, namepairA_B, asta, bsta, distpairA_B = get_bins_pairs(p,\
                            staTarget1, staTarget2, comp1, comp2)
    
    pool = mp.Pool(processes=p['nthreads'])
    threader = [pool.apply_async(corr_c3_pair, args=(p, [bin1[ic], bin2[ic]],\
            [gdist(asta, csta), gdist(bsta, csta)], [asta, bsta, csta], ic))\
            for ic, csta in enumerate(p['stations'])]
    for t in threader:
        try:
            res = t.get()
            if res[2] > p['distThreshold']:
                print '>> DistThreshold exceeded: %s'%res[2]
                continue
            if np.all(res[0]==0.) or np.all(res[1]==0.):
                print 'At least one trace is empty'
                continue
            else:
                datac3 = np.vstack((datac3, res[0]))
                datac3 = np.vstack((datac3, res[1]))
        except Exception as e: 
            print e

    print "STACK"
    corrc3 = stack(datac3, stack_method=p['stack_method'])
    
    try:
        os.makedirs('c3/%s/%s.%s/'%(staTarget1, comp1, comp2))
   #     os.makedirs('c1/%s/'%staTarget1)
    except: pass
    t = Trace()
    t.stats.station = 'c3'
    t.stats.sampling_rate = p['df']
    t.data= np.array(corrc3[::-1])
    t.stats.starttime -= ( len(corrc3) / 2 ) / p['df']
    t.write('c3/%s/%s.%s/BO.c3.%s.%s.%s.mseed'%(staTarget1,\
        comp1, comp2, namepairA_B, comp1, comp2), format='MSEED')
    # t2 = Trace()
    # t2.stats.station = 'c1'
    # t2.stats.sampling_rate = df
    # t2.data= np.array(aa)
    # t2.stats.starttime -= ( len(aa) / 2 ) / df
    # t2.write('c1/%s/c1.%s.mseed'%(staTarget1, pair), format='MSEED')

if __name__=='__main__':
    #staTarget1 = '235713'
    #staTarget2 = '236977'
    #depth = 0
    t = time.time()
    c3_from_bin(sys.argv[1:])
    print 'It took :%s s'%(time.time()-t)


#EOF

import os, glob, time
import sys, getopt
import numpy as np
from pyrocko import model, orthodrome,io
from obspy.core import read, Trace, Stats
try:from obspy.signal.invsim import cosine_taper
except:from obspy.signal.invsim import cosTaper as cosine_taper
from obspy.signal.filter import *
from params import *
from zeispy.obin import *
from zeispy.c3 import *
from zeispy.mycorr import nextpow2
from zeispy.stack import stack


def c3_from_bin_fun(s1, s2, d, binfileA, binfileB, time):
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
    stations = p['stations']
    sources = p['sources']
    receivers = p['receivers']
    #bindir = p['bindir'] 
    df = p['df']          
    coda_len = p['coda_len'] #coda lenght (s)
    appvel=p['appvel']       #apparent velocity (km/s)
    Xcoda=p['Xcoda']         # = 2 times the coda teoretical time
    shortWin = int(p['shortWin']*df)  # short win if want a moving window. If shortwin=coda_len -> no moving window
    overlap= int(p['overlap']*df)     # overlap of moving window if used   
    onebit = p['onebit']
    distThreshold = p['distThreshold']
    stack_method = p['stack_method']
    freqmin = p['freqmin']
    freqmax = p['freqmax']
    staTarget1 = s1
    staTarget2 = s2
    depth = d

    c3 = np.zeros( shortWin * 2 - 1 ) # len of output
    datac3 = np.zeros(len(c3))
    for ia, asta in enumerate(sources):
        if asta.station == staTarget1: 
            print asta.station 
        else: 
            continue
        if asta.depth != depth: 
            #Important condition about the depth of the receivers at Groningen!
            print 'Passing depth %s'%depth 
            continue
        for ib, bsta in enumerate(receivers):#[ia+1:]):
            if bsta.station == staTarget2: 
                print bsta.station
                nb = ib
                na = ia
            else: 
                continue
            if bsta.depth != depth: 
                print 'Passing depth %s'%depth
                continue#

            A = asta.station
            B = bsta.station
            pair = A+'.'+B
            # get dist between pair
            dd = orthodrome.distance_accurate50m(asta,bsta)/1000.
            print  '>> Dist = %s km'%dd, 'for pair %s'%pair
            # get respective binfile
            bin1 = binfileA#glob.glob('%s/*%s*'%(bindir, A))#+1 becaus ia starts at 0 but stalist at 1
            bin2 = binfileB#glob.glob('%s/*%s*'%(bindir, B))
            # if both files exists open
            if bin1 !=[] and bin1!=[]:
                bin1 = raw_open(bin1[0], dt=1./df)#, shape=shape,)
                bin2 = raw_open(bin2[0], dt=1./df)#, shape=shape,)
            else: 
                print '>> No files in common...'
                continue

            # Now we must only correlate C1 at a same depth becauz binfile contains all the C1 at all depths.
            for ic, csta in enumerate(stations):
                if csta.station in [staTarget1, staTarget2]: # meaning auto-corr
                    continue
                if csta.depth != depth: 
                    print '>> Not a good depth: %s'%csta.depth
                    continue 
                # for each sta in list we get the respective data
                data1 = bin1[ic]#[::-1]
                data2 = bin2[ic]#[::-1]
                if ic==0:
                    aa = bin1[nb] / np.max(bin1[nb])
                print '>> We look for ', pair, ' with ', csta.station, '-%s- '%time, '| ',ic 
                # some can be empty
                if np.all(data1==0.) or np.all(data2==0.):
                    print 'At least one trace is empty'
                    continue
                # filtering
                data1 = bandpass(data1, freqmin, freqmax, df, zerophase=True)    
                data2 = bandpass(data2, freqmin, freqmax, df, zerophase=True)    
                datapair = [data1, data2]
                # get dists betweem S-A and S-B
                d1 = int(orthodrome.distance_accurate50m(asta, csta) / 1000.)
                d2 = int(orthodrome.distance_accurate50m(bsta, csta) / 1000.)
                if d1 > distThreshold or d2 > distThreshold: 
                    '>> DistThreshold exceeded: %s'%(np.max([d1, d2]))
                    continue
                distpair = [d1, d2]

                l1, r1, l2, r2 = getCodas(datapair, distpair, df, appvel, Xcoda,\
                        coda_len, detrend=True, onebit=onebit )

                l1 -= np.mean(l1)
                r1 -= np.mean(r1)
                l2 -= np.mean(l2)
                r2 -= np.mean(r2)

                l1 = bandpass(l1, freqmin, freqmax, df, zerophase=True)    
                l2 = bandpass(l2, freqmin, freqmax, df, zerophase=True)    
                r1 = bandpass(r1, freqmin, freqmax, df, zerophase=True)    
                r2 = bandpass(r2, freqmin, freqmax, df, zerophase=True)    
                c3l = corrCodas(l1, l2, shortWin, overlap, stack_method=stack_method)
                c3r = corrCodas(r1, r2, shortWin , overlap, stack_method=stack_method)

                datac3 = np.vstack((datac3, c3l))
                datac3 = np.vstack((datac3, c3r))

            print "STACK"
            corrc3 = stack(datac3, stack_method=stack_method)

            c3 = np.array(corrc3[::-1])

    return c3




def prestackc3(argv):
   
    p = get_params()
    df = p['df']  # short win if want a moving window. If shortwin=coda_len -> no moving window
    shortWin = int(p['shortWin']*df)  # short win if want a moving window. If shortwin=coda_len -> no moving window
    dir_prestack = '../c1/bins/corr/rotation/'
    sorted = True#False

    try:
        opts, args = getopt.getopt(argv,"ha:b:d:v",["help", "staTarget1=", "staTarget2=", "depth=", "verbose"])
    except getopt.GetoptError,  e:
        print e
        print "binc3.py -a <staTarget1> -b <staTarget2> -d <depth> -v <verbose> -h <help>"
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print "binc2.py -a <staTarget1> -b <staTarget2> -d <depth> -v <verbose> -h <help>"
            sys.exit(2)
        elif opt in ('-a','--staTarget1'):
            s1= arg
        elif opt in ('-b','--staTarget2'):
            s2 = arg
        elif opt in ('-d','--depth'):
            d = int(arg)
    
    c3 = np.zeros( shortWin * 2 - 1 ) # len of output
    datac3 = np.zeros(len(c3))
    d=0
    
    ga = glob.glob(dir_prestack+s1+'/*ZZ*')
    gb = glob.glob(dir_prestack+s2+'/*ZZ*')
    
    am = np.argmax([len(ga), len(gb)])
    if am == 0: 
        gmajor = ga
        smajor = s1
        gminor = gb
        sminor = s2
    else: 
        gmajor = gb
        smajor = s2
        gminor = ga
        sminor = s1

    if sorted:
        for ia, a in enumerate(gmajor):
#            if ia >50:continue
            time = os.path.basename(a)[6:14]
            time_name =  '%s%s/*%s*ZZ*' %( dir_prestack, sminor, time)
            gtmp = glob.glob(time_name)#dir_prestack+smajor+'TA*')
            binfileA = [a]
            binfileB = gtmp
            if gtmp != []:
                c3 = c3_from_bin_fun(s1, s2, d, binfileA, binfileB, time)
                datac3 = np.vstack((datac3, c3))
            else: continue
            print np.shape(datac3) 
   
    else:
        for ia, a in enumerate(gminor):
 #           if ia >50:continue
            binfileA = [a]
            binfileB = [gmajor[ia]]
            c3 = c3_from_bin_fun(s1, s2, d, binfileA, binfileB)
            datac3 = np.vstack((datac3, c3))
            print np.shape(datac3) 
            
        
    corrc3 = stack(datac3, stack_method='linear')

    t = Trace()
    t.stats.station = 'pc3'
    t.stats.sampling_rate = df
    t.data= np.array(corrc3[::-1])
    t.stats.starttime -= ( len(corrc3) / 2 ) / df
    if sorted:
        t.write('pc3/prestackc3.%s.%s.mseed'%(s1,s2), format='MSEED')
    else:
        t.write('pc3/random_prestackc3.%s.%s.mseed'%(s1,s2), format='MSEED')






if __name__=='__main__':
    #staTarget1 = '235713'
    #staTarget2 = '236977'
    #depth = 0.
    prestackc3(sys.argv[1:])


#EOF

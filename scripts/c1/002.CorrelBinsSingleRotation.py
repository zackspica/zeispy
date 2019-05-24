import sys, getopt
import numpy as np
from zeispy.mycorr import *
from params import get_params
from obspy.signal.invsim import cosine_taper, simulate_seismometer 
import os, glob, time
from scipy.signal import detrend
import multiprocessing as mp
from obspy.geodetics import gps2dist_azimuth


def azimuth(lat1, lon1, lat2, lon2):
    dist, azim, bazim = gps2dist_azimuth(lat1, lon1, lat2, lon2)
    return azim


def rm_nans(bin):
    for il, l in enumerate(bin):
        try:
            if np.mean(l) > 10**6: 
                b[il] = np.zeros(len(l))
            if np.isinf(l[0]):
                b[il] = np.zeros(len(l))
            if np.isnan(l[0]):
                b[il] = np.zeros(len(l))
        except:continue


def get_slices(params):
    """Create slice for iterations"""
    slices = int(params['corr_duration'] * params['df'] / params['npts'])
    begins = []
    ends = []
    i = 0
    while i <=  (params['bin_duration'] - params['npts']/params['df']):
        begins.append(int(i * params['df']))
        ends.append(int(i * params['df'] + params['npts']))
        i += int(params['npts']/params['df'] * (1.0-params['overlap']))
    slices = len(begins)

    return slices, begins, ends


def mkdir(WORKDIR, sta1):
    """Create a dir"""
    try:
        os.makedirs('%s/bins/corr/%s'%(WORKDIR,sta1))
    except :
        pass


def test_empty_sta(a, verbose, sta1, outname, comp):
    if np.all(a==0):
        if verbose: 
            print 'Station %s is empty for %s'%(sta1, outname)
            print 'Going out...'
        os.system('echo "Empty station... %s %s" >> o.txt'%(sta1, outname))
        sys.exit(2)
    else:
        if verbose: 
            print ' >> Starting with sta1 = %s for %s'%(sta1, comp)        


def simulate(data, params, remove_sensitivity=True ):
    """Remove instrumental response"""
    data -= np.mean(data)
    data = simulate_seismometer(data, params['df'], paz_remove=params['paz'], \
           pre_filt=params['pre_filt'], remove_sensitivity=remove_sensitivity)
    return data


def looping(ib, b, a, cp, params, begins, ends, corrType, stop, verbose):
    """Main function for cross correlation of binaryfiles"""
    sta2 = params['stations'][ib].station
    s2 = params['stations'][ib]
    if verbose: print  'sta2 =', sta2, ib
    i = 1
    cc = np.zeros( ((params['corr_duration']*params['df'])*2 -1 ))
    if np.all(b==0.): 
        if verbose:
            print 'sta2 = %s is empty... passing'%sta2
        os.system('echo ">> sta2 = %s is empty... passing" >> o.txt'%sta2)
    if np.isnan(a).any(): print 'NANs';  a = np.zeros(len(a))
    if np.isnan(b).any(): print 'NANs'; b = np.zeros(len(b))

    else:
        ###Respons###Hard coded.In case two different sensor response can be correcte here too. 
        #a = simulate(a, params)
        if params['temp_norm'] == 1 or params['temp_norm'] == 0:
            rmsa = np.std(np.abs(a))
            rmsb = np.std(np.abs(b))
        if params['temp_norm'] == -1:
            a = np.sign(a); b = np.sign(b)

        for islice, (begin, end) in enumerate(zip(begins,ends)):
            vb, va = b[begin:end], a[begin:end]
            if params['temp_norm'] == 1:
                indexes = np.where(np.abs(vb) > (params['clipping'] * rmsb))[0]
                vb[indexes] = (vb[indexes] / np.abs(vb[indexes])) * params['clipping'] * rmsb
                indexes = np.where(np.abs(va) > (params['clipping'] * rmsa))[0]
                va[indexes] = (va[indexes] / np.abs(va[indexes])) * params['clipping'] * rmsa
            if params['temp_norm'] == 0:
                if np.std(np.abs(vb)) >= rmsb * params['clipping']:
                    #print '>> Eq detected: removing window'
                    continue
                if np.std(np.abs(va)) >= rmsa * params['clipping']:
                    #print '>> Eq detected: removing window'
                    continue

            vb, va = vb[:params['npts']], va[:params['npts']]
            if np.all(va==0.) or np.all(vb==0.): continue
            else:
                """hard coded!!!!"""
                va = detrend(va, type='linear'); vb = detrend(vb, type='linear')
                va *= cp; vb *= cp
                if corrType == 'crosscoherence':
                    corr = crosscoherence(vb, va)#, npt=15.,  eps=1e-7)
                if corrType == 'convolution':
                    corr = convolution(vb, va)#, npt=15.,  eps=1e-7)
                cc += corr
                i+=1.
                del corr, va, vb 
     
    return cc/i, ib-stop



def corrbin(argv):
    """
    Function to correlate one source with all the receivers. 
    Usage: python CorrelBinsSingleRotation.py --help to display help message"
    """

    t0 = time.time()
    binfile = ''
    verbose = False
    params = get_params()
    nthreads = params['nthreads']
    corrType = params['corrType']

    try:
        opts, args = getopt.getopt(argv,"hp:d:s:c:t:v",["help", "path2bin=", "date=", "station=", "components", "threads=", "verbose"])
    except getopt.GetoptError,  e:
        print e
        print "correlbins.py -p <path2bin> -d <date> -s <sourceName> -c <components> -t <nthreads> -v <verbose> -h <help>"
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print "correlbins.py -p <path2bin> -d <nate> -s <sourceName> -c <components> -t <nthreads> -v <verbose> -h <help>"
            sys.exit(2)
        elif opt in ('-p','--path'):
            path = arg
        elif opt in ('-d','--date'):
            date = arg
        elif opt in ('-s','--source'):
            staToCorr = arg
        elif opt in ('-c','--component'):
            comps = arg
        elif opt in ('-t','--threads'):
            nthreads = int(arg)
        elif opt in ('-v', '--verbose'):
            verbose = True
    
    if path.endswith('/'): pass
    else: path = path+'/'
    #print path, date
    files = sorted(glob.glob(path+'*%s*'%(date))) #ENZ
    print files
    stop = params['stop'] 
    sizeout = params['sizeout'] #size of the out matrix
    slices, begins, ends = get_slices(params)
    cp = cosine_taper((ends[0]-begins[0]), 0.05)
    binout = np.zeros( ( sizeout, ((params['corr_duration']*params['df'])*2-1) ))
    pool = mp.Pool(processes=nthreads)

    for ia, sta in enumerate(params['stations']):
        if sta.station==staToCorr:
            s1 = params['stations'][ia]
            sta1 = s1.station
            break
            
    outname = '%s.%s.%s.%s.%s'%(params['outname'],date, sta1, comps, corrType)
    
    print ' >> Starting %s'% outname

    if comps == 'ZZ':
        bin = np.load(files[2])
        a = bin[ia]
        test_empty_sta(a, verbose, sta1, outname, comps)
        
        results = [pool.apply_async(looping, args=(ib+stop, b, a, cp, params, begins, ends, corrType, stop, verbose))\
                for ib, b in enumerate(bin[stop:])]
        outputs = []
        for p in results:
            try:
                tmp = p.get()
                outputs.append( tmp )
            except Exception as e: 
                outputs.append((tmp[0]*0, tmp[1]+1))
        #outputs = [p.get() for p in results]
    
    else:
        trace_len = int(params['bin_duration']*params['df'])
        if verbose: print "Creating frame_a of size :%s for %s"%(trace_len, comps[0])
        allAz=[]
        for is2,  s2 in enumerate(params['stations'][stop:]):
            allAz.append(np.deg2rad(azimuth(s1.lat, s1.lon, s2.lat, s2.lon)))
        
        if comps[0]=='Z':
            a = np.load(files[2])[ia]
            test_empty_sta(a, verbose, sta1, outname, comps[0])
            frame_a = np.vstack([a]*binout.shape[0])
        else:
            frame_a = np.zeros((binout.shape[0], trace_len)) 
            n = np.load(files[1])[ia]
            e = np.load(files[0])[ia]
            test_empty_sta(n, verbose, sta1, outname, comps[0]); test_empty_sta(e, verbose, sta1, outname, comps[0])
            for iazi, azi in enumerate(allAz):
                if comps[0]=='R':
                    if azi != 0 :
                        frame_a[iazi] = n * np.cos(azi) + e * np.sin(azi)
                    else: 
                        frame_a[iazi] = e 
                elif comps[0]=='T':
                    if azi != 0 :
                        frame_a[iazi] = n * np.sin(azi) - e * np.cos(azi)
                    else: 
                        frame_a[iazi] = n 
        if verbose: print "Creating frame_a of size :%s for %s"%(trace_len, comps[1])
        if comps[1]=='Z':
            frame_b = np.load(files[2])[stop:]
        else:
            frame_b = np.zeros((binout.shape[0], trace_len))
            N = np.load(files[1])[stop:]
            E = np.load(files[0])[stop:]
            for iazi, (azi, n, e) in enumerate(zip(allAz, N, E)):
                if comps[1]=='R':
                    if azi != 0 :
                        frame_b[iazi] = n * np.cos(azi) + e * np.sin(azi)
                    else: 
                        frame_b[iazi] = e 
                elif comps[1]=='T':
                    if azi != 0 :
                        frame_b[iazi] = n * np.sin(azi) - e * np.cos(azi)
                    else: 
                        frame_b[iazi] = n 
        #rm_nans(frame_a)
        #rm_nans(frame_b)
        #pool = mp.Pool(processes=nthreads)
        results = [pool.apply_async(looping, args=(ib+stop, b, a, cp, params, begins, ends, corrType, stop, verbose))\
                for ib, (b, a) in enumerate(zip(frame_b, frame_a))]
        outputs = []
        for p in results:
            try:
                tmp = p.get()
                outputs.append( tmp )
            except Exception as e: 
                outputs.append((tmp[0]*0, tmp[1]+1))
        #outputs = [p.get() for p in results]
        
    for out in outputs:
        if len(out[0]) == np.shape(binout)[1]:
            binout[out[1]] = out[0]
    if np.all(binout==0):
        if verbose: 
            print ' >> Empty bin... %s %s'%(sta1, outname)
            print 'Going out...'
        os.system('echo "Empty bin... %s %s" >> o.txt'%(sta1, outname))
        sys.exit(2)
    else:
        mkdir(params['WORKDIR'], sta1)
        filename = '%s/bins/corr/%s/%s' % (params['WORKDIR'],sta1, outname)
        if verbose: 
            os.system('echo ">>> Saving bin for %s %s %s" >> o.txt'%(sta1, outname, comps))
            os.system('echo ">>> It took: %s min for %s with all others" >> o.txt'%((time.time() - t0)/60., sta1))
        binout.astype(np.float32)
        np.save(filename, binout)



if __name__=='__main__':
    t0 = time.time() 
    corrbin(sys.argv[1:])
    print ' >>> It took: ', (time.time() - t0), 'sec'

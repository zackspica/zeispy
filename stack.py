#!/usr/bin/env python
import numpy as np
import scipy as sp
import os, time, glob
from params import get_params
#np.set_printoptions(threshold=np.nan)
import matplotlib.pyplot as plt
#from multiprocessing import Pool



def simple_stack_bins(folder, compo='ZZ', stack_method='linear', df=20. ):
    """ Stacking function to stack bins together
    The output file will be saved in bins/cc_average/
    and sorted by stations.

    :param folder: The folder where are the cc bins (dayly?)
    :type folder: str
    :param compo: The component to compute
    :type compo: str
    :param stack_method: 'linear' or 'pws'
    :type stack_method: str
    :param df: sampling rate
    :type df: int, float
    """

    params = get_params()
    if folder.endswith('/'): pass
    else: folder+='/'
    sourcesta = os.path.dirname(folder).split('/')[-1]
    fileout = params['corrType']+'.'+sourcesta+'.'+compo
    folderout = '%s/bins/cc_average'%params['WORKDIR']
    out = os.path.join(folderout,fileout)
    print (folder)
    try:os.makedirs(folderout)
    except: pass

    bins = sorted(glob.glob(folder+'*%s*'%compo))
    print (len(bins),' CC for component %s for station %s'%(compo, sourcesta))

    for ibin, bin in enumerate(bins):
        day = os.path.basename(bin)
        try:matrix = np.load(bin)
        except:continue
        if ibin == 0.:
            divider = np.ones(np.shape(matrix)[0])
            if glob.glob(out+'_*.npy') != []:
                print ('>> One stack file already exist: continue.')
                stacker = np.load(bin, allow_pickle=True)
            else:
                print ('>> New stacking...')
                stacker = matrix
            continue
        else:
            for icc, cc in enumerate(matrix):
                if np.all(cc==0.):
                    continue
                else:
                    try:
                        stacker[icc] = stack(np.vstack((stacker[icc], cc)), \
                            stack_method=stack_method, df=df)
                       # divider[icc]+=1.
                    except Exception as e:
                        print (e)
        print (sourcesta, day)
    #for i in np.arange(len(stacker)):
    #    stacker[i] /= divider[i]
    
    np.save(out+'.%s'%ibin, stacker)
    
    del matrix, stacker


def stack(data, stack_method='linear', pws_timegate=10, pws_power=2, df=20.):
    """Original stack function from MSNoise

    :param data: matrix in which every line is a cc to be stacked.
    :type data: ndarray
    :param stack_method: If ``stack_method`` is 'linear', then a simple mean CFF 
            is returned. If ``stack_method`` is 'pws', then the Phase Weighted 
            Stack (PWS) is computed. The PWS is done in two steps: first the mean 
            coherence between the instataneous phases of all windows is calculated, 
            and eventually serves a weighting factor on the mean.                 
    :type stack_method: 'linear' or 'pws'
    :param pws_timegate: The smoothness of this weighting array.     
    :type pws_timegate: int
    :param pws_power: The weighting array is the power of the mean coherence array. 
            If ``pws_power`` is equal to 0, a linear stack is done (it's faster to
            set ``method``='linear'). Usual value is 2.
    :type pws_power: int
    :param df: sampling rate
    :type df: int, float
    """

    data = data[~np.isnan(data).any(axis=1)]
    
    for i in range(data.shape[0]):
        data[i] -= data[i].mean()

    if stack_method == "linear":
        corr = data.mean(axis=0)

    elif stack_method == "pws":
        corr = np.zeros(data.shape[1], dtype='f8')
        phasestack = np.zeros(data.shape[1], dtype='c8')
        for i in range(data.shape[0]):
            data[i] -= data[i].mean()
        for c in data:
            phase = np.angle(sp.signal.hilbert(c))
            phasestack.real += np.cos(phase)
            phasestack.imag += np.sin(phase)
        coh = 1. / data.shape[0] * np.abs(phasestack)

        timegate_samples = int(pws_timegate * df)
        coh = np.convolve(sp.signal.boxcar(timegate_samples) /
                timegate_samples, coh, 'same')
        for c in data:
            corr += c * np.power(coh, pws_power)
        corr /= data.shape[0]

    return corr





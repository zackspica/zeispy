import commands
from popen2 import popen2
import time, glob, os, sys
from pyrocko import model
import subprocess
import numpy as np

def run(pairs):
    k=0
    for p in pairs:
      for d in np.r_[356:367,1:24]:
        k+=1
        if k < 1678302: continue #written on 30/01/2019
        
        staTarget1 = p[0]
        staTarget2 = p[1]
        
#        target = sta.station
        print k , staTarget1, staTarget2, k
        #continue
        # Open a pipe to the qsub command.
        output, input = popen2('qsub')
        # Customize your options here
        job_name = p[0]+'.'+p[1]+'.%03d.pc3'%d
        nnodes = 1
        nthreads = 1
#        path = '/data/lawrence2/zspica/Groningen/TA/bins/split/'
        processors = "nodes=1:ppn=%s"%nnodes
        command = "python prestack_perday_binc3.py -a %s -b %s -d %s" %(staTarget1, staTarget2, d )
        node = 'default'
        #node = 'beroza'
        #if k % 100.==0.:node= 'Beroza'
        #if k % 20.==0.:node= 'beroza'
    #    if ibin % 100.==0.:node= 'Q2a'

        job_string = """
        #!/bin/bash\n\
        #PBS -N %s\n\
        #PBS -q %s\n\
        #PBS -l %s\n\
        #PBS -o out/%s.out\n\
        #PBS -e out/%s.err\n\
        cd $PBS_O_WORKDIR\n\
        %s""" % (job_name, node, processors, job_name, job_name, command)
        
        # Send job_string to qsub
        input.write(job_string)
        input.close()
        
        # Print your job and the system response to the screen as it's submitted
        print job_string
        print output.read()
        time.sleep(0.8000)
        
       # os.system('echo "%s %s %s" >> o.txt'%(bin, target, k ))
        try:
            njob = int(commands.getoutput('qstat -u zspica | wc -l'))
        except:
            pass
        while njob >= 199:
            try:
                njob = int(commands.getoutput('qstat -u zspica | wc -l'))
            except: 
                time.sleep(2)
                continue
            time.sleep(10)


if __name__=='__main__':
    try: os.makedirs('out/')
    except: pass
    stations = model.load_stations('stationsBO.txt')
    # Loop over your jobs
    pairs = []
    staTarget = [sta.station for sta in stations]#
    for staTarget1 in staTarget:
        list = np.ndarray.tolist(np.load('liststa.npy'))
        print list
        if staTarget1 in list:
            continue
        for staTarget2 in stations:
            if staTarget2.station == staTarget1: 
                continue
            pairs.append([staTarget1, staTarget2.station])
    print pairs, len(pairs)
    run(pairs)
    

import commands
from popen2 import popen2
import time, glob, os, sys
from pyrocko import model
import subprocess

def run(pairs):
    k=0
    for p in pairs:
        
        k+=1
        #if k < 376300: continue
        
        staTarget1 = p[0]
        staTarget2 = p[1]
#        target = sta.station
        print staTarget1, staTarget2, k
        #continue
        # Open a pipe to the qsub command.
        output, input = popen2('qsub')
        # Customize your options here
        job_name = p[0]+'.'+p[1]+'.pc3'
        nnodes = 1
        nthreads = 1
#        path = '/data/lawrence2/zspica/Groningen/TA/bins/split/'
        processors = "nodes=1:ppn=%s"%nnodes
        command = "/home/zspica/local/bin/python2.7 prestack_binc3.py -a %s -b %s -d 0" %(staTarget1, staTarget2)
        node = 'default'
        node = 'beroza'
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
        time.sleep(0.1000)
        
       # os.system('echo "%s %s %s" >> o.txt'%(bin, target, k ))
        njob = int(commands.getoutput('qstat -u zspica | wc -l'))
        while njob >= 206:
            njob = int(commands.getoutput('qstat -u zspica | wc -l'))
            time.sleep(30)


if __name__=='__main__':
    try: os.makedir('out/')
    except: pass
    stations = model.load_stations('stations.txt')
    # Loop over your jobs
    pairs = []
    staTarget1 = '236367'
    for staTarget2 in stations:
        if staTarget2.station == staTarget1: 
            continue
        pairs.append([staTarget1, staTarget2.station])
    print pairs
    run(pairs)
    

import commands
from popen2 import popen2
import time, glob, os, sys
from pyrocko import model
import subprocess



#print binfiles
#time.sleep(99)
def run():#binfiles, stations):
  k=0
  try: os.mkdir('out/')
  except: pass
  #for ibin, bin in enumerate(binfiles):
  #  for ista, sta in enumerate(stations[2665:2680]):
  dates = []
  for day in range(356, 367,1):
    dates.append(['%03d'%day, '2016'])
  for day in range(1,24,1):
    dates.append(['%03d'%day, '2017'])

          #  for day in range(306,330,1):
  print dates
  for comp in ['Z','E','N']:#,'N']:
    for day in dates:
        print day[0],day[1], comp
        #continue
        # Open a pipe to the qsub command.
        output, input = popen2('qsub')
        # Customize your options here
        job_name = 'BO.bin.%s.%s'%(day[0],comp)
        nnodes = 8
        processors = "nodes=1:ppn=%s"%nnodes
        command = "python2.7 001.makebins.py %s %s %s %s" %(day[0], day[0], comp, day[1])
        node = 'default'
        #node = 'beroza'
        #if k % 100.==0.:node= 'Q2a'
    #    if ista % 100.==0.:node= 'Q2a'
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
        time.sleep(0.1)
        
        #os.system('echo "%s %s %s" >> o.txt'%(bin, target, k ))
        #njob = int(commands.getoutput('qstat -u zspica | wc -l'))
        #while njob >= 15:
        #    njob = int(commands.getoutput('qstat -u zspica | wc -l'))
        #    time.sleep(60)


if __name__=='__main__':
    #try: os.makedir('out/')
    #except: pass
    #stations = model.load_stations('stations.txt')
    # Loop over your jobs
    #binfiles = sorted(glob.glob('bins/full/*'))#start at day all_253
    run()#binfiles,stations)
    #binfiles = dubble_check(binfiles)
    #run(binfiles,stations)
    #binfiles = error_check()
    #run(binfiles,stations)
    

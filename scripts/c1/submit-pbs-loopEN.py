import commands
from popen2 import popen2
import time, glob, os, sys
from pyrocko import model
import subprocess

def dubble_check(binfiles):
    tmp = glob.glob('corr/*/*')
    tt = []
    for t in tmp:
        n = os.path.basename(t).split('.')[-2]
        tt.append(n)
    bins= []
    for b in binfiles:
        n = os.path.basename(b).split('.')[-2]
        if n in tt:
            print 'exist'
        else : 
            print 'add'
            bins.append(b)
    print bins, len(bins)
    binfiles=bins
    return binfiles

#def error_check():
#    g = glob.glob('out/*err')
#    binfiles=[]
#    for f in g:
#        if os.path.getsize(f) != 0:
#            #check if wall time
#            lines = open(f, 'r+').readlines()
#            if 'walltime' in lines[0]:
#                num = f.split('_')[-1].split('.')[0]
#                bname = 'bins/all_%s.bin'%num
#                binfiles.append(bname)
#    return binfiles



#time.sleep(99)
#binfiles = sorted(glob.glob('bins/full/*238*'))#start at day all_253

#print binfiles
#time.sleep(99)
def run(stations, components, dates):
  k=0
  for comps in components:
    for date in dates:
      for ista, sta in enumerate(stations):
#    for ibin, bin in enumerate(binfiles):
#        if sta.station not in ['004','007']: continue
        k+=1
        #if k < 91643: continue 
        target = sta.station
        print date, comps,  sta.station, k
        #continue
        # Open a pipe to the qsub command.
        output, input = popen2('qsub')
        # Customize your options here
        job_name = date+'.'+target+'.'+comps
        nnodes = 6
        nthreads = 6 
        path = '/data/beroza/zspica/NL-TA/bins/split/'
        processors = "nodes=1:ppn=%s"%nnodes
        command = "/home/zspica/local/bin/python2.7 002.CorrelBinsSingleRotation.py -p %s -d %s -s %s -c %s -t %s -v" %(path, date, target, comps,  nthreads)
        node = 'default'
       # if k % 500.==0.:node= 'beroza'
        #if k % 75.==0.:node= 'jfl'
    #    if ibin % 100.==0.:node= 'Q2a'

        job_string = """
        #!/bin/bash\n\
        #PBS -N %s\n\
        #PBS -q %s\n\
        #PBS -l %s\n\
        #PBS -o out2/%s.out\n\
        #PBS -e out2/%s.err\n\
        md $PBS_O_WORKDIR\n\
        %s""" % (job_name, node, processors, job_name, job_name, command)
        
        # Send job_string to qsub
        input.write(job_string)
        input.close()
        
        # Print your job and the system response to the screen as it's submitted
        print job_string
        print output.read()
        time.sleep(0.1000)
        
        os.system('echo "%s %s %s" >> o.txt'%(bin, target, k ))
        try:
            njob = int(commands.getoutput('qstat -u zspica | wc -l'))
        except:
            pass
        while njob >= 306:
            try:
                njob = int(commands.getoutput('qstat -u zspica | wc -l'))
            except: 
                time.sleep(2)
                continue
            time.sleep(20)


if __name__=='__main__':
    try: os.makedir('out/')
    except: pass
    stations = model.load_stations('stationsNL.txt')
    # Loop over your jobs
    components = ['TT','RR']#'ZR','ZT','RZ','RT','TZ','TR']
    dates = []
    for n in xrange(284,330,1):
      for ix in xrange(1,25,1):
        dates.append('2016.%03d.%02d'%(n,ix))
    print dates
    run(stations, components, dates)
#    binfiles = dubble_check(binfiles)
#    print binfiles
    #run(binfiles,stations)
   # binfiles = error_check()
   # run(binfiles,stations)
    

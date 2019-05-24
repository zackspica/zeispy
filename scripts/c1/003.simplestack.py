import os
from zeispy.stack import simple_stack_bins
import glob

if __name__=='__main__':
    path = 'bins/corr/no_rotation/'
    g = glob.glob('%s*'%path)
    stations = []
    for s in g:
        stations.append(os.path.basename(s))
#    stations = ['235555',  '235564' , '235573' , '235589',  '235708' , '235714' , '235758',  '235844'  ,'235855'  ,'235879',  '235902' , '235910','235916','236065'  ,'236341',  '236393' ,'236410' ,'236978','237010','237029','235559',  '235568','235583','235596','235713','235736','235838','235846',  '235876',  '235899',  '235904', '235914' , '235918',  '236111',  '236386' , '236399'  ,'236412'  ,'237005',  '237015',  '237061']
    #compos = ['ZZ','RR','TT']#,'ZR','ZT','RZ','RT','TZ','TR']
    print stations
    #compos = ['ZZ']#,'RZ','RT','TZ','TR']
    #compos = ['ZZ', 'TT', 'RR','ZR','ZT','RZ','RT','TZ','TR']
    compos = ['EE', 'NN','ZR','ZN','EZ','EN','NZ','NE']
    #compos = ['TZ','TR']
    for comp in compos:
        for sta in stations:
            folder = '%s%s'%(path,sta)
            try:
                simple_stack_bins(folder, compo=comp)
            except:
                pass

import obspy, os, math, sys, datetime

# ------------------------------------------------------------
def MergeTraces(st):
  count=0
  if len(st) >= 2:
    stime0 = st[0].stats.starttime
    etime0 = st[0].stats.endtime
    for itr in range(len(st)-1):
      stime1 = st[1+count].stats.starttime
      etime1 = st[1+count].stats.endtime
      if (stime0 > etime1) or (stime1 > etime0):
        st[0] = st[0].__add__(st[1+count],fill_value=0)
        st.remove(st[1+count])
      else:
        count+=1

# ------------------------------------------------------------
def Detrend(tr):
  for inum in range(0,len(tr.data)):
    if math.isnan(tr.data[inum]):
      tr.data[inum]=0.0

  tr.detrend(type='constant')
  tr.detrend(type='linear')


# ------------------------------------------------------------
def CorrectStartSample(tr):
  dt=obspy.core.UTCDateTime(tr.stats.starttime)
  sr=tr.stats.sampling_rate
  npts = tr.stats.npts
  msec = (int)(math.ceil(dt.microsecond * sr / 1000000)*1000000/sr)
  if msec >= 1000000:
    msec=msec-1000000
    dt = dt+1
    npts = npts-1

  dt.microsecond=msec

  try:
    tr.interpolate(sampling_rate=tr.stats.sampling_rate,starttime=dt,npts=npts)
  except:
    try:
      tr.interpolate(sampling_rate=tr.stats.sampling_rate,starttime=dt,npts=npts-1)
    except:
      pass

# ------------------------------------------------------------
def StartMidnight(tr):
  dt = obspy.core.UTCDateTime(tr.stats.starttime)
  dt.hour        = 0
  dt.minute      = 0
  dt.second      = 0
  dt.microsecond = 0

  tr.stats.sampling_rate = round(tr.stats.sampling_rate*1000)/1000
  tr.trim(dt,dt+86400-1/tr.stats.sampling_rate,pad=True,fill_value=0)

# ------------------------------------------------------------
#def ResampDeconResp_RESPfile(st, inv, f_prefilt, resamp):
def ResampDeconResp_RESPfile(tr, paz, f_prefilt, resamp):
  #for tr in st:
  #  tr.resample(resamp)
  #  tr.stats.location=''
  #  if tr.stats.channel == 'HG1':
  #      tr.stats.channel = 'HG2'
  #  elif tr.stats.channel == 'HH1':
  #      tr.stats.channel = 'HH2'
   # print tr.id
  #seedresp={'filename': respfile, 'units': "VEL" }
    tr.simulate(paz_remove=paz, pre_filt=f_prefilt, taper=True, taper_fraction=0.05, water_level=60.0)
    #tr.simulate(paz_remove=paz, pre_filt=f_prefilt, seedresp=seedresp, taper=True, taper_fraction=0.05, water_level=600.0)
  
  # the routine automatically picks the correct response for each trace

  #st.remove_response(inventory=inv, pre_filt=f_prefilt, output='VEL',taper=True, taper_fraction=0.05, water_level=60.0)  

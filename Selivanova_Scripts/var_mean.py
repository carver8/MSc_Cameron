# mean steric level; can be used for any variable 


import numpy as np
import netCDF4 as nc
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
#import cartopy
import matplotlib.colors as clrs
import xarray as xr
import matplotlib.path as mpath
import numpy.ma as ma
import os.path
import glob
from scipy import stats
from scipy import signal
from scipy.stats import linregress
from netCDF4 import Dataset

start=1993
end=2010
timelen=end-start+1
mons={0:'jan',1:'feb',2:'mar',3:'apr',4:'may',5:'jun',6:'jul',7:'aug',8:'sep',9:'oct',10:'nov',11:'dec'}
exp='CG8E7'

lev='700'


prod=['cice','cice_init','hadisst']
sth_y=np.zeros((timelen,1050,1440))

stt_y=np.zeros((timelen,1050,1440))
darea='/work/cmcc/aspect/CESM2/inputdata/STATIC_DA_NEMO42/mesh_mask.nc'
area=xr.open_dataset(darea,engine='netcdf4')
dx=area.e1t
dy=area.e2t
s=dx*dy
s=s[0]
tmask=area.tmask[0,0]
weights =s 
weights.name = "weights"
weights=weights.fillna(0)







'/data/cmcc/ac28919/R025L75/ORCA025-DATA/v3.6/mesh_mask.nc'
darea='/data/cmcc/ac28919/R025L75/ORCA025-DATA/v3.6/mesh_mask.nc'
area=xr.open_dataset(darea,engine='netcdf4')
dx=area.e1t
dy=area.e2t
s=dx*dy
s=s[0]
tmask2=area.tmask[0,0]
weights2 =s
weights2.name = "weights"
weights2=weights2.fillna(0)





#cice=xr.open_dataset('INIT_MB0.cice.h.1974-12-31.nc')
'''ste=xr.open_dataset('/work/cmcc/c-glors/steric_level/steric_0-300_'+exp+'_1993-'+str(end)+'.nc')

stt=ste.thermosteric.where(tmask==1)
stt_weighted = stt.weighted(weights)
stt_weighted=stt_weighted.mean(dim=['x','y']).values



#sth=ste.halosteric.where(tmask==1).mean(dim=['x','y']).values
sth=ste.halosteric.where(tmask==1)
sth_weighted = sth.weighted(weights)
sth_weighted=sth_weighted.mean(dim=['x','y']).values
'''
ste2=xr.open_dataset('/work/cmcc/c-glors/steric_level/steric_0-'+lev+'_'+exp+'_1993-'+str(end)+'.nc')
stt2=ste2.thermosteric.where(tmask==1)
stt_weighted2 = stt2.weighted(weights)
stt_weighted2=stt_weighted2.mean(dim=['x','y']).values

#sth=ste.halosteric.where(tmask==1).mean(dim=['x','y']).values
sth2=ste2.halosteric.where(tmask==1)
sth_weighted2 = sth2.weighted(weights)
sth_weighted2=sth_weighted2.mean(dim=['x','y']).values

ste=xr.open_dataset('/work/cmcc/c-glors/steric_level/CGLORS7/steric_0-'+lev+'_CGLORS7_1993-'+str(end)+'.nc')
stt=ste.thermosteric.where(tmask==1)
stt_weighted = stt.weighted(weights)
stt_weighted=stt_weighted.mean(dim=['x','y']).values

#sth=ste.halosteric.where(tmask==1).mean(dim=['x','y']).values
sth=ste.halosteric.where(tmask==1)
sth_weighted = sth.weighted(weights)
sth_weighted=sth_weighted.mean(dim=['x','y']).values



'''ste3=xr.open_dataset('/work/cmcc/c-glors/steric_level/steric_0-6000_'+exp+'_1993-'+str(end)+'.nc')

stt3=ste3.thermosteric.where(tmask==1)
stt_weighted3 = stt3.weighted(weights)
stt_weighted3=stt_weighted3.mean(dim=['x','y']).values



#sth=ste.halosteric.where(tmask==1).mean(dim=['x','y']).values
sth3=ste3.halosteric.where(tmask==1)
sth_weighted3 = sth3.weighted(weights)
sth_weighted3=sth_weighted3.mean(dim=['x','y']).values
'''

cice2=xr.open_dataset('/work/cmcc/aspect/CESM2/rea_archive/SLAMB0/DAILY/SLAMB0.cice.h.2000-08-16.nc')
lat=cice2.TLAT
lon=cice2.TLON


print(stt_weighted)

fig = plt.figure(figsize=(8,4))
ax = plt.gca()
ax.grid(color='black', linestyle='dotted', axis='x',linewidth=0.8,alpha=0.2)
ax.grid(color='black', linestyle='dotted', axis='y',linewidth=0.8,alpha=0.2)
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
ax.set_xlim(0,timelen*12)
ax.set_xticks(np.arange(0,timelen*12,12))
ax.set_xticklabels(np.arange(start,end+1,1),rotation=300)
ax.set_ylabel('mm',fontsize=14,labelpad=1)
ax.plot(np.arange(0,timelen*12,1),sth_weighted*1000,linewidth=0.8,color='#de2e39',label='CGLORSv7')
ax.plot(np.arange(0,timelen*12,1),sth_weighted2*1000,linewidth=0.8,color='#0089ff',label='CG8E7')
ax.legend()
#ax.plot(np.arange(0,timelen*12,1),sth_weighted,linewidth=0.8,color='blue',label='Halosteric')
plt.savefig(exp+'_halosteric_timeseries_'+lev+'.png', format='png',  bbox_inches='tight', dpi=400)


fig = plt.figure(figsize=(8,4))
ax = plt.gca()
ax.grid(color='black', linestyle='dotted', axis='x',linewidth=0.8,alpha=0.2)
ax.grid(color='black', linestyle='dotted', axis='y',linewidth=0.8,alpha=0.2)
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
ax.set_xlim(0,timelen*12)
#ax.set_ylim(0,10)
ax.set_xticks(np.arange(0,timelen*12,12))
ax.set_xticklabels(np.arange(start,end+1,1),rotation=300)
ax.set_ylabel('mm',fontsize=14,labelpad=1)
ax.plot(np.arange(0,timelen*12,1),sth_weighted*1000+stt_weighted*1000,linewidth=0.8,color='#de2e39',label='CGLORSv7')
ax.plot(np.arange(0,timelen*12,1),1000*sth_weighted2+stt_weighted2*1000,linewidth=0.8,color='#0089ff',label='CG8E7')
ax.legend()
#ax.plot(np.arange(0,timelen*12,1),sth_weighted,linewidth=0.8,color='blue',label='Halosteric')
plt.savefig(exp+'_totalsteric_timeseries_'+lev+'.png', format='png',  bbox_inches='tight', dpi=400)


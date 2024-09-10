# sea ice area/extent/volume calculation.


import netCDF4 as nc
from netCDF4 import Dataset
import xarray as xr
import numpy as np
import cartopy
import matplotlib.pyplot as plt
import pandas as pd
import glob 
import os.path
from pathlib import Path
import re
#d='/work/csp/aspect/CESM2/rea_archive/MB0/DAILY/'
#d='/work/oda/is14120/input/SLAMB1/'
d='/work/cmcc/c-glors/CESM2/rea_archive/CG8E8/MONTHLY/'
#d='/data/csp/aspect/CESM2/rea_archive/MB0/MONTHLY_ICE/'
#d='/data/csp/aspect/CESM2/rea_archive/INIT_MB0/MONTHLY_ICE/'
#d='/work/oda/ac28919/C-GLORSv7/ICE/'
#dfile='SLAMB1.cice.h.1993_2020.nc'
#dfile='INIT_MB0_1m.cice.*'
#dfile='MB0_1m.cice.*'
dfile='CG8E8.cice.h.*'
#dfile='NEMO_1d_*'
darea='/work/cmcc/aspect/CESM2/rea_archive/SLAMB0/DAILY/SLAMB0.cice.h.2000-08-16.nc'
#filename="/work/oda/ac28919/CESM2/inputdata/STATIC_DA_NEMO42/mesh_mask.nc"
#da1=Dataset(filename,'r')


ll=[]
'''list=Path(d).glob(dfile+'.nc')
for path in list:
  #  y_str=re.findall('[0-9]+', str(path))
    ll.append(str(path)[-9:-5])
for path in sorted(ll): #Path(d).glob('nemo*_sie_nh.dat'):
    print(path)    
print(min(ll),max(ll))
start=int(min(ll))
end=int(max(ll))'''
start=1992
end=1996
timelen=end-start+1
hem='sh'
var='sie'
def tofile(data_list,name, var=var,hem=hem):
   # mylist=np.reshape(mylist,(len(data_list),len(data_list[0])))
    outfile ='/users_home/cmcc/is14120/dat/'+name+'_'+var+'_'+hem+'.dat'
    outf = open((outfile),'w')
    outf.write(var)
    outf.write('\n')
    for t in (data_list):
        outf.write(str(t)+' ')
        outf.write('\n')
    outf.close()
    return outf
# CICE
ds = xr.open_mfdataset(d+dfile,combine='by_coords',engine='netcdf4')
area=xr.open_dataset(darea,engine='netcdf4')
s=area.tarea/1.e12
sie_n=s.where((ds.aice>0.15)&(ds.TLAT>0)).sum(dim=['ni','nj'])
sia_n=(s.where((ds.aice>0.15)&(ds.TLAT>0))*ds.aice).sum(dim=['ni','nj'])
vol_n=(s.where((ds.aice>0.15)&(ds.TLAT>0))*ds.hi).sum(dim=['ni','nj'])


sie_s=s.where((ds.aice>0.15)&(ds.TLAT<0)).sum(dim=['ni','nj'])
sia_s=(s.where((ds.aice>0.15)&(ds.TLAT<0))*ds.aice).sum(dim=['ni','nj'])
vol_s=(s.where((ds.aice>0.15)&(ds.TLAT<0))*ds.hi).sum(dim=['ni','nj'])
#print(sie_n.values)
name='CG8E8_'+str(start)+'_'+str(end)
tofile(sie_n.values,name,var='sie',hem='nh')
tofile(sia_n.values,name,var='sia',hem='nh')
tofile(sie_s.values,name,var='sie',hem='sh')
tofile(sia_s.values,name,var='sia',hem='sh')
tofile(vol_n.values,name,var='vol',hem='nh')
tofile(vol_s.values,name,var='vol',hem='sh')

start=2010
end=2023
timelen=end-start+1
#CRYO
name='cryo'
ds = xr.open_mfdataset('/work/cmcc/c-glors/CESM2/inputdata/SIT/CryoSat_SMOS_merged/cryo_merged.nc',engine='netcdf4')
sit=ds.analysis_sea_ice_thickness
sic=ds.sea_ice_concentration
sit2=ds.cryosat_sea_ice_thickness
area=(25.*25.)/1.e6 #/1.e12
latmask=ds.lat.where(ds.lat>60,0)
latmask=latmask.where(latmask==0,1)


print(latmask)
area_a=np.full((432, 432), area)
#vol_s=area*(sit.where((sic>15)&(sit>0)&(ds.lat<0))*sic.where((sic>15)&(sit>0)&(ds.lat<0))).sum(dim=['xc','yc'])
vol_n=(area_a*sit.where((sic>15)&(sit>0))*sic.where((sic>15)&(sit>0))).sum(dim=['xc','yc'])/100
vol_n2=(area_a*sit2.where((sic>15)&(sit2>0))*sic.where((sic>15)&(sit2>0))).sum(dim=['xc','yc'])/100
#sit_n=sit.where((sic>15)&(sit>0)).mean(dim=['xc','yc'])
#sit_n2=sit.where((sic>15)&(sit2>0)).mean(dim=['xc','yc'])
#si_n=(area_a*sic.where((sic>1)&(sit>0))).sum(dim=['xc','yc'])/100
#si_n2=(area_a*sic.where((sic>1)&(sit2>0))).sum(dim=['xc','yc'])/100
#vol_n=area*sit.where((sic>15)&(sit>0)&(ds.lat>0)).sum(dim=['xc','yc'])
area_n=area*latmask.sum(dim=['xc','yc'])
#vol_s=(area.where((sic>15)&(sit>0)&(ds.lat<0))*sit).sum(dim=['yc','xc'])
#vol_n=(area.where((sic>15)&(sit>0)&(ds.lat>0))*sit).sum(dim=['yc','xc'])


print(area_n.values)
kk=0
new_vol=[]
for nn in range(timelen-1):
    for mm in range(12):
        kk=kk+1
        if mm in [1,2,3,4,5]:
            new_vol.append(vol_n.values[nn*7+mm])
        else:
            new_vol.append(np.nan)

        



#print(new_vol)
tofile(new_vol,name,var='vol',hem='nh')

#tofile(vol_s.values/100,name,var='vol',hem='sh')
#OSI
'''dfile='osisaf_nh.nc'
ds = xr.open_mfdataset(d+dfile,engine='netcdf4')
s = 25.*25./1.e6 
sic=ds.ice_conc.where(ds.ice_conc>15)
sic=sic/sic
sie_n_osi=sic.sum(dim=['xc','yc'])*s
sia_n_osi=s*ds.ice_conc.where(ds.ice_conc>15).sum(dim=['xc','yc'])/1.e2
tofile(sia_n_osi.values,'osisaf',var='sia',hem='nh')
tofile(sie_n_osi.values,'osisaf',var='sie',hem='nh')
dfile='osisaf_sh.nc'
ds = xr.open_mfdataset(d+dfile,engine='netcdf4')
s = 25.*25./1.e6
sic=ds.ice_conc.where(ds.ice_conc>15)
sic=sic/sic
sie_s_osi=sic.sum(dim=['xc','yc'])*s
sia_s_osi=s*ds.ice_conc.where(ds.ice_conc>15).sum(dim=['xc','yc'])/1.e2
tofile(sia_s_osi.values,'osisaf',var='sia',hem='sh')
tofile(sie_s_osi.values,'osisaf',var='sie',hem='sh')
#HADISST
#d='/work/csp/aspect/CESM2/inputdata/SIC/HADISST/'
dfile='HadISST_ice.nc'
ds = xr.open_mfdataset(d+dfile,engine='netcdf4')
s=(111.1**2)*np.cos(np.deg2rad(ds.latitude))/1.e6
sie_n_had=s.where((ds.sic>0.15)&(ds.latitude>0)).sum(dim=['latitude','longitude'])   #.sel(time=slice('1979','1985'))
sia_n_had=(s.where((ds.sic>0.15)&(ds.latitude>0))*ds.sic).sum(dim=['latitude','longitude'])    #.sel(time=slice('1979','1985'))
sie_s_had=s.where((ds.sic>0.15)&(ds.latitude<0)).sum(dim=['latitude','longitude'])   #.sel(time=slice('1979','1985'))
sia_s_had=(s.where((ds.sic>0.15)&(ds.latitude<0))*ds.sic).sum(dim=['latitude','longitude'])    #.sel(time=slice('1979','1985'))

tofile(sia_n_had.values,'hadisst',var='sia',hem='nh')
tofile(sie_n_had.values,'hadisst',var='sie',hem='nh')

tofile(sia_s_had.values,'hadisst',var='sia',hem='sh')
tofile(sie_s_had.values,'hadisst',var='sie',hem='sh')
'''



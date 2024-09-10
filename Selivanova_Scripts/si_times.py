


# sea ice area times calc (volume is calculated with sivol variable)
# this script is messy because of various configurations/meshmasks etc. to see only the formulas look into timeseries.py

import netCDF4 as nc
from netCDF4 import Dataset
import xarray as xr
import numpy as np
import cartopy
import matplotlib.pyplot as plt
import pandas as pd
import os.path
import glob
import operator


d='/data/oda/is14120/highresmip/work/'
mdir='/data/oda/is14120/grids/'
out='/work/oda/is14120/dat/highresmip/'
start=1979
end=2014
timelen=end-start+1
hem='NH'


def hemi(lati,hem=hem):
    cond=lati<0 if hem=='SH' else lati>0 
    return cond






#sim='highres-future'    
sim='hist-1950'

vv=1 
#var='siconc'
vars=['siconc','sivol','uo']
var_names=['sia','vol','uo']

sectors={'total':[[0,360],[0,90]],'CA':[[0,360],[85,90]], 'B-K':[[10,100],[65,85]],'LV':[[100,145],[65,85]],'ESS':[[145,185],[65,85]],'B-C':[[185,235],[65,85]],'GD':[[-125,10],[65,85]]}
# lon==lon if key_s=='GD' else lon=np.where(lon<0,lon+360,lon)
#if sim=='hist-1950':
prod={ 0:['ECMWF','ECMWF','IFS',['LR','HR'],['orca_292x362_.nc','meshmask50_ORCA025_.nc']],   
      1:['EC-Earth-Consortium','EC','Earth3P',['','HR'],['orca_292x362_.nc','orca025.nc']],\
          2:['CNRM-CERFACS','CNRM','CM6-1',['','HR'],['orca_292x362_.nc','orca025.nc']],\
          3:['MOHC','HadGEM3','GC31',['LL','MM','HM'],['areacella_HadGEM3-GC31-LL.nc','areacella_HadGEM3-GC31-MM.nc','areacella_HadGEM3-GC31-HM.nc']],
          4:['CMCC','CMCC','CM2',['HR4','VHR4'],['orca025.nc']*2],\
          5:['MPI-M','MPI','ESM1-2',['HR','XR'],['MPI_grid.nc','MPI_grid.nc']]}
         
#else:      
# highres-future
 #   prod={1:['EC-Earth-Consortium','EC','Earth3P',['','HR'],['r1i1p2f1']*2,['v20190421','v20190802'],['orca_292x362.nc','orca_1050x1442.nc']],\
 #         2:['CNRM-CERFACS','CNRM','CM6-1',['','HR'],['r1i1p1f2']*2,['v20190314','v20190920'],['mesh_mask_CNRM-CM6_eORCA1L75.nc','mesh_mask_CNRM-CM6_eORCA025L75.nc']],\
 #         3:['MOHC','HadGEM3','GC31',['LL','MM','HM'],['r1i1p1f1']*3,['v20181206','v20190222','v20190301']],\
 #         4:['CMCC','CMCC','CM2',['HR4','VHR4'],['r1i1p1f1']*2,['v20190509','v20190509'],['orca_1050x1442.nc']*2],\
 #         5:['MPI-M','MPI','ESM1-2',['HR','XR'],['r1i1p1f1']*2,['v20190517','v20190517']]}

for vv in [1]:
    thr=15 if vv==0 else 0.02. # SIC threshold for sea ice extent calculation
    a=1.e14 if vv==0 else 1.e12
    for key_s in sectors:
       
        ### write the output to dat file
        def tofile(data_list,key_pr, var=var_names[vv],hem=hem):
            mylist=[]
            for t in range(len(data_list)):
                m2=np.asarray(m1[t])
                mylist.append(m2)
            mylist=np.reshape(mylist,(len(data_list),len(data_list[0])))
            suff='' if key_s=='' else '_'
            outfile = out+prod[key_pr][1]+'_'+prod[key_pr][2]+'_'+var_names[vv]+'_'+hem+'_'+sim+suff+key_s+'.dat'
            outf = open((outfile),'w')  
            outf.write(" ".join(['.' if x=='' else x for x in prod[key_pr][3]]))
            outf.write('\n')
            for n in range(len(data_list[0])):
                for t in range(len(data_list)):
                    outf.write(str(mylist[t,n])+' ')
                outf.write('\n') 
            outf.close()      
            return outf


        nu=-1 
        for key_pr in [0]:
            var=vars[vv]
            m1=[]
            for n, k in enumerate(prod[key_pr][3]):
                nu=nu+1  
                suf='' if (key_pr==1 and n==0) or (key_pr==2 and n==0) else '-'
                if key_pr in [0,1]:
                    mfdataDIR_si ='/data/oda/is14120/highresmip/work/'+var+'_'+prod[key_pr][1]+'-'+prod[key_pr][2]+suf+k+'_'+str(start)+'_'+str(end)+'.nc'
                
                else: 
                    mfdataDIR_si ='/data/oda/is14120/highresmip/work/'+var+'_'+prod[key_pr][1]+'-'+prod[key_pr][2]+suf+k+'_'+str(start)+'_'+str(end)+'_2.nc'
                mfdataDIR_a ='/data/oda/is14120/grids/'+prod[key_pr][4][n]
               
                ds = xr.open_mfdataset(mfdataDIR_si)
                ds_a=xr.open_dataset(mfdataDIR_a)
               # if key_pr in [5]:
                #    lat=ds_a.lat if 'lat' in list(ds_a.coords) else ds_a.latitude
                 #   lon=ds_a.lon if 'lon' in list(ds_a.coords) else ds_a.longitude
               # else:

                lat=ds.latitude  #if 'lat' in list(ds.coords) else ds.latitude
                lon=ds.longitude # if 'lon' in list(ds.coords) else ds.longitude
              #  print(list(ds_a.coords))
                lon2=ds_a.longitude
                lat2=ds_a.latitude  
                if key_s == 'GD':   
                   lon=lon.where(lon<180,lon-360,drop=False) # lon from sea ice file
                   lon2=lon2.where(lon2<180,lon2-360,drop=False) # lon from meshmask file
                else:

                    lon=lon.where(lon>0,lon+360,drop=False) 
                    lon2=lon2.where(lon2>0,lon2+360,drop=False)
                if len(np.shape(lat))==1 :
                    lon,lat=np.meshgrid(lon,lat)
        #            lon=xr.DataArray(lon, dims=['i','j'], attrs=ds.attrs)
        #            lat=xr.DataArray(lat, dims=['i','j'], attrs=ds.attrs)
                if key_pr in [0,1,2,4]: 
                    ###area calculation where the meshmasks present
                    dx=ds_a.e1t
                    dy=ds_a.e2t
                    s=dx*dy
                    s=s[0][0] if key_pr==0 and n in [1,2] else s[0]
             #       print(ds_a.dims,s)
                elif key_pr==3:
                    s=ds_a.areacella
                elif key_pr==2 and n==1:
                    s=ds.area
                else:
                    s=ds_a.areacello
                    #s=s[0]
               # if n in [0,2]:
                #print(s.dims)
               # if 'i' not in list(s.dims):
               #     if key_pr in [0]:
               #         s = s.rename({'x': 'j','y': 'i'})
               #     elif key_pr in [3]:
               #         s = s.rename({'lat': 'i','lon': 'j'})
               #     elif key_pr in [1,4]:
               #         s = s.rename({'x': 'i','y': 'j'})
                #    s = s.rename({'x': 'i','y': 'j'})
                s=s[:,2:] if key_pr in [4] else s # all these conditional statements because of different meshmask files
                siconc=ds.siconc if vv==0 else ds.sivol
                siconc=siconc.where((siconc >= 15) & (siconc <=100)) if vv==0 else siconc
                 # if key_pr in [5]:
              #      siconc = siconc.rename({'x': 'i','y': 'j'})
                #if key_pr in [3]:
                 #   siconc = siconc.rename({'lat': 'i','lon': 'j'})
                siconc=siconc[:,2:,:] if (key_pr==2 and n==0) else siconc ### slice only because of the meshmask that have two rows more
                lon=lon[2:,:] if key_pr in [2] and n==0 else lon
                lat=lat[2:,:] if key_pr in [2] and n==0 else lat
                
                lons=sectors[key_s][0]
                lats=sectors[key_s][1]
                cond=(lon>=lons[0])&(lon<lons[1])&(lat>=lats[0])&(lat<=lats[1]) # diffenent sectors
                cond2=(lon2>=lons[0])&(lon2<lons[1])&(lat2>=lats[0])&(lat2<=lats[1]) # the same conditions using meshmask lon and lat just to check 

                if key_pr in [3]:
                    #sia=s.where((cond2)&(siconc>=thr)).sum(dim=['lat','lon']). # sea ice extent (sum of grid cells where sic>15)
                    sia=(siconc.where((siconc>=thr)&cond)*s).sum(dim=['lat','lon'])  # sea ice area or sea ice volume (if variable is sivol) 
                else: 
                    #sia=s.where((cond2)&(siconc>=thr)).sum(dim=['i','j']) 
                    sia=(siconc.where((siconc>=thr)&cond)*s).sum(dim=['i','j'])   
                m1.append(tuple(sia.values/a))
               # print(sia.values/1e14, key_pr, k, key_s)    
                #print(tuple(sia.values/1.e14)) 
           # m1=np.reshape(m1,(len(prod[key_pr][3]),1))    
            m1=np.reshape(m1,(len(prod[key_pr][3]),timelen*12))
            tofile(m1,key_pr) # write to dat file

                  
    

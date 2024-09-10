


# seasonal cycle of zonal mean wind

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

d='/gws/nopw/j04/primavera5/stream1/CMIP6/HighResMIP/'
mdir='/gws/nopw/j04/primavera1/masks/mesh_masks/'
wpath='/work/scratch-nopw/is1420/'
start=1979
end=2014
timelen=end-start+1
hem='SH'

### func for hemisphere
def hemi(lati,hem=hem):
    cond=lati<0 if hem=='SH' else lati>0 
    return cond





### present or future run
#sim='highres-future'    
sim='hist-1950'

### choice of variable
vv=2 # 0 FOR SIT; 1 for SIV 2 for SIA
#var='siconc'
vars=['tos','sos','ua','rlds','tas','rsds','rsus','rlus','uo']
var_names=['tos','sos','ua, (m/s)','rlds','tas','rsds','rsus','rlus','uo']

sectors={'total':[[0,360],[0,90]]}#,'CA':[[0,360],[85,90]], 'B-K':[[10,100],[65,85]],'LV':[[100,145],[65,85]],'ESS':[[145,185],[65,85]],'B-C':[[185,235],[65,85]],'GD':[[-125,10],[65,85]]}
# lon==lon if key_s=='GD' else lon=np.where(lon<0,lon+360,lon)
#if sim=='hist-1950':
###dictionary with the models and their parameters
prod={0:['ECMWF','ECMWF','IFS',['LR','MR','HR'],['r1i1p1f1']*3,['v20180221','v20181119','v20170915'],[1,1,1],['#2c5f0f','#69bb1f','#69bb1f'],['solid','solid','dashed']],\
      1:['EC-Earth-Consortium','EC','Earth3P',['','HR'],['r1i1p2f1']*2,['v20190314','v20181212'],[1,1],['#6b6d69','#6b6d69'],['solid','dashed']],\
      2:['CNRM-CERFACS','CNRM','CM6-1',['','HR'],['r1i1p1f2']*2,['v20190401','v20190221'],[1.5,0.5],['#ff55c9','#ff55c9'],['solid','dashed']],\
       #3:['MOHC','HadGEM3','GC31',['LL'],['r1i1p1f1'],['v20170921']],\
      3:['MOHC','HadGEM3','GC31',['LL','MM','HM'],['r1i1p1f1']*3,['v20170921','v20170928','v20180730'],[1.8,1.5,1.5],['#910000','#ff1111','#ff1111'],['solid','solid','dashed']],\
      4:['CMCC','CMCC','CM2',['HR4','VHR4'],['r1i1p1f1']*2,['v20200917','v20200917'],[1,1],['#24b0ff','#24b0ff'],['solid','dashed']],
      5:['MPI-M','MPI','ESM1-2',['HR','XR'],['r1i1p1f1']*2,['v20180606','v20180606'],[1,1],['#ff922e','#ff922e'],['solid','dashed']]}


shmons=['J','F','M',"A",'M','J','J','A','S','O','N','D']

#colors = distinctipy.get_colors(14)
grep_conf={'NH':'Arctic_Ocean','SH':'Southern_Ocean'}
hem='SH'
mon=5
### geogr sectors
sectors={'P':[150,290],'I':[20,150],'A':[290,20],'total':[0,360]}
sectors1={'P':[150,-70],'I':[20,150],'A':[-70,20]}


#else:      
# highres-future
#    prod={1:['EC-Earth-Consortium','EC','Earth3P',['','HR'],['r1i1p2f1']*2,['v20190421','v20190802'],['orca_292x362.nc','orca_1050x1442.nc']],\
  #        2:['CNRM-CERFACS','CNRM','CM6-1',['','HR'],['r1i1p1f2']*2,['v20190314','v20190920'],['mesh_mask_CNRM-CM6_eORCA1L75.nc','mesh_mask_CNRM-CM6_eORCA025L75.nc']],\
    #      3:['MOHC','HadGEM3','GC31',['LL','MM','HM'],['r1i1p1f1']*3,['v20181206','v20190222','v20190301']],\
    #      4:['CMCC','CMCC','CM2',['HR4','VHR4'],['r1i1p1f1']*2,['v20190509','v20190509'],['orca_1050x1442.nc']*2],\
      #    5:['MPI-M','MPI','ESM1-2',['HR','XR'],['r1i1p1f1']*2,['v20190517','v20190517']]}

for vv in [2]:
    
   

        

    for key_s in sectors:



        ### func write output to dat
        def tofile(data_list,key_pr, var=var_names[vv],hem=hem):
            mylist=[]
            for t in range(len(data_list)):
                m2=np.asarray(m1[t])
                mylist.append(m2)
            mylist=np.reshape(mylist,(len(data_list),len(data_list[0])))
            #suff='_+sic' if key_pr==4 else ''
            suff='' if key_s=='' else '_'
            outfile =  prod[key_pr][1]+'_'+prod[key_pr][2]+'_'+var+'_lat_'+hem+'_'+sim+suff+key_s+'.dat' #!!!!!
            outf = open((outfile),'w')  
            outf.write(" ".join(['.' if x=='' else x for x in prod[key_pr][3]]))
            outf.write('\n')
            for n in range(len(data_list[0])):
                for t in range(len(data_list)):
                    outf.write(str(mylist[t,n])+' ')
                outf.write('\n') 
            outf.close()      
            return outf
        ### figure
        fig = plt.figure(figsize=(8,8))
        fig.set_facecolor("white")

        ax = plt.gca()
        ax.grid(color='black', linestyle='dotted', axis='x',linewidth=0.8,alpha=0.2)
        ax.grid(color='black', linestyle='dotted', axis='y',linewidth=0.8,alpha=0.2)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        ax.set_xlim(0,17)
        ax.set_xticks(np.arange(0,18,2))
        ax.set_xticklabels(np.arange(-76,-41,4))
        #ax.set_ylabel(variables[key_var][1],fontsize=16,labelpad=1)
        #ax.set_ylabel('MIZF (%)',fontsize=16,labelpad=1)
      #  ax.set_ylabel('MIZ area ($\mathregular{10^6}$ $\mathregular{km^2}$)',fontsize=16,labelpad=1)
        ax.set_ylabel(var_names[vv],fontsize=16,labelpad=1)
        if key_s in ['A','I','P']:
            ax.set_ylim(-8,15)




        ### model loop    
        for key_pr in prod:
            var=vars[vv]
            
            
            m1=[]
            ### model config loop
            for n, k in enumerate(prod[key_pr][3]):  
                suf='' if (key_pr==1 and n==0) or (key_pr==2 and n==0) else '-'
                
                #data=nc.Dataset(glob.glob('/work/scratch-nopw/is14120/'+var_lib+'/'+var+'_'+prod[key_pr][1]+'-'+prod[key_pr][2]+suf+k+'_'+str(start)+'_'+str(end)+'_monclim.nc')[0])
                
               # mfdataDIR='/work/scratch-nopw/is14120/'+var_lib+'/'+var+'_'+prod[key_pr][1]+'-'+prod[key_pr][2]+suf+k+'_'+str(start)+'_'+str(end)+'_monclim.nc'
                
                ### suffix in the file names
                suf='' if (key_pr==1 and n==0) or (key_pr==2 and n==0) else '-'
               
                mfdataDIR ='/work/scratch-nopw/is14120/work/'+var+'_'+prod[key_pr][1]+'-'+prod[key_pr][2]+suf+k+'_'+str(start)+'_'+str(end)+'.nc'
                ds = xr.open_mfdataset(mfdataDIR)
                lat=ds.lat if 'lat' in list(ds.coords) else ds.latitude
                lon=ds.lon if 'lon' in list(ds.coords) else ds.longitude

                ua=ds.data_vars[vars[vv]]
                ua=ua[:,1,:]
                ua=ua[mon::12]
                ua=ua.mean(dim='time')
                a=lon.where((lon>340)).values
               
                if key_s in ['P']:
                    if sum(a[np.isnan(a)==False])>0:
                        mask = ((lon >= sectors[key_s][0]) & (lon < sectors[key_s][1])& (lat<-40))
                    else:    
                        mask = ((lon >= sectors1[key_s][0]) | (lon < sectors1[key_s][1])& (lat<-40))
                elif key_s in ['I']:
                    mask = ((lon >= sectors[key_s][0]) & (lon < sectors[key_s][1])& (lat<-40))
                elif key_s in ['A']:    
                    if sum(a[np.isnan(a)==False])>0:
                        mask = ((lon >= sectors[key_s][0]) | (lon < sectors[key_s][1])& (lat<-40))
                    else:    
                        mask = ((lon >= sectors1[key_s][0]) & (lon < sectors1[key_s][1])& (lat<-40))
                elif key_s in ['total']:
                    mask= (lat<-40)
               # ua=ua.where(lat<=-39)
                ### mask the region
                ua=ua.where(mask) 
  # define two-degree wide latitude bins
  
                ### zonal mean calculation (should be area weighted if the grid is not regular)
                grad=prod[key_pr][6][n]
                lat_bins = np.arange(-77, -40.75, grad)
                d='stacked_i_j' if key_pr in [0] else 'stacked_j_i' if key_pr in [1,3,4] else 'stacked_y_2_x_2' if key_pr==5 else'stacked_y_x'
                #zonal = ua.groupby_bins(lat, lat_bins).mean(dim=[d],skipna=True)
                zonal = ua.groupby_bins(lat, lat_bins).mean(dim=['lat','lon'],skipna=True)
                print(len(zonal))
                zonal=np.array(zonal)
             #   print(grad)
                ax.plot(np.arange(0,18,grad/2),zonal,color=prod[key_pr][7][n],linestyle=prod[key_pr][8][n],linewidth=1.5,label=prod[key_pr][1]+' '+k)
              
                
        ###era5 dataset
        mfdataDIR ='/work/scratch-nopw/is14120/work/era5_1979_2014.nc'
        ds = xr.open_mfdataset(mfdataDIR)
        lat=ds.latitude 
        lon=ds.longitude
        u=ds.u.sel(time=slice('1979-01','2014-12'))
        u=u[mon::12]
        u=u.mean(dim='time')
        u=u.where(lat<=-39)
        a=lon.where((lon>340)).values
      #  sector=[150,290] 
       # sector1=[150,-70]
                
        if key_s in ['P']:
            if sum(a[np.isnan(a)==False])>0:
                mask = ((lon >= sectors[key_s][0]) & (lon < sectors[key_s][1])& (lat<-40))
            else:    
                mask = ((lon >= sectors1[key_s][0]) | (lon < sector1s[key_s][1])& (lat<-40))
        elif key_s in ['I']:
            mask = ((lon >= sectors[key_s][0]) & (lon < sectors[key_s][1])& (lat<-40))
        elif key_s in ['A']:    
            if sum(a[np.isnan(a)==False])>0:
                mask = ((lon >= sectors[key_s][0]) | (lon < sectors[key_s][1])& (lat<-40))
            else:    
                mask = ((lon >= sectors1[key_s][0]) & (lon < sectors1[key_s][1])& (lat<-40)) 
        elif key_s in ['total']:
            mask= (lat<-40)        
        u=u.where(mask)      

  ### define latitude bins
        lat_bins = np.arange(-77, -40, 1)
        zonal = u.groupby_bins(lat, lat_bins).mean(dim=['latitude','longitude'],skipna=True)
        print(len(zonal))
        zonal=np.array(zonal)
        ax.plot(np.arange(0,18,0.5),zonal,color='black',linestyle='solid',linewidth=1.5,label='ERA5')
 
        ### jra dataset
        mfdataDIR ='/work/scratch-nopw/is14120/work/jra_1979_2014.nc'
        ds = xr.open_mfdataset(mfdataDIR)
        lat=ds.latitude 
        lon=ds.longitude
        u=ds.u.sel(time=slice('1979-01','2014-12'))
        u=u[mon::12]
        u=u.mean(dim='time')
        u=u.where(lat<=-39)
        #sector=[150,290] 
        #sector1=[150,-70]
        a=lon.where((lon>340)).values
        if key_s in ['P']:
            if sum(a[np.isnan(a)==False])>0:
                mask = ((lon >= sectors[key_s][0]) & (lon < sectors[key_s][1])& (lat<-40))
            else:    
                mask = ((lon >= sectors1[key_s][0]) | (lon < sectors1[key_s][1])& (lat<-40))
        elif key_s in ['I']:
            mask = ((lon >= sectors[key_s][0]) & (lon < sectors[key_s][1])& (lat<-40))
        elif key_s in ['A']:    
            if sum(a[np.isnan(a)==False])>0:
                mask = ((lon >= sectors[key_s][0]) | (lon < sectors[key_s][1])& (lat<-40))
            else:    
                mask = ((lon >= sectors1[key_s][0]) & (lon < sectors1[key_s][1])& (lat<-40))
        elif key_s in ['total']:
            mask= (lat<-40)        
        u=u.where(mask)      

        
  ### define latitude bins
        lat_bins = np.arange(-77, -40, 2)
        zonal = u.groupby_bins(lat, lat_bins).mean(dim=['latitude','longitude'],skipna=True)
        print(len(zonal))
        zonal=np.array(zonal)
        ax.plot(np.arange(0,18,1),zonal,color='black',linestyle='--',linewidth=1.5,label='JRA55')        
        ax.legend(fontsize=7,loc=2,ncol=3)             
        plt.savefig(vars[vv]+'_'+hem+'_seas_'+str(start)+'_'+str(end)+'_'+str(mon)+key_s+'.png', format='png',  bbox_inches='tight', dpi=400)

                
      
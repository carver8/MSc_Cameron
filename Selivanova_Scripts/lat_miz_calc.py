# lat miz calculation

import matplotlib as mpl
#from scipy.ndimage.filters import gaussian_filter1d
#from scipy.stats import skewnorm
from matplotlib.patches import Polygon
#import cmocean as cmocean
import matplotlib.pyplot as plt
#import scipy.stats
#from scipy import stats
import numpy as np  #numerical python tool
import pandas as pd
import netCDF4 as nc
from netCDF4 import Dataset
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#from scipy.stats import pearsonr 
import os.path
import numpy.ma as ma
reans=['CGLORS','GLORYS2V4','FOAM','ORAS5','ENSEMBLE_MEAN'] 
reans1=['CGLORS','GLORYS2V4','FOAM','ORAS5','ENS']    
sat_datasets=['cdr','osi','ifr']
sat_colors=['#8da0cb','#fc8d62','#66c2a5']

#colors=['blue','green','darkviolet','orange','red']
#colors=['#008dff','#59972c','#591071', '#fbb74e', '#e24390']
colors=['#008dff','#59972c','#591071', '#fbb74e', '#e24390']

shmons=['J','F','M',"A",'M','J','J','A','S','O','N','D']
shnames=['S Ocean','Indian', 'West P','Ross', 'A&B','Weddell']
satst=['solid', '--','dotted']
bounds=[[0,360],[20,90],[90,160],[160,230],[230,300],[-60,20]] # lon0, the last - lon1
names=['Southern Ocean','Indian Ocean', 'Western Pacific','Ross Sea', 'A&B','Weddell Sea']
names1=['Southern_Ocean','Indian_Ocean', 'Western_Pacific','Ross_Sea', 'A&B','Weddell_Sea']
path='/work/data/ORA-IP/'
path_sat='/work/data/sat_ice/nsidc/'
path_home='/work2/iuliiaselivanova/'
#area 
mons = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
#sat_datasets=['cdr']

#satst=['solid']
start=2020
end=2020
timelen=end-start+1

mons=['02','05','09','12']
k=0

'''# cdr
sat_lat=np.zeros((4,timelen))
for n in range(4):
    all_sic=np.zeros((timelen,1021,1442))
    for yr in range(start,end+1):
        data_sat=nc.Dataset(path_sat+'orca_'+str(yr)+mons[n]+'.nc')
        sic_cdr=data_sat.variables['seaice_conc_monthly_cdr'][0,:] 
        sic_cdr=np.where(np.logical_and(0.01 < sic_cdr, sic_cdr <= 1),sic_cdr*100,np.nan)
        all_sic[yr-start,:]=sic_cdr
    mean_lat3=[]

    for yr in range(timelen):
        mylist1=[]
        for kk in range(-180,180,4): 
            cond=lat[np.where((lat<-50) & (abs(lon-kk)<0.3)& (all_sic[yr]<=80)&(all_sic[yr]>15))]

            mylist1.append(np.median(cond))
        mean_lat3.append(np.nanmedian(np.where(np.asarray(mylist1)==0,np.nan,np.asarray(mylist1))))
    sat_lat[n]=mean_lat3
outfile =  path_home+'dat/MIZ_lat_'+'cdr'+'.txt'
outf = open((outfile),'w')  
for n in range(timelen):
    for t in range(4):
        outf.write(str(sat_lat[t,n])+' ')
    outf.write('\n') 
outf.close() '''  


# reans
sat_lat=np.zeros((4,timelen))
start=2020

for n in range(4):

    all_sic=np.zeros((timelen,1050,1442))
   # all_sic=np.zeros((timelen,1021,1442))

    for yr in range(start,end+1):
        if os.path.isfile(path+reans[k]+"/"+reans[k]+"_"+str(yr)+mons[n]+"_icemod.nc"):
            data=nc.Dataset(path+reans[k]+"/"+reans[k]+"_"+str(yr)+mons[n]+"_icemod.nc")
            lat=data.variables['nav_lat'][:]
            lon=data.variables['nav_lon'][:]

        
        if k==0 and yr==2020:
            sic=data.variables['soicecov'][0,:]
        else:
            sic=data.variables['ileadfra'][0,:]
        sic=np.where(sic>0.01,sic,0)
        sic=np.where(sic<1,sic*100,np.nan)     
        all_sic[yr-start,:]=sic


    mean_lat3=[]

    for yr in range(start-start,end-start+1):
        mylist1=[]
        for kk in range(-180,180,4): 
            cond=lat[np.where((lat<-50) & (abs(lon-kk)<0.3)& (all_sic[yr]<=80)&(all_sic[yr]>15))]

            mylist1.append(np.median(cond))
        print(mylist1)
        mean_lat3.append(np.nanmedian(np.where(np.asarray(mylist1)==0,np.nan,np.asarray(mylist1))))
    sat_lat[n]=mean_lat3

print(sat_lat[0,:])



outfile =  path_home+'dat/MIZ_lat_'+reans[k]+'.txt'
outf = open((outfile),'a')  
for n in range(timelen):
    for t in range(4):
        outf.write(str(sat_lat[t,n])+' ')
    outf.write('\n') 
outf.close() 
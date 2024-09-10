# 2 subplots seasonal cycle
# here ensemble members of cnrm model



import netCDF4 as nc
from netCDF4 import Dataset
import xarray as xr
import numpy as np
import cartopy
import matplotlib.pyplot as plt
import pandas as pd
#from distinctipy import distinctipy
import os.path
from scipy.stats import iqr


start=1979 # obs
end=2014
timelen=end-start+1


start1=1979 # models
timelen1=end-start1+1
st=start1-1950
sims=['hist-1950','highres-future']
#mydict={0:['HR','black'],1:['MR','darkgreen'],2:['LR','darkblue']}
#prod={0:['ECMWF','ECMWF','IFS',['HR','MR','LR'],'v20170915']}
prod={0:['ECMWF','ECMWF','IFS',['LR','HR'],['r1i1p1f1']*2,['v20180221','v20181119','v20170915'],['#2c5f0f','#69bb1f','#69bb1f'],['solid','dashed']],\
      1:['EC-Earth-Consortium','EC','Earth3P',['','HR'],['r1i1p2f1']*2,['v20190314','v20181212'],['#6b6d69','#6b6d69'],['solid','dashed']],\
      2:['CNRM-CERFACS','CNRM','CM6-1',['','HR'],['r1i1p1f2']*2,['v20190401','v20190221'],['#ff55c9','#ff55c9'],['solid','dashed']],\
      3:['MOHC','HadGEM3','GC31',['LL','MM','HM'],['r1i1p1f1']*3,['v20170921','v20170928','v20180730'],['#882255','#cc6677','#cc6677'],['solid','solid','dashed']],\
       # 3:['MOHC','HadGEM3','GC31',['LL','HM'],['r1i1p1f1']*2,['v20170921','v20180730'],\
#['#910000','#ff1111'],['solid','dashed']],\
      4:['CMCC','CMCC','CM2',['HR4','VHR4'],['r1i1p1f1']*2,['v20200917','v20200917'],['#24b0ff','#24b0ff'],['solid','dashed']],
      5:['MPI-M','MPI','ESM1-2',['HR','XR'],['r1i1p1f1']*2,['v20180606','v20180606'],['#ff922e','#ff922e'],['solid','dashed']]}
shmons=['J','F','M',"A",'M','J','J','A','S','O','N','D']

#colors = distinctipy.get_colors(14)
grep_conf={'NH':'Arctic_Ocean','SH':'Southern_Ocean'}
hem='NH'

variables={0:['siconc','SIC (%)'],1:['sithick','SIT (m)'],2:['sia','SIA ($\mathregular{10^6}$ $\mathregular{km^2}$)'],3:['vol','SIV ($\mathregular{10^3}$ $\mathregular{km^3}$)'],4:['sisnthick','snow depth (m)'],5:['tas','T'+ u'\u2103'], 6:['tos','SST'+ u'\u2103)'],7:['sos','SSS, (psu)'],8:['ua','m/s'],
9:['rlds','W/$\mathregular{m^2}$'],10:['rsds','W/$\mathregular{m^2}$'],11:['rlus','W/$\mathregular{m^2}$'],12:['rsus','W/$\mathregular{m^2}$'],13:['clt', 'cloud fraction (%)']}
# for ua max add manually 'max'
#key_var=2
#var=variables[key_var][0]
ice_class={0:'',1:'_pack',2:'_miz',3:'_mizf'}


def len_fix(data,orig_start,orig_end,start=start,end=end):
    if end >=orig_end:
        data=np.concatenate((data,[np.nan]*12*(end-orig_end)))
    else:
        data=data[:12*(end-orig_start+1)]   
    if start<orig_start:
        data=np.concatenate(([np.nan]*12*(orig_start-start),data))
    else: 
        data=data[12*(start-orig_start):]    
    return data
def tofile(data_list, var='sia',hem=hem):
   # mylist=np.reshape(mylist,(len(data_list),len(data_list[0])))
    outfile =  'CNRM_CM6-1'+'_'+var+'HR_'+hem+'_ens.dat'
    outf = open((outfile),'w')
    outf.write(var)
    outf.write('\n')
    for t in (data_list):
        outf.write(str(t)+' ')
        outf.write('\n')
    outf.close()
    return outf
d='/home/users/is14120/'







fig, (ax1, ax2) = plt.subplots(1,2, figsize=(8.27,3.5))
#fig.subplots_adjust(hspace=0.25,wspace=0.3)
key_var=2
ii=0
key_i=ice_class[ii]
var=variables[key_var][0]


if key_var==2 and end <2022:
       # key_i='_miz'
        
    cdr=pd.read_csv(d+'dat_jasmin/cdr'+'_'+var+'_'+hem+'_total'+key_i+'.dat',delim_whitespace=True,header=0)[:]
    cdr=np.asarray(cdr) 
       
    cdr=len_fix(cdr[:,0],1979,2021)
    cdr=np.where(cdr==0,np.nan,cdr)  
    cdr=np.nanmean(np.reshape(cdr,(12,timelen),order='F'),axis=1)   
    '''cdr1=pd.read_csv('cdr'+'_'+var+'_'+hem+'.dat',delim_whitespace=True,header=0)[:]
    cdr1=np.asarray(cdr1)     
    cdr1=len_fix(cdr1[:,0],1979,2021)
    cdr1=np.where(cdr1==0,np.nan,cdr1)
    cdr1=np.nanmean(np.reshape(cdr1,(12,timelen),order='F'),axis=1)   '''
       # cdr=100*cdr/cdr1
        
      
    osi=pd.read_csv(d+'dat_jasmin/osisaf'+'_'+var+'_'+hem+''+key_i+'.dat',delim_whitespace=True,header=0)[:]
    osi=np.asarray(osi) 
    osi=len_fix(osi[:,0],1979,2021)
    osi=np.where(osi==0,np.nan,osi) 
    osi=np.nanmean(np.reshape(osi,(12,timelen),order='F'),axis=1)
        
    osi1=pd.read_csv(d+'dat_jasmin/osisaf'+'_'+var+'_'+hem+'.dat',delim_whitespace=True,header=0)[:]
    osi1=np.asarray(osi1) 
    osi1=len_fix(osi1[:,0],1979,2021)
    osi1=np.where(osi1==0,np.nan,osi1)
    osi1=np.nanmean(np.reshape(osi1,(12,timelen),order='F'),axis=1)       
      #  osi=100*osi/osi1


   # fig.set_facecolor("white")
#ax1 = plt.gca()
ax1.grid(color='black', linestyle='dotted', axis='x',linewidth=0.8,alpha=0.2)
ax1.grid(color='black', linestyle='dotted', axis='y',linewidth=0.8,alpha=0.2)
ax1.tick_params(axis='x', labelsize=8)
ax1.tick_params(axis='y', labelsize=8)
ax1.set_xlim(0,11)
ax1.set_ylim(0,24)

ax1.set_xticks(np.arange(0,12,1))
ax1.set_xticklabels(shmons)

ax1.set_ylabel(variables[key_var][1],fontsize=10,labelpad=1)

nu=0

multi=[]
key_pr=2
for n,k in enumerate(prod[key_pr][3]):   
    ens=[]
    for m in ['_total','_2','_3']:
        if key_pr<8: 
            #try:
              #  df1=pd.read_csv(d+'dat_jasmin/'+prod[key_pr][1]+'_'+prod[key_pr][2]+'_'+var+'_'+hem+'_'+sims[0]+'_total'+key_i+m+'.dat',delim_whitespace=True,header=0)[:]
               
            
           # except:
            df1=pd.read_csv(d+'dat_jasmin/'+prod[key_pr][1]+'_'+prod[key_pr][2]+'_'+var+'_'+hem+'_'+sims[0]+''+key_i+m+'.dat',delim_whitespace=True,header=0)[:]

            data=np.asarray(df1) 
            if  len(data)>500:
                data=data[(start-1950)*12:,:]     #!!!! 
            data=data[:(end-start+1)*12,:]
            ens.append(data[:,n])
        print(var,m,np.nanmean(np.reshape(data[:,n],(12,timelen1),order='F'),axis=1)) 
        #if m=='_ensmean':
         #   ax1.plot(np.arange(0,12,1),np.nanmean(np.reshape(data[:,n],(12,timelen1),order='F'),axis=1),color='#2c5f0f',linestyle=prod[key_pr][7][n],linewidth=1.5,label=k+' EM')      
      #  else:    
        ax1.plot(np.arange(0,12,1),np.nanmean(np.reshape(data[:,n],(12,timelen1),order='F'),axis=1),color='#69bb1f',linestyle=prod[key_pr][7][n],linewidth=0.8)
        nu += 1
    ens=np.nanmean(ens,axis=0)      
    if n==1:
        tofile(ens)
    #print(ens)
    ax1.plot(np.arange(0,12,1),np.nanmean(np.reshape(ens,(12,timelen1),order='F'),axis=1),color='#2c5f0f',linestyle=prod[key_pr][7][n],linewidth=1.5,label=k+' EM')
if key_var in [1,3] and end<2022: 
    if hem=='NH':
        ax1.plot(np.arange(0,12,1),pi,linewidth=1.5,label='PIOMAS',color='black')
    elif hem=='SH':
        ax1.plot(np.arange(0,12,1),pi,linewidth=1.5,label='GIOMAS',color='black')
        ax1.plot(np.arange(0,12,1),grep,linewidth=1.5,label='GREP',color='cyan')
     
       # ax.plot(np.arange(0,12,1),gi,linewidth=2,label='GIOMAS',color='black',linestyle='dashed')
       
elif key_var==2 and end<2022:
       
    ax1.plot(np.arange(0,12,1),cdr,linewidth=1.5,label='CDR',color='black')      
    ax1.plot(np.arange(0,12,1),osi,linewidth=1.5,label='OSISAF',color='black',linestyle='dashed')

ax1.legend(fontsize=4.5,loc=2,ncol=2) 
        
 
ii=0
key_i=ice_class[ii]
key_var=2
var=variables[key_var][0]

if key_var in [1,3] and end<2022 and hem=='NH':
    pi=pd.read_csv(d+'dat_jasmin/piomas'+'_'+var+'_'+hem+'.dat',delim_whitespace=True,header=0)[:]
    pi=np.asarray(pi) 
    pi=len_fix(pi,1979,2020)
    pi=np.nanmean(np.reshape(pi,(12,timelen),order='F'),axis=1)
elif key_var in [3] and end<2022 and hem=='SH':
    pi=pd.read_csv('giomas'+'_'+var+'_'+hem+'_total.dat',delim_whitespace=True,header=0)[:]
    pi=np.asarray(pi) 
    pi=len_fix(pi,1979,2020)
    pi=np.nanmean(np.reshape(pi,(12,timelen),order='F'),axis=1)
    grep=pd.read_csv(grep_conf[hem]+'_ENSEMBLE_MEAN_timeseries_'+variables[key_var][0]+'.dat',delim_whitespace=True,header=0)[:]
    grep=np.asarray(grep) 
    grep=len_fix(grep[:,ii],1993,2019)
    grep=np.nanmean(np.reshape(grep,(12,timelen),order='F'),axis=1)

   # fig.set_facecolor("white")
#ax1 = plt.gca()
ax2.grid(color='black', linestyle='dotted', axis='x',linewidth=0.8,alpha=0.2)
ax2.grid(color='black', linestyle='dotted', axis='y',linewidth=0.8,alpha=0.2)
ax2.tick_params(axis='x', labelsize=8)
ax2.tick_params(axis='y', labelsize=8)
ax2.set_xlim(0,11)
#ax2.set_ylim(0,100)
ax2.set_ylim(0,95)   # NH
#ax2.set_ylim(0,30)

ax2.set_xticks(np.arange(0,12,1))
ax2.set_xticklabels(shmons)
ax2.set_ylabel(variables[key_var][1],fontsize=10,labelpad=1)

nu=0
multi=[]
key_pr=2
for n,k in enumerate(prod[key_pr][3]):   
    ens=[]

    for m in ['_total','_2','_3']:
        print(m)
        if key_pr<8: 
            try:
                df1=pd.read_csv(d+'dat_jasmin/'+prod[key_pr][1]+'_'+prod[key_pr][2]+'_'+var+'_'+hem+'_'+sims[0]+'_total'+key_i+m+'.dat',delim_whitespace=True,header=0)[:]
               
            
            except:
                df1=pd.read_csv(d+'dat_jasmin/'+prod[key_pr][1]+'_'+prod[key_pr][2]+'_'+var+'_'+hem+'_'+sims[0]+''+key_i+m+'.dat',delim_whitespace=True,header=0)[:]


            data=np.asarray(df1) 
          
            if  len(data)>500:
                data=data[(start-1950)*12:,:]     #!!!! 
            data=data[:(end-start+1)*12,:]
            ens.append(data[:,n])
  
        ax2.plot(np.arange(0,12,1),np.nanmean(np.reshape(data[:,n],(12,timelen1),order='F'),axis=1),color='#69bb1f',linestyle=prod[key_pr][7][n],linewidth=0.8)
        
        nu += 1
    print(np.shape(ens))   
    ens=np.nanmean(ens,axis=0)      
    
    ax2.plot(np.arange(0,12,1),np.nanmean(np.reshape(ens,(12,timelen1),order='F'),axis=1),color='#2c5f0f',linestyle=prod[key_pr][7][n],linewidth=1.5,label=k+' EM')
if key_var in [1,3] and end<2022: 
    if hem=='NH':
        ax2.plot(np.arange(0,12,1),pi,linewidth=1.5,label='PIOMAS',color='black')
    elif hem=='SH':
        ax2.plot(np.arange(0,12,1),pi,linewidth=1.5,label='GIOMAS',color='black')
        ax2.plot(np.arange(0,12,1),grep,linewidth=1.5,label='GREP',color='cyan')
     
       # ax.plot(np.arange(0,12,1),gi,linewidth=2,label='GIOMAS',color='black',linestyle='dashed')
       
elif key_var==2 and end<2022:
       
    ax2.plot(np.arange(0,12,1),cdr,linewidth=1.5,label='CDR',color='black')      
    ax2.plot(np.arange(0,12,1),osi,linewidth=1.5,label='OSISAF',color='black',linestyle='dashed')
    ax.plot(np.arange(0,12,1),grep,linewidth=2,label='GREP',color='cyan')

        
 
ax2.legend(fontsize=4.5,loc=1,ncol=2) 

#ax2.legend(fontsize=6,loc=2,ncol=1) 
plt.savefig(var+'_'+hem+'_seas_'+str(start)[2:]+str(end)[2:]+key_i+'_2subs_cnrm.png', format='png',  bbox_inches='tight', dpi=500)



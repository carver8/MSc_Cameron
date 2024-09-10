# times plot and trends 
from scipy import stats
import netCDF4 as nc
from netCDF4 import Dataset
import xarray as xr
import numpy as np
import cartopy
import matplotlib.pyplot as plt
import pandas as pd
import os.path
from scipy.stats import pearsonr
from scipy import signal

d='/home/users/is14120/dat_jasmin/'


start=1979
end=2014
timelen=end-start+1
sectors={'CA':[[0,360],[85,90]], 'B-K':[[10,100],[65,85]],'LV':[[100,145],[65,85]],'ESS':[[145,185],[65,85]],'B-C':[[185,235],[65,85]],'GD':[[-125,10],[65,85]]}#,'':[[0,360],[65,90]]}
sectors={'':[[0,360],[0,90]]}
prod={0:['ECMWF','ECMWF','IFS',['LR','MR','HR'],['r1i1p1f1']*3,['v20180221','v20181119','v20170915'],['#2c5f0f','#69bb1f','#69bb1f'],['solid','solid','dashed']],\
      1:['EC-Earth-Consortium','EC','Earth3P',['','HR'],['r1i1p2f1']*2,['v20190314','v20181212'],['#6b6d69','#6b6d69'],['solid','dashed']],\
      2:['CNRM-CERFACS','CNRM','CM6-1',['','HR'],['r1i1p1f2']*2,['v20190401','v20190221'],['#ff55c9','#ff55c9'],['solid','dashed']],\
      3:['MOHC','HadGEM3','GC31',['LL','MM','HM'],['r1i1p1f1']*3,['v20170921','v20170928','v20180730'], ['#882255','#cc6677','#cc6677'],['solid','solid','dashed']],\
      4:['CMCC','CMCC','CM2',['HR4','VHR4'],['r1i1p1f1']*2,['v20200917','v20200917'],['#24b0ff','#24b0ff'],['solid','dashed']],
      5:['MPI-M','MPI','ESM1-2',['HR','XR'],['r1i1p1f1']*2,['v20180606','v20180606'],['#ff922e','#ff922e'],['solid','dashed']]}

hem='NH'
sims=['hist-1950','highres-future']
variables={0:['siconc','%'],1:['sithick','SIT (m)'],2:['sia','SIA ($\mathregular{10^6}$ $\mathregular{km^2}$)'],3:['vol','SIV ($\mathregular{10^3}$ $\mathregular{km^3}$)'],4:['tas','T'+ u'\u2103'],7:['sos','SSS, (psu)'],8:['tos','SST'+ u'\u2103)'],9:['rsds','W/$\mathregular{m^2}$'],10:['ua','m/s']}
key_var=2
var=variables[key_var][0]
ice_class={0:'',1:'_miz',2:'_pack',3:'_mizf'}
key_i=ice_class[0]

#substract seasonal cycle
def an_seas(v):
    seas=np.nanmean(np.reshape(v,(12,timelen),order='F'),axis=1)
    #seas=[np.nanmean(v[n::12]) for n in range(12)]
   
    new_v=[v[n+12*k]-seas[n] for k in range(len(v)//12) for n in range(12)]
    return new_v

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
for key_s in sectors:
    if key_var==2 and end<=2021:
       
        cdr=pd.read_csv(d+'cdr'+'_'+var+'_'+hem+key_s+'_total.dat',delim_whitespace=True,header=0)[:]
        cdr=np.asarray(cdr) 

        cdr=len_fix(cdr[:,0],1979,2021)
        cdr=np.where(cdr==0,np.nan,cdr)
       # print(np.std(signal.detrend(cdr[8::12])))
        cdr=an_seas(cdr)
       # mask = ~np.isnan(cdr)
       
       # cdr=np.array(cdr)[mask]
        osi=pd.read_csv(d+'osisaf'+'_'+var+'_'+hem+key_s+key_i+'.dat',delim_whitespace=True,header=0)[:]
        print('osisaf'+'_'+var+'_'+hem+key_s+key_i+'.dat')
        osi=np.asarray(osi)
        osi=len_fix(osi[:,0],1979,2021)
        osi=np.where(osi==0,np.nan,osi)
        osi=an_seas(osi)
        osi=np.asarray(osi)
        osi=np.where(osi>2,np.nan,osi)
   
    if key_var in [1,3] and hem=='SH':
    
        gi=pd.read_csv(d+'giomas'+'_'+var+'_'+hem+key_s+'_total.dat',delim_whitespace=True,header=0)[:]
        gi=np.asarray(gi) 
        gi=gi[:,0]
    
        gi=len_fix(gi,1979,2020)
        
    if key_var in [1,3] and hem=='NH' and end<=2021:
        suff='' if key_s=='' else '_'
        pi=pd.read_csv(d+'piomas'+'_'+var+'_'+hem+suff+key_s+'.dat',delim_whitespace=True,header=0)[:]
        pi=np.asarray(pi)
        pi=pi[:,0]
        print(pi)
        pi=len_fix(pi,1979,2020)




   
  #  suff='' if key_s=='' else '_'
    suff=''
    outfile =  'trends_'+var+'_'+str(start)+'-'+str(end)+'_'+hem+suff+key_s+key_i+'.dat'
    outf = open((outfile),'w') # open dat file to write the trends
    #outf.write('TREND 10^3 km^3/year')
    outf.write('TREND . cm/year')
    outf.write('\n') 
    outf.write('product'+' '+'slope'+' '+'p_value' + ' ' + 'conf_lev_95')
    outf.write('\n') 

    fig = plt.figure(figsize=(12,6))
    ax = plt.gca()
    ax.grid(color='black', linestyle='dotted', axis='x',linewidth=0.8,alpha=0.2)
    ax.grid(color='black', linestyle='dotted', axis='y',linewidth=0.8,alpha=0.2)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    ax.set_xlim(0,timelen*12)
    ax.set_xticks(np.arange(0,timelen*12,24))
    ax.set_xticklabels(np.arange(start,end+1,2),rotation=300)
    ax.set_ylabel(variables[key_var][1],fontsize=14,labelpad=1)
    for key_pr in prod:
        suff='' if key_s=='' else '_'
        df1=pd.read_csv(d+prod[key_pr][1]+'_'+prod[key_pr][2]+'_'+var+'_'+hem+'_'+sims[0]+suff+key_s+key_i+'.dat',delim_whitespace=True,header=0)[:]
        print(prod[key_pr][1]+'_'+prod[key_pr][2]+'_'+var+'_'+hem+'_'+sims[1]+suff+key_s+key_i+'.dat')
        if os.path.isfile(d+prod[key_pr][1]+'_'+prod[key_pr][2]+'_'+var+'_'+hem+'_'+sims[1]+suff+key_s+key_i+'.dat'):
            print(prod[key_pr][1]+'_'+prod[key_pr][2]+'_'+var+'_'+hem+'_'+sims[1]+suff+key_s+key_i+'.dat')
            df2=pd.read_csv(d+prod[key_pr][1]+'_'+prod[key_pr][2]+'_'+var+'_'+hem+'_'+sims[1]+suff+key_s+key_i+'.dat',delim_whitespace=True,header=0)[:]
            data = df1._append(df2) 
        else:
            data=df1._append(pd.DataFrame(np.full([(2050-2014)*12, len(prod[key_pr][3])], np.nan),columns=prod[key_pr][3])) 

        data=np.asarray(data) 
        if len(data)>500:
         
            data=data[(start-1950)*12:,:] # !!!  
        data=data[:(end-start+1)*12,:]
        print(data)
        for n,k in enumerate(prod[key_pr][3]):
            slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(0,timelen,1), np.nanmean(np.reshape(data[:,n],(timelen,12)),axis=1))
            a=1000 if key_var in [2,3] and key_i!='_mizf' else 1 
            outf.write(prod[key_pr][1]+'_'+k+' '+f"{a*slope:.3f}"+ ' ' + f"{p_value:.5f}" + ' ' + f"{a*std_err*1.96:.3f}")
            outf.write('\n')    
            mask = ~np.isnan(data[:,n]) 
            ax.plot(np.arange(0,timelen*12,1),an_seas(data[:,n]),linewidth=1,color=prod[key_pr][6][n],linestyle=prod[key_pr][7][n],label=prod[key_pr][1]+' '+k)
            mask = ~np.isnan(cdr)
    if key_var in [1,3] : 
        if hem=='NH':   # northern hemisphere
           
            ax.plot(np.arange(0,timelen*12,1),an_seas(pi),linewidth=1.5,label='PIOMAS',color='black')
            slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(0,timelen,1), np.nanmean(np.reshape(1000*pi[:],(timelen,12)),axis=1))
            outf.write('piomas'+' '+f"{slope:.3f}"+ ' ' + f"{p_value:.3f}" + ' ' + f"{std_err*1.96:.3f}")
            outf.write('\n')
        
        else:
            ax.plot(np.arange(0,timelen*12,1),an_seas(gi),linewidth=1.5,label='GIOMAS',color='black',linestyle='solid')
            slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(0,timelen,1), np.nanmean(np.reshape(1000*gi[:],(timelen,12)),axis=1))
            outf.write('giomas'+' '+f"{slope:.3f}"+ ' ' + f"{p_value:.3f}" + ' ' + f"{std_err*1.96:.3f}")
            outf.write('\n')
            outf.write('grep'+' '+f"{1000*slope:.3f}"+ ' ' + f"{p_value:.3f}" + ' ' + f"{1000*std_err*1.96:.3f}")
            outf.write('\n')        
    if key_var==2 and end<=2021:  
        ax.plot(np.arange(0,timelen*12,1),(cdr),color='black',linewidth=1.5,label='CDR')
        ax.plot(np.arange(0,timelen*12,1),(osi),color='black',linewidth=1.5,label='OSISAF',linestyle='dashed')
        slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(0,timelen,1), np.nanmean(np.reshape(1000*np.asarray(cdr),(timelen,12)),axis=1))
        outf.write('cdr'+' '+f"{slope:.3f}"+ ' ' + f"{p_value:.3f}" + ' ' + f"{std_err*1.96:.3f}")
        outf.write('\n')
        slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(0,timelen,1), np.nanmean(np.reshape(1000*np.asarray(osi),(timelen,12)),axis=1))
        outf.write('osisaf'+' '+f"{slope:.3f}"+ ' ' + f"{p_value:.3f}" + ' ' + f"{std_err*1.96:.3f}")
        outf.write('\n')


      #  slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(0,timelen,1), np.nanmean(np.reshape(grep,(timelen,12)),axis=1))
       # outf.write('grep'+' '+f"{1000*slope:.3f}"+ ' ' + f"{p_value:.3f}" + ' ' + f"{1000*std_err*1.96:.3f}")
        #outf.write('\n')
        outf.close()

    ax.legend(fontsize=8,loc=1,ncol=8)

        #ax.text(5,75,hem,fontsize=12)
        #plt.savefig(prod[key_pr][1]+'_'+prod[key_pr][2]+'_'+var+'_'+hem+'_'+sim+'.png', format='png',  bbox_inches='tight', dpi=400)
    suff='' if key_s=='' else '_'
    plt.savefig(var+'_'+hem+'_an_'+str(start)[2:]+str(end)[2:]+key_i+suff+key_s+'.png', format='png',  bbox_inches='tight', dpi=400)



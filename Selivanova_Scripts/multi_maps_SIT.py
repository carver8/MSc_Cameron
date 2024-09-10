# multi maps
import numpy as np
import netCDF4 as nc
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy
import matplotlib.colors as clrs
import xarray as xr
import matplotlib.path as mpath
import cmocean as cmocean
from cartopy.util import add_cyclic_point
import glob
import matplotlib as mpl
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from scipy import stats
import numpy.ma as ma
from mpl_toolkits.basemap import addcyclic
#d='/gws/nopw/j04/primavera5/stream1/CMIP6/HighResMIP/'
## sivol cmcc not working
path='/work/scratch-nopw/is14120/'
path='nc/'
start=1979
end=2014
timelen=end-start+1
domain={0:['NH',ccrs.NorthPolarStereo(),[-180, 180, 48, 90]],1:['SH',ccrs.SouthPolarStereo(),[-180, 180, -90, -48]]}
key_dom=1
hem=domain[key_dom][0]
def hemi(lati,key_dom=key_dom):
    cond=lati<0 if key_dom==1 else lati>0 
    return cond
    
path='/home/users/is14120/nc/'
d1='/home/users/is14120/fig/'    
#for thickness
colors=[(255/256,189/256,171/256),(255/256,224/256,112/256),(135/256,191/256,0/256),(85/256,170/256,255/256),(0,119/256,196/256)]
cmap = mpl.colors.ListedColormap(colors) 
cmap_name = 'my_list'
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=256)

###  model dictionary
prod={0:['ECMWF','ECMWF','IFS',['LR','MR','HR'],['r1i1p1f1']*3,['v20180221','v20181119','v20170915']],\
        1:['MOHC','HadGEM3','GC31',['LL','MM','HM'],['r1i1p1f1']*3,['v20170921','v20170928','v20180730']], 
      2:['EC-Earth-Consortium','EC','Earth3P',['','HR',''],['r1i1p2f1']*3,['v20190314','v20181212']],
      3:['CNRM-CERFACS','CNRM','CM6-1',['','HR',''],['r1i1p1f2']*3,['v20190401','v20190221']],
       4:['AWI','AWI','CM-1-1',['LR','HR',''],['r1i1p1f2']*2],
      5:['CESM','CESM1','CAM5-SE',['LR','HR',''],['r1i1p1f1']*2],    
      6:['CMCC','CMCC','CM2',['HR4','VHR4'],['r1i1p1f1']*2,['v20200917','v20200917']],
      7:['MPI-M','MPI','ESM1-2',['HR','XR'],['r1i1p1f1']*2,['v20180606','v20180606']],
      8:['Envisat'],
      9:['GIOMAS','lat_scaler','lon_scaler']}
      #9:['PIOMAS','lat','lon'],
 #     10:['era','latitude','longitude'],11:['merra','latitude','longitude']}

### variables dictionary
variables={0:['siconc','RdYlBu_r',0,100,'SIC, (%)','mean_siconc'],
           1:['sithick','inferno',0.05,4,'SIT (m)','mean_sithick'],
           2:['sisnthick','RdYlBu_r',0,0.6,'snow depth (m)'],
           3:['sivol','inferno_r',0,3.2,'SIT (m)','mean_sivol'],
           4:['siconc','bwr',-3,3,'SIC trend (%/year)','trend_siconc'],
           5:['sithick','bwr',-0.6,0.6,'SIT trend (m/decade)'],
           6:['sisnthick','bwr',-0.6,0.6,'Snow depth trend (m/decade)','trend_sisnthick'],
           7:['tos','RdYlBu_r',-2,12,'SST','mean_tos'],
          # 8:['tos','bwr',-0.2,0.2,'SST (u'\u2103'/year)','trend_tos'],
           9:['tos','bwr',-0.02,0.02,'SSS (psu/year)','trend_sos'],
           10:['sos',cmocean.cm.haline,33,35.5,'SSS','mean_sos'],
          11:['ua','RdYlBu_r',-5,18, 'm/s','mean_ua'],
          12:['rlds',cmocean.cm.thermal,100,350,'W/$\mathregular{m^2}$','mean_rlds'],
        13:['tas','RdYlBu_r',-20,20,'T','mean_tos'],
           14:['clt','RdYlBu_r',0,100,'cloud fraction (%)','mean_tos']}
key_var=3

var=variables[key_var][0]
var_in=variables[key_var][5]
col='#373b96' if key_var in [0,2,7] else 'black' if key_var in [1,3] else 'white'
mons={0:'jan',1:'feb',2:'mar',3:'apr',4:'may',5:'jun',6:'jul',7:'aug',8:'sep',9:'oct',10:'nov',11:'dec'}
mon=8


def adjust_longitude(dataset: xr.Dataset) -> xr.Dataset:
        """Swaps longitude coordinates from range (0, 360) to (-180, 180)
        Args:
            dataset (xr.Dataset): xarray Dataset
        Returns:
            xr.Dataset: xarray Dataset with swapped longitude dimensions
        """
        lon_name = "longitude"  # whatever name is in the data

        # Adjust lon values to make sure they are within (-180, 180)
        dataset["_longitude_adjusted"] = xr.where(
            dataset[lon_name] > 180, dataset[lon_name] - 360, dataset[lon_name]
        )
        dataset = (
            dataset.swap_dims({lon_name: "_longitude_adjusted"})
            .sel(**{"_longitude_adjusted": sorted(dataset._longitude_adjusted)})
            .drop(lon_name)
        )

        dataset = dataset.rename({"_longitude_adjusted": lon_name})
        return dataset

cropped_dataset = adjust_longitude(cropped_dataset)

cropped_dataset = cropped_dataset.sel(
        latitude=slice(lat_max, lat_min), 
        longitude=slice(lon_min, lon_max)
    )

plot_dataset(cropped_dataset)


####func to plot contours       
def z_masked_overlap(axe, X, Y, Z, source_projection=None):
    """
    for data in projection axe.projection
    find and mask the overlaps (more 1/2 the axe.projection range)

    X, Y either the coordinates in axe.projection or longitudes latitudes
    Z the data
    operation one of 'pcorlor', 'pcolormesh', 'countour', 'countourf'

    if source_projection is a geodetic CRS data is in geodetic coordinates
    and should first be projected in axe.projection

    X, Y are 2D same dimension as Z for contour and contourf
    same dimension as Z or with an extra row and column for pcolor
    and pcolormesh

    return ptx, pty, Z
    """
    if not hasattr(axe, 'projection'):
        return Z
    if not isinstance(axe.projection, ccrs.Projection):
        return Z

    if len(X.shape) != 2 or len(Y.shape) != 2:
        return Z

    if (source_projection is not None and
            isinstance(source_projection, ccrs.Geodetic)):
        transformed_pts = axe.projection.transform_points(
            source_projection, X, Y)
        ptx, pty = transformed_pts[..., 0], transformed_pts[..., 1]
    else:
        ptx, pty = X, Y


    with np.errstate(invalid='ignore'):
        # diagonals have one less row and one less columns
        diagonal0_lengths = np.hypot(
            ptx[1:, 1:] - ptx[:-1, :-1],
            pty[1:, 1:] - pty[:-1, :-1]
        )
        diagonal1_lengths = np.hypot(
            ptx[1:, :-1] - ptx[:-1, 1:],
            pty[1:, :-1] - pty[:-1, 1:]
        )
        to_mask = (
            (diagonal0_lengths > (
                abs(axe.projection.x_limits[1]
                    - axe.projection.x_limits[0])) / 2) |
            np.isnan(diagonal0_lengths) |
            (diagonal1_lengths > (
                abs(axe.projection.x_limits[1]
                    - axe.projection.x_limits[0])) / 2) |
            np.isnan(diagonal1_lengths)
        )

        # TODO check if we need to do something about surrounding vertices

        # add one extra colum and row for contour and contourf
        if (to_mask.shape[0] == Z.shape[0] - 1 and
                to_mask.shape[1] == Z.shape[1] - 1):
            to_mask_extended = np.zeros(Z.shape, dtype=bool)
            to_mask_extended[:-1, :-1] = to_mask
            to_mask_extended[-1, :] = to_mask_extended[-2, :]
            to_mask_extended[:, -1] = to_mask_extended[:, -2]
            to_mask = to_mask_extended
        if np.any(to_mask):

            Z_mask = getattr(Z, 'mask', None)
            to_mask = to_mask if Z_mask is None else to_mask | Z_mask

            Z = ma.masked_where(to_mask, Z)

        return ptx, pty, Z
####



### envisat sit
env50=nc.Dataset(path+'SIT_SH_2002_2012_ENV_SnowAMSR.ease2_12500_smth50000.nc')
lat_l=env50.variables['latitude'][:]
lon_l=env50.variables['longitude'][:]  
env50_sit=env50.variables['SIT'][1:,:]
env50_sit=np.where(np.logical_and(env50_sit>=0.1,env50_sit<=20),env50_sit,np.nan)
env50_sit=np.nanmean(env50_sit[4::6,:],axis=0) ### there are only 6 months in the year - from may to october (4 is september)

### figure
fig, axs = plt.subplots(nrows=8,ncols=3,
                        subplot_kw={'projection':domain[key_dom][1]},
                        figsize=(10,28),
                        constrained_layout = True)
fig.set_facecolor("white")
#fig.subplots_adjust(hspace=0.1,wspace=0.1)

axs=axs.flatten()
nu=-1
#Loop over all of the models

for i,key_pr in enumerate(prod):
    if key_pr<8:
        for n, k in enumerate(prod[key_pr][3]):
            nu=nu+1
            if nu in [8,11,14,17]:
                axs[nu].axis('off')
          #  elif nu==20:
          #      exit
            else:
                nu=21 if (key_pr==7 and n==0) else nu    
                nu=22 if (key_pr==7 and n==1) else nu    
                print(i,key_pr, nu)
                axs[nu].set_title(prod[key_pr][1]+' '+prod[key_pr][3][n],fontsize=10) 
                suf='' if (key_pr==3 and n==0) or (key_pr==2 and n==0) else '-'    
        
                ### sit   
                data=nc.Dataset(glob.glob(path+'/'+var+'_'+prod[key_pr][1]+'-'+prod[key_pr][2]+suf+k+'_'+str(start)+'_'+str(end)+'_monclim.nc')[0])        
                all_sic_mean=data.variables[var][mon,:]  #!!!!
               # print(np.nanmax(all_sic_mean),key_pr)
        
                all_sic_mean=np.where((all_sic_mean>0.05)&(all_sic_mean<10),all_sic_mean, np.nan) ### mask 
                 
                   # if key_pr in [0,1,2,3,4,5]:
                if 'lat' in data.variables:
                    lat=data.variables['lat'][:]
                    lon=data.variables['lon'][:]  
                    
                else:    
                    lat=data.variables['latitude'][:]
                    lon=data.variables['longitude'][:] 
                if len(np.shape(lat))==1 :
                    lon,lat=np.meshgrid(lon,lat)        
           
                ### sic
                data1=nc.Dataset(glob.glob(path+'siconc_'+prod[key_pr][1]+'-'+prod[key_pr][2]+suf+k+'_'+str(start)+'_'+str(end)+'_monclim.nc')[0])   
                sic=data1.variables['siconc'][mon,:]
                sic=np.where((sic>0.01) & (sic<=100),sic, np.nan)
                if 'lat' in data1.variables:
                    lat1=data1.variables['lat'][:]
                    lon1=data1.variables['lon'][:]  
                    
                else:    
                    lat1=data1.variables['latitude'][:]
                    lon1=data1.variables['longitude'][:] 
                if len(np.shape(lat1))==1 :
                    lon1,lat1=np.meshgrid(lon1,lat1)                       
                        
                        
              
              #all_sic_mean=all_sic_mean*100 if key_pr==5 and key_var in [0,4] else all_sic_mean
                #print(all_sic_mean[all_sic_mean<2])
           
            
                axs[nu].set_extent(domain[key_dom][2], ccrs.PlateCarree())
                axs[nu].coastlines(color='k', linewidth=0.5, zorder=3)
                axs[nu].add_feature(cartopy.feature.LAND, zorder=3, color='#b1b8b8')
                axs[nu].add_feature(cartopy.feature.OCEAN, color='white')
                theta = np.linspace(0, 2*np.pi, 50)
                center, radius = [0.5, 0.5], 0.5
                verts = np.vstack([np.sin(theta), np.cos(theta)]).T
                circle = mpath.Path(verts * radius + center)
                axs[nu].set_boundary(circle, transform=axs[nu].transAxes)
                gl=axs[nu].gridlines(linewidth=0.4,alpha=0.5, linestyle='--', draw_labels=False)
            
                    
                    
                   # p=axs[nu].pcolormesh(lon,lat,all_sic_mean,cmap=variables[key_var][1], transform=ccrs.PlateCarree(),shading='flat',vmin=variables[key_var][2],vmax=variables[key_var][3])#,ax=ax)
                all_sic_mean, lon = addcyclic(all_sic_mean, lon)
                lat = addcyclic(lat)
                sic=sic*100 if key_pr==7 else sic #!!!
        
        
                x,y, masked = z_masked_overlap(axs[nu], lon,lat,all_sic_mean,source_projection=ccrs.Geodetic())   
                p=axs[nu].contourf(x,y,masked,50,cmap=variables[key_var][1],vmin=variables[key_var][2],vmax=variables[key_var][3])#,ax=ax)
                sic1=np.where(sic<5,sic,np.nan)
                x,y, masked = z_masked_overlap(axs[nu], lon1,lat1,sic1,source_projection=ccrs.Geodetic())   
                axs[nu].contourf(x,y,masked,level=[15],cmap='binary',vmin=0,vmax=1000)#,ax=ax)
        
        
                    
              #  levels=np.arange(variables[key_var][2],variables[key_var][3]+0.1,0.2)
                x,y, masked_sic = z_masked_overlap(axs[nu], lon1, lat1, sic,source_projection=ccrs.Geodetic())   
                axs[nu].contour(x,y,masked_sic,colors='black',linestyles='solid',levels=[15],linewidths=0.9)                     
                axs[nu].contour(x,y,masked_sic,colors='black',linestyles='dashed',levels=[80],linewidths=0.9)
     #   nu=nu+1
    #elif key_pr in [8]:
 
        #x,y, masked_sic = z_masked_overlap(axs[nu], lon1, lat1, sic,source_projection=ccrs.Geodetic())    
        #axs[nu].contour(x,y,masked_sic,colors='black',linestyles='solid',levels=[15],linewidths=0.9)
        #axs[nu].contour(x,y,masked_sic,colors='black',linestyles='dashed',levels=[80],linewidths=0.9)        

         
                                   
    elif key_pr in [9] and key_var in [3]:
        nu=23
        data_gi=nc.Dataset('/home/users/is14120/'+'giomas/heff_'+str(start)+'_'+str(end)+'_monclim.nc')
        lat=data_gi.variables[prod[key_pr][1]][:]
        lon=data_gi.variables[prod[key_pr][2]][:]  
        all_sic_mean=data_gi.variables['heff'][mon,:,:]
        print(np.shape(all_sic_mean))
        all_sic_mean=np.where(np.logical_and(all_sic_mean>0,all_sic_mean<10),all_sic_mean,np.nan)
        
        
        
        ### plot giomas        
        axs[nu].set_title('GIOMAS',fontsize=8)
        axs[nu].set_extent(domain[key_dom][2], ccrs.PlateCarree())
        axs[nu].coastlines(color='k', linewidth=0.5, zorder=3)
        axs[nu].add_feature(cartopy.feature.LAND, zorder=3, color='#b1b8b8')
        axs[nu].add_feature(cartopy.feature.OCEAN, color='white')   
        theta = np.linspace(0, 2*np.pi, 50)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        axs[nu].set_boundary(circle, transform=axs[nu].transAxes)
        gl=axs[nu].gridlines(linewidth=0.4,alpha=0.5, linestyle='--', draw_labels=False)
        all_sic_mean, lon = addcyclic(all_sic_mean, lon)
        lat = addcyclic(lat)
        x,y, masked = z_masked_overlap(axs[nu], lon, lat, all_sic_mean,source_projection=ccrs.Geodetic())   
        masked = np.ma.masked_where(np.isnan(masked), masked)
        p=axs[nu].contourf(x,y,masked,levels=np.arange(variables[key_var][2],variables[key_var][3]+0.1,0.2),cmap=variables[key_var][1],vmin=variables[key_var][2],vmax=variables[key_var][3],extend='max')#,ax=ax)
        data1=nc.Dataset('/home/users/is14120/cdr/seaice_conc_mon_sh_'+str(start)+'_'+str(end)+'_monclim.nc')  
            #data=nc.Dataset('work/'+prod[key_pr].lower()+'_map_trend_siconc_sep_1979_2014_NH.nc')      
        lat1=data1.variables['latitude'][:]
        lon1=data1.variables['longitude'][:] 
        
        sic=data1.variables['cdr_seaice_conc'][mon,:] 
        
        
        sic1=np.where(sic<5,sic,np.nan)
        x,y, masked = z_masked_overlap(axs[nu], lon1,lat1,sic1,source_projection=ccrs.Geodetic())   
        axs[nu].contourf(x,y,masked,level=[15],cmap='binary',vmin=0,vmax=1000)#,ax=ax)
        x,y, masked_sic = z_masked_overlap(axs[nu], lon1, lat1, sic,source_projection=ccrs.Geodetic())    
        axs[nu].contour(x,y,masked_sic,colors='black',linestyles='solid',levels=[15],linewidths=0.9)
        axs[nu].contour(x,y,masked_sic,colors='black',linestyles='dashed',levels=[80],linewidths=0.9)

    elif key_pr in [8] and mon==8:
        nu=20
        axs[nu].set_title('Envisat',fontsize=8)
        axs[nu].set_extent(domain[key_dom][2], ccrs.PlateCarree())
        axs[nu].coastlines(color='k', linewidth=0.5, zorder=3)
        axs[nu].add_feature(cartopy.feature.LAND, zorder=3, color='#b1b8b8')
        axs[nu].add_feature(cartopy.feature.OCEAN, color='white')   
        theta = np.linspace(0, 2*np.pi, 50)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        axs[nu].set_boundary(circle, transform=axs[nu].transAxes)
        gl=axs[nu].gridlines(linewidth=0.4,alpha=0.5, linestyle='--', draw_labels=False)
        all_sic_mean, lon = addcyclic(all_sic_mean, lon)
        lat = addcyclic(lat)
        x,y, masked = z_masked_overlap(axs[nu], lon_l, lat_l, env50_sit,source_projection=ccrs.Geodetic())   
        masked = np.ma.masked_where(np.isnan(masked), masked)
        p=axs[nu].contourf(x,y,masked,levels=np.arange(variables[key_var][2],variables[key_var][3]+0.1,0.2),cmap=variables[key_var][1],vmin=variables[key_var][2],vmax=variables[key_var][3],extend='max')#,ax=ax)
        data1=nc.Dataset('/home/users/is14120/cdr/seaice_conc_mon_sh_'+str(start)+'_'+str(end)+'_monclim.nc')  
                    #data=nc.Dataset('work/'+prod[key_pr].lower()+'_map_trend_siconc_sep_1979_2014_NH.nc')      
        lat1=data1.variables['latitude'][:]
        lon1=data1.variables['longitude'][:] 
                
        sic=data1.variables['cdr_seaice_conc'][mon,:] 
                
                
        sic1=np.where(sic<5,sic,np.nan)
        x,y, masked = z_masked_overlap(axs[nu], lon1,lat1,sic1,source_projection=ccrs.Geodetic())   
        axs[nu].contourf(x,y,masked,level=[15],cmap='binary',vmin=0,vmax=1000)#,ax=ax)
        data1=nc.Dataset('/home/users/is14120/cdr/seaice_conc_mon_sh_'+str(start)+'_'+str(end)+'_monclim.nc')  
            #data=nc.Dataset('work/'+prod[key_pr].lower()+'_map_trend_siconc_sep_1979_2014_NH.nc')      
        lat1=data1.variables['latitude'][:]
        lon1=data1.variables['longitude'][:] 
        
        sic=data1.variables['cdr_seaice_conc'][mon,:] 
        
        
        sic1=np.where(sic<5,sic,np.nan)
        x,y, masked = z_masked_overlap(axs[nu], lon1,lat1,sic1,source_projection=ccrs.Geodetic())   
        axs[nu].contourf(x,y,masked,level=[15],cmap='binary',vmin=0,vmax=1000)#,ax=ax)
        x,y, masked_sic = z_masked_overlap(axs[nu], lon1, lat1, sic,source_projection=ccrs.Geodetic())    
        axs[nu].contour(x,y,masked_sic,colors='black',linestyles='solid',levels=[15],linewidths=0.9)
        axs[nu].contour(x,y,masked_sic,colors='black',linestyles='dashed',levels=[80],linewidths=0.9)
        



           # p=axs[nu].pcolormesh(lon,lat,all_sic_mean,cmap=variables[key_var][1], transform=ccrs.PlateCarree(),shading='gouraud',vmin=variables[key_var][2],vmax=variables[key_var][3])#,ax=ax) 
          #  nu=nu+1
            
            
            
            ### contours


      #  axs[nu].contour(lon1,lat1,sic,levels=[15],color='white',transform=ccrs.PlateCarree())
       # axs[nu].contour(lon1,lat1,sic,levels=[80],color='white',transform=ccrs.PlateCarree())'''
             
if mon==1:
    axs[20].axis('off')    
        
        
#for nu in [14]:
 #   fig.delaxes(axs[nu])
#fig.subplots_adjust(bottom=0.3, tp=0.9, left=0.1, right=0.9,
     #               wspace=0.02, hspace=0.3)
fig.subplots_adjust(wspace=0.2, hspace=0.15)
# Add a colorbar axis at the bottom of the graph
#cbar_ax = fig.add_axes([0.26, 0.06, 0.5, 0.005])  #[0.26, 0.24, 0.5, 0.02]
# Draw the colorbar
cbar_ax = fig.add_axes([0.774, 0.505, 0.03, 0.18]) 

cbar=fig.colorbar(p, cax=cbar_ax,orientation='vertical',extend='both')
cbar.set_label(variables[key_var][4],fontsize=12,rotation=270,labelpad=15)

plt.savefig(d1+variables[key_var][0]+'_'+domain[key_dom][0]+'_map_'+mons[mon]+'_'+str(start)+'_'+str(end)+'_cesm_.png', format='png',  bbox_inches='tight', dpi=500)

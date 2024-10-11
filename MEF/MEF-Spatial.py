#### REPROJECTION OF MODEL DATA ONTO NSIDC GRID ####

# %% Import Relevant Libraries
import os; os.chdir('/Users/crcarver/Desktop/AOS_THESIS/MSc_Cameron/Data/')
import xarray as xr
import calendar
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.basemap import Basemap
# %%% Import Observation Data
access_obs_data='OBSERVATIONS/SICONC/seaice_conc_monthly_sh_n07_v04r00_all.nc' #access combined monthly observational mean data from NSIDC
obs_ds = xr.open_dataset(access_obs_data)                   # Open dataset
obs_ds = obs_ds.sel(tdim=slice('1979-01-01','2009-12-31'))  # Time slice of data
obs_da = obs_ds.cdr_seaice_conc_monthly                     # Select desired variable
obs_da = obs_da.where(obs_da <= 2.51, other=float('nan'))   # Set all flagged values to nan
# %%% Import Model Data 
# access_model_data = 'MODELS/SICONC/siconc_SImon_CMCC-CM2-SR5_omip1_r1i1p1f1_gn_163801-200912_reproj.nc' #CMCC-SR5-OMIP1
# access_model_data = 'MODELS/SICONC/siconc_SImon_CMCC-CM2-SR5_omip2_r1i1p1f1_gn_165301-201812_reproj.nc' #CMCC-SR5-OMIP2
access_model_data = 'MODELS/SICONC/siconc_SImon_MRI-ESM2-0_omip1_r1i1p1f1_gn_170001-200912_reproj.nc' #MRI-ESM2-OMIP1
# access_model_data = 'MODELS/SICONC/siconc_SImon_GFDL-CM4_omip1_r1i1p1f1_gn_190801-200712_reproj.nc' # GFDL-CM4-OMIP1 ##!! MUST ADJUST TIME TO END 2007

model_ds = xr.open_dataset(access_model_data)                   # Open dataset
model_ds = model_ds.sel(time=slice('1979-01-01','2009-12-31'))  #Time slice of data

model_ds = model_ds.assign_coords(time=obs_ds['time'])  # Set time coords of model ds to be the same as the obs ds
model_ds = model_ds.swap_dims({'tdim': 'time'})
model_ds = model_ds.swap_dims({'time': 'tdim'})

model_da = model_ds.siconc / 100    # Select desired variable and scale to match obs scale of 0-1
model_name, exp_id = model_ds.model_name, model_ds.exp_id   # Call model name from attributes
# %% Plotting Details
# Import reference grid from NSIDC anciliarry file
access_ref_data = 'OBSERVATIONS/SICONC/G02202-cdr-ancillary-sh.nc'
ref_ds = xr.open_dataset(access_ref_data)
# maskv = ref_ds.valid_ice_mask.values
levels=[ 0, 0.2, 0.5, 0.65, 1.0]
colors = ['peru', 'yellowgreen', 'deepskyblue', 'royalblue']
cmap = mcolors.ListedColormap(colors)
cmap.set_under('darkred')

m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')   # Create a Basemap instance with South Polar Stereographic projection
x, y = m(ref_ds.longitude.values, ref_ds.latitude.values)       # Convert latitude and longitude to map projection coordinates

# %% MEF - Cummulative Mean
MEF_num = (obs_da-model_da) ** 2       # cell-wise RMSD
MEF_num = MEF_num.sum(dim=['tdim'])    # Temperol sum of cell-wise RMSD

o_mean = obs_da.mean(dim='tdim')       # Temporal mean of observational product
MEF_den = (obs_da-o_mean) ** 2         # cell-wise SD
MEF_den = MEF_den.sum(dim='tdim')      # Temporal sum of cell-wise SD

MEF = (1 - (MEF_num/MEF_den))
# %%% Plot Cummulative Mean
plt.figure(figsize=(8, 10))
m.contourf(x, y, MEF.values, cmap=cmap, levels=levels, extend='min')
m.drawcoastlines()
m.drawparallels(np.arange(-90., 0., 10.), labels=[1,0,0,0])
m.drawmeridians(np.arange(-180., 181., 30.), labels=[0,0,0,1])
plt.colorbar(label='MEF', orientation='horizontal', pad=0.05,aspect=40)
plt.title(f'{model_name}-{exp_id} Annual Mean MEF ', fontsize=18)
# %% New Metrics - Cummulative Monthly WORKING!!
# MEF T-S-T
MEF_num = (obs_da-model_da) ** 2        # cell-wise RMSD
MEF_num = MEF_num.groupby('time.month').sum(dim=['tdim'])    # Temperol mean of cell-wise RMSD

o_mean = obs_da.mean(dim='tdim')        # Temporal mean of observational product
MEF_den = (obs_da-o_mean) ** 2          # cell-wise SD
MEF_den = MEF_den.groupby('time.month').sum(dim='tdim')      # Temporal mean of cell-wise SD

MEF = (1 - (MEF_num/MEF_den))
vmask = xr.where(MEF == 1, 0, 1)    # Mask out values of 1
MEF_m = MEF*vmask
# %%% Plot Monthly Mean - 4x3
months = list(calendar.month_abbr)[1:]
fig, axes = plt.subplots(4, 3, figsize=(12, 16))

for i in range(12):
    ax = axes[i // 3, i % 3]
    m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l', ax=ax)
    x, y = m(ref_ds.longitude.values, ref_ds.latitude.values)
    contour = m.contourf(x, y, MEF_m[i,:,:].values, levels=levels, cmap=cmap, extend='min')
    m.drawcoastlines()
    m.drawparallels(np.arange(-90., 0., 10.), labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180., 181., 30.), labels=[0,0,0,1])
    ax.set_title(f'{months[i]}')

fig.suptitle(f'{model_name}-{exp_id} Monthly Mean MEF', fontsize=24, y=0.92)
cbar = fig.colorbar(contour, ax=axes.ravel().tolist(), pad=0.05, orientation='horizontal', aspect=50)
cbar.set_label('MEF')
# %%%% Plot Individual Monthly Mean
# months = list(calendar.month_abbr)[1:]
# for i in range(0,12,1):
#     # plt.figure(figsize=(8, 8))
#     m.contourf(x, y, MEF_m[i,:,:].values, cmap=cmap, levels=levels, extend='min')
#     m.drawcoastlines()
#     m.drawparallels(np.arange(-90., 0., 10.), labels=[1,0,0,0])
#     m.drawmeridians(np.arange(-180., 181., 30.), labels=[0,0,0,1])
#     plt.colorbar(label='MEF')
#     plt.title(f'MEF for Month {months[i]}')

# %% MEF - Cummulative Seasonal
MEF_num = (obs_da-model_da) ** 2        # cell-wise RMSD
MEF_num = MEF_num.groupby('time.season').sum(dim=['tdim'])    # Temperol mean of cell-wise RMSD

o_mean = obs_da.mean(dim='tdim')        # Temporal mean of observational product
MEF_den = (obs_da-o_mean) ** 2          # cell-wise SD
MEF_den = MEF_den.groupby('time.season').sum(dim='tdim')      # Temporal mean of cell-wise SD

MEF = (1 - (MEF_num/MEF_den))
vmask = xr.where(MEF == 1, 0, 1)    # Mask out values of 1
MEF_s = MEF*vmask
# %%% Plot Seasonal Mean - 2x2
season = MEF_s['season'].values
fig, axes = plt.subplots(2, 2, figsize=(10, 12))

for i in range(4):
    ax = axes[i // 2, i % 2]
    m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l', ax=ax)
    x, y = m(ref_ds.longitude.values, ref_ds.latitude.values)
    contour = m.contourf(x, y, MEF_s[i,:,:].values, levels=levels, cmap=cmap, extend='min')
    m.drawcoastlines()
    m.drawparallels(np.arange(-90., 0., 10.), labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180., 181., 30.), labels=[0,0,0,1])
    ax.set_title(f'{season[i]}')

fig.suptitle(f'{model_name}-{exp_id} Seasonal Mean MEF', fontsize=24, y=0.94)
cbar = fig.colorbar(contour, ax=axes.ravel().tolist(), pad=0.05, orientation='horizontal', aspect=50)
cbar.set_label('MEF')
# %%%% Plot Individual Seasonal Mean
# season = MEF_s['season'].values
# for i in range(4):
#     plt.figure(figsize=(8, 8))
#     m.contourf(x, y, MEF_s[i,:,:].values, cmap=cmap, levels=levels, extend='min')
#     m.drawcoastlines()
#     m.drawparallels(np.arange(-90., 0., 10.), labels=[1,0,0,0])
#     m.drawmeridians(np.arange(-180., 181., 30.), labels=[0,0,0,1])
#     plt.colorbar(label='MEF')
#     plt.title(f'MEF for Season {season[i]}')

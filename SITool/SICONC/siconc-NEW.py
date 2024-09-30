#### SITool Sea Ice Concentration HIGH RESOLUTION COMPATIBLE ####
# Cameron Carver Sept - 2024
# University of Cape Town - Applied Ocean Science MSc

# %% Import relevant packages and data
# %%% Import packages
import os
os.chdir('/Users/crcarver/Desktop/AOS_THESIS/MSc_Cameron/Data/') # Set base directory
import xarray as xr

# %%% Import Model Data 
access_model_data = 'MODELS/SICONC/siconc_SImon_CMCC-CM2-SR5_omip1_r1i1p1f1_gn_163801-200912_reproj.nc' #access reprojected model sea ice dataset
model_ds = xr.open_dataset(access_model_data) # Open dataset
model_ds = model_ds.sel(time=slice('1979-01-01','2009-12-31')) #Time slice of data
# print(model_ds)

# %%% Import Observation Data
access_obs_data='OBSERVATIONS/SICONC/seaice_conc_monthly_sh_n07_v04r00_all.nc' #access combined monthly observational mean data from NSIDC
obs_ds = xr.open_dataset(access_obs_data) # Open dataset
obs_ds = obs_ds.sel(tdim=slice('1979-01-01','2009-12-31')) #Time slice of data
# print(obs_ds.attrs)

# %% Sea Ice Concentration Metrics

# INPUT model_ds & obs_ds with 'siconc" variable
# OUTPUT error_mean. error_std, & error_trend

# %%% Monthly Mean Cycle
    #1. Evaluation of the mean seasonal cycle
print('             ')
print('MEAN CYCLE')
print('=============')
## Monthly Mean
    # Model
model_mean_mon = model_ds.siconc.groupby('time.month').mean('time')  
conc1 = model_mean_mon   

    # Observational
obs_mean_mon = obs_ds.siconc.groupby('time.month').mean('tdim')
conc2 = obs_mean_mon

## Error of Mean Cycle
error_mean_conc = abs(conc1-conc2) ## Mean Absolute Error ##

## Mask
    # Create Mask for NaN and for when both datasets == 0
mask = xr.where(((conc1 == 0) & (conc2 == 0)), 0, 1) # Where both datasets = 0
mask = mask.where(conc1.notnull() | conc2.notnull(), 0, 1) # Where model OR obs is NaN

## Global Error
    # Calculate the global error by hemisphere
error_mean_monthly = ((error_mean_conc*mask).sum(dim=["x","y"]))/((mask).sum(dim=["x","y"])) 

ndpm=[31,28.25,31,30,31,30,31,31,30,31,30,31]; #Number of days per month
error_mean = ((error_mean_monthly*ndpm).sum())/365.25
print(error_mean.values)

# %%% Anomaly Variance
    #2.Evaluate the variance of anomalies. The anomalies are defined as the signal minus the mean seasonal cycle
print('                ')
print('ANOMALY VARIANCE')
print('================')   
    
## Anomolies
monthly_mean_aligned = model_ds.siconc.groupby('time.month') - model_mean_mon
ano_conc1 = model_ds.siconc.groupby('time') - monthly_mean_aligned

monthly_mean_aligned = obs_ds.siconc.groupby('time.month') - obs_mean_mon
ano_conc2 = obs_ds.siconc.groupby('tdim') - monthly_mean_aligned

## Std Dev of Anomolies
std_ano_conc1 = ano_conc1.std(dim=["time"])
std_ano_conc2 = ano_conc2.std(dim=["tdim"])

## Errors
error_std_conc = abs(std_ano_conc1-std_ano_conc2) ## Mean absolute anomoly deviation ##

## Mask
# Create Mask for NaN and for when both std devs == 0
mask_std = xr.where(((std_ano_conc1 == 0) & (std_ano_conc2 == 0)), 0, 1) # Where both datasets = 0
mask_std = mask_std.where(std_ano_conc1.notnull() | std_ano_conc2.notnull(), 0, 1) # Where model OR obs is NaN

## Monthly variabiliy ##

## Global Error
# Compute global error on std
error_std = ((error_std_conc*mask_std).sum(dim=["x","y"]))/((mask_std).sum(dim=["x","y"])) 
print(error_std.values)

# %%% Trend
    #3.Evaluate the skill in each cell to reproduce the trend 
print('     ')
print('TREND')
print('=====')   
## Trend of each cell
trend_ano_conc1 = model_ds.siconc.polyfit(dim='time', deg=1)
trend_ano_conc2 = obs_ds.siconc.polyfit(dim='tdim', deg=1)

## Error of Trend
error_trend_conc=12*10*abs(trend_ano_conc1-trend_ano_conc2)#/decade

## Mask
mask_trend_conc = xr.where(((trend_ano_conc1 == 0) & (trend_ano_conc2 == 0)), 0, 1) # Where both datasets = 0
mask_trend_conc = mask_trend_conc.where(trend_ano_conc1.notnull() | trend_ano_conc2.notnull(), 0, 1) # Where model OR obs is NaN

## Global Error
error_trend = ((error_trend_conc*mask_trend_conc).sum(dim=["x","y"]))/((mask_trend_conc).sum(dim=["x","y"])) 
print(error_trend['polyfit_coefficients'].values)

# %% Cluster Dendrogram

# Assuming the datasets have the same variables and dimensions

model_mean_mon = model_ds.siconc.resample(time='ME').mean()  
obs_mean_mon = obs_ds.siconc.resample(time='ME').mean()

import cftime

# Convert da1 to Gregorian calendar
model_mean_mon['time'] = model_mean_mon.indexes['time'].to_datetime()

# Convert da2 to Gregorian calendar
obs_mean_mon['time'] = obs_mean_mon.indexes['time'].to_datetimeindex()


model_mean_mon

combined = xr.concat([model_mean_mon, obs_mean_mon], dim='dataset')
import numpy as np

# Flatten the data
data1 = conc1.values.tolist()
data2 = conc2.values.tolist()

data = [data1, data2]

from scipy.cluster.hierarchy import dendrogram, linkage

# Perform hierarchical/agglomerative clustering
Z = linkage(data, 'ward')


# %%
import matplotlib.pyplot as plt
# Plot the dendrogram
plt.figure(figsize=(10, 7))
dendrogram(Z)
plt.title('Cluster Dendrogram')
plt.xlabel('Sample index')
plt.ylabel('Distance')
plt.show()

# %% Visualisation
import calendar
import matplotlib.pyplot as plt

months = list(calendar.month_abbr)[1:]

for i in range(0,12,1):
    diff = (error_mean_conc[i,:,:])
    diff.plot(cmap="magma", vmin=0, vmax=85)
    plt.title(f'Model to Obs Difference for {months[i]}')
    plt.xlabel("x (meters)"); plt.ylabel("y (meters)");
    plt.show()

# for i in range(0,12,1):
#     error_mean.plot(cmap="Blues_r")
#     plt.title(f'Mean Err for {months[i]}')
#     plt.xlabel("x (meters)"); plt.ylabel("y (meters)");
#     plt.show()
# print(diff)
# import xmovie
# import ffmpeg
# print(error_mean_conc)
# movie = xmovie.Movie(error_mean_conc, framedim="month", vmin=0, vmax=85, cmap='magma')
# movie.save('mean_siconc.gif', overwrite_existing=True, gif_framerate=12)
 


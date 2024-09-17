#### SITool Sea Ice Concentration HIGH RESOLUTION COMPATIBLE ####
# Cameron Carver Sept - 2024
# University of Cape Town - Applied Ocean Science MSc

# %% Import relevant packages and data
# %%% Import packages
import os
os.chdir('/Users/crcarver/Desktop/AOS_THESIS/MSc_Cameron/Data/') # Set base directory
import xarray as xr

# %%% Import Model Data 
access_model_data = 'siconc_SImon_CMCC-CM2-SR5_omip1_r1i1p1f1_gn_163801-200912_reproj.nc' #access reprojected model sea ice dataset
model_ds = xr.open_dataset(access_model_data) # Open dataset
model_ds = model_ds.sel(time=slice('1979-01-01','2009-12-31')) #Time slice of data
# print(model_ds)

# %%% Import Observation Data
access_obs_data='seaice_conc_monthly_sh_n07_v04r00_all.nc' #access combined monthly observational mean data from NSIDC
obs_ds = xr.open_dataset(access_obs_data) # Open dataset
obs_ds = obs_ds.sel(tdim=slice('1979-01-01','2009-12-31')) #Time slice of data
# print(obs_ds)

# %% Sea Ice Concentration Metrics

# INPUT model_ds & obs_ds with 'siconc" variable
# OUTPUT error_mean. error_std, & error_trend

# %%% Monthly Mean Cycle
    #1. Evaluation of the mean seasonal cycle
print('             ')
print('1. MEAN CYCLE')
print('=============')
## Monthly Mean
    # Model
model_mean_mon = model_ds.siconc.groupby('time.month').mean('time')  
conc1 = model_mean_mon   
# print(model_mean_mon)

    # Observational
obs_mean_mon = obs_ds.siconc.groupby('time.month').mean('tdim')
conc2 = obs_mean_mon
# print(obs_mean_mon)

## Error of Mean Cycle
error_mean_conc = abs(conc1-conc2)

## Mask
    # Create Mask for NaN and for when both datasets == 0
mask = xr.where(((conc1 == 0) & (conc2 == 0)), 0, 1) # Where both datasets = 0
mask = mask.where(conc1.notnull(), 0, 1) #Where model is NaN
mask = mask.where(conc2.notnull(), 0, 1) #Where obs is NaN
mask_ant = xr.where(conc2 >= 99, 0, 1) # Mask out Antarctica
mask = mask*mask_ant
# print(mask)
# mask.sum(dim=["x","y"])

## Global Error
    # Calculate the global error by hemisphere
error_mean_monthly = ((error_mean_conc*mask).sum(dim=["x","y"]))/((mask).sum(dim=["x","y"])) 
        #Removed cell area as it is a constant in both num and denom
# print(error_mean_monthly)

ndpm=[31,28,31,30,31,30,31,31,30,31,30,31]; #Number of days per month
error_mean = ((error_mean_monthly*ndpm).sum())/365
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
error_std_conc  =abs(std_ano_conc1-std_ano_conc2)

## Mask
# Create Mask for NaN and for when both std devs == 0
mask_std = xr.where(((std_ano_conc1 == 0) & (std_ano_conc2 == 0)), 0, 1) # Where both datasets = 0
mask_std = mask_std.where(std_ano_conc1.notnull(), 0, 1)
mask_std = mask_std.where(std_ano_conc2.notnull(), 0, 1)
mask_ant = xr.where(conc2[1,:,:] >= 99, 0, 1) # Mask out Antarctica
mask_std = mask_std*mask_ant
# print(mask_std)

## Global Error
# Compute global error on std
error_std = ((error_std_conc*mask_std).sum(dim=["x","y"]))/((mask_std).sum(dim=["x","y"])) 
        #Removed cell area as it is a constant in both num and denom
print(error_std.values)

# %%% Trend
    #3.Evaluate the skill in each cell to reproduce the trend 
print('     ')
print('TREND')
print('=====')   
## Trend of each cell
trend_ano_conc1 = model_ds.siconc.polyfit(dim='time', deg=1)
trend_ano_conc2 = obs_ds.siconc.polyfit(dim='tdim', deg=1)

# trend_coeff_conc1 = trend_ano_conc1['polyfit_coefficients']
# trend_coeff_conc2 = trend_ano_conc2['polyfit_coefficients']

## Error of Trend
error_trend_conc=12*10*abs(trend_ano_conc1-trend_ano_conc2)#/decade
# print(error_trend_conc)

## Mask
mask_trend_conc = xr.where(((trend_ano_conc1 == 0) & (trend_ano_conc2 == 0)), 0, 1) # Where both datasets = 0
mask_trend_conc = mask_trend_conc.where(trend_ano_conc1.notnull(), 0, 1)
mask_trend_conc = mask_trend_conc.where(trend_ano_conc2.notnull(), 0, 1)
mask_ant = xr.where(conc2[1,:,:] >= 99, 0, 1) # Mask out Antarctica
mask_trend_conc = mask_trend_conc*mask_ant
# print(mask_trend_conc)
# mask_trend_conc.dims

## Global Error
error_trend = ((error_trend_conc*mask_trend_conc).sum(dim=["x","y"]))/((mask_trend_conc).sum(dim=["x","y"])) 
print(error_trend['polyfit_coefficients'].values)

# %% Heatmap function  |
import matplotlib.pyplot as plt

mask_ant = xr.where(obs_ds.siconc[1,:,:] >= 99, float('nan'), 1)
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
for i in range(0,12,1):
    diff = (model_mean_mon[i,:,:] - obs_mean_mon[i,:,:])*mask_ant
    diff.plot(cmap="coolwarm")
    plt.title(f'Model to Obs Difference for {months[i]}')
    plt.xlabel("x (meters)"); plt.ylabel("y (meters)");
    plt.show()

# for i in range(0,12,1):
#     err = error_mean_conc[i,:,:]*mask_ant
#     err.plot(cmap="OrRd")
#     plt.title(f'Model to Obs Difference for {months[i]}')
#     plt.xlabel("x (meters)"); plt.ylabel("y (meters)");
#     plt.show()





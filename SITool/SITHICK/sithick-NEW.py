#### SITool Sea Ice Thickness HIGH RESOLUTION COMPATIBLE ####
# Cameron Carver Sept - 2024
# University of Cape Town - Applied Ocean Science MSc
# 2) A function that computes sea ice thickness errors between two datasets
#    on the mean state in order to get the metrics;
#    The input data is a numpy array but a NetCDF file 
# 6) A script calls the function 2) and then computes the ice thickness metrics



# %% Import relevant packages and data
# %%% Import packages
import os
os.chdir('/Users/crcarver/Desktop/AOS_THESIS/MSc_Cameron/Data/') # Set base directory
import xarray as xr

# %%% Import Model Data 
access_model_data = 'sithick_SImon_CMCC-CM2-SR5_omip1_r1i1p1f1_gn_163801-200912_reproj.nc' #access reprojected model sea ice dataset
model_ds = xr.open_dataset(access_model_data) # Open dataset
model_ds = model_ds.sel(time=slice('2003-01-01','2007-12-31')) #Time slice of data
print(model_ds)

# %%% Import Observation Data
access_obs_data='SITHICK/SIT_sh_200205_201203_env.nc' #access combined monthly observational mean data from NSIDC
obs_ds = xr.open_dataset(access_obs_data) # Open dataset
obs_ds = obs_ds.sel(time=slice('2003-01-01','2007-12-31')) #Time slice of data
obs_ds = obs_ds.rename({'sea_ice_thickness':'sithick'})
obs_ds = obs_ds['sithick'].to_dataset()
print(obs_ds)

# %% Sea Ice Concentration Metrics

# INPUT model_ds & obs_ds with 'sithick" variable
# OUTPUT error_mean
# %%% Monthly Mean Cycle
    #1. Evaluation of the mean seasonal cycle
print('             ')
print('MEAN CYCLE')
print('=============')
## Monthly Mean
    # Model
model_mean_mon = model_ds.sithick.groupby('time.month').mean('time')  
conc1 = model_mean_mon   

    # Observational
obs_mean_mon = obs_ds.sithick.groupby('time.month').mean('time')
conc2 = obs_mean_mon

## Error of Mean Cycle
error_mean_conc = abs(conc1-conc2) ## Mean Absolute Error ##
error_mean_conc
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



# %% Visualisation
import calendar
import matplotlib.pyplot as plt

months = list(calendar.month_abbr)[1:]

for i in range(0,12,1):
    diff = (obs_mean_mon[i,:,:])
    diff.plot(cmap="magma", vmin=0, vmax=3)
    plt.title(f'Model to Obs Difference for {months[i]}')
    plt.xlabel("x (meters)"); plt.ylabel("y (meters)");
    plt.show()

#### SITool Sea Ice Edge Extent HIGH RESOLUTION COMPATIBLE ####
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
# print(obs_ds.attrs)

# %% Sea Ice Edge Extent Metrics

# INPUT model_ds & obs_ds with 'siconc" variable
# OUTPUT IIEE, AEE, ME, error_mean

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
conc1mask = xr.where((conc1 >= 15.0), 1, 0) # Where conc1 >= 15
conc2mask = xr.where((conc2 >= 15.0) & (conc2 <= 95), 1, 0) # Where conc2 >= 15

conc2mask[1,:,:].plot()
cellarea = 1

IIEE = ((abs(conc1mask-conc2mask))*cellarea).sum(dim=["x","y"])
AEE = abs((((conc1mask-conc2mask))*cellarea).sum(dim=["x","y"]))
ME=IIEE-AEE

ME

ndpm=[31,28.25,31,30,31,30,31,31,30,31,30,31]; #Number of days per month
error_mean = ((IIEE*ndpm).sum())/365.25
print(error_mean.values)




#### SITool Sea Ice Concentration HIGH RESOLUTION COMPATIBLE ####
# Cameron Carver Sept - 2024
# University of Cape Town - Applied Ocean Science MSc

# %% Import relevant packages
import os
os.chdir('/Users/crcarver/Desktop/AOS_THESIS/MSc_Cameron/Data/') # Set base directory
import xarray as xr

# %% Import Model Data 
access_model_data = 'siconc_SImon_CMCC-CM2-SR5_omip1_r1i1p1f1_gn_163801-200912_reproj.nc' #access reprojected model sea ice dataset
model_ds = xr.open_dataset(access_model_data) # Open dataset
model_ds = model_ds.sel(time=slice('1979-01-01','2009-12-31')) #Time slice of data
# print(model_ds)

# %% Import Observation Data
access_obs_data='seaice_conc_monthly_sh_n07_v04r00_all.nc' #access combined monthly observational mean data from NSIDC
obs_ds = xr.open_dataset(access_obs_data) # Open dataset
obs_ds = obs_ds.sel(tdim=slice('1979-01-01','2009-12-31')) #Time slice of data
print(obs_ds)

# %% Sea Ice Concentration Metrics
# %%% Monthly Mean Cycle

## Monthly Mean
# Model
model_mean_mon = model_ds.siconc.groupby('time.month').mean()     
model_mean_ds = xr.Dataset({"siconc": model_mean_mon})
# print(model_mean_ds)

#Observational
obs_mean_mon = obs_ds.cdr_seaice_conc_monthly.groupby('time.month').mean()
obs_mean_ds = xr.Dataset({"cdr_seaice_conc_monthly": obs_mean_mon})
siconc = obs_mean_ds['cdr_seaice_conc_monthly']/2.55*100 # = obs_mean_ds.cdr_seaice_conc_monthly/2.55*100
obs_mean_ds['siconc'] = siconc
# print(obs_mean_ds)

## Error of Mean
error_mean_conc = abs(model_mean_ds["siconc"] - obs_mean_ds["siconc"])
error_mean_ds= xr.Dataset({"siconc": error_mean_conc})
# print(error_mean_ds)

# %% Plotter
error_mean_ds.siconc.isel(month=2).plot()

#### REPROJECTION OF MODEL DATA ONTO NSIDC GRID ####

# %% Import Relevant Libraries
import os
os.chdir('/Users/crcarver/Desktop/AOS_THESIS/MSc_Cameron/Data/')
import xarray as xr
import rioxarray
import netCDF4 as nc


# %% Import Model Data (SOURCE)
access_model_data = '/Users/crcarver/Desktop/AOS_THESIS/MSc_Cameron/SITool/SICONC/Model-data/siconc_SImon_CMCC-CM2-SR5_omip1_r1i1p1f1_gn_163801-200912.nc' #access sea ice dataset
model_ds = xr.open_dataset(access_model_data)
# model_ds_reset = model_ds.reset_cooref_ds(["latitude","longitude"])
model_ds = model_ds.sel(time=slice('1979-01-01', '2017-12-31')) # i=slice(0, 360), j=slice(0,270)
print(model_ds)
model_ds.variables
print(model_ds['type'])
# model_ds_reset.siconc.isel(time=1).plot()

model_lat = model_ds['latitude']
model_lon = model_ds['longitude']
print(model_lat)

# %% Import Reference Data (TARGET)
access_ref_data = 'G02202-cdr-ancillary-sh.nc'
ref_ds = xr.open_dataset(access_ref_data)
print(ref_ds)
ref_ds.variables

ref_lon=ref_ds['longitude']
ref_lat=ref_ds['latitude']


# %% Observation Data
access_obs_data='seaice_conc_monthly_sh_n07_197811_v03r01.nc'
obs_ds = xr.open_dataset(access_obs_data)
print(obs_ds)
obs_ds.cdr_seaice_conc_monthly.isel(tdim=0).plot()


# %% Reproject
model_ds = model_ds.rio.write_crs(ref_ds.rio.crs)

import rasterio
from rasterio.warp import reproject, Resampling
import numpy as np

# Define the source and target coordinate reference systems (CRS)
model_ds = model_ds.rio.write_crs("EPSG:4326")
ref_ds = ref_ds.rio.write_crs("EPSG:3976")
src_crs = model_ds.rio.crs
dst_crs = ref_ds.rio.crs

# Reproject the data
model_ds = model_ds.rename({'i':'x','j':'y'})
model_ds = model_ds.rio.set_spatial_dims(x_dim='x',y_dim='y', inplace=True)
ref_ds = ref_ds.rio.set_spatial_dims(x_dim='x',y_dim='y', inplace=True)

reprojected_ds = model_ds.rio.reproject_match(ref_ds)


reproject(
    source=model_ds['siconc'].values,
    destination=model_ds['siconc'].values,
    src_transform=model_ds.rio.transform(),
    src_crs=src_crs,
    dst_transform=ref_ds.rio.transform(),
    dst_crs=dst_crs,
    resampling=Resampling.nearest
)
print(model_ds)
model_ds.siconc.isel(time=1).plot()






#### NetCDF Monthly Dataset Combiner ####
# Cameron Carver Sept - 2024
# University of Cape Town - Applied Ocean Science MSc

import glob
import xarray as xr
# Define path to directory that contains all netCDF files to be combined
data_dir = glob.glob('/Users/crcarver/Desktop/AOS_THESIS/MSc_Cameron/Data-Obs/G02202_V4_south_monthly/*.nc')

datasets = []
for file_path in data_dir:
    obs_ds = xr.open_dataset(file_path)
    datasets.append(obs_ds)
    
combined = xr.concat(datasets, dim='tdim')
combined = combined.assign_coords(tdim=combined['time'], x=combined['xgrid'], y=combined['ygrid'])
combined = combined.sortby('tdim')
siconc = combined['cdr_seaice_conc_monthly']/2.55*100 # Convert siconc to a percentage
combined['siconc'] = siconc

# Define desired name of combined file
combined.to_netcdf('seaice_conc_monthly_sh_n07_v04r00_all.nc')
print(combined)

# %% Prompts

## Data directory
# # Prompt user for directory path
# user_path = input("Please enter the directory path: /Users/.../xxx/ -> ")
# directory_end = "*.nc"
# full_path = "%s%s" % (user_path, directory_end)
# data_dir = glob.glob(full_path)

## Filename prompt
# # Prompt user for destination directory path
# dest_path = input("Destination directory path: /Users/.../xxx/ -> ")
# import os
# os.chdir(dest_path)
# file_name = input("Desired name of combined file: ")
# file_end = ".nc"
# file_name = "%s%s" % (file_name, file_end)
# combined.to_netcdf(file_name)
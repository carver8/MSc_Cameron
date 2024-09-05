#### REPROJECTION OF MODEL DATA ONTO NSIDC GRID ####

# %% Import Relevant Libraries
import os
os.chdir('/Users/crcarver/Desktop/AOS_THESIS/MSc_Cameron/Data/')
import copy
import json
import warnings

import cf_xarray as cfxr
import matplotlib.pyplot as plt
import shapely
import xarray as xr
import xesmf as xe
from clisops.core.subset import subset_bbox  # For subsetting
# from xclim.testing import open_dataset  # For opening xclim's test data

cmap = copy.copy(plt.cm.get_cmap("viridis"))
cmap.set_bad("lightgray")
# %% Import Model Data (SOURCE)
access_model_data = 'siconc_SImon_CMCC-CM2-SR5_omip1_r1i1p1f1_gn_163801-200912.nc' #access sea ice dataset
model_ds = xr.open_dataset(access_model_data)

# %% Import Reference Data (TARGET)
access_ref_data = 'G02202-cdr-ancillary-sh.nc'
ref_ds = xr.open_dataset(access_ref_data)

# %% Observation Data
access_obs_data='seaice_conc_monthly_sh_f13_200001_v03r01.nc'
obs_ds = xr.open_dataset(access_obs_data)

# %% INPUT - MODEL
# Let's look at the grid shape itself and the data for one time step
fig, axs = plt.subplots(ncols=2, figsize=(12, 4))

axs[0].scatter(x=model_ds.longitude.values, y=model_ds.latitude.values, s=0.1)
axs[0].set_title(
    "The input horizontal grid points as seen on a lat/lon map.\nOnly the northern hemisphere is shown."
)
axs[0].set_ylim(-90, -50)
axs[0].set_ylabel(f"latitude [{model_ds.latitude.units}]")
axs[0].set_xlabel(f"longitude [{model_ds.longitude.units}]")
model_ds.siconc.isel(time=0).plot(ax=axs[1], cmap=cmap)
axs[1].set_title("Sea Ice concentration - CMCC-CM2-SR5")
axs[1].set_ylim(0, 90)
fig.tight_layout()

# %% TARGET - REFERENCE
# For this example, we're not interested in the observation data, only its underlying grid, so we'll select a single time step.

# Subset over the Hudson Bay and the Labrador Sea for the example
bbox = dict(lon_bnds=[-180, 180], lat_bnds=[-80, -50])
ds_tgt = subset_bbox(ref_ds, **bbox)
ds_tgt

ds_tgt.cf.plot.scatter(x="longitude", y="latitude", s=0.1)
plt.title("Target regular grid");

model_ds.cf.describe()

reg_bil = xe.Regridder(model_ds, ds_tgt, "bilinear")
reg_bil  # Show information about the regridding

# NBVAL_IGNORE_OUTPUT

# xesmf/frontend.py:476: FutureWarning: ``output_sizes`` should be given in the ``dask_gufunc_kwargs`` parameter. It will be removed as direct parameter in a future version.
warnings.filterwarnings("ignore", category=FutureWarning)

# Apply the regridding weights to the input sea ice concentration data
sic_bil = reg_bil(model_ds.siconc)

# %% Plot the results

# for i in range(0,1,1):
#     fig = plt.figure()
sic_bil.isel(time=1).plot(cmap=cmap)
plt.title("Regridded siconc data");





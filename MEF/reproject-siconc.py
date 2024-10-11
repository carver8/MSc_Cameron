#### REPROJECTION OF MODEL DATA ONTO NSIDC GRID ####
## Adapted from: https://pavics-sdi.readthedocs.io/en/latest/notebooks/regridding.html
# %% Import Relevant Libraries
import os
os.chdir('/Users/crcarver/Desktop/AOS_THESIS/MSc_Cameron/Data/')
import copy
import warnings
import matplotlib.pyplot as plt
import xarray as xr
import xesmf as xe
from clisops.core.subset import subset_bbox  # For subsetting

cmap = copy.copy(plt.get_cmap("magma"))
cmap.set_bad("lightgrey")

# %% Import Reference Data (TARGET)
access_ref_data = 'OBSERVATIONS/SICONC/G02202-cdr-ancillary-sh.nc'
ref_ds = xr.open_dataset(access_ref_data)

# %% Import Model Data (SOURCE)
access_model_data = 'MODELS/SICONC/siconc_SImon_CMCC-CM2-SR5_omip1_r1i1p1f1_gn_163801-200912.nc' #CMCC-SR5-OMIP1
# access_model_data = 'MODELS/SICONC/siconc_SImon_CMCC-CM2-SR5_omip2_r1i1p1f1_gn_165301-201812.nc' #CMCC-SR5-OMIP2
# access_model_data = 'MODELS/SICONC/siconc_SImon_MRI-ESM2-0_omip1_r1i1p1f1_gn_170001-200912.nc' #MRI-ESM2-OMIP1
# access_model_data = 'MODELS/SICONC/siconc_SImon_GFDL-CM4_omip1_r1i1p1f1_gn_190801-200712.nc' # GFDL-CM4-OMIP1 ##!! MUST ADJUST GEOLAT/LON
model_ds = xr.open_dataset(access_model_data)
model_name = model_ds.source_id; exp_id = model_ds.experiment_id
# %% INPUT - MODEL

#GFDL ADJUSTMENT - also need to change scatterplot to lon and lat
# model_ds = model_ds.assign_coords(lat=model_ds['GEOLAT'])
# model_ds = model_ds.assign_coords(lon=model_ds['GEOLON'])

# Let's look at the grid shape itself and the data for one time step
fig, axs = plt.subplots(ncols=2, figsize=(12, 4))

axs[0].scatter(x=model_ds.longitude.values, y=model_ds.latitude.values, s=0.1)
axs[0].set_title("The input horizontal grid points as seen on a lat/lon map.")
axs[0].set_ylim(-90, 0)
# axs[0].set_ylabel(f"latitude [{model_ds.latitude.units}]")
# axs[0].set_xlabel(f"longitude [{model_ds.longitude.units}]")
model_ds.siconc.isel(time=8).plot(ax=axs[1], cmap=cmap)
axs[1].set_title("Sea Ice concentration - CMCC-CM2-SR5")
axs[1].set_ylim(-90, 0)
fig.tight_layout()

# %% TARGET - REFERENCE

# Visualization of input data and its corresponding grid
bbox = dict(lon_bnds=[-180, 180], lat_bnds=[-80, -50])
ds_tgt = subset_bbox(ref_ds, **bbox)

ds_tgt.cf.plot.scatter(x="x", y="y", s=0.1)
plt.title("Target regular grid");

reg_bil = xe.Regridder(model_ds, ds_tgt, "bilinear", periodic=True,)
reg_bil  # Show information about the regridding

# xesmf/frontend.py:476: FutureWarning: ``output_sizes`` should be given in the ``dask_gufunc_kwargs`` parameter. It will be removed as direct parameter in a future version.
warnings.filterwarnings("ignore", category=FutureWarning)

# Apply the regridding weights to the input sea ice concentration data
sic_bil = reg_bil(model_ds.siconc)

model_rp = xr.Dataset({"siconc": sic_bil})

# Overlay land mask from observational ancillary file
mask = xr.where(ref_ds.landmask > 0, float('nan'), 1)
tct = model_rp.coords['time'].size
for i in range(0,tct,1):
    model_rp['siconc'][i,:,:] = model_rp['siconc'][i,:,:]*mask

model_rp = model_rp.assign_attrs(model_name=model_name, exp_id=exp_id)
model_rp.to_netcdf(f'{access_model_data[:-3]}_reproj.nc')

# %% Visualisation

# for i in range(0,11,1):
    # fig = plt.figure()
model_rp.siconc.isel(time=8).plot(cmap=cmap)
plt.title("Regridded siconc data");

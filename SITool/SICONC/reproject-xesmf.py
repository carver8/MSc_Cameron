#### REPROJECTION OF MODEL DATA ONTO NSIDC GRID ####

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
access_model_data = 'MODELS/SICONC/siconc_SImon_MIROC6_omip2_r1i1p1f1_gn_195301-201812.nc' #access sea ice dataset
model_ds = xr.open_dataset(access_model_data)
model_ds.siconc[1,:,:].plot()

# %% INPUT - MODEL
# Let's look at the grid shape itself and the data for one time step
fig, axs = plt.subplots(ncols=2, figsize=(12, 4))

axs[0].scatter(x=model_ds.longitude.values, y=model_ds.latitude.values, s=0.1)
axs[0].set_title(
    "The input horizontal grid points as seen on a lat/lon map."
)
axs[0].set_ylim(-90, -50)
axs[0].set_ylabel(f"latitude [{model_ds.latitude.units}]")
axs[0].set_xlabel(f"longitude [{model_ds.longitude.units}]")
model_ds.siconc.isel(time=8).plot(ax=axs[1], cmap=cmap)
axs[1].set_title("Sea Ice concentration - CMCC-CM2-SR5")
axs[1].set_ylim(0, 90)
fig.tight_layout()

# %% TARGET - REFERENCE

# Visualization of input data and its corresponding grid
bbox = dict(lon_bnds=[-180, 180], lat_bnds=[-80, -50])
ds_tgt = subset_bbox(ref_ds, **bbox)
# print(ds_tgt)

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

model_rp.to_netcdf('siconc_SImon_CMCC-CM2-SR5_omip1_r1i1p1f1_gn_163801-200912_reproj.nc')

# %% Visualisation

# for i in range(0,11,1):
#     fig = plt.figure()
model_rp.siconc.isel(time=10).plot(cmap=cmap)
plt.title("Regridded siconc data");

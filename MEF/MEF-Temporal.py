#### MEF Metric ####
# Cameron Carver Oct - 2024
# University of Cape Town - Applied Ocean Science MSc

# %% Import relevant packages and data
# %%% Import packages
import os; os.chdir('/Users/crcarver/Desktop/AOS_THESIS/MSc_Cameron/Data/') # Set base directory
import xarray as xr
import pandas as pd
import numpy as np
import calendar
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# %% Plot Details
bins = [-np.inf, 0, 0.2, 0.5, 0.65, np.inf]
bins2 = [-np.inf, 1, 2, 3, 4, np.inf]
labels = ['Bad', 'Poor', 'Good', 'Very Good', 'Excellent',]
colors = {
    'Bad': 'firebrick',
    'Poor': 'sandybrown',
    'Good': 'yellowgreen',
    'Very Good': 'skyblue',
    'Excellent': 'royalblue',
}
cmap = ListedColormap(colors)
# %%% Import Observation Data
access_obs_data='OBSERVATIONS/SICONC/seaice_conc_monthly_sh_n07_v04r00_all.nc' #access combined monthly observational mean data from NSIDC
obs_ds = xr.open_dataset(access_obs_data) # Open dataset
obs_ds = obs_ds.sel(tdim=slice('1979-01-01','2009-12-31')) #Time slice of data
obs_da = obs_ds.cdr_seaice_conc_monthly
obs_da = obs_da.where(obs_da <= 2.51, other=float('nan'))

# %%% Import Model Data 
# access_model_data = 'MODELS/SICONC/siconc_SImon_CMCC-CM2-SR5_omip1_r1i1p1f1_gn_163801-200912_reproj.nc' #CMCC-SR5-OMIP1
# access_model_data = 'MODELS/SICONC/siconc_SImon_CMCC-CM2-SR5_omip2_r1i1p1f1_gn_165301-201812_reproj.nc' #CMCC-SR5-OMIP2
access_model_data = 'MODELS/SICONC/siconc_SImon_MRI-ESM2-0_omip1_r1i1p1f1_gn_170001-200912_reproj.nc' #MRI-ESM2-OMIP1
# access_model_data = 'MODELS/SICONC/siconc_SImon_GFDL-CM4_omip1_r1i1p1f1_gn_190801-200712_reproj.nc' # GFDL-CM4-OMIP1 ##!! MUST ADJUST TIME TO END 2007

model_ds = xr.open_dataset(access_model_data) # Open dataset
model_ds = model_ds.sel(time=slice('1979-01-01','2009-12-31')) #Time slice of data

model_ds = model_ds.assign_coords(time=obs_ds['time'])  # Set time coords of model ds to be the same as the obs ds
model_ds = model_ds.swap_dims({'tdim': 'time'})
model_ds = model_ds.swap_dims({'time': 'tdim'})
model_da = model_ds.siconc / 100

model_name, exp_id = model_ds.model_name, model_ds.exp_id   # Call model name from attributes
# %% New Metrics
# MEF T-S-T
MEF_num = (obs_da-model_da) ** 2#* mask
MEF_num = MEF_num.mean(dim=['x','y'])

o_mean = obs_da.mean(dim='tdim')
MEF_den = (obs_da-o_mean) ** 2 #* mask
MEF_den = MEF_den.mean(dim=['x','y'])

vmask = xr.where(MEF_den != 0, 1, float('nan'))
MEF = (1 - (MEF_num/MEF_den))#*vmask
MEFm = MEF*vmask

plt.figure(figsize=(6, 4))
MEFm.plot()
plt.xlabel('Year')
plt.ylabel('MEF')
plt.title(f'{model_name}-{exp_id}')

# %%%  Line Plot with each month as own line
months = list(calendar.month_abbr)[1:]
colors2 = [
    'darkred', 'darkred', 'blue', 'blue', 'blue', 'peru', 
    'peru', 'peru', 'green', 'green', 'green', 'darkred'
]
plt.figure(figsize=(8, 5))
for month in range(1, 13):
    monthly_data = MEFm.sel(tdim=MEFm['time.month'] == month)
    plt.plot(monthly_data['tdim'], monthly_data, label=f'{months[month-1]}', color=colors2[month-1])
plt.xlabel('Year'); plt.ylabel('MEF'); plt.title('Time Series of MEF for Each Month')
plt.legend(title='Month', bbox_to_anchor=(1.05, 1), loc='best', borderaxespad=0)
plt.show()
# %% Stacked Bar Graph - monthly
monthly_data = MEFm.resample(tdim='MS').mean()

def bin_data(array, bins):
    return np.digitize(array, bins)

binned_data = xr.apply_ufunc(
    bin_data, 
    monthly_data, 
    input_core_dims=[[]],  # No core dimensions
    kwargs={'bins': bins}
)

df = binned_data.to_dataframe(name='value').reset_index()

# Count occurrences in each bin
df['bin'] = pd.cut(df['value'], bins=bins2, labels=labels)
bin_counts = df.groupby([df['tdim'].dt.month, 'bin'], observed=False).size().unstack(fill_value=0)
bin_counts.index = pd.to_datetime(bin_counts.index, format='%m').month_name()

# Plot stacked bar graph
ax = bin_counts.plot(kind='bar', stacked=True, color=[colors[label] for label in bin_counts.columns])

ax.set_xticklabels(bin_counts.index, rotation=45)

plt.xlabel('Month'); plt.ylabel('Count'); plt.title(f'{model_name}-{exp_id}')
plt.legend(title='MEF Category', bbox_to_anchor=(1.05, 1), loc='best', borderaxespad=0)
plt.show()

# %% Stacked Bar Graph - seasonaly
seasonal_data = MEFm.resample(tdim='QS-DEC').mean()

def bin_data(array, bins):
    return np.digitize(array, bins)

binned_data = xr.apply_ufunc(
    bin_data, 
    seasonal_data, 
    input_core_dims=[[]],  # No core dimensions
    kwargs={'bins': bins}
)

df = binned_data.to_dataframe(name='value').reset_index()

# Count occurrences in each bin
df['bin'] = pd.cut(df['value'], bins=bins2, labels=labels)
bin_counts1 = df.groupby([df['tdim'].dt.month, 'bin'], observed=False).size().unstack(fill_value=0)
bin_counts1.index = ['MAM', 'JJA', 'SON', 'DJF']

# Plot stacked bar graph
ax = bin_counts1.plot(kind='bar', stacked=True, color=[colors[label] for label in bin_counts1.columns])

plt.xlabel('Season'); plt.ylabel('Count'); plt.title(f'{model_name}-{exp_id}')
plt.legend(title='MEF Category', bbox_to_anchor=(1.05, 1), loc='best', borderaxespad=0)
plt.show()
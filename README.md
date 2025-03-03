# ColdSnapTracking

This package allows you to track cold spell events, by 
- Deseasonalising the data
- Applying a persistence and a percentile based threshold
- Finding contours of these extreme temperature anomalies
- Tracking these contours in space and time

Please contact Weronika Osmolska (eewo@leeds.ac.uk) if you have any questions or require more information.


## Installation

The dependencies for this package are listed in requirements.txt and can be downloaded through mamba as

```Bash
$ mamba install --file requirements.txt
```

Or conda as 

```Bash
$ conda install -c conda-forge --file requirements.txt 
```

Once complete, the package can then install with

```Bash
$ pip install git+https://github.com/Zun-1-1/ColdSnapTracking
```

Note: this package is only compatible with Python 3.9.

## Example

The following shows how to run the code. These cells need to be run through slurm using a python script
on a supercomputer. In some cases in many parallel runs.

Firstly, we deseasonalise the data

```Python

import xarray as xr
import coldsnap as cs
import numpy as np

#______________INPUTS______________
#Directories
path = 'path_to_the_temperature/temperature_sf*.nc'
out_file = 'dir_to_save_deseasonalised_t_/deseason_linear_fit' 
tmpdev_package = 'path_to_tmpdev_package/'
#Time range
year_start = 1940
year_end = 2023
#Choose a start and end for the latitude, to allow more parallel runs at lower ram
start=0
end=10
#______________INPUTS______________



#Open datasets using xarrays.
year_list = list(range(year_start,year_end,1))

temperature_surface = xr.open_mfdataset(path,parallel=True,chunks={'time': 2})
temperature_surface = temperature_surface.sel(time=temperature_surface.time.dt.year.isin(year_list))
#Run deseasonalisation
d_ave = cs.deseason(temperature_surface,start,end,year_start,year_end,poly_fit=1)
#Save
comp = dict(zlib=True, complevel=9)
job = d_ave.load().astype(dtype='float64').to_netcdf(out_file + str(sys.argv[1]) + '_' + '.nc', compute=False, encoding= {'t2m' : comp})
job.compute()
```

Then we only want to keep only the extended winter season

```Python
#______________INPUTS______________
#Deseasonalised and taking just the extended winter once we complete all the runs:
extnd_winter_temps = 'location_of_deseasonalised_and_extended_winter_temperature_anomalies'
#______________INPUTS______________


dave = xr.open_mfdataset(out_file+'*.nc')
daily_average = dave.resample(time='1D').mean(dim='time')
daily_average_EXT_winter = daily_average.sel(time=daily_average.time.dt.month.isin([1, 2, 3, 11, 12])).chunk(dict(time=-1))
comp = dict(zlib=True, complevel=9)
job = daily_average_EXT_winter['t2m'][:,:,:].load().astype(dtype='float32').to_netcdf(extnd_winter_temps + '_' + '.nc', compute=False, encoding= {'t2m' : comp})
job.compute()
```


Once you have the full dataset of deseasonalised data we need to find data that meets persistence and 5th percentile threshold


```Python
import xarray as xr
import numpy as np
import coldspell as cs
#______________INPUTS______________
#Place to save data of t2m and t2m-anomaly when persistance and temperature thresholds are met
outfile_data = 'dir_1'
outfile_sta = 'dir_2'
#______________INPUTS______________


temperature_surface = xr.open_mfdataset('dir_to_resampled_extended_winter_temperature') 
#Follow above points for directions how to resample

dave = xr.open_mfdataset(extnd_winter_temps)


#Run the code
al = cs.cao_diag_st(temperature_surface['t2m'],dave['t2m'],3,0.05)

comp = dict(zlib=True, complevel=5)
job = al['data'][:,:,:].load().astype(dtype='float32').to_netcdf(outfile_data + '_' + '.nc', compute=False, encoding= {'data' : comp})
job.compute()


comp = dict(zlib=True, complevel=5)
job = al['sta'][:,:,:].load().astype(dtype='float32').to_netcdf(outfile_sta + '_' + '.nc', compute=False, encoding= {'sta' : comp})
job.compute()

```

Now we need to find contours of extreme temperature anomalies

```Python
import xarray as xr
import pickle 
import lzma
import coldsnap as cs


#The deseasonalised+persistent+extreme data
STA = xr.open_mfdataset(outfile_sta)['sta'][start:end,:,:].load()

#The weighted average wind
uwind = xr.open_mfdataset('path_to_u_winds/uwind/*')
uwind = uwind.sel(time=uwind.time.dt.month.isin([1, 2, 3, 11, 12]))
uwind = uwind['u'][start:end,:,:].load()
vwind = xr.open_mfdataset('path_to_v_winds/vwind/*')
vwind = vwind.sel(time=vwind.time.dt.month.isin([1, 2, 3, 11, 12]))
vwind = vwind['v'][start:end,:,:].load()

#Get the area matrix
dA = cs.degrees_into_area(STA['latitude'],STA['longitude'])

#CETA datasets
dictionary_ = cs.area_constraint(STA,dA,uwind,vwind)
file = out_file+str(start)+'.xz'

del dictionary_['cntr']
with lzma.open(file, "wb") as f:
    pickle.dump(dictionary_, f)
```


Now we need to track all CETAs once we found the contours for all time periods


```Python 
import coldsnap as cs

#______________INPUTS______________
save_tracks = 'where_to_save_the_tracks'
year = 1940 #Run parallel with different years
speed = 20 # speed threshold for advection
advection = 'Total' # this is advection
min_ar = 70000 # minimum area
#______________INPUTS______________



components = ["time","COM","area","uwind","vwind","contour_no", "gridpoints_inside","uwind_cntr","vwind_cntr","contour_pnts","uwind_points","vwind_points"]
#This loads back the contour data in the same increments you saved it
data = cs.load_by_year(year,end-start,components,file_name= out_file)


    
filename = save_tracks+'_'+advection+'_year_'+str(year)+'_areathresh_'+str(min_ar)+'_maxspeed_'+str(speed)+'.npy'
track = cs.trajectories(data,filename,resolution=0.25,area_thres=min_ar,overlap=False,advection=advection,max_speed=speed)

```

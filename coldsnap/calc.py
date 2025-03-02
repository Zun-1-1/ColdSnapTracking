import xarray as xr
import numpy as np
from . import persistance as ps
from . import area as ar
#from . import deseason as de
import scipy.stats as stats

import pandas as pd 
def cao_diag_st(org_data,surface_temperature,N=3,percentile=0.05,temperature_th=None):
    """
    Calculating Surface Temperature Anomalies
    -------------------------------------------------------------------------------------
    Notes:
    Calculates the surface temperature anomalies based on a X percentile that we are
    looking below, and Nth consequetive day that this occurs over. Extra area constraint
    is also present. 

    NOTE: org_data and surface_temperature(deseasonalised) need to be time sorted in the
    same order.
    -------------------------------------------------------------------------------------

    Parameters (IN):

    surface_temperature: xarray, dim: {time,lat,lon}, units [K]
                        Surface temperature throughout the whole dataset (deseasonalised).

    N: int, dim{1}, dimensionless
        Number of consecutive days we want to consider.

    percentile: float, dim{1}, dimensionless
                Percentile of temperatures that we want to consider.

    temperature_th : int, dim{1}, [degree/kelvin change]
                Temperature above zero that we want to remove from data. Optional, i.e if inputted 5 
                then anything below 5 degree will be considered as a CAO

    Parameters (OUT):

    xarray.Dataset:
        -data: Only the temperatures of regions where conditions are met are
        specified otherwise the value is nan.
        -percentile: The percentiles of these regions.
        -sta: The surface temperature anomalies of these regions
    """
    #Dimension variables
    time, lat, lon = surface_temperature.shape
    #Resample original data:
    #org_data = resample(org_data)
    #Average surface temperature over each day.
    daily_average = surface_temperature.resample(time='1D').mean(dim=surface_temperature.dims[0])
    daily_average_EXT_winter = daily_average.sel(time=daily_average.time.dt.month.isin([1, 2, 3, 11, 12])).chunk(dict(time=-1))
    #Find index where we are time discontinous
    index = index_finder(daily_average_EXT_winter)

    #Value below which the data has to be below, to be in the percentile.
    quantile_value = daily_average_EXT_winter.chunk(dict(time=-1)).quantile(percentile,dim=daily_average_EXT_winter.dims[0],keep_attrs=True)
    #Show only data where this we are below this percentile
    percentile_data =  np.array(daily_average_EXT_winter.where(daily_average_EXT_winter<quantile_value))
    #Cold air outbreaks are only accounted when the temperature is 0 or less.
    orginal_2m = np.array(org_data)
    #If we want a temperature threshold
    if temperature_th is not None:
        percentile_data[orginal_2m>(273+temperature_th)] = np.NaN
    #ATTEMPT TO SAVE RAM (alternative to above)
    # for i in range(time):
    #     print(i)
    #     #Original 2m temperature dataset in a numpy array
    #     orginal_2m = np.array(org_data[i])
    #     percentile_data[i][orginal_2m>(273+temperature_th)] = np.NaN
    #     del orginal_2m
    
    #Has to be consequtive for 3 days at least
    percentile_data = ps.consecutiveness(percentile_data,index,N)
    
    #Convert data from detrended surface temperature anomalies to original temperature for the CAO datapoints
    percentile_data_cp = percentile_data.copy()
    percentile_data_cp[~np.isnan(percentile_data_cp)] = orginal_2m[~np.isnan(percentile_data_cp)] 
    #ATTEMPT TO SAVE RAM (alternative option to above)
    # for i in range(time):
    #     mesh = ~np.isnan(percentile_data[i])
    #     orginal_2m = np.array(org_data[i])
    #     percentile_data[i][mesh] = orginal_2m[mesh] 
    #     del orginal_2m

    
    #Find the percentiles of all data points
#    distribution = rel_percentile(np.array(surface_temperature),percentile_data)

    # #Area constraint
    #dA = ar.degrees_into_area(surface_temperature[lat],surface_temperature[lon],degree=True,equal_spacing=True)
    # print("done")
    #filled_areas,day = ar.area_constraint(percentile_data,area,dA)
    # print("done")

    #Find the percentile data without the removed season and determine its deviation away from
    #the percentile.
    # daily_average_np = np.array(surface_temperature)
    # daily_average_np[np.isnan(percentile_data)] = np.NaN

    #Percentage deviation away from the percentile
    #distribution =  ((-daily_average_np + np.array(quantile_value)) *100)/np.array(quantile_value)
    #Coodinates for the xarray
    coordinates = {
                **{ 
                    key: value
                    for key, value in surface_temperature.coords.items()},
                    'ConsecutiveDays' : np.array(N),
                    'Quantile' : np.array(percentile)
                }
    #Return xarray.
    return xr.Dataset({'data': (surface_temperature.dims,percentile_data_cp,{'Description': 'Data filtered by quantile and consequtive days [K]'})
                      # ,'percentiles': (surface_temperature.dims,distribution,{'Description': 'Percentiles [%]'})
                       ,'sta': (surface_temperature.dims,percentile_data,{'Description': 'Surface temperature anomalies [K]'})
                },coords=coordinates)#,day


def rel_percentile(surface_temperature_,percentile_data_):
    """
    Relative percentiles of the remaining temperatures.
    -----------------------------------------------------------------

    Parameters (IN):

    surface_temperature : numpy.ndarray, dim = {time,lat,lon}, units [K]
                Deseasonalised surface temperature.
    
    percentile_data : numpy.ndarray, dim = {time,lat,lon}, units [K]
                The daily average deseasonalised temperature for grid points where
                temperature is in the xth percentile. 

    Parameters (OUT):

    distribution : np.ndarray, dim = {time,lat,lon}, units [%]
                The temperatures where the percentiles are below (?)%,
                will be written as the percentile.
                
    """
    percentile_data = np.array(percentile_data_.to_array())[0]
    surface_temperature = np.array(surface_temperature_.to_array())[0]
    #Find the percentiles of all data points
    time, lat, lon = surface_temperature.shape
    #Empty Array
    distribution = np.zeros((time,lat,lon))
    for i in range(lat):
        print(i)
        for j in range(lon):
            #Need to use pandas for this
            df = pd.DataFrame(surface_temperature[:,i,j], columns=['t2m'])
            #Find percentiles
            perc = df["t2m"].apply(lambda x: stats.percentileofscore(df["t2m"],x, kind = 'weak',nan_policy = 'propagate'))
            #mask all values that are nan in the CAO data
            mask = ~np.isnan(percentile_data[:,i,j])
            #Copy the CAO data
            new_array = np.copy(percentile_data[:,i,j])
            #Anywhere where non nan CAO data, convert data to the percentile data
            new_array[mask] = np.array(perc)[mask]
            #Add to an array
            distribution[:,i,j] = new_array
    coordinates = {
                **{ 
                    key: value
                    for key, value in surface_temperature_.coords.items()}
                }
    #Return xarray.
    return xr.Dataset(#{'data': (surface_temperature.dims,percentile_data,{'Description': 'Data filtered by quantile and consequtive days [K]'})
                       {'dev': (surface_temperature_[list(surface_temperature_.keys())].dims,distribution,{'Description': 'Percentiles [%]'})
                       #,'area': (surface_temperature.dims,filled_areas,{'Description': 'Areas above the threshold'})
                },coords=coordinates)#,day


def index_finder(temperature_surface):
    """
    Finds the index of the first month of October, in order to prevent
    the consecutiveness(...) function, from finding consecutiveness between
    March and October.
    ------------------------------------------------------------------------
    Parameters (IN):
    temperature_surface : xarray, dim: {time,lat,lon}, units [K]
                        Surface temperature throughout the whole dataset, for the extended winter.
    Parameters (OUT):
    """
    Discontinuity = temperature_surface.where(((temperature_surface['time.day'] == 1) & (temperature_surface['time.month'] == 11)), drop=True)
    print(Discontinuity)
    Discontinuity = Discontinuity.where((Discontinuity['time.hour'] == 0), drop=True)
    print(Discontinuity,np.array(Discontinuity['time']))
    index = []
    for i in np.array(Discontinuity['time']):
        index.append(np.where(np.array(temperature_surface['time'])==i)[0][0])
    #index = np.where(np.array(temperature_surface['time'])==np.array(Discontinuity['time']))
    #print(index)
    return np.array(index)

def resample(func):
    """
    Resamples the function so that it is now in daily time,
    with times only in extended winter period.
    ------------------------------------------------------------------------
    Parameters (IN):
    func : xarray, dim: {Contain time, lat and lon}
            
    Parameters (OUT):
    func : xarray, dim: {as before but with altered time dimension}
    """ 
    func = func.sortby(['time'],ascending=True)

    func = func.resample(time='1D').mean(dim='time')
    func = func.sel(time=func.time.dt.month.isin([1, 2, 3, 11, 12])).chunk(dict(time=-1))

    return func
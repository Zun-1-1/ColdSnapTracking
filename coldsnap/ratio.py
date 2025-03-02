import numpy as np
import xarray as xr

def ratio(dp,cao_,start=None,end=None):
    """
    Overlap between CAM days and CAO days per gridcell.
    ---------------------------------------------------------
    Parameters (IN):

        dp : xarray, dim={time,lat,lon}, units:[hPa]
            Cold air mass depth (resampled)

        cao_ : xarray, dim={time,lat,lon}, unitless:
            Cold air outbreaks days

    Parameters (OUT):

        percentages : xarray, dim={lat,lon}
                Percentage of time cold air outbreak days and cold air mass occur at the same time
                per all cold air outbreak days.
    """
    #Sort by time, take each day into groups, average each group, then select extended winter
    dp = dp.sortby("time").resample(time='1D').mean(dim="time")
    dp = dp.sel(time=dp.time.dt.month.isin([11,12,1,2,3]))
    if start is not None and end is not None:
        data = cao_[:,start:end,:]
        dp = dp[:,start:end,:]
    else:
        dp = dp[:,:,:]
    time,lat,lon = data.shape
    percentages = np.zeros((lat,lon))
    for i in range(lat):
        dp1 = dp[:,i,:].load()
        data1_ = cao_[:,i,:].load()
        #Turn data into boolean
        #Starting with dp, where we remove dp=0 or dp=np.NaN
        zero_dp = dp1 < 10
        nans_dp = np.isnan(dp1)
        zero_dp |= nans_dp
        dp2 = ~zero_dp
        data2 = ~np.isnan(data1_)

                #Collapse the data to when both CAM and CAO occurred at the same time.
        cross_sec = np.logical_and(dp2,data2)
                #Sum of all CAO days / CAM days in the dataset
        cs_tot = np.array(cross_sec.sum(dim="time"))
        #dp_tot = np.array(dp2.sum(dim="time"))
        data_tot = np.array(data2.sum(dim="time"))

                #Input percentage into array
        percentages[i,:] = np.array((cs_tot * 100)/data_tot)
        percentages[i][data_tot==0] = -999

    #Make the dataset
    new_coords = {
   # 'longitude': xr.DataArray(
    #    data1_['longitude']),
    #'latitude': np.array(xr.DataArray(
       # data1_['latitude']))
       # }   
     **{
            key: value
            for key, value in data.drop("time").coords.items()
        }
    }
    percentages = np.transpose(percentages)

    #Return numpy.
    print(percentages.shape,lon,lat)
    return xr.Dataset({'overlap': ({'longitude': lon, 'latitude': lat},percentages,{'Description': 'Percentage overlap between cold air mass and cold air outbreaks as a percentage of CAO days'})
                       #,'dev': (surface_temperature.dims,distribution,{'Description': 'Percentiles [%]'})
                       #,'area': (surface_temperature.dims,filled_areas,{'Description': 'Areas above the threshold'})
                },coords=new_coords)#,day

            

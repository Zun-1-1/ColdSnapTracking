import numpy as np
import xarray as xr
import numpy as np
from numpy.polynomial import polynomial as P
from scipy.fft import fft, ifft, fftfreq
import cmath

def deseason(temperature_surface,start,end,year_start,year_end,poly_fit=2,var='t2m'):
    """
    Function to deseasonalise by taking first 3 Fourier Frequencies 
    and removing the trend in the data
    -------------------------------------------------------------------------
    NOTE: For RAM and speed resons I have added a start and end latitude,
    so that we can deseasonalise data in chuncks of latitudes and save the files
    seperately, then join them up via e.g xarray.open_mfdataset(...)
    -------------------------------------------------------------------------
    temperature_surface : xarray, dim:{time,lat,lon}, units = [K]
        The temperature at the surface for all recorded time.

    start, end: int
        This is the start and end of the latitude elements -- used to save ram. See example

    year_start,year_end: int
        This is the year range the dataset is running through e.h. 1940 to 2023

    poly_fit : int
        The fitting of the global increase in temperatures . 2 would be quadratic
        1 would be linear

    start, end : int
        The latitudes we want to loop through, can be for the whole dataset latitude
        e.f -90 to 90, or can do it in chunks (see notes)

        """
    #Deseasonalisation
    #Load the 2m temperature and find the daily mean
    daily_average = temperature_surface.resample(valid_time='1D').mean(dim=temperature_surface[var].dims[0])#.load()
    #Use only a subset of latitudes of the 2m temperature daily mean
    d_ave = daily_average[var][:,start:end,:].load()
    #Shape time,lat,lon
    _,lat,lon = daily_average[var].shape
    #Months in a year
    months = [31,29,31,30,31,30,31,31,30,31,30,31]

    #Between the choosen latitudes
    for no,i in enumerate(range(start,end)):
        #For all longitudes
        for j in range(lon):
            #Empty array for the mean surface temperature across day 0 to 365 of all the years
            mean_surf_t = np.zeros((np.sum(months)))
            #Start a counter for each day 
            count=0
            #Daily temperature average in a specific gridpoint
            T_ij = d_ave[:,no,j]
            tester = d_ave[:,no,j]
            #Select each day in the year, (by taking each day in the month) and then averaging the day over the years.
            for mth in range(0,12):
                #Choose the days in the month
                moth = months[mth]
                
                #For each day in the month
                for day in range(1,moth+1):
                    #Select each day in the year
                    month_ = T_ij.sel(valid_time=T_ij.valid_time.dt.month.isin([mth+1]))
                    day_ = month_.sel(valid_time=month_.valid_time.dt.day.isin([day]))
                    #Put the day into the daily temperature average
                    mean_surf_t[count] = day_.mean()
                    count+=1
                
            #Fourier transform the daily temperature average over the years
            yf = fft(mean_surf_t)
            #Steps of 1 since its daily
            xf = fftfreq(len(mean_surf_t),1)
            #Take the first 6 frequencies (both positive and negative count as one)
            bound = np.sort(np.abs(xf))[10]
            #Remove anything other than that
            yf[np.logical_or(xf>=bound,xf<=-bound)] = 0+0j
            #Inverse Fourier Transform
            #Only interested in the magnitude not the phase, we take the real part
            filtered = ifft(yf).real
            # #Can plot:
            # plt.plot(filtered,label="filtered")
            # plt.plot(mean_surf_t[:,0,0],label="Averaged T")
            # plt.xlabel("Time")
            # plt.ylabel("Temperature")
            #Select each year and remove the filtered data
            non_leap = np.concatenate([filtered[:60],filtered[61:]])
            for t in range(year_start,year_end):
            #for t in range(1940,1979):
    
                #Select each 365 or 366 days
                times = T_ij.sel(valid_time=str(t))
                #If the year is a leap year, then remove the whole filtered data (filtered consists of 366 days as it
                # encompasses the mean from each day from the dataset)
                if len(times)==366:
                    T_ij.loc[times['valid_time']] = T_ij.loc[times['valid_time']] - filtered
                #If the year is not a leap year, ignore the last day
                else:
                    T_ij.loc[times['valid_time']]= T_ij.loc[times['valid_time']] - non_leap
    
            #THIS PART REMOVES THE LONG TIME TREND
            #X is the range of the time, in steps of 1.
            X = np.array(list(range(len(T_ij['valid_time']))))
            # Fit a quadratic to the filtered data 
            fit = P.polyfit(X,T_ij,poly_fit)
            #Define the polynomial variable that will be fitted
            poly=np.zeros(len(X))
            # Loop over the quadratic variables to form a array we can plot.
            for q in range(len(fit)):
                poly = poly +  fit[q] * X**q
            #Remove this polynomial offset from the dataset.
            dave =  T_ij - poly
            d_ave[:,no,j] = T_ij - poly

    return d_ave
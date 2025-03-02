import numpy as np


def consecutiveness(percentile_data,index,N):
    """
    PERSISTENCE TEST OF ANOMOLOUS COLD TEMPERATURES
    ---------------------------------------------------------
    Parameters (IN):

        percentile_data : xarray, dim={time,lat,lon}, units:[K]
            Data that is filtered by a percentile.

        index : numpy.ndarray, dim={no of years}, unitless:
            Points to where discontinuities in time occur in the dataset, 
            since we are only looking at extended winter period.

        N : int, dim=1
            Number of days we would like anomolous temperatures to 
            prevail for at a single gridpoint. 

    
    Parameters (OUT):

        percentile_data : xarray, dim={time,lat,lon}, units:[K]
                Data that is filtered by a percentile and that 
                only datapoints are non zero, when the gridpoint
                is below the given percentile for at least N days.
    """
    index = np.array(index).flatten()
    #Shape of matrix
    time,lat,lon = percentile_data.shape
    #Empty False mask
    mask = np.zeros((time,lat,lon),dtype=bool)
    
    #Loop for all gridpoints
    for i in range(lat):
        for j in range(lon):
            #Take the time slice across each gridpoint
            data = percentile_data[:,i,j].astype(float)
            #Boolean matrix of NaN = True, else False
            mat = np.array(np.isnan(data))
            #Make sure it does not cumulatively sums between the discontinuity of March to October
            #By creating an additional value that we will remove later.
            #Insert Extra NaNs between these discontinuities 
            mat = np.insert(mat,index,[True]*len(index))
            #Cumulatively sum the booleans, True = 1, False = 0
            sum = mat.cumsum()
            #Mask all NaNs using sum[~mat]
            #Find all the unique values and the number of unique values
            val, no = np.unique(sum[~mat], return_counts=True)
            #print(val,no)
            #Remove unqiue values, where the number of unique values
            #is less than N
            val1 = val[no>=N]
            #Empty False Boolean
            boolean = np.array([False] * len(sum))
            for v in val1:
                #Turn each unqiue value in the sum to True, otherwise 
                #Value is False
                boolean |= (sum == v)
            for ix in index:
                #Delete this new values now.
                boolean = np.delete(np.array(boolean),[ix])
            #Save the boolean into mask and loop across all gridpoints
            mask[:,i,j] = ~boolean
    #NaN all values where value != unique value
    percentile_data[mask] = np.nan

    return percentile_data
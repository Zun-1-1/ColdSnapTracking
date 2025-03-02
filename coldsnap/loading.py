import xarray as xr
import numpy as np
from datetime import datetime
from datetime import date
import calendar
import pickle 
import lzma

def load_by_year(year,step,lists,file_name= '/home/users/vera_os/canari/era5_data/CAO_ST_diag/COM_new_2/',start_year=1940):
    """
    NOTE: DOES NOT WORK BECAUSE FILE PATH IS BOUND TO INDIVIDUAL.
    Load a file if you know the step size, the year you are considering, assuming ERA5 data and you want extended winter period.

    NOTE2: Need to fix the year 1940 -- i.e the first year you are putting in, since here we start in janurary.
    ------------------------------------------------------------------------------------------------------------------------------
    Parameter (IN):

        year: int()
            Year of data
        step: int()
            Step size data was saved in
        list: list of parameter output wanted in the dictionary e.g ["time","com","area_","dpave","uave","vave","countour_no", "points"]
        file_name: string
            Should be directory leading to the file, without the number.npz part
        start_year: int
            The year at which your data start, so we only select the Janurary dates to March here.
    Parameter (OUT):
        data: dict

    """
    days=0
    #Account for initial offset:
    #Time strings
    time_string_1 = '01/01/'+str(start_year)
    time_string_2 = '01/04/' + str(start_year)
    #Time display
    date_format = "%d/%m/%Y"
    #
    a = datetime.strptime(time_string_1, date_format)
    b = datetime.strptime(time_string_2, date_format)
    delta = b - a
    days += delta.days
    #If this is the year we are considering take a note of the start day
    # and the end day
    if year == start_year:
        initial_no = 0
        final_no = days
    for i in range(1941,year+1,1):
        #Otherwise we are looking at the extended winter
        #Time string of the year
        time_string_1 = '01/11/'+str(i-1)
        time_string_2 = '01/04/'+str(i)   
        a = datetime.strptime(time_string_1, date_format)
        b = datetime.strptime(time_string_2, date_format)
        delta = b - a
        #If this is the year we are considering take a note of the start day
        # and the end day
        if i ==year:
            initial_no = days
            #Make sure you account for leap years!
            final_no = days + delta.days #+ calendar.isleap(i+1)
        days += delta.days #+ calendar.isleap(i+1)

    file_1 = int(initial_no/step)
    file_n = int(final_no/step)
    #List of files we need to load and concanate
    file_list = list(range(file_1*step,file_n*step+step,step))
    #What are the indicies for the arrays we load and concatanate
    index_1 = initial_no - file_1 * step
    index_2 = index_1 + (final_no - initial_no)
    
    components = lists
    data = {c:[] for c in components}
    print(file_list)
    for file in file_list:
        print(file)
        #con = np.load(file_name+str(file)+'.npz', allow_pickle=True)
        with lzma.open(file_name+str(file)+'.xz', 'rb') as f:
            con = pickle.load(f)
        for c in components:
            if type(con[c]) is np.ndarray:
                data[c].extend(con[c].tolist())
            else:    
                data[c].extend(con[c])
        print("done")
    for c in components:
       
        
        if c == 'NHC_generation':
            print(len(data[c][10]),len(data['NHC_adv'][10]))
        data[c] = data[c][index_1:index_2]
        if c == 'NHC_generation':
            print(len(data[c][10]),len(data['NHC_adv'][10]))    
    return data
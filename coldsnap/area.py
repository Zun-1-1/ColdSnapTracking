import warnings
import numpy as np
import cv2 as cv
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy.ndimage import gaussian_filter
from shapely import geometry
from shapely.geometry import Point, Polygon
from pyproj import Geod
import xarray as xr
from . import constants as cs


def degrees_into_area(latitude,longitude,degree=True,equal_spacing=True,Npole=True,chunk=False,res_lon=0.25,res_lat=0.25):
    """
    Calculating the area occupied by each point
    ---------------------------------------------------------
    Parameters (IN):
        latitude : xarray or numpy.ndarray, dim={lat}, units degrees or radians
                latitude coordinates
        longitude : xarray or numpy.ndarray, dim={lon}, units degrees or radians
                longitude coordinates
        degrees : bool, optional, default=True
                If not in degrees, will not convert to radians for calculating sinusoidal functions
        equal_spacing : bool, default = True,
                If not True, the results might be incorrect. Will show Warning!
        pole : optional
            If false, the north pole is excluded
        chunk : optional
            If true, then we are only looking at the area of a segment
    Parameters (OUT):
        dA : numpy.ndarray, dim{lat,lon}, [m^2]
            Area occupied by a gridpoint
            
    """

    shape_lat = len(latitude)
    shape_lon = len(longitude)
        
    #Flag if we do not have the complete dataset but we start at pole.
    flag1 = 0
    if latitude[0] == -90:
        flag1 = 1

    #Needs equal spacing
    if equal_spacing == False:
        warnings.warn('Dataset needs to have equal spacing')

    #Needs to be degrees
    if degree == False:
            warnings.warn('Code needs to be in degrees')

    theta = np.array(latitude)
    if chunk:
        phi = longitude
    else:
        phi = np.append(longitude,longitude[0] + 360)

    dphi = phi[1:] - phi[:-1]
    dtheta = theta[:-1] - theta[1:]
    if Npole == False:
        dphi = res_lon
        dtheta = res_lat
    #Upper and lower integration values for non-south pole valies
    if flag1 == 1:
        upper_b = (theta[1:-1] + dtheta[:-1]/2) * cs.degree_to_rad
        lower_b = (theta[1:-1] - dtheta[1:]/2) * cs.degree_to_rad
    #If we dont have neither poles
    elif Npole == False:
        upper_b = (theta + dtheta/2) * cs.degree_to_rad
        lower_b = (theta - dtheta/2) * cs.degree_to_rad        
    #If we dont have south pole
    else:
        upper_b = (theta[1:] + dtheta/2) * cs.degree_to_rad
        lower_b = (theta[1:] - dtheta/2) * cs.degree_to_rad

    #Upper and lower values at the pole
    sin_CD = np.sin(upper_b)
    sin_AB = np.sin(lower_b)


    #Shape of the array
    sin_CD_ = np.zeros((shape_lat,shape_lon))
    sin_AB_ = np.zeros((shape_lat,shape_lon))

    if flag1 == 1:
        for lon in range(shape_lon):
            sin_CD_[1:-1,lon] = sin_CD
            sin_AB_[1:-1,lon] = sin_AB
    if flag1 ==0 :
        if Npole:
            for lon in range(shape_lon):
                sin_CD_[1:,lon] = sin_CD
                sin_AB_[1:,lon] = sin_AB
        else:
            for lon in range(shape_lon):
                sin_CD_[:,lon] = sin_CD
                sin_AB_[:,lon] = sin_AB

    #Not pole
    dA =  cs.radius_e**2 * (sin_CD_ - sin_AB_) * dphi * cs.degree_to_rad
    if Npole == True:
        #Pole
        #print(np.sin(np.pi/2) - sin_CD[0],sin_CD_ - sin_AB_,np.sin(np.pi/2),sin_CD[0])
        N_pole = np.sum( cs.radius_e**2 * (np.sin(np.pi/2) - sin_CD[0]) * dphi * cs.degree_to_rad)
        if flag1 == 1:
            S_pole = np.sum( cs.radius_e**2 * (sin_AB[-1] - np.sin(-np.pi/2)) * dphi * cs.degree_to_rad)
            #Stitch
            #South Pole
            dA[-1,:] = S_pole
    
        #Stitch
        #North Pole
        dA[0,:] = N_pole
    return dA
def area_constraint(variable,dA,uwind,vwind,*args,resolution=0.25,COM=False,COM_var=0,min_val=True):
    """
    AREA CONSTRAINT
    ---------------------------------------------------------
    NOTES:
    The area over which anomolous temperatures occurs. This
    is calculated in pixels in the code, but knowing
    the resolution of the gridpoints allows us to 
    determine the spatial area.

    NOTE 2 : Make sure that uwind, vwind and dp are resampled to
    daily data and only contain extended winter season, 
    if this is the format of percentile data. Use calc.resample(...)

    NOTE 3: The xarrays should put through into the function
    in data_array format, i.e dataset['key']

    NOTE 4: The code takes an average excluding the zero -- 
    since in the case I worked with, it means data is unavailable
    and hence it can skew the average
    ---------------------------------------------------------
    Parameters (IN):

        variable: xarray, dim={time,lat,lon}
            Data you would like to cluster
            
        dA : numpy.ndarray, dim={lat,lon}, [m]
            The area each gridpoint represents

        uwind: xarray, dim={time,lat,lon}, units:[ms**-2]
            Zonal wind at 10m.

        vwind: xarray, dim={time,lat,lon}, units:[ms**-2]
            Meridional wind at 10m.
        
        *args: xarray, dim={time,lat,lon}
            Arguments that you would like to save in the catalogue.

        resolution: float
            Resolution of the dataset

        COM: boolean
            True: We find center of mass of a cluster 
            False: It's center of mass is the position of STA minimum,
            the advantage of this is that it is always inside the cluster
        COM_var: int
            ONLY IF COM = FALSE
            Position of the COM variable in the arguments of the definiton 
            (including the compulsary ones)
        min_val: bool
            When true we find the minimum value of the COM_var, otherwise 
            we find the max
    Parameters (OUT):
        filled_areas : numpy_ndarray, dim={time,contour,lat,lon}, units:[m]
            Contours of each CAO.
    """
    geod = Geod(ellps="WGS84")

   # kernel = np.int8([ 
   # [-1, -1, -1],
   # [-1, +1, -1],
   # [-1, -1, -1],
   # ])

    
    #Create a dictionary
    dictionary = {
    "contour_pnts": [],
    "contour_no": [],
    "cntr": [],
    "time": [],
    "COM": [] ,
    "gridpoints_inside": [],
    "uwind_cntr": [],
    "vwind_cntr": [],
    "uwind_points": [],
    "vwind_points": [],
    "variable": [],
    "area": [],
    "uwind": [],
    "vwind": [],
    }
#    result = 0
#    # Iterating over the Python args tuple
    for x in args:
        if type(x)==xr.core.dataarray.DataArray:
            x_ = x.to_dataset()
            #Get the last key
            key_of_arg = list(x_.keys())[-1]

        dictionary[str(key_of_arg)] = []
    
    #Change to type int8
    #uint8 does not like nans or invalid values so make them zero
    variable_data = np.array(variable).copy()
    variable_data[~np.isfinite(np.array(variable))] = 0
    #Since uint8 of 256 is actually 0 which would lead to 2 clusters
    variable_data[np.isfinite(np.array(variable))] = 255

    variable_data_ = np.array(variable_data,dtype='uint8')
    #Change data into binary
    variable_data_[variable_data_!=0] = 1
    #time
    time = np.array(uwind['time'].dt.strftime("%m/%d/%Y"))
    #Shape
    time_,lat,lon = variable.shape
    #List of all the required argumenst
    req_args = [np.array(variable),np.array(dA),np.array(uwind),np.array(vwind)]
    #Append the args
    for _ in args:
        req_args.append(np.array(_))
    for i in range(time_):
        print(i,"/",time_)
        flag = 0
        area_points,cntr_pts = [],[]
        uwind_points,vwind_points = [],[]
        uwind_cont,vwind_cont = [] ,[]
        #Displacing the datasets.
        percentile_data1 = variable_data_[i]
        #Count the iterations of displacement
        error_ct = 0

        #Case where the pole is anomalous.
        #Displacement based on this won't work unless
        #We discount the pole
        loop_contour = False

        
        Boundary_err = False
        #This tries to account for cyclicly.
        if (percentile_data1[:,-1] == 1).any():
            #We have a discontinuity
            Boundary_err = True
            #The position we are discontinous at:
            position = np.where(percentile_data1[:,-1] == 1)
            percentile_data1 = np.tile(percentile_data1,2)
            #percentile_data1[percentile_data1 == 0] = np.nan
            #Repeat the new array twice along the longitude
            #dA_ = np.concatenate((dA,dA),axis=1)
            #Concatenate dp,uwind, surface temperature and vwind to get the average 
            #STA,percentile,surf_temp,dA,uwind,vwind
            data_concatenated = []
            for counts,arguments in enumerate(req_args):
                #We have already done this one!
                if counts==1:
                    data_concatenated.append(np.tile(dA,2))
                else:
                    ar = arguments[i]
                    data_concatenated.append(np.tile(ar,2))
                
        
        #Convert to RGB colour notation
        image_8bit = percentile_data1 * 255
        #Set threshold
        ret, thresh = cv.threshold(image_8bit, 127, 255, cv.THRESH_BINARY)
        #Find Contours
        contours, _ = cv.findContours(thresh, cv.RETR_EXTERNAL, cv.CHAIN_APPROX_NONE)

        #Depending on whether we have a contour at the edge or not
        if Boundary_err == True:
            filled_areas = np.zeros((len(contours),lat,2*lon))
            #We need a counter for the actual number of contours we are using, since now we probably have about double (due to repeats)
            counter = 0
        else:
            filled_areas = np.zeros((len(contours),lat,lon))
        #If we have no contours break the loop but also append to the list that we have zeros in these areas.
            
        #Making empty arrays
        #For average u_wind and v_wind and dp
        average_matricies = np.zeros((len(req_args),len(contours)))
        #Loop over to find areas of contours
        #We need a new variable to make sure that we are accounting of the 
        #area correctly of each gridpoint.
        area_square = np.zeros((lat,lon))
        #For the COM coordinates
        com_coo = np.zeros((len(contours),2))
        

        for j in range(len(contours)):
            #Where there is discontinuity at the boundary.
            if Boundary_err == True:
                print("Contour at edge")
                #The image now is duplicated so we need to change the area:
                area_square = np.zeros((lat,2*lon))
                img_pl = np.zeros((lat,2*lon,3))
                #For ease define the lat and lon points given from the cv package
                #Find where in longitude it happens and latitude.
                arr = np.array(cv.drawContours(img_pl, contours, j, color=(1,0,0),thickness=-1))
                #When arr has a value, set as True, else False
                mesh = arr[:,:,0] > 0
                #To account that CV has hierarchy and that we might have holes that we need to account for
                mesh &= np.array(percentile_data1,dtype=bool)
                #Account for contours inside holes that we need to remove -- to stop double count 
     #           if np.sum(mesh)>1:
     #               #we have more than just one pixel otherwise this doesnt work
     #               mesh_2 = np.array(mesh,dtype='uint8')
     #               #Convert from binary to 0 255
     # 
     #               mesh_2[mesh_2>0] = 255
     #               neighbors_all_zero = cv.morphologyEx(src=mesh_2, op=cv.MORPH_HITMISS, kernel=kernel)
     #               
     #               mesh = np.array(mesh_2 & ~neighbors_all_zero,dtype=bool)
     #               #Convert back to binary
     #               #mesh[mesh>0] = 1
                #Set a new variable to be dA, and mesh it with the
                #filled contours
                #The area
                area_square[:,:] = data_concatenated[1]
                area_square[~mesh] = 0

                #Where do we have values?
                lat_pos,lon_pos  = np.where(area_square>0)
                for latval in position[0]:
                    #Find where the latitude is equal to the correct value, and the longitude value is at the border
                    coo_bounds_correct = np.logical_and([latval] in lat_pos,[lon-1] in lon_pos)
                    #When the contour is not on the edges of the image and not in the non-duplicated image
                    # or the correct contour that was discontinuous is now present but in the duplicated version.
                    #find the contour and the area.
                    if lon_pos.min() != 0 and ( lon_pos.max() < lon or coo_bounds_correct==True) :
                        #When the value here == lon, then we are at the correct point
                        #Find the grid as an area grid
                        filled_areas[counter,:,:] = area_square
                        #Bug check step, we shouldnt find the same contour twice.
                       # if np.max(np.sum(np.array(filled_areas,dtype=bool),axis=0))>1:
                            #return filled_areas
                       #     raise ValueError("At day: ", i, "And counter",counter," We have accounted for a cold spell twice")
                        
                        outside_lon,outside_lat,u_int,v_int = wind_At_contour(contours[j],uwind[i],vwind[i],resolution)
                        uwind_cont.append(u_int)
                        vwind_cont.append(v_int)
                        #Then find nan mean
                        for number_,args_ in enumerate(data_concatenated):
                            #COM coordinates
                            if number_==COM_var:
                                if COM:
                                    com_coo[counter,0] = ndimage.center_of_mass(filled_areas[counter])[0]
                                    com_coo[counter,1] = ndimage.center_of_mass(filled_areas[counter])[1]
                                else:
                                    
                                    cluster_in_sta = np.where(mesh,args_,np.NaN)
                                    if min_val:
                                        minimum_val = np.where(cluster_in_sta == np.nanmin(cluster_in_sta))
                                    else:
                                        minimum_val = np.where(cluster_in_sta == np.nanmax(cluster_in_sta))
                                    com_coo[counter,0] = minimum_val[0][0]
                                    com_coo[counter,1] = minimum_val[1][0]%(int(360/resolution))
                            if number_==1:
                                #Sum the values up in the square, and convert to km
                                average_matricies[number_,counter] = np.sum(area_square)/(1000*1000)
                    
                            else:
                                
                                cluster_arg = args_[mesh]
                                average_matricies[number_,counter] = np.nanmean(cluster_arg[cluster_arg!=0])

                        
                        #Find the points that are inside the contour.
                        lat_pos_,lon_pos_  = np.where(filled_areas[counter,:,:]>0)
                        coordinates_for_sat = zip(lat_pos_,lon_pos_)
                        lon_pos_ = lon_pos_%(int(360/resolution))
                        coordinates_for_sat = zip(lat_pos_,lon_pos_)
                        #Find uwinds_coods
                        u___ =[]
                        v___ =[]
                        for coordinates_ in coordinates_for_sat:
                            u___.extend([data_concatenated[2][int(coordinates_[0])][int(coordinates_[1])]])
                            v___.extend([data_concatenated[3][int(coordinates_[0])][int(coordinates_[1])]])
                        uwind_points.append(u___)
                        vwind_points.append(v___)
                        area_points.append([lat_pos_.tolist(),lon_pos_.tolist()])
                        cntr_pts.append([outside_lat,outside_lon])
                        if COM==False:
                            check_pos = list(zip(lat_pos_,lon_pos_))
                            #if counter ==73:
                            #    print("here1")
                            #    return minimum_val,mesh,minimum_val[0][0],minimum_val[1][0],check_pos
                            if (minimum_val[0][0],minimum_val[1][0]%(int(360/resolution))) not in check_pos:
                                raise ValueError("COM not in contour at day",i,"and contour",counter)
                        counter +=1
                        #Break we have found the contour this time, with the correct edge discontinouity
                        break
                    #In this case it just overlaps the whole duplicated image so we divide by two the area.
                    elif lon_pos.min() == 0 and lon_pos.max() >= lon-1:

                        #This is a special case that needs extra care. Since this loops all the way round
                        #To avoid any contouring errors we need to fold it back into normal size
                        mesh[:,:lon] |= mesh[:,lon:]
                        #Similarly with area_square, lets repeat it
                        area_square[:,:] = data_concatenated[1]
                        area_square[~mesh] = 0

                        #When the value here == lon, then we are at the correct point
                        #Find the grid as an area gri
                        filled_areas[counter,:lat,:lon] = area_square[:lat,:lon]
                        #Bug check step, we shouldnt find the same contour twice.
                        #if np.max(np.sum(np.array(filled_areas,dtype=bool),axis=0))>1:
                            #if np.sum(np.array(filled_areas,dtype=bool),axis=0)
                        #    raise ValueError("At day: ", i, "And counter",counter," We have accounted for a cold spell twice")
                        outside_lon,outside_lat,u_int,v_int = wind_At_contour(contours[j],uwind[i],vwind[i],resolution)
                        uwind_cont.append(u_int)
                        vwind_cont.append(v_int)

                        #Then find nan mean
                        for number_,args_ in enumerate(data_concatenated):

                            #COM coordinates
                            if number_==COM_var:
                                if COM:
                                    com_coo[counter,0] = ndimage.center_of_mass(filled_areas[counter])[0]
                                    com_coo[counter,1] = ndimage.center_of_mass(filled_areas[counter])[1]
                                else:
                                    
                                    cluster_in_sta = np.where(mesh,args_,np.NaN)
                                    if min_val:
                                        minimum_val = np.where(cluster_in_sta == np.nanmin(cluster_in_sta))
                                    else:
                                        minimum_val = np.where(cluster_in_sta == np.nanmax(cluster_in_sta))
                                    com_coo[counter,0] = minimum_val[0][0]
                                    com_coo[counter,1] = minimum_val[1][0]%(int(360/resolution))
                            if number_==1:
                                #Sum the values up in the square, and convert to km
                                average_matricies[number_,counter] = np.sum(area_square[:lat,:lon])/(1000*1000)
                    
                            else:
                                cluster_arg = args_[mesh]
                                average_matricies[number_,counter] = np.nanmean(cluster_arg[cluster_arg!=0])

                        #Find the points that are inside the contour.
                        lat_pos_,lon_pos_  = np.where(filled_areas[counter,:,:]>0)
                        lon_pos_ = lon_pos_%(int(360/resolution))
                        coordinates_for_sat = zip(lat_pos_,lon_pos_)
                        #Find uwinds_coods
                        u___ =[]
                        v___ =[]
                        for coordinates_ in coordinates_for_sat:
                            u___.extend([data_concatenated[2][int(coordinates_[0])][int(coordinates_[1])]])
                            v___.extend([data_concatenated[3][int(coordinates_[0])][int(coordinates_[1])]])
                        uwind_points.append(u___)
                        vwind_points.append(v___)
                        area_points.append([lat_pos_.tolist(),lon_pos_.tolist()])
                        cntr_pts.append([outside_lat,outside_lon])
                        if COM ==False:
                            check_pos = list(zip(lat_pos_,lon_pos_))
                            if (minimum_val[0][0],minimum_val[1][0]%(int(360/resolution))) not in check_pos:
                                raise ValueError("COM not in contour at day",i,"and contour",counter)
                        counter += 1
                        break
                    #For the case where the contour just intecepts the index = [0] line but is not actually discontinous
                    elif lon_pos.min() == lon:
                        #When the value here == lon, then we are at the correct point

                        #Sum the values up in the square, and convert to km
                        #Find the grid as an area grid
                        filled_areas[counter,:,:] = area_square
                        #Bug check step, we shouldnt find the same contour twice.
                        #if np.max(np.sum(np.array(filled_areas,dtype=bool),axis=0))>1:
                        #    raise ValueError("At day: ", i, "And counter",counter," We have accounted for a cold spell twice")
                        outside_lon,outside_lat,u_int,v_int = wind_At_contour(contours[j],uwind[i],vwind[i],resolution)
                        uwind_cont.append(u_int)
                        vwind_cont.append(v_int)
                        #First make sure that any zeros are nan
                         #Then find nan mean
                        for number_,args_ in enumerate(data_concatenated):
                            #COM coordinates
                            if number_==COM_var:
                                if COM:
                                    com_coo[counter,0] = ndimage.center_of_mass(filled_areas[counter])[0]
                                    com_coo[counter,1] = ndimage.center_of_mass(filled_areas[counter])[1]
                                else:
                                    
                                    cluster_in_sta = np.where(mesh,args_,np.NaN)
                                    if min_val:
                                        minimum_val = np.where(cluster_in_sta == np.nanmin(cluster_in_sta))
                                    else:
                                        minimum_val = np.where(cluster_in_sta == np.nanmax(cluster_in_sta))
                                    com_coo[counter,0] = minimum_val[0][0]
                                    com_coo[counter,1] = minimum_val[1][0]%(int(360/resolution))
                            if number_==1:
                                #Sum the values up in the square, and convert to km
                                average_matricies[number_,counter] = np.sum(area_square)/(1000*1000)
                            else:
                                cluster_arg = args_[mesh]
                                average_matricies[number_,counter] = np.nanmean(cluster_arg[cluster_arg!=0])

                    

                        #Find the points that are inside the contour.
                        lat_pos_,lon_pos_  = np.where(filled_areas[counter,:,:]>0)
                        lon_pos_ = lon_pos_%(int(360/resolution))
                        coordinates_for_sat = zip(lat_pos_,lon_pos_)
                        #Find uwinds_coods
                        u___ =[]
                        v___ =[]
                        for coordinates_ in coordinates_for_sat:
                            u___.extend([data_concatenated[2][int(coordinates_[0])][int(coordinates_[1])]])
                            v___.extend([data_concatenated[3][int(coordinates_[0])][int(coordinates_[1])]])
                        uwind_points.append(u___)
                        vwind_points.append(v___)
                        #Bug check again to see if min_STA points match with the points inside a contour
                        if COM==False:
                            check_pos = list(zip(lat_pos_,lon_pos_))
                            if (minimum_val[0][0],minimum_val[1][0]%(int(360/resolution))) not in check_pos:
                                raise ValueError("COM not in contour at day",i,"and contour",counter)
                        area_points.append([lat_pos_.tolist(),lon_pos_.tolist()])
                        cntr_pts.append([outside_lat,outside_lon])
                        counter+=1
                        break
                    #Otherwise dont find contour.
                    else:
                        pass
                    #Delete the values so it can be overwritten
                    del coo_bounds_correct

                #If there is  no discontinuity at boundary find area normally
            else:
                #Empty array with dim {lat,lon,3}, 3 accounts for 3 RGB variables
                img_pl = np.zeros((lat,lon,3))
                #Areas that are filled by contours of area greater than threshold
                
                #All values are in the RED matrix (RGB)
                arr = np.array(cv.drawContours(img_pl, contours, j, color=(1,0,0),thickness=-1))
                #When arr has a value, set as True, else False
                mesh = arr[:,:,0] > 0
                mesh &= np.array(percentile_data1,dtype=bool)
                #Account for contours inside holes that we need to remove -- to stop double count 
             #   if np.sum(mesh)>1:
             #       #we have more than just one pixel otherwise this doesnt work
             #       mesh_2 = np.array(mesh,dtype='uint8')
             #       #Convert from binary to 0 255
    
   #                 mesh_2[mesh_2>0] = 255
    #                neighbors_all_zero = cv.morphologyEx(src=mesh_2, op=cv.MORPH_HITMISS, kernel=kernel)
     #               
      #              mesh = np.array(mesh_2 & ~neighbors_all_zero,dtype=bool)
       #             #Convert back to binary
        #            #mesh[mesh>0] = 1
                #Set a new variable to be dA, and mesh it with the
                #filled contours
                area_square[:,:] = dA
                area_square[~mesh] = 0

                outside_lon,outside_lat,u_int,v_int = wind_At_contour(contours[j],uwind[i],vwind[i],resolution)
                uwind_cont.append(u_int)
                vwind_cont.append(v_int)
                #Contour in full display
                filled_areas[j,:,:] = area_square
                if np.max(np.sum(np.array(filled_areas,dtype=bool),axis=0))>1:
                    raise ValueError("At day: ", i, "And contour",j," We have accounted for a cold spell twice")

                #Change zeros to nan
                for number_,args_ in enumerate(req_args):
                    #COM coordinates
                    if number_==COM_var:
                        if COM:
                            com_coo[j,0] = ndimage.center_of_mass(filled_areas[j])[0]
                            com_coo[j,1] = ndimage.center_of_mass(filled_areas[j])[1]
                        else:                                  
                            
                            cluster_in_sta = np.where(mesh,args_[i],np.NaN)
                            if min_val:
                                minimum_val = np.where(cluster_in_sta == np.nanmin(cluster_in_sta))
                            else:
                                minimum_val = np.where(cluster_in_sta == np.nanmax(cluster_in_sta))
                            com_coo[j,0] = minimum_val[0][0]
                            com_coo[j,1] = minimum_val[1][0]%(int(360/resolution))
                    if number_==1:
                        #Sum the values up in the square, and convert to km
                        average_matricies[number_,j] = np.sum(area_square)/(1000*1000)
                    else:
                        cluster_arg = args_[i][mesh]
                        average_matricies[number_,j] = np.nanmean(cluster_arg[cluster_arg!=0])
            
                #Find the points that are inside the contour.
                lat_pos_,lon_pos_ = np.where(filled_areas[j,:,:]>0)
                #Make it cyclical
                lon_pos_ = lon_pos_ %(int(360/resolution))
                coordinates_for_sat = zip(lat_pos_,lon_pos_)
                #Find uwinds_coods
                print(len(req_args[2]),len(req_args[2][0][0]))
                #Find uwinds_coods
                u___ =[]
                v___ =[]
                for coordinates_ in coordinates_for_sat:
                    u___.extend([req_args[2][i][int(coordinates_[0])][int(coordinates_[1])]])
                    v___.extend([req_args[3][i][int(coordinates_[0])][int(coordinates_[1])]])
                uwind_points.append(u___)
                vwind_points.append(v___)
                if COM==False:
                    check_pos = list(zip(lat_pos_,lon_pos_))
                    if (minimum_val[0][0],minimum_val[1][0]%(int(360/resolution))) not in check_pos:
                        raise ValueError("COM not in contour at day",i,"and contour",counter)
                area_points.append([lat_pos_.tolist(),lon_pos_.tolist()])
                cntr_pts.append([outside_lat,outside_lon])
        #Save each day


        dictionary['contour_pnts'].append(cntr_pts)
        if Boundary_err == True:
            dictionary['contour_no'].append(counter)
        else:
            dictionary['contour_no'].append(len(contours))
        dictionary['cntr'].append(filled_areas.tolist())
        dictionary['COM'].append(com_coo.tolist())
        dictionary['gridpoints_inside'].append(area_points)
        dictionary['uwind_cntr'].append(uwind_cont)
        dictionary['vwind_cntr'].append(vwind_cont)
        dictionary['uwind_points'].append(uwind_points)
        dictionary['vwind_points'].append(vwind_points)
        counting_ = 0
        for key_no,keys in enumerate(dictionary.keys()):
            if key_no>9:
                dictionary[keys].append(average_matricies[counting_].tolist())
                counting_+=1
        
        #Check that we didn't double count a contour again
        if Boundary_err:
            all = np.array(filled_areas,dtype=bool)
            check_contours = np.sum(all,axis=0)
            check_contours[:,:lon] += check_contours[:,lon:]
            if np.max(check_contours)>1 and np.max(check_contours)<3:
                
                #This can happen when a point is near the border on the left but when it is duplicated that point is inside
                #a closed space, which therefore accounts for this point, but we also account for the point near the border
                #since without duplication there is no reason why to doubt it is another contour.
                
                all = all[:,:,:lon] + all[:,:,lon:]
                loc_x_,loc_y_ = np.where(np.sum(all,axis=0)>1)

                location = list(zip(loc_x_.tolist(),loc_y_.tolist()))

                to_pop = {k : set() for k in dictionary}
                
                for loc_x,loc_y in location:
                    contour_loc = np.where(all[:,loc_x,loc_y]>0)[0]
                    contour_to_get_rid_of = []
                    if np.sum(all[contour_loc[0]])<np.sum(all[contour_loc[1]]):
                        larger_contour = 1
                    else:
                        larger_contour =0
                        
                
                    #Remove this key
                    for keys in dictionary.keys():
                       # print(keys,len(dictionary[keys]),dictionary[keys],contour_loc[(larger_contour+1)%2])
                       # dictionary[keys][0].pop(contour_loc[(larger_contour+1)%2])
                        if keys == 'contour_no':
                            pass
                            #dictionary[keys][i] = dictionary[keys][i]  -1 
                        elif keys == "time":
                            pass
                        else:
                            to_pop[keys].add(contour_loc[(larger_contour+1)%2])
                       # dictionary[keys][i].pop(contour_loc[(larger_contour+1)%2])
                    check_contours[loc_x,loc_y]=1

                for keys, pops in to_pop.items():
                    for index in reversed(sorted(pops)):
                        dictionary[keys][i].pop(index)
                    dictionary['contour_no'][i] = dictionary['contour_no'][i] - len(pops)
                        

            
        #Bug check all counter were accounted for
            if np.sum(np.array(percentile_data1[:,:lon],dtype=bool) ^ check_contours[:,:lon])>0:
                raise ValueError("Not all contours accounted for")
        elif np.sum(np.array(percentile_data1,dtype=bool) ^ np.sum(np.array(filled_areas,dtype=bool),axis=0))>0:
            raise ValueError("Not all contours accounted for")
        
    dictionary['time'].extend(time.tolist())
    
    return dictionary


def xr_add_cyclic_points(da):
    """
    Turning values into cyclic ones from -180 to 180
    NOTE: ASSUMES THAT THE ARRAY IS IN LATITUDE,LONGITUDE POSITION
    """
    
    da_ = np.zeros((len(da['latitude']),len(da['longitude'])+1))
    da_[:,:len(da['longitude'])] = np.array(da)
    da_[:,-1] = np.array(da.sel(longitude =0))
    lon_vals = list(da['longitude'].values) + [da['longitude'].values[0]+360]
    # Generate output DataArray 
    return xr.DataArray(data=da_, 
                           coords=dict(latitude=da.latitude, longitude=lon_vals), 
                           dims=da.dims, 
                           attrs=da.attrs)

def wind_At_contour(contour,uwind,vwind,resolution,NH=True):
    """
    WIND DIRECTION AT THE EDGES OF THE CONTOUR 
    ---------------------------------------------------------
    The function interpolates the wind onto the 4 lines of a gridcell
    and then retains only the coordinates on the outer edge of the contour
             - - - -| - - - - 
            |               |
            |               |
            -       X       -
            |               |
            |               |
             - - - -| - - - - 
                 (x_intr,y_intr)
    ---------------------------------------------------------
    Parameters (IN):
    
    contour: openCV format
        Coordinates of the boundary of the contour
        
    uwind: numpy.ndarray, uniste = [ms-1]
        The zonal wind direction

    vwind: numpy.ndarray, uniste = [ms-1]
        The meridional wind direction

    resolution: float
        Resolution of the dataset

    NH: bool
        If we are only considering part of the hemisphere, coordinates
        stop being cyclical
        
    Parameters (OUT):
    
    outside_lon: list
        The longitude of the interpolated edge of contour
    outside_lat: list
        The latitude of the interpolated edge of contour
    u_int: numpy.ndarray
        The uwind at the corresponding latitude and longitude
    v_int:
        The vwind at the corresponding latitude and longitude 
    """

    #Turn contour points into 
    lat_pt = (90 - np.array(contour[:,:,1].flatten())*resolution).astype(float)
    lon_pt = (-180 + (np.array(contour[:,:,0].flatten()))*resolution).astype(float)
    
    #Finding edge points
    #Top,right,bottom,left
    edge_pts_x = (lon_pt).tolist() +  \
                (lon_pt+resolution/2).tolist() +  \
                     (lon_pt).tolist() +  \
                     (lon_pt-resolution/2).tolist() 
                    

    edge_pts_y = (lat_pt+resolution/2).tolist()  + \
                     (lat_pt).tolist()  +\
                     (lat_pt-resolution/2).tolist()  +\
                     (lat_pt).tolist()
                    
    #print(edge_pts_x,edge_pts_y)
    #Concatenate the points
    xy = list(zip(edge_pts_x,edge_pts_y))
    #The second term is to make the polygon cyclical
    original_xy = list(zip(lon_pt,lat_pt)) + [(lon_pt[0],lat_pt[0])]
    if len(original_xy)>3:
        #Create a polygon
        poly = geometry.Polygon(original_xy)
        #Find a point and determine whether it is within the polygon or not
        contains = np.vectorize(lambda p: poly.contains(Point(p)), signature='(n)->()')
        arr = contains(xy)
        #Pick out the coordinates that are the outer bounds
        outside = np.array(xy)[~arr]
    else:
        #If it is less than 3 gridpoints, then all corners are outside
        outside = xy
    #Split it into x and y instead of xy coordinates
    outside_lon,outside_lat= zip(*outside)
    #Make sure the boundary is accounted for
    outside_lon = np.array(outside_lon)
    outside_lat = np.array(outside_lat)
    #Longitude boundary is cyclical
    outside_lon = np.where(outside_lon>180,outside_lon-360,outside_lon)
    outside_lon = np.where(outside_lon<-180,outside_lon+360,outside_lon)
    #Latitude boundary is somewhat cyclical, since we are only doing NH
    outside_lat = np.where(outside_lat>90,-outside_lat +180 , outside_lat)
    if NH:
        outside_lat = np.where(outside_lat<0, 0, outside_lat)
    else:
        outside_lat = np.where(outside_lat<-90,-outside_lat-180,outside_lat)
    #Make it cyclic at the boundary
    uwind_ = xr_add_cyclic_points(uwind)
    vwind_ = xr_add_cyclic_points(vwind)
    #Need the uwind and vwind of these coordinates
    u_int = uwind_.interp(latitude=list(outside_lat),longitude=list(outside_lon))
    v_int = vwind_.interp(latitude=list(outside_lat),longitude=list(outside_lon))
    uwind_interp_list = []
    vwind_interp_list = []
    #We need to extract the correct points now
    for q in range(len(outside_lat)):
        uwind_interp_list.extend([np.array(u_int[q,q]).tolist()])
        vwind_interp_list.extend([np.array(v_int[q,q]).tolist()])
    return list(outside_lon),list(outside_lat),uwind_interp_list,vwind_interp_list




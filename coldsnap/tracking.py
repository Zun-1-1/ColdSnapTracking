from blist import blist
import numpy as np
from cartopy.geodesic import Geodesic
import pyproj
from shapely import geometry
from shapely.geometry import Point, Polygon
import numpy as np
from datetime import datetime
from shapely.ops import unary_union
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from shapely import geometry,affinity
from shapely.geometry import Point, Polygon,LineString
from shapely.ops import split
from time import perf_counter


geo = Geodesic()
radius_e = 6.371*10**6 #Radius of the Earth average in [m]
g = pyproj.Geod(ellps='WGS84')


def create_polygon(latitude_pos,longitude_pos,u_ave,v_ave,resolution=0.25,advection="Total",fraction=1,max_speed=None):
    """
    Creates a polygon from geodesic equations 
    -------------------------------------------------------------
    Geodesic equations involve calculating distance from
    average speed in the cold spell area and the azimuth angle
    of wind vector.
    -------------------------------------------------------------

    Parameters (IN):

        latitude_pos: float
            Latitudinal position of the cold spells contour
        longitude_pos: float
            Longitudinal position of the cold spells contour
        u_ave: float
            Average zonal wind speed in the area of the cold spell
        v_ave: float
            Average meridional wind speed in the area of the cold spell
        correction: float
            Will be +- to the azimuth angle to account for divergences
            and uncertainties in the advection to the azimuthal angle
        resolution: float
            Resolution of your dataset
        advection: str
            "Zero" : No advection is present t0(t=0)
            "Extreme" : Advection using the position of the cluster at t0(t+24h)
            "Total" : Advection including all the points from t0(t=t0) to t0(t+24h) 
        Fraction: float
            What fraction of "TOTAL" advection. This would be multiplied by the windspeed to either increase or decrease
            advection. Optional.
        max_speed: float, dim:[ms-1]
            The maximum speed the cold spell can move at
    Parameters (OUT):

        poly: shapely function
            Polygon shape that satisfies the distance that can be advected by 
            wind speed and the azimuthal angle +- correction.

        distance: float
            Distance that we expect the cold spell to traverse in a day.

    """
    v_mag = np.sqrt((np.array(u_ave))**2+(np.array(v_ave))**2) #[ms^-1]
    if max_speed is not None:
        v_mag[v_mag>max_speed] = max_speed
    #Forming a geodesic polygon, depending on advection type
    if advection == "Zero":
        distance = v_mag * 0 #1*10**(-10) #[m]
        poly = None
        return poly,distance.max()
    elif advection == "Extreme":
        distance = v_mag * 60*60*24 #[m]
    elif advection=="Total":
        distance = v_mag * 60*60*24*fraction #[m]
    else:
        raise ValueError("Incorrect advection parameter, needs to be: 'Zero', 'Extreme' or 'Total', was given: ", advection)
    
        
    #print(distance/1000,"distance")
    azimuth_ = 180/np.pi * np.arctan2(np.array(u_ave),np.array(v_ave)) #[degrees]


    new_pos_x,new_pos_y = [],[]
    if advection == "Extreme":
        for coo_number in range(len(v_mag)):
            #Advect by geodesics
            lon,lat,azimth = geo.direct([longitude_pos[coo_number],latitude_pos[coo_number]],azimuth_[coo_number],distance[coo_number])[0]
            #Append the new coordinates to a list.
            new_pos_x.append(lon)
            new_pos_y.append(lat)
    if advection == "Total":
        for coo_number in range(len(v_mag)):
            #Advect by geodesics
            lon,lat,azimth = geo.direct([longitude_pos[coo_number],latitude_pos[coo_number]],azimuth_[coo_number],(distance[coo_number]))[0]
            new_pos_x.extend([lon])
            new_pos_y.extend([lat])
            lonlats = g.npts(longitude_pos[coo_number], latitude_pos[coo_number], lon, lat,
                1 + int(distance[coo_number] / 10000))
            x_,y_ = map(list,zip(*lonlats))
            #Append the coordinates joining the two geodesic points
            new_pos_x.extend(x_)
            new_pos_y.extend(y_)

   # print(latitude_pos,longitude_pos,azimuth[0],distance[0])
   # print(new_pos_x,new_pos_y)
    if max(new_pos_x)*min(new_pos_x)<0 and max(new_pos_x) - min(new_pos_x)>180:
       # print("before", new_pos_x)
        new_pos_x.extend([new_pos_x[0].tolist()])
        new_pos_y.extend([new_pos_y[0].tolist()])
      #  print("after",new_pos_x)
        new_pos_x = np.array(new_pos_x)
        new_pos_x[new_pos_x<0] += 360
        #Concatenate everything
        poly_coo = list(zip(new_pos_x,new_pos_y))
        #Make it cyclical
        poly_coo = poly_coo + [poly_coo[0]]
        #Convex hull problem: create the biggest polygon
        hull = ConvexHull(poly_coo)
        #Export the verticies to coordinate system and then form the large polygon
        
        x_vertex = np.array(new_pos_x)[hull.vertices]
        y_vertex = np.array(new_pos_y)[hull.vertices]
        
    
        poly_coo2 = list(zip(x_vertex,y_vertex))
        poly_coo2 = poly_coo2 + [poly_coo2[0]]
        poly_ = geometry.Polygon(poly_coo2)
        line = geometry.LineString([(180,-90), (180,90)])
        # line_ = line.buffer(0.000001) 
        # multi_pgon = poly.difference(line_)
        # multi_pgon.geoms[1]
        result = split(poly_, line)
        if max(list(result.geoms[1].exterior.coords))[0]>180:
            shifted = affinity.translate(result.geoms[1],xoff=-360)
            poly = result.geoms[0].union(shifted)
        
        elif max(list(result.geoms[0].exterior.coords))[0]>180:
            shifted = affinity.translate(result.geoms[0],xoff=-360)
            poly = result.geoms[1].union(shifted)
        else:
            raise ValueError("None are at greater than 180 degrees")

    else:
       # print("before", new_pos_x,new_pos_x[0],new_pos_x[0].tolist(),type(new_pos_x[0]),type(new_pos_x))
        new_pos_x.extend([new_pos_x[0].tolist()])
        new_pos_y.extend([new_pos_y[0].tolist()])
       # print("after",new_pos_x)
        poly_coo = list(zip(new_pos_x,new_pos_y))
        #Make it cyclical
        #poly_coo = poly_coo + [poly_coo[0]]
        
        #Convex hull problem: create the biggest polygon
        hull = ConvexHull(poly_coo)
        #Export the verticies to coordinate system and then form the large polygon
       # print(len(new_pos_x),len(new_pos_y),hull.vertices)
        x_vertex = np.array(new_pos_x)[hull.vertices]
        y_vertex = np.array(new_pos_y)[hull.vertices]
        poly_ = list(zip(x_vertex,y_vertex))
        poly = geometry.Polygon(poly_)
    return poly,distance.max()

def append_extend(cao_dict,data,day,l,append=True,indx=None,name=None,
                  resolution=0.25):
    """
    Appending and inserting to the dictionary, values such as COM, in the correct 
    position
    ---------------------------------------------------------------------------------

    Parameters (IN):

        cao_dict: dict
            Contains all the trajectories
        data: numpy.ndarray 
            Stores all the data about each cold spell
        day: int
            The day of the extended winter period we are in
        l: int
            The contour we are considering in a specific day
        append: bool
            True = Appends to the end of the dictionary. Usually 
            when there is no parent cold spell, and this is a new trajectory
            False = Appends to the correct position in the dictionary using index,
            this will be to the parent index + , i.e appends to the right of it.
        name: str
            Name of the cold spell that we will append/insert into
        indx: (optional) int
            Index as mentioned in append (False) that we insert data into
        resolution: float
            Resolution of the grid being used.
    
    Parameters (OUT):
        cao_dict: dict
            Dictionary with now newly appended values

    """
    #Appending or inserting into structures
    
    input_list = [
                    data['time'][day+1],
                    [90 - resolution*data['COM'][day+1][l][0],
                     -180 + resolution*(data['COM'][day+1][l][1]%(360/resolution))],
                    [round(data['uwind'][day+1][l],3),round(data['vwind'][day+1][l],3)],
                    [day+1,l],
                    round(data['area'][day+1][l]),
                    [data['gridpoints_inside'][day+1][l]],
                    [data['uwind_points'][day+1][l]],
                    [data['vwind_points'][day+1][l]]
                ]

    if [input_list[1][0],input_list[1][1]] in input_list[5]:
        raise ValueError("COM and gridpoints don't match, check day: ", day+1,"and contour: ",l) 
    
    keys_in = ['time','COM','uwind','vwind','contour_no','area','gridpoints_inside','uwind_points','vwind_points']
    not_needed_keys = ['contour_pnts','uwind_cntr','vwind_cntr','Ordinary Duplicates','Track Merge']
    extra_args = set(list(data.keys())) ^ set(keys_in)


    for arguments in extra_args:
        #print(arguments,day+1,l)
        if arguments not in not_needed_keys:
            input_list.extend([round(data[arguments][day+1][l],6)])
            err_check = arguments
    count = 0
    
    if append == True:
        for key in cao_dict[name][name+'.1'].keys():
            if key in not_needed_keys:
                pass
            else:
                cao_dict[name][name+'.1'][str(key)].append(input_list[count])
                count +=1
        if len(cao_dict[name][name+'.1'][err_check]) != len(cao_dict[name][name+'.1']['Date']):
            raise ValueError("Lengths of key mismatch")

    else:
        if indx is not None:
            #Extract the trajectory number
            name_first = name.split('.')[0]
            for key in cao_dict[name_first][name].keys():
                if key in not_needed_keys:
                    pass
                else:
                    cao_dict[name_first][name][str(key)].insert(indx,input_list[count])
                    count+=1
            if len(cao_dict[name_first][name][err_check]) != len(cao_dict[name_first][name]['Date']):
                raise ValueError("Lengths of key mismatch",len(cao_dict[name_first][name][err_check]),len(cao_dict[name_first][name]['Date']),name,name_first,err_check)

        else:
            raise NameError("Did not supply with index position or name position")


    return cao_dict


def contour_overlap_in_next_the_day(cao_dict,data,poly,l,day,latitude_COM,
                                    longitude_COM,distance,parent,name,overlap_day_one,uave,
                                    index=None,area_thres=0,resolution=0.25,overlap=False,advection="Total"):
    
    """
    Finds whether cold spell in the next day corresponds
    to the cold spell in the previous
    ------------------------------------------------------------
    Checks if the contour of cold spell in day+1
    matches with the cold spell in the previous say
    - if so append data to this trajectory
    - Else check another contour
    ------------------------------------------------------------

    Parameters (IN):

        cao_dict: dict
            Contains all the trajectories
        data: numpy.ndarray 
            Stores all the data about each cold spell
        day: int
            The day of the extended winter period we are in
        poly: shapely func
            Polygon covering the area where we expect
            the cold spell to be in day + 1
        l: int
            The contour we are considering in a specific day        
        latitude_COM: float
            Latitudinal position of the cold spells contours
        longitude_COM: float
            Longitudinal position of the cold spells contours
        distance: float
            Distance that we expect the cold spell to traverse in a day.
        parent: bool
            True = The cold spell has a parent cold spell
            False = The cold spell is a new trajectory
        name: str
            Name of the cold spell that we will append/insert into
        uave: float
            Used to calculate distance
        overlap_day_one: only used if overlap=True, and this tracks also clusters
            that are overlaping between two days
        indx: (optional) int
            Index as mentioned in append (False) that we insert data into
        area_thres: (optional) float
            The area threshold each cold spell has to be in order to 
            be added to the dict
        resolution: float
            Resolution of the grid being used.
        advection: str
            Type of advection used, see previous definitions
    
    Parameters (OUT):
        
        cao_dict: dict
            Dictionary with now newly appended values 

    """
    #Area threshold
    if data['area'][day+1][l]<area_thres:
        pass
    else:
        #Latitude and longitude points of COM
        longitude_COM_d = -180 + resolution*(data['COM'][day+1][l][1]%(360/resolution))
        #Make sure the longitude is within -180 and 180 degrees
        if longitude_COM_d>=180:
            raise ValueError("Longitude out of range")
       # while longitude_COM_d>=180:
        #    longitude_COM_d = longitude_COM_d-360
         #   if longitude_COM_d<-180:
          #      raise ValueError("Longitude out of range")
    #Bug test that lat points and lon points coincide with COM
    
#    check_test = list(zip(data['gridpoints_inside'][day+1][l][0],data['gridpoints_inside'][day+1][l][1]))
#    if [int(data['COM'][day+1][l][0]),int(data['COM'][day+1][l][1])] in check_test:
#        raise ValueError("COM and gridpoints don't match, check day: ", day+1,"and contour: ",l) 

    
    #Turn latitude and longitude points inside the contours into degrees
    lat_pt = 90 - np.array(data['gridpoints_inside'][day+1][l][0])*resolution
    lon_pt = -180 + (np.array(data['gridpoints_inside'][day+1][l][1]))*resolution
        
    #Check to speed up the runs, find the smallest distance between a point on the CAO and the COM (day-1), if the distance is much bigger than
    #Expected we can break the loop.
    #longitude_COM_,lon_points = fix_coordinates(longitude_COM,lon_pt,uave)
    #Using Haversine formula
    dist = [] 
    for pts_ in range(len(lat_pt)):
        dist.extend(Haversine_formula(lat_pt[pts_],latitude_COM,lon_pt[pts_],longitude_COM))
    #At zero advection just find overlap
    if advection == "Zero":
        if max(lon_pt) > 180:
            raise ValueError("Longitude greater than 180")
        overlap_day_two = list(zip(data['gridpoints_inside'][day+1][l][0],data['gridpoints_inside'][day+1][l][1]))
        if set(overlap_day_two).intersection(overlap_day_one):
            #Find a point and determine whether it is within the polygon or not
            #If there isn't a parent append it to the newly created one
            if parent == False:
                append_extend(cao_dict,data,day,l,append=True,name=name)
            #If there is a parent append it to the parent
            #(Duplicates are accounted for see aboce)
            if parent == True:
                append_extend(cao_dict,data,day,l,append=False,indx=index,name=name)
    #Break the loop if the distance is too big, as this CAO can never be according to our definition a CAO
    elif min(dist) <= distance:
        #For all the contour points in day + 1 find whether they overlap with day 1 polygon.
        #lon_pt[lon_pt>=180] = lon_pt[lon_pt>=180] -360
        if max(lon_pt) > 180:
            raise ValueError("Longitude greater than 180")
        
        test = [[x,y] for x,y in zip(lon_pt,lat_pt)]
        #Find a point and determine whether it is within the polygon or not
        contains = np.vectorize(lambda p: poly.contains(Point(p)), signature='(n)->()')
        arr = contains(np.array([test]))
        overlap_day_two = list(zip(data['gridpoints_inside'][day+1][l][0],data['gridpoints_inside'][day+1][l][1]))
        if overlap:
            ov = set(overlap_day_two).intersection(overlap_day_one)
        else:
            ov = False
        if arr.any() == True or ov:
            #Find a point and determine whether it is within the polygon or not
            #If there isn't a parent append it to the newly created one
            if parent == False:
                append_extend(cao_dict,data,day,l,append=True,name=name)
            #If there is a parent append it to the parent
            #(Duplicates are accounted for see aboce)
            if parent == True:
                append_extend(cao_dict,data,day,l,append=False,indx=index,name=name)
    return cao_dict

def delete_duplicates(cao_dict,first_name,second_name,index):
    """
    Deletes the index which contains the cold spell present 
    in another pathway, so we don't track it twice
    ----------------------------------------------------------
    
    Parameters (IN):

        cao_dict: dict
            Contains all the trajectories, with duplicates
        first_name: str
            Unique name of the trajectory
        second_name: str
            Unique name of the pathway in
            the trajectory.
        index: int
            Index of the cold spell we are deleting from 
            other pathways and tracks
    
    Parameters (OUT):

            cao_dict: dict
            Contains all the trajectories, without
            duplicates
        
    """
    #Delete the index from all the variables
    for key in cao_dict[first_name][second_name].keys():
        if key == 'Ordinary Duplicates' or key == 'Track Merge':
            pass
        else:
            del cao_dict[first_name][second_name][str(key)][index]
    return cao_dict

def previous_position_finder(cao_dict,first_name,second_name,index):
    """
    When we have a duplicated cold spell in different paths/
    trajectories, we can use the index to find the index of a
    cold spell connecting to this duplicated one. So that once
    the duplicated index is deleted we can remember from which
    index the paths merged into the duplicated one.
    ------------------------------------------------------------
    
    Parameters (IN):

        cao_dict: dict
            Contains all the trajectories, with duplicates
        first_name: str
            Unique name of the trajectory
        second_name: str
            Unique name of the pathway in
            the trajectory.
        index: int
            Index of the cold spell we are deleting from 
            other pathways and tracks
    
    Parameters (OUT):

        position: int
            Index position of the cold spell leading
            to the duplicate
    
    """
    date_format = "%m/%d/%Y"
    one = datetime.strptime(cao_dict[first_name][second_name]['Date'][index], date_format)
    #Pick anything from start to the index
    #Reverse the order
    for position,date_ in enumerate(cao_dict[first_name][second_name]['Date'][:index][::-1]):
        two = datetime.strptime(date_, date_format)
        delta = one - two
        if delta.days == 1:
            #Append the one that goes into the merging point
            position  = len(cao_dict[first_name][second_name]['Date'][:index][::-1]) - position - 1
            #Final check that this is in fact correct
            if (one - datetime.strptime(cao_dict[first_name][second_name]['Date'][position], date_format)).days == 1:
                #print("correct,one,two",one, datetime.strptime(cao_dict[first_name][name]['Date'][position], date_format))
                return position
            else:
                raise NameError("Previous position not found")
    
    raise NameError("Previous position not found")

def Haversine_formula(lat1,lat2,lon1,lon2,radius_e=6.371*10**6):
    """
    Calculates distance between two points on a sphere
    ------------------------------------------------------------
        
    Parameters (IN):

        lat1: float
            Latitudinal position of the first coordinate.
        lat2: float
            Latitudinal position of the second coordinate.
        lon1: float
            Longitudinal position of the first coordinate.
        lon2: float
            Longitudinal position of second coordinate
        radius_e: float, (optional), dim: [m]
            Average radius of the earth 
    
    Parameters (OUT):

        dist: float, dim:[m]
            Distance between two points on a sphere.
    
    """
    #Using Haversine formula

    try:
        dist = 2*radius_e*np.arcsin((np.sin((lat1-lat2)*(np.pi/180)/2)**2 \
                                                    + np.cos(lat1*(np.pi/180)) \
                                                        *np.cos(lat2*(np.pi/180)) \
                                                            *np.sin((lon1-lon2) \
                                                                    *(np.pi/180)/2)**2)**(1/2))
    except RuntimeError:
        print(lat1,lat2,lon1,lon2)
    
    return dist

def name_index_finder(cao_dict,day,i,parent):
    """
    Finds the position and index, of where
    the cold spell on the specific day and of contour
    position i, occurs in the dictionary
    --------------------------------------------------------
        
    Parameters (IN):

        cao_dict: dict
            Contains all the trajectories
        day: int
            The day of the extended winter period we are in
        i: int
            The contour we are considering in a specific day
        parent: bool
            True = The cold spell has a parent cold spell
            False = The cold spell is a new trajectory
    
    Parameters (OUT):

        second_name: str
            The name of the pathway in the trajectory
        index: int
            The index at which to find this specific
            cold spell in the dictionary
        first_name: str
            The name of the trajectory
        parent: bool
            True = The cold spell has a parent cold spell
            False = The cold spell is a new trajectory
    """
    second_name = []
    index = []
    first_name = []
    for m in cao_dict:
    #Find all the different contours
        #If this len > 0 there is at least one parent contour
        for o in range(1, len(cao_dict[m])+1,1):
            if len(cao_dict[m][m+'.'+str(o)]['Contours']) > 0:
                counts = 1
                for n in cao_dict[m][m+'.'+ str(o)]['Contours']:
                    #Check if they match the contour number present
                    if day == n[0] and i == n[1]:
                        #A parent CAO already exits
                        parent = True
                        #Remember which parents it belongs to and the index
                        first_name.append(m)
                        second_name.append(m+'.'+ str(o))
                        index.append(counts)
                    counts += 1
    return second_name,index,first_name,parent

def merging_intra_events(cao_dict,index,name_new,name_og,index_new,day,i):
    """
    Merges trajectories together, takes note of which index/pathway
    merges into what in the 'Track Merge' key. Deletes the old 
    trajectory.
    ---------------------------------------------------------------
    
    Parameters (IN):
    
        cao_dict: dict
            Contains all the trajectories, with duplicates
        index: int
            Index of the cold spell that is merging into another
            trajectory
        name_new: str
            The main trajectory branch name
        name_og: str
            The trajectory branch name that is being merged into the
            main
        index_new: int
            The index (cold spell) of the main trajectory that 
            coalesces with the trajectory that will be merged.
        day: int
            The day of the cold spell
        i: int
            The contour number of cold spell

    Parameter (OUT):
    
        cao_dict: dict
            Contains all the trajectories, with the two 
            trajectories merged

    """

    #When merging we don't need the index of the one following but of the correct one,
    #Hence it is index -1
    index = index - 1
    #Find the COM that connect the two trajectories

    first_name_new = name_new.split('.')[0]
    first_name_og = name_og.split('.')[0]
     #Date during the merge, now find the index that connects to this.
    position = previous_position_finder(cao_dict,first_name_og,name_og,index)
    

    #Bug check
    if position >= len(cao_dict[first_name_og][name_og]['Area']):
        raise ValueError("Position is incorrect")

    #Append the merging one
    cao_dict[first_name_og][name_og]['Track Merge'].append([[cao_dict[first_name_og][name_og]['COM'][position],\
                                                                    cao_dict[first_name_og][name_og]['COM'][index]],[
                                                                    cao_dict[first_name_og][name_og]['Date'][position],
                                                                    cao_dict[first_name_og][name_og]['Date'][index]
                                                                    ],[name_og,name_new], [cao_dict[first_name_og][name_og]['Contours'][position],cao_dict[first_name_new][name_new]['Contours'][index_new-1]]])   # Puts in COM, DATE, the name of the joining pathway, index of the 
                                                                    #position of join in merging branch and the position of the point on the original branch
    #Bug test, the COM of the original on the main branch and the merging trajectory cold spell COM should be equal
    if cao_dict[first_name_og][name_og]['COM'][index] != cao_dict[first_name_new][name_new]['COM'][index_new-1]:
        print("Incorrect",cao_dict[first_name_og][name_og]['COM'][index],cao_dict[first_name_new][name_new]['COM'][index_new-1])
        raise ValueError("Something went wrong...")

    #Check that the main branch name is of lower order than the merging one
    if int(first_name_new) > int(first_name_og):
        #print(first_name_new,first_name_og,int(first_name_new),int(first_name_og))
        raise ValueError("Incorrect first_name")
    #Record the time as well

    contour_no = [day,i]
    if contour_no != cao_dict[first_name_og][name_og]['Contours'][index]:
        raise ValueError("Wrong index")
    #Delete the index from all the variables
    cao_dict = delete_duplicates(cao_dict,first_name_og,name_og,index)
    #Check that we did this correctly, by checking that this index is no longer in the subset.
    if contour_no in cao_dict[first_name_og][name_og]['Contours']:
        raise ValueError("Index deleting wrong elements")

    #Take out all the string merging numbers that need to be put in
    into_new_dict = list(range(len(cao_dict[first_name_og]) + len(cao_dict[first_name_new]) + 1))[len(cao_dict[first_name_new])+1:] 
    #plus one accounts for pythons start at 0

    #Merges the two events
    for no in range(len(cao_dict[first_name_og])):
        #Change the name of the pathways for it to be .n where n-1 is the what the last pathway of the main branch is
        cao_dict[first_name_og][first_name_new+'.'+str(into_new_dict[no])] = cao_dict[first_name_og][first_name_og+'.'+str(no+1)]
        del cao_dict[first_name_og][first_name_og+'.'+str(no+1)] 
    #Merge
    cao_dict[first_name_new] = cao_dict[first_name_new] | cao_dict[first_name_og]
    #Delete all the elements corresponding to the duplicate
  #  for paths in cao_dict[first_name_new].keys():
  #      for number,grdpt in enumerate(cao_dict[first_name_new][paths]['Gridpoints inside']):
 #           #Turn latitude and longitude points inside the contours into degrees
  #          resolution=0.25
  #          lat_pt = 90 - np.array(grdpt[0][0])*resolution
  #          lon_pt = -180 + (np.array(grdpt[0][1]))*resolution
  #          check_test = list(zip(lat_pt,lon_pt))
  #          if (cao_dict[first_name_new][paths]['COM'][number][0],cao_dict[first_name_new][paths]['COM'][number][1]) not in check_test:
  #              print(cao_dict[first_name_new][paths]['COM'][number][0],cao_dict[first_name_new][paths]['COM'][number][1])
  #              print(lat_pt,lon_pt)
                #print(int(cao_dict[first_name_new][paths]['COM'][number][0]),int(cao_dict[first_name_new][paths]['COM'][number][1])
  #              raise ValueError("COM and gridpoints don't match, check day: ", name_new,name_og) 
        
    del cao_dict[first_name_og]
    return cao_dict

def append_merging_to_dict(cao_dict,first_name,second_name,alt_q,initial_pos,index,
                           counter___,track_merge,
                           intra):
    """
    Appends the index that is merging and its previous position, as a key
    so that once the duplicated index is deleted, we can trace back where 
    the two paths/trajectories merged.
    ----------------------------------------------------------------------------------

    Parameters (IN):

        cao_dict: dict
            Contains all the trajectories, with duplicates  
        first_name: str
            Name of the trajectory
        second_name: str
            Name of the pathway
        initial_pos: int
            See previous_position_finder.
        index: list
            Index positions of all the duplicated cold spells that are hence
            that are merging.
        alt_q: int
            Describes which position we are looking at in the index list.
        counter___: int
            Accounts for the fact that if we remove an cold spell from a 
            pathway, then the data in the cold spell at index[alt_q] would
            be shifted by -1.
        track_merge: bool
            True: Trajectories are merging, therefore we append the initial 
            and final values before and after the merge into the key 
            'Tracks Merge'
            False: Pathways are merging, therefore we append the initial
            and final values (which we delete) of the merge to 'Ordinary 
            Duplicates'.
            These keys allow us to understand where a merge was made, 
            but prevent having duplicated data.
        intra: int
            This is required when we have pathways merging not 
            trajectories. We need to know what index represents the main 
            pathway, by the name intra.
    
    Parameters (OUT):

        cao_dict: dict
            Contains all the trajectories, with duplicates, but now the
            merges are documented under the 'Track Merges' and 
            'Ordinary Duplicates' keys.

    """
    #Time at the two points
    int_final_t  = [cao_dict[first_name][second_name[alt_q]]['Date'][initial_pos],\
                    cao_dict[first_name][second_name[alt_q]]['Date'][index[alt_q]+counter___]]
    com_in_final = [cao_dict[first_name][second_name[alt_q]]['COM'][initial_pos],\
                    cao_dict[first_name][second_name[alt_q]]['COM'][index[alt_q]+counter___]]

    #Bug check
    if initial_pos >= len(cao_dict[first_name][second_name[alt_q]]['Area']):
        raise ValueError("Indices are incorrect")

    #Minus one because we want to be at the correct value not one above
    #If we have merging tracks then we need to append the name of the track origin, otherwise
    
    if track_merge:
        key_app = 'Track Merge'
        cao_dict[first_name][second_name[alt_q]][key_app].append([com_in_final, #COM initial and final because of sorting
                                int_final_t #DATE initial and final because it is already sorted
                                ,[second_name[alt_q], #Name of original and name of merge
                                second_name[0]
                                ]
                                ,[cao_dict[first_name][second_name[alt_q]]['Contours'][initial_pos],
                                cao_dict[second_name[0].split('.')[0]][second_name[0]]['Contours'][index[0]-1]] #INDEX NUMBER final at the main branch AND the one joining to it
                                ])

    else:
        key_app = 'Ordinary Duplicates'
        cao_dict[first_name][second_name[alt_q]][key_app].append([com_in_final,#COM initial and final because of sorting
                                int_final_t #DATE initial and final because it is already sorted
                                ,[
                                cao_dict[first_name][second_name[alt_q]]['Contours'][initial_pos],
                                cao_dict[first_name][second_name[intra]]['Contours'][index[intra]-1]] #INDEX NUMBER FINAL in that pathway AND INITIAL what joins onto it index
                                ])
    return cao_dict

def reduction(name,index):
    """
    Reduce the name list and index list
    to an string and integer
    ----------------------------------------------------
    Once we have taken into consideration all
    the duplicates, through mergings of pathways and 
    trajectories, we should be left with one main
    trajectory under name[0] and the index of the
    cold spell should be under index[0] of the 
    dictionary[name]
    ----------------------------------------------------
    
    Parameters (IN):

        name: list
            List of pathway names
        index: list
            List of indices where the specific cold 
            spell appears in the dictionary under name
    
    Parameter (OUT):

        name: int
            Main pathways name
        index: int
            Main pathways index 
    """
    name = name[0]
    index = index[0]
    return name,index

def removing_duplicates(first_name,cao_dict,second_name,index,day,i):
    """
    Removes duplicated cold spells that appear in different
    pathways and trajectories, by taking note of the position
    they occur in the key variables, merging the dictionaries
    of different trajectories, if required, and then deleting
    the duplicated cold spell, but using the index number.
    ------------------------------------------------------------

    Parameters (IN):

        first_name: lists
            List of dictionary trajectories that contain this specific
            cold spell event. I.e if more than one then there is 
            trajectories merging
        cao_dict: dict
            Contains all the trajectories, with duplicates  
        second_name: list
            List of dictionary pathways that contain this specific
            cold spell event. I.e if more than one then there is 
            pathways merging
        index: list
            List of dictionary indexes that point to this specific
            cold spell event in different trajectories, see first_name
            list.
        day: int
            Day number
        i: int
            Contour Number

    Parameters (OUT):

        second_name: str
            Points to the main pathway containing this cold spell
        index: int
            Points to the position of this cold spell in the 
            trajectory, and the pathway of second_name.
        cao_dict: dict
            Contains all the trajectories, with no duplicates  
    """
    if len(first_name) != len(second_name):
        raise ValueError("First name and second name lengths do not match")
    #If len(name) > 1 then there are mergings.
    if len(first_name)>1:
        #If this statement is equal then all the new mergings are from different seperate trajectories
        if len(set(first_name)) != len(first_name):
            if len(set(first_name))>1:
                #MIXED MERGING OF PATHWAYS AND TRAJECTORIES
                track_merge = True                
            else:
                #PURE PATHWAY MERGING
                track_merge = False
            #Each intra_indices is a list of indices inside the same trajectory, but pathway merging
            for intra in range(len(set(first_name))):
                position_ = []
                name_pos = []
                #This finds different seperate trajectories merging
                for pos in set(first_name):
                    position_.append(first_name.index(pos))
                #This finds seperate pathway and trajectory positions
                for nam_pos in set(second_name):
                    name_pos.append(second_name.index(nam_pos))
                #Sort the indices
                position = sorted(position_)
                #from pos index [1] to pos index [2] (exclusive) it is a single trajectory merging. 
                #[2] represents a new trajectory merging/
                #We want to exclude the first point as we will be continuing the trajectory on that one
                #But the other points that are merging to it will be removed.
                #print("intra",intra,"len(position)",len(position))
                if intra == len(position) -1 :
                    #Looping between the final position of merging tracks inside same pathway,
                    #to the end position
                    #We do not need to do len(name)+1, since len(name)+1 would give us an index
                    #we cannot access.
                    intra_indicies = list(range(position[intra]+1,len(second_name)))
                    #print(Counter(name[position[intra]+1:len(name)]))
                else:
                    #Looping between tracks inside the same pathway
                    intra_indicies = list(range(position[intra]+1,position[intra+1]))
                #intra is only correct position if we account for mergings in the index.
                #q will loop around all names that are the same pathway
                #Need to change the index everytime we delete a duplicate
                counter___ = 0
                #Counter for removing name and index from list
                counter_for_index_removal = 0
                for q in intra_indicies:
                    
                    alt_q = q+counter_for_index_removal
                    first_name_ = second_name[alt_q].split('.')[0]                                     
                    #This q doesn't change as name_pos indicies do not get removed
                    #q only changes if there is a pathway i.e 1956_71.1, 1956_71.1 twice
                    #As this would mean there is one less index to look at
                    if q in name_pos:
                        counter___ = -1
                    else:
                        counter___ -= 1
                    #Reset the counter when there is a new name, i.e two pathways that are already 
                    #merged, but there is another merging                                        
                    #Add position of the mergings to the Normal duplicates key values
                    #Find the one connecting to it...
                    initial_pos = previous_position_finder(cao_dict,first_name_,second_name[alt_q],index[alt_q]+counter___)
                    #Check that the indices match from main pathway/trajectory to the one we are merging
                    if cao_dict[first_name_][second_name[alt_q]]['COM'][index[alt_q]+counter___] != cao_dict[first_name[0]][second_name[0]]['COM'][index[0]-1]:
                        raise ValueError("Indices went wrong")
                    #Check that the index never goes to -1
                    if index[0]-1 == -1:
                        raise ValueError("Index values went wrong")
                    #Check that the index never goes to zero
                    
                    #Bug test:
                    #for paths in cao_dict[first_name_].keys():
                    #    for number,grdpt in enumerate(cao_dict[first_name_][paths]['Gridpoints inside']):
                    #        #Turn latitude and longitude points inside the contours into degrees
                    #        resolution=0.25
                    #        lat_pt = 90 - np.array(grdpt[0][0])*resolution
                    #        lon_pt = -180 + (np.array(grdpt[0][1]))*resolution
                    #        check_test = list(zip(lat_pt,lon_pt))
                    #        if (cao_dict[first_name_][paths]['COM'][number][0],cao_dict[first_name_][paths]['COM'][number][1]) not in check_test:
                    #            print(cao_dict[first_name_][paths]['COM'][number][0],cao_dict[first_name_][paths]['COM'][number][1])
                    #            print(lat_pt,lon_pt)
                    #            print("before")
                    #            #print(int(cao_dict[first_name_new][paths]['COM'][number][0]),int(cao_dict[first_name_new][paths]['COM'][number][1])
                    #            raise ValueError("COM and gridpoints don't match, check day: ", first_name_,paths) 
                        
                    cao_dict = append_merging_to_dict(cao_dict,first_name_,second_name,alt_q,initial_pos,index,counter___,track_merge,intra)
                    #Check that we are deleting the correct index, find the unique contour_no
                    contour_no = [day,i]
                    #Check if this matches the index
                    if contour_no != cao_dict[first_name_][second_name[alt_q]]['Contours'][index[alt_q]+counter___]:
                        raise ValueError("We are looking at the wrong day")
                    #Remove name after each time we use it so we should be left with only one
                    cao_dict = delete_duplicates(cao_dict,first_name_,second_name[alt_q],index[alt_q]+counter___)
                    #Check that we are deleting the correct index
                    if contour_no in cao_dict[first_name_][second_name[alt_q]]['Contours']:
                        #This is flagging up because we have this main pathway which will have one contour
                        #Do we have another index after this one? 
                        if alt_q+1 > len(second_name)-1:
                            next_val = False
                        else: 
                            next_val = True
                        #If so is the next index the same as the previous one, because then it is fine 
                        #We still have the same contour number
                        if position[intra] == alt_q-1 and cao_dict[first_name_][second_name[alt_q]]['Contours'].count(contour_no)==1:
                            pass
                        #If there is more than one in the main branch then it is a problem
                        elif next_val:
                            if second_name[alt_q] != second_name[alt_q+1]:
                                raise ValueError("Index is deleting wrong elements")
                            else:
                                pass   
                        else:
                            raise ValueError("Index is deleting wrong elements")

                    del second_name[alt_q], index[alt_q],first_name[alt_q]
                    counter_for_index_removal -=1
                    
            #Seperate merging trajectories 
            if len(set(first_name))>1:
                #len(position) should hold the number of unique trajectories
                #And we pick each one after the first to merge, which is why we -1
                for intra_ in range(len(position)-1):
                    #No need to delete indicies as any deleted indicies would be larger/ further into the list
                    #So should not interfere with any indexing
                    cao_dict = merging_intra_events(cao_dict,index[intra_+1],second_name[0],second_name[intra_+1],index[0],day,i)
                #Now only need to look at the first name and index, after the merge
                second_name = second_name[:1]
                index = index[:1]
                pass
        #If there are just two different trajectories merging
        else:
            #PURE TRAJECTORY MERGES
            for posit,inter in enumerate(second_name[1:]):
                cao_dict = merging_intra_events(cao_dict,index[posit+1],second_name[0],inter,index[0],day,i)
            second_name = second_name[:1]
            #Comment out later
            index = index[:1]
            #print(cao_dict[first_name_])#[str(name[0])]['Contours'][index[0]-1])
            #return print("Stop")
    #Take it out of the list so no need for indexing later
    if len(second_name) == 1:
        second_name,index = reduction(second_name,index)
    else:
        raise ValueError("Name contains more than one value. Name contains:", second_name)
    return second_name,index,cao_dict
    


def naming_structure(names, data,latitude_COM,longitude_COM,u_ave,v_ave,day,i):
    """
    Naming each new track:
    """
    all_keys = list(data.keys())
    keys_in = ['time','COM','area','gridpoints_inside','contour_no','uwind','vwind','uwind_points','vwind_points']
    extra_args = set(all_keys) ^ set(keys_in)
    dict = {names+'.1': {
    #The key ones:
    'Date': blist([data['time'][day]]),
    'COM': blist([[latitude_COM,longitude_COM]]),
    'Average Wind Speed': blist([[round(u_ave,3),round(v_ave,3)]]),
    'Contours': blist([[day,i]]),
    'Area': blist([round(data['area'][day][i])]),
    'Gridpoints inside': blist([[data['gridpoints_inside'][day][i]]]),
    'Uwind points': blist([[data['uwind_points'][day][i]]]),
    'Vwind points': blist([[data['vwind_points'][day][i]]]),
    'Ordinary Duplicates': blist([]), #When two points go into the same point from the track
    'Track Merge': blist([])}} # two different tracks merge
    not_needed_keys = ['contour_pnts','uwind_cntr','vwind_cntr']
    for add_keys in extra_args:
        if add_keys not in not_needed_keys: 
            dict[names+'.1'][add_keys] = blist([round(data[add_keys][day][i],6)])
    #The extra aguments:
#      'NHC generation': blist([round(data['nhcgenave'][day][i],6)]),
#      'NHC flux divergence': blist([round(data['nhcadvave'][day][i],6)]),
#      'Average Surface temperature': blist([round(data['STave'][day][i],2)]),
#      'DP average at 270K': blist([round(data['dp_lave'][day][i],3)]),
#      'DP average at 280K': blist([round(dp_ave,3)]),
            
    return dict

    
    
def track_at_each_contour(cao_year,cao_dict,data,day,i,counter,area_thres=0,resolution=0.25,overlap=False,advection="Total",max_lat=None,fraction=1,max_speed=None):
    """
    Main block of code that tracks all the contours
    ----------------------------------------------------------

    Parameters (IN):

    cao_year: str
        The year we are currently looking at
    cao_dict: dict
        The current dictionary
    data: dict
        The data of contours
    day: int
        Current day
    i: int
        Current contour
    resolution: float
        The resolution of the database
    counter: int
        Counts the trajectories found
    overlap: bool
        True: include overlap as a criteria
        False: Do not include overlap as criteria
    Advection: str
        Type of advection, see create_polygon(...) definition
  #  Full_globe: bool
  #      True: we have the full globe of data
  #      False: we have partial data --> needs max_lat
  #  max_lat: float
  #      The position of the cut off
    Fraction: float
        What fraction of advection. This would be multiplied by the windspeed to either increase or decrease
        advection. Optional.
    max_speed: float, dim:[ms-1]
        The maximum speed the cold spell can move at

    Paramets (OUT):

    cao_dict: dict
        The new altered dictionary with the new contour.
    """
            
    #print("Day:",day,"/", len(data['time'])," Contour:",i,"/",int(data['contour_no'][day]))
    #Find data relevant to the specific contour
    u_ave = data['uwind'][day][i]#[ms^-1]
    v_ave = data['vwind'][day][i]#[ms^-1]

    #Latitude and longitude points of COM
    latitude_COM = 90 - resolution*data['COM'][day][i][0] #[degrees]
    longitude_COM = -180 + resolution*(data['COM'][day][i][1]%(360/resolution))#[degrees]
    latitude_pos = np.array(data['contour_pnts'][day][i][0])
    longitude_pos = np.array(data['contour_pnts'][day][i][1])%(360/resolution)
    #print(longitude_pos,latitude_pos,i)
    #Find the polygon a STA(day+1) should be in and the advected distance
    poly,distance = create_polygon(latitude_pos,longitude_pos,data['uwind_cntr'][day][i],data['vwind_cntr'][day][i],advection=advection,fraction=fraction,max_speed=max_speed)
    overlap_day_one = list(zip(data['gridpoints_inside'][day][i][0],data['gridpoints_inside'][day][i][1]))
    #Currently we don't know if the CAO has a parent
    parent = False

    #Make sure the longitude is within -180 and 180 degrees,
    #This is because when we found "data", we duplicated the image
    #To understand the contours at the interface of the rectangularly
    #projected STA
    if longitude_COM>=180:
        raise ValueError("Longitude out of range")
   # while longitude_COM>=180:
    #    longitude_COM = longitude_COM-360
     #   if longitude_COM<-180:
      #      raise ValueError("Longitude out of range")
    #If we have nothing in the dictionary form a dictionary.
    if len(cao_dict)==0:
        #Give the new CAO a namemyList = [round(x) for x in myList]
        counter += 1
        names = str(data['time'][day][6:]) + '_' + str(counter)
        #print(data['time'][day])
        cao_dict[names] = naming_structure(names, data,latitude_COM,longitude_COM,u_ave,v_ave,day,i)
        
    #Check if we have this contour already.
    else:
        #Find under which names does this contour appear under
        name,index,first_name,parent = name_index_finder(cao_dict,day,i,parent)
        #test this should be equal
        if len(name) != len(first_name):
            raise ValueError("Different number of first names and second names")
        #Delete the duplicate names, there is either pure mergings from the same trajectory
        #pure mergings of different trajectories
        #or there are mergings from same trajectory and from different trajectories into one
        #This should all be put as one trajectory
        if parent:
            name,index,cao_dict = removing_duplicates(first_name,cao_dict,name,index,day,i)
            #Lets just check that we have removed all the duplicated by find the name indices again
            if len(name_index_finder(cao_dict,day,i,parent)[0])>1:
                raise ValueError("Duplicates Found")

        #If there is no parent create a new CAO trajectory.
        else:
            #Give the new CAO a name
            counter += 1
            names = str(data['time'][day][6:]) + '_' + str(counter)
            cao_dict[names] = naming_structure(names, data,latitude_COM,longitude_COM,u_ave,v_ave,day,i)
    #Delete this later
    # if counter > 500:
    #     return cao_dict
    
    #Break the loop at the last day since we have no information about the next day
   # if day==11:  
    if len(data['time']) -1 == day: 
        return cao_dict,counter #np.save(cao_year+'.npy',  cao_dict) 
    else:
        #Loop over all the contours in the next day
        for l in range(0,int(data['contour_no'][day+1]),1):
            if data['area'][day+1][l]>area_thres:
               # if max_lat not in (90 - resolution*np.array(data['gridpoints_inside'][day+1][l][0])):
    
                if parent == True:
                    cao_dict = contour_overlap_in_next_the_day(cao_dict,data,poly,l,day,latitude_pos,
                                longitude_pos,distance,parent,name,overlap_day_one,u_ave, index,area_thres=area_thres,overlap=overlap,advection=advection)
                if parent == False:
                    cao_dict = contour_overlap_in_next_the_day(cao_dict,data,poly,l,day,latitude_pos,
                            longitude_pos,distance,parent,names,overlap_day_one,u_ave,area_thres=area_thres,overlap=overlap,advection=advection)
    return cao_dict,counter



def trajectories(data,filename=None,resolution=0.25,area_thres=0,overlap=False,advection="Total",fraction=1,max_speed = None): #,Full_globe=False,max_lat = None
    """
    Calculating the trajectories of CAO
    ------------------------------------------------------------------------------------------
    NOTE: data must contain the keys formed from contour package for one whole year.
    Other notes:
    The resulting function will give a dictionary of all CAO in the year, with the area of
    impact, dp and wind speed average, contour indexes corresponding to the CAOs, the center
    of mass of the CAO according to the area, and the time of the CAO.
    -----------------------------------------------------------------------------------------

    Parameters (IN):

        data : dict or numpy.ndarray
            A object type numpy array which has a length of the number of days in a year, 
            and then it contains the 'time', 'dpave', 'contours', 'uave','vave','area_',
            'COM'.
        filename: str
            Location and name of the file containing the catalogue
        resolution: float
            The resolution of the grid
        area_thres: float, units: [km**2]
            The threshold area that cold spells have be in order to be placed 
            into the dictionary.
        overlap: bool
            Do you want the clusters that are overlaping in time to also be part of the trajectory?
       # Full_globe: bool
       #     Are you considering the entire globe, including the SH and NH (True) or just part of it (False)?
        Fraction: float
            What fraction of advection. This would be multiplied by the windspeed to either increase or decrease
            advection. Optional.
        max_speed: float, dim:[ms-1]
            The maximum speed the cold spell can move at
            
    Parameters (OUT):

        cao_dict: dict
            Contains all the trajectories of CAO labelled from 1 ... n. If the CAOs are 
            interconnected they are under the parent label.
    
    """
   # if Full_globe ==False and max_lat == None:
   #     raise ValueError("Need max_lat variable, i.e the latitude at which your data cuts off")
    counter = 0
    #Do we already have a dictionary we can append to? 
    cao_dict = dict()
    cao_year =  str(data['time'][-1][6:])
    print("cao_year: ", cao_year)
    if filename==None:
        filename='Catalogue_Advection_'+advection+'_year_'+cao_year+'_areathresh_'+str(area_thres)+'_maxspeed_'+str(max_speed)+'.npy'

    #For loop all contours for all days and all its contours
    for day in range(0,len(data['time']),1):
        print("Day:",day,"/", len(data['time']))
        for i in range(0,int(data['contour_no'][day]),1):
            #Area Threshold:
            if data['area'][day][i]>area_thres:
                #If we are not considering the full globe we need to remove the
                #spells at the edge of the map -- since they are not complete
              #  if not Full_globe:
              #      if max_lat in (90 - resolution*np.array(data['gridpoints_inside'][day][i][0])):

                        #We are at the edge, ignore this contour
              #          pass
              #      else:
                        #Find the trajectory of this contour
                cao_dict,counter = track_at_each_contour(cao_year,cao_dict,data,day,i,counter,area_thres=area_thres,resolution=resolution,overlap=overlap,advection=advection,fraction=fraction,max_speed=max_speed)
                        #If we have the full globe then just find every track directio
               # else:
#                    cao_dict,counter = track_at_each_contour(cao_year,cao_dict,data,day,i,counter,area_thres=area_thres,resolution=resolution,overlap=overlap,advection=advection,max_lat=max_lat,fraction=fraction)
        if len(data['time']) -1 == day: 

            return np.save(filename,  cao_dict)   
    return np.save(filename,  cao_dict)   

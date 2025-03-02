import xarray as xr
import dask
import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

def poly_mask(lons_list,lats_list,resolution):
    """
    Boolean mask that masks anything not inside a area of polygon.
    ------------------------------------------------------------------
    Parameters (IN):

    lons_list: list
            The longitude coordinates of the polygon that describes the area
    lats_list: list
            The latitude coordinates of the polygon that describes tha area
    resolution: float
            The resolution of the dataset

    Parameters (OUT):
    lon,lat: list, {len(lon)} and {len(lat)},
            List of the coordinates inside the area
        
    """
    if max(lons_list)>180:
        base = -180
        lons_list = ((np.array(lons_list) - (base)) % 360) + (base)
    # we do not account for intersections so we manually increase the box by the resolution
    lons_list = np.where(lons_list == np.max(lons_list),lons_list+resolution,lons_list)
    lons_list = np.where(lons_list == np.min(lons_list),lons_list-resolution,lons_list)
    lats_list = np.where(lats_list == np.max(lats_list),lats_list+resolution,lats_list)
    lats_list = np.where(lats_list == np.min(lats_list),lats_list-resolution,lats_list)
    #Create a polygon and check if inside or outside
    XY = list(zip(lons_list,lats_list))
    polygon = Polygon(XY)

    contains = np.vectorize(lambda p: polygon.contains(Point(p)), signature='(n)->()')
    #Longitude and latitude values that are in our array
    lon_vals = np.arange(-180,180,resolution)
    lat_vals = np.arange(90,0-resolution,-resolution)
    lat_lon = np.meshgrid(lon_vals,lat_vals)
    #List of lats and lons that are inside the polygon
    correct_coo = []
    for i in range(len(lat_lon[0])):
        la = list(zip(lat_lon[0][i],lat_lon[1][i]))
        arr = contains(np.array(la))
        if len(np.where(arr==True)[0])>1:
            for pnt in np.where(arr==True)[0]:
                correct_coo.append(la[pnt])

    #long and latitude values
    #booly = np.full((len(lat_vals),len(lon_vals)),False)
    lon = []
    lat = []
    for point in correct_coo:
        lon.append(int((point[0] + 180)/resolution))
        lat.append(int((-point[1]+90)/resolution))
    return lon,lat


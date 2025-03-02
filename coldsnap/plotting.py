from datetime import datetime
import matplotlib.pyplot as plt
from  matplotlib.pyplot import cm
import cartopy.crs as ccrs
import pyproj
import numpy as np
radius_e = 6.371*10**6 #Radius of the Earth average in [m]
"""

"""
def plot_all(cao_dict,year):

    g = pyproj.Geod(ellps='WGS84')

    fig = plt.figure(figsize=(10,10))

    #ax = fig.add_subplot(1,1,1, projection=crs.Orthographic(central_longitude=0,central_latitude=90))
    myProj = ccrs.Orthographic(central_longitude=0, central_latitude=90)
    myProj._threshold = myProj._threshold/50. 

    ax = plt.axes(projection = myProj)

    ax.set_global()
    ax.gridlines()
    ax.set_global()
    resolution=0.25
    #ax.add_feature(edgecolor="tomato")
    #ax.gridlines(edgecolor="tomato")
    ax.coastlines(resolution='50m')
    #plt.show()
    #ax.set_extent([-180, 180, 90, 20], crs=crs.PlateCarree())
    date_format = "%m/%d/%Y"
    #How much colors do we need? We need as much as CAOs we have

    color = cm.rainbow(np.linspace(0, 1, len(cao_dict)))
    start = datetime.strptime(cao_dict[str(year)+'_1']['1']['Date'][0], date_format)

    list_of_areas = []
    #Find maximum area
    for key,trajectory in list(cao_dict.items()):
        for key_pth,path in list(trajectory.items()):
            list_of_areas.append(max(path['Area']))

    #Allowed interval for now
    #allowed = 0

    counter = 0
    number_cm = 0 
    for key,trajectory in list(cao_dict.items()):
        for key_pth,path in list(trajectory.items()):
            #Sorting out the scattering so it doesn't loop around the pole.
            #glossary = '1941_201'
            lat, lon = map(list, zip(*path['COM']))

            #print(path['COM'],key_pth)
            if np.max(lon)-np.min(lon) > 300 and np.min(lon)<-100:
                base = 0
            elif np.max(lon)>180:
                base = 360
            else:
                base = -180
            #check from here
            #Algorithm to connect the dots in the way they are saved
            for i,date_ in enumerate(np.flip(path["Date"])):
                one = datetime.strptime(date_, date_format)
                #print(date_)
                #print(np.flip(path["Date"][:-(i+1)]))
                for no,date_try in enumerate(np.flip(path["Date"][:-(i+1)])):
                    two = datetime.strptime(date_try, date_format)
                    delta = one - two
                    if delta.days == 1:
                        counter +=1


                        com_last = np.flip(path["COM"],axis=0)[i].copy()
                        com_next = np.flip(path["COM"][:-(i+1)],axis=0)[no].copy()

                        com_last[1] = ((com_last[1] - (base)) % 360) + (base)

                        com_next[1] = ((com_next[1] - (base)) % 360) + (base)

                        #Wanting the area on the initial one
                        plt.scatter([com_last[1],com_next[1]], [com_last[0],com_next[0]],
                            transform=ccrs.PlateCarree(),color=color[number_cm],
                            s=np.flip(path["Area"][:-(i+1)],axis=0)[no]/max(list_of_areas)*500)
                        plt.plot([com_last[1],com_next[1]], [com_last[0],com_next[0]],
                            transform=ccrs.PlateCarree(),color=color[number_cm], linewidth=1)

            for merges_pth in path['Ordinary Duplicates']:
                COM_merge_init = merges_pth[0][0]
                COM_merge_lst = merges_pth[0][1]
                #check for the largest lon value
                lon = np.array([COM_merge_init[1],COM_merge_lst[1]])
                #print(path['COM'],key_pth)
                if np.max(lon)-np.min(lon) > 300 and np.min(lon)<-100:
                    base = 0
                elif np.max(lon)>180:
                    base = 360
                else:
                    base = -180
                #Fix the longitude value

                COM_merge_init[1] = ((COM_merge_init[1] - (base)) % 360) + (base)
                COM_merge_lst[1] = ((COM_merge_lst[1] - (base)) % 360) + (base)


                plt.scatter([COM_merge_init[1],COM_merge_lst[1]], [COM_merge_init[0],COM_merge_lst[0]],
                            transform=ccrs.PlateCarree(),color=color[number_cm],
                            s=path['Area'][int(merges_pth[3][0])]/max(list_of_areas)*500)
                plt.plot([COM_merge_init[1],COM_merge_lst[1]], [COM_merge_init[0],COM_merge_lst[0]],
                            transform=ccrs.PlateCarree(),color=color[number_cm], linewidth=1)
            for merges_traj in path['Track Merge']:
                COM_merge_init = merges_traj[0][0]
                COM_merge_lst = merges_traj[0][1]
                #check for the largest lon value
                lon = np.array([COM_merge_init[1],COM_merge_lst[1]])
                #print(path['COM'],key_pth)
                if np.max(lon)-np.min(lon) > 300 and np.min(lon) < -100:
                    base = 0
                elif np.max(lon)>180:
                    base = 360
                else:
                    base = -180
                #Fix the longitude value

                COM_merge_init[1] = ((COM_merge_init[1] - (base)) % 360) + (base)
                COM_merge_lst[1] = ((COM_merge_lst[1] - (base)) % 360) + (base)
                #print(merges_traj[3],merges_traj[3][0],key_pth,merges_traj[3],path['Area'])
                plt.scatter([COM_merge_init[1],COM_merge_lst[1]], [COM_merge_init[0],COM_merge_lst[0]],
                            transform=ccrs.PlateCarree(),color=color[number_cm],
                            s=path['Area'][merges_traj[3][0]]/max(list_of_areas)*500)

                plt.plot([COM_merge_init[1],COM_merge_lst[1]], [COM_merge_init[0],COM_merge_lst[0]],
                            transform=ccrs.PlateCarree(),color=color[number_cm], linewidth=1)
        number_cm +=1
    return
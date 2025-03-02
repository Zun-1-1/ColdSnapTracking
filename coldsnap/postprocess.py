import numpy as np
from datetime import datetime
from collections import Counter
from . import area as ar

def fixing_keys(data,year):
    """
    Change the naming system so we don't have missing 
    trajectory numbers after the merges
    --------------------------------------------------
    Notes: It also changes the 'Merge Trajectory'
    and 'Ordinary Duplicates' keys so that we have the
    initial and final center of mass, but also the
    renamed pathways in the correct naming format
    and the indicies in the path instead of the 
    unique contour and day numbers.
    It is:
    COM1 --> COM2, Date1 --> Date2, 'Name of connecting
    pathway', Index in this pathway and index in the 
    connecting pathway
    --------------------------------------------------

    Parameters (IN):

    data: dict
        Contains all the different trajectories.
    year: int
        The year that this data corresponds to.

    Parameters (OUT):

    data: dict
        Contains all the different trajectory
        consecutively labelled.
    """

    #Get the values out of the dictionary
    kv = list(data.items())
    #Clear the dictionary
    data.clear()
    counter = 1
    #We need to note the path changes.
    previous_path = []
    new_path = []
    for key, value in kv:
        #Rename the keys
        kv_path = list(value.items())
        counter_path = 1
        value2 = {}
        for key_path,value_path in kv_path:
            value2[str(counter_path)] = value_path
            previous_path.append(key_path)#.split('.')[0])
            new_path.append(str(counter_path))
            counter_path +=1
        data[str(year)+'_'+str(counter)] = value2
        counter +=1

    #This is for checks that we are looking at the correct data
    date_format = "%m/%d/%Y"
    #Change namings of pathways in the merges keys now that
    #We changed the pathway names in the keys
    # for key,traj in list(data.items()):
    #     path_count = 0
    #     for k,path in list(traj.items())[::-1]:
    #         for pathies in path['Track Merge']:
            
    #             #Replace Track Merge pathway numbers, if we know what the previous path is
    #             print(pathies[2])
    #             print(previous_path)
    #             if pathies[2][1] in previous_path:
    #                 pathies[2][1] = new_path[previous_path.index(pathies[2][1])]
    #                 pathies[2].pop(0)
    #             else:
    #                 #If its more than 2 different trajectories involved find in the pathways
    #                 # of the trajectory 
    #                 Break = False
    #                 for keys,paths in list(traj.items())[::-1][path_count:]:
    #                     if Break:
    #                         break
    #                     for merges in paths['Track Merge']:

    #                         #If there is that pathway that links to another main pathway
    #                         #See supplementary notes, we can find the original pathway from that.
    #                         if pathies[2][1] in merges[2][0]:
    #                             pathies[2][1] = keys
    #                             pathies[2].pop(0)
    #                             Break = True
    #                             break
    #                 #Otherwise check for the trajectory merge
    #                 if not Break:
    #                     for keys,paths in list(traj.items())[::-1][path_count:]:
    #                         for merges in paths['Track Merge']:
    #                             if pathies[2][1].split('.')[0] in merges[2][0].split('.')[0]:
    #                                 #Find the correct pathway by looking at the pathway numbers of other trajectories
    #                                 #with the same trajectory number
    #                                 difference_in_paths = int(merges[2][0].split('.')[1])- int(pathies[2][1].split('.')[1])
    #                                 key_of_merged = int(keys)
    #                                 correct_path = key_of_merged-difference_in_paths
    #                                 pathies[2][1] = str(correct_path)
    #                                 pathies[2].pop(0)
    #                                 Break = True
    #                                 break
    #                         if Break:
    #                             break
    #                 if not Break:
    #                     raise ValueError("Not found correct Track Merge",pathies)
    #        path_count +=1

    #Change the [contour number,day] (unique variables of where each track
    # merges) to index position and the key before to the new key
    for key,traj in list(data.items()):
        path_count = 0
        for k,path in list(traj.items()):
            #Track Merge
            for pathies in path['Track Merge']:
                initial_day = pathies[3][0]
                final_day = pathies[3][1]
                start_index = path['Contours'].index(initial_day)
                #There should be only one unique contour number 
                #Anywhere in the pathways so we can just
                #Loop over and find it
                contour_found = False
                check_no = 0
                for keys,paths in list(traj.items()):
                    if final_day in paths['Contours']:
                        pathies[3][0] = start_index
                        pathies[3][1] = paths['Contours'].index(final_day)
                        pathies[2][1] = str(keys)
                        pathies[2].pop(0)
                        contour_found = True
                        check_no +=1
                        #Check that everything is still correct with the days being one
                        #Apart
                        if (datetime.strptime(paths['Date'][pathies[3][1]], date_format) - \
                            datetime.strptime(path['Date'][start_index], date_format)).days != 1:
                            raise ValueError("Not the correct pairing")
                #Bug checks
                if check_no >1:
                    raise ValueError("Contour number not unique")
                if not contour_found:
                    raise ValueError("Missing Contour numbers")

    #Do the same for the ordinary duplicates
    for key,traj in list(data.items()):
        path_count = 0
        for k,path in list(traj.items()):
            #Track Merge
            for pathies in path['Ordinary Duplicates']:
                initial_day = pathies[2][0]
                final_day = pathies[2][1]
                start_index = path['Contours'].index(initial_day)
                #There should be only one unique contour number 
                #Anywhere in the pathways so we can just
                #Loop over and find it
                contour_found = False
                check_no = 0
                for keys,paths in list(traj.items()):
                    if final_day in paths['Contours']:
                        pathies[2][0] = start_index
                        pathies[2][1] = paths['Contours'].index(final_day)                        
                        #Insert the path containing this key
                        pathies.insert(2,[keys])
                        #Check that everything is still correct with the days being one
                        #Apart
                        if (datetime.strptime(paths['Date'][pathies[3][1]], date_format) - \
                            datetime.strptime(path['Date'][start_index], date_format)).days != 1:
                            raise ValueError("Not the correct pairing")
                        contour_found = True
                        check_no +=1
                #Bug checks
                if check_no >1:
                    raise ValueError("Contour number not unique")
                if not contour_found:
                    raise ValueError("Missing Contour numbers")


    counter = 1
    #Check that all the keys are now corectly labelled
    for key in data.keys():
        if str(year)+'_'+str(counter) != key:
            #print(key)
            raise ValueError
        counter +=1
    return data


def threshold_process(data,year,area_thresh=0.8*10**6):
    """
    Place thresholds on the cold spells that occur
    --------------------------------------------------

    Parameters (IN):

    data: dict
        Contains all the different trajectories.
    year: int
        The year that this data corresponds to.
    area_thresh: float, dims: [km**2]
        The area threshold
    temp_thresh: float, dims: [K]
        The temperature threshold

    Parameters (OUT):

    data: dict
        Contains all the different trajectory
        consecutively labelled.
    """
    #Get the values out of the dictionary
    kv = list(data.items())
    #Clear the dictionary
    data.clear()
    counter = 1 
    for key,value in kv:
        Area_thres_met = False
        for k_path,v_path in list(value.items()):
            #The area needs to be bigger than what is defined and the temperature needs to be
            #lower than what is defined
            if max(v_path['Area']) > area_thresh: #and min(v_path['Average Surface temperature']) < temp_thresh:
                Area_thres_met = True
                break
        #If those conditions were satisfied in any of the pathways then let it remain in the dictionary
        if Area_thres_met == True:
            data[str(year)+'_'+str(counter)] = value
            counter +=1
    #Return the new reprocessed dictionary
    return data



def paths_to_keep(data,key,no,area,lsm,quasi_thresh =30,keep=True,resolution=0.25):
    """
    Finds the paths that contain quasi-stationary behaviour.
    ------------------------------------------------------------------
    This is done by checking how many times we overlap the same 
    gridcell in a single pathway by using the quasi_thresh metric.
    ------------------------------------------------------------------
    Parameters (IN):
    
    data: dict
        Contains all the different trajectories.
    key: str
        The key we are considering
    no: str
        The path we are considering
    area: float, dims:[km**2]
        The area threshold
    quasi_thresh: (optional) int
        The number of times we allow overlap over same gridcell
        in same pathway until we discard that pathway for being
        quasi-stationary.

    Parameters (OUT):

    list(set(paths_kept)): list
        Altered paths that should be kept
    list(set(paths_to_remove)):
        Altered paths that should be removed
    data: dict
        Contains all the different trajectories, not in
        Consecutive order

    quasi_stationary: bool
        quasi_stationary = True
            ==> We found quasi-stationary behaviour 
        quasi_stationary = False 
             ==> No quasi-stationary behaviour
        
    """
    lat = np.arange(90,-0.25,-resolution)
    lon = np.arange(-180,180,resolution)
    all_coords = []
    quasi_stationary = False
    paths_kept = []
    #Find coordinates
    for q,each_coo in enumerate(data[key][no]['Gridpoints inside']):
        all_coords.extend(list(zip(*each_coo[0])))
    Passing = False
    if max(Counter(all_coords).values())>quasi_thresh:
        matrix = np.zeros((int(90/resolution+1),int(360/resolution)))
        #Check the land to sea ratio
        for coo in all_coords:
            matrix[int(coo[0]),int(coo[1])] +=1
        area_globe = ar.degrees_into_area(lat,lon)
        #Sea area:
        sea_ = np.where(np.array(lsm)==0,matrix,0)
        #Land
        land_ = np.where(np.array(lsm)>0,matrix,0)
        #Total area_globe:
        
        if np.sum((sea_*area_globe))/np.sum((matrix*area_globe)) > 0.8:
            #print("PASS",key,np.sum((sea_*area_globe))/np.sum((matrix*area_globe)))
            Passing = True

    #See if the coordinates repeat in a path
    if Passing:
        quasi_stationary = True
        #If they do delete that path and any that relates to it that does not fit the area criteria
        #1. check which pathways meet area criteria / or change this to a criteria of your choice
        #Change 'area' to the str of choice.
        for other_paths in range(1,len(data[key])+1):    
            if max(data[key][str(other_paths)]['Area'])> area and str(other_paths) != str(no):
                paths_kept.extend([str(other_paths)])
                #Check what other pathways it is connected to as they can be kept as well
                for i in data[key][str(other_paths)]['Track Merge']:
                    if i[2][0] != no:
                        paths_kept.extend([str(i[2][0])])
                for i in data[key][str(other_paths)]['Ordinary Duplicates']:
                    if i[2][0] != no:
                        paths_kept.extend([str(i[2][0])])
        #2. Check if any of the paths are part of paths_kept, or merge to it.
        #We repeat this until we get all the paths present
        number = len(set(paths_kept)) - 1
        while len(set(paths_kept)) != number:
            number = len(set(paths_kept))
            for other_paths in range(1,len(data[key])+1):    
                #We only consider the paths not already kept
                if str(other_paths) not in paths_kept:
                    for i in data[key][str(other_paths)]['Track Merge']:
                        #If it merges to one of the paths we keep then place it into paths_kept
                        if i[2][0] in paths_kept and str(other_paths) != no:
                            #And take all the paths from this key
                            paths_kept.extend([str(other_paths)])
                    for i in data[key][str(other_paths)]['Ordinary Duplicates']:
                        #If it merges to one of the paths we keep then place it into paths_kept
                        if i[2][0] in paths_kept and str(other_paths) != no:
                            #And take all the paths from this key
                            paths_kept.extend([str(other_paths)])
                #Make sure any newly taken in paths are checked for pathways
                else:
                    for i in data[key][str(other_paths)]['Ordinary Duplicates']:
                        if i[2][0] != no and i[2][0] not in paths_kept:
                            paths_kept.extend([str(i[2][0])])    
                    for i in data[key][str(other_paths)]['Track Merge']:
                        if i[2][0] != no and i[2][0] not in paths_kept:
                            paths_kept.extend([str(i[2][0])])
                           # print(str(i[2][0]),"here4",other_paths)
        #3.From the paths we keep, delete the pathway we are not considering anymore
        
        
        for other_paths in range(1,len(data[key])+1):  
            deleting = []
            for i,path__ in enumerate(data[key][str(other_paths)]['Track Merge']):
                if path__[2][0] == no:
                    deleting.extend([i])
            #Now remove the bad track_merges
            for each_no in sorted(deleting,reverse=True):
                del data[key][str(other_paths)]['Track Merge'][each_no]
            deleting = []
            for i,path__ in enumerate(data[key][str(other_paths)]['Ordinary Duplicates']):
                if path__[2][0] == no:
                    deleting.extend([i])
                    
            #Now remove the bad track_merges
            for each_no in sorted(deleting,reverse=True):
                del data[key][str(other_paths)]['Ordinary Duplicates'][each_no]
        #Delete all track merges and ordinary duplicates from the special pathway
        if keep:
            deleting = []
            for i,path__ in enumerate(data[key][str(no)]['Track Merge']):
                deleting.extend([i])
            #Now remove the bad track_merges
            for each_no in sorted(deleting,reverse=True):
                del data[key][str(no)]['Track Merge'][each_no]
            deleting = []
            for i,path__ in enumerate(data[key][str(no)]['Ordinary Duplicates']):
                deleting.extend([i])
            for each_no in sorted(deleting,reverse=True):
                del data[key][str(no)]['Ordinary Duplicates'][each_no]         
    return list(set(paths_kept)),data,quasi_stationary,no


def filter_quasi_stat(data,key,no,area,lsm,quasi_thresh=30,keep=True):
    """
    Filtering process of Quasi stationary behaviour
    ------------------------------------------------------------------
    Now that we know the paths to remove and paths to keep we need to
    alter the dictionary.
    ------------------------------------------------------------------

    Parameters (IN):
    
    data: dict
        Contains all the different trajectories.
    key: str
        The key we are considering
    no: str
        The path we are considering
    area: float, dims:[km**2]
        The area threshold
    quasi_thresh: (optional) int
        The number of times we allow overlap over same gridcell
        in same pathway until we discard that pathway for being
        quasi-stationary.
    keep: (optional)
        True: Keeps the special pathways as a special trajectory, with a flag that 
        it is special. These key trajectories will have an 's' in front of their key number
        False: Removes these pathways completely.

    Parameters (OUT):

    data: dict
        The dictionary now seperated the paths that contained
        quasi-stationary.
    quasi_stationary: bool
        quasi_stationary = True
            ==> We found quasi-stationary behaviour 
        quasi_stationary = False 
             ==> No quasi-stationary behaviour
    list(new_key.keys()) : 
        quasi_stationary = False
             ==> 0
        quasi_stationary = True :
            ==> Keys of new altered trajectories when
        
    """
    year = key[:4]
    
    paths_kept,new_data,quasi_stationary,no = paths_to_keep(data,key,no,area,lsm,quasi_thresh =quasi_thresh,keep=keep)
    
    if quasi_stationary:
        print("Quasi stationary cold spells found on key", key)
        if len(paths_kept)>0:    
            key_no = int(list(data.keys())[-1].split("_")[1])
            new_key = dict()
            new_key[str(year)+'_'+str(key_no+1)] = { '1' : []}
            int_paths = [eval(i) for i in paths_kept]
            connection_pth = []
            connection_org = []
            for path in sorted(int_paths):
            
                check = ['Track Merge', 'Ordinary Duplicates']
                #connection = []
                for key_merges in check:
                    for merge in new_data[key][str(path)][str(key_merges)]:
                        if path != int(merge[2][0]):
                            #If the path it merges to had a higher number we should switch the connections for the
                            #Algorithm to still work
                            if path<int(merge[2][0]):
                                connection_pth.extend([int(merge[2][0])])
                                connection_org.extend([path])
                            else:
                                connection_pth.extend([path])
                                connection_org.extend([int(merge[2][0])])
                        
            count = 1
            con = connections(connection_pth,connection_org,int_paths)
            groups_done = []
            #We need to know which path is going to which key_no
            key_no_counter = []
            path_count_og = []
            path_count_new = []
            
            for path in sorted(int_paths):
                if count == 1:
                    new_key[str(year)+'_'+str(key_no+1)]['1'] = new_data[key][str(path)]
                    #But now we need to change the merges and duplicate paths so that they match the newly labeled path
                    groups_done.extend([con[path]])
                    key_no_counter.extend([str(year)+'_'+str(key_no+1)])
                    path_count_og.extend([str(path)])  
                    path_count_new.extend(['1']) 
                    count +=1
                else:
                    if con[path] in groups_done:
                        index_no = groups_done.index(con[path])
                        key_on = key_no_counter[index_no]
                        number_of_paths = len(new_key[key_on].keys())
                        new_key[key_on][str(number_of_paths+1)] = new_data[key][str(path)]
                        #print(path)
                        groups_done.extend([con[path]])
                        key_no_counter.extend([str(year)+'_'+str(key_no+1)])
                        path_count_og.extend([str(path)])  
                        path_count_new.extend([str(number_of_paths+1)]) 
                    else:
                        
                        key_no +=1
                        new_key[str(year)+'_'+str(key_no+1)] = { '1' : new_data[key][str(path)]}
                        #But now we need to change the merges and duplicate paths so that they match the newly labeled path
                        groups_done.extend([con[path]])
                        key_no_counter.extend([str(year)+'_'+str(key_no+1)])
                        path_count_og.extend([str(path)]) 
                        path_count_new.extend(['1'])  




            #Now we need to delete the previous key and append it to the original datafile
            if keep and max(data[key][str(no)]['Area'])> area:
                key_no +=1  
                new_key[str(year)+'_'+str(key_no+1)+'_s'] = { '1' : new_data[key][str(no)]}
            elif keep:
                print("QS too small")
            del data[key]

            for all_new_keys in new_key:
                data[str(all_new_keys)] = new_key[str(all_new_keys)] 
           #     all_others = np.zeros((361,1440))
           #     for paths1 in new_key[str(all_new_keys)]:
           #         for number,points in enumerate(new_key[str(all_new_keys)][str(paths1)]['Gridpoints inside']):
           #             coo = list(zip(points[0][1],points[0][0]))
                        #The unique gridpoints in a dictionary
                        
           #             for unique_gridpoints in set(coo ):
           #                 all_others[int(unique_gridpoints[1]),int(unique_gridpoints[0])] +=1 
           #     plt.imshow(all_others)
           #     plt.show()
           #     print(all_new_keys)
            new_keys_list = list(new_key.keys())
            #We need to account for all the new key_merges and unmerges
            for keys_new in new_key:
                for paths_new in new_key[str(keys_new)]:
                    for key_merges in check:
                        for merge_no,merge in enumerate(new_key[str(keys_new)][str(paths_new)][key_merges]):

                            index_old = path_count_og.index(str(merge[2][0]))
                            new_key[str(keys_new)][str(paths_new)][key_merges][merge_no][2] = [str(path_count_new[index_old])]
                            
        else:
            if keep and max(data[key][str(no)]['Area'])> area:
                key_no = int(list(data.keys())[-1].split("_")[1])

                
                data[str(year)+'_'+str(key_no+1)+'_s'] = { '1' : new_data[key][str(no)]}
                print(str(year)+'_'+str(key_no+1)+'_s')
            elif keep:
                print("QS too small")
            del data[key]
            quasi_stationary = True
            new_keys_list = [] 
        
        return data,quasi_stationary,new_keys_list
    else:
        return data,quasi_stationary,0
        
def QS_filter_full(data,area,lsm,year,quasi_thresh=30,keep=True):
    """
    This is to filter all data instead of just a specific key and path
    ------------------------------------------------------------------------
    Parameters (IN):
    
    data: dict
        Unfiltered data
    area: float, dims:[km**2]
        The area threshold
    quasi_thresh: (optional) int
        The number of times we allow overlap over same gridcell
        in same pathway until we discard that pathway for being
        quasi-stationary.
    keep: (optional)
        True: Keeps the special pathways as a special trajectory, with a flag that 
        it is special. These key trajectories will have an 's' in front of their key number
        False: Removes these pathways completely.
    lsm:  dims:[lat,lon]
        The land sea mask, used to evaluate the relative position of the cold spell over
        sea and over land. Quasi stationary that we are looking for are over sea.
        If 95% of the area is over the sea, then consider it QS.
        
    Parameters (OUT):

    data_new: dict
        Filtered data
    """
    if lsm is None:
        raise ValueError("Input a landmask")
    all_values = list(data.keys())
    for keys in all_values:
        if keys[-1] == 's':
            pass
        else:
            for paths in data[str(keys)]:
                data,quasi_stationary,no_of_new_keys = filter_quasi_stat(data,keys,paths,area,lsm,quasi_thresh=quasi_thresh,keep=keep)
                if quasi_stationary == True:
                    all_values.extend(no_of_new_keys)
                    break
    #year = keys[:4]
    data_new = {}
    total_number_of_keys = len(data.keys())
    counter =1
    for keys in data:
        #print("From ", keys, "To ",str(year)+'_'+str(counter))
        if keys[-1] == 's':
            data_new[str(year)+'_'+str(counter)+'_s'] = data[str(keys)]
        else:
            data_new[str(year)+'_'+str(counter)] = data[str(keys)] = data[str(keys)]
       # 
       # all_others = np.zeros((361,1440))
       # for paths1 in data_new[str(year)+'_'+str(counter)]:
       #     for number,points in enumerate(data_new[str(year)+'_'+str(counter)][str(paths1)]['Gridpoints inside']):
       #         coo = list(zip(points[0][1],points[0][0]))
                #The unique gridpoints in a dictionary
                
       #         for unique_gridpoints in set(coo ):
       #             all_others[int(unique_gridpoints[1]),int(unique_gridpoints[0])] +=1 
       # plt.imshow(all_others)
       # plt.show()
        counter+=1
    return data_new
    
def connections(keys_from,keys_to,keys_ind):
    connected = set()
    for k1, k2 in zip(keys_from+keys_ind, keys_to+keys_ind):
        to_add = tuple(sorted([k1, k2]))
    
        matching_cons = []
        for con in connected:
            if to_add[0] in con or to_add[1] in con:
                matching_cons.append(con)
    
        if not matching_cons:
            connected.add(to_add)
        else:
            combined = []
            for con in matching_cons:
                connected.remove(con)
                combined += list(con)        
            new_con = tuple(set(sorted([*to_add, *combined])))
            connected.add(new_con)


    mapping = {}
    for i, con in enumerate(connected):
        for c in con:
            mapping[c]=i
    return mapping

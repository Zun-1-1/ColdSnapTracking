from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import numpy as np 

def inside_pts(data):
    """
    Check the number of COMs that are inside the cold air anomaly and the number that are outside
    --------------------------------------------------------------------------------------------------------------------------------------------
    Parameters (IN):
        data: numpy.ndarray, object type
            Contains dictionary of contour_no, area and contour_pnts.
    
    Parameters (OUT):
        inside: float
            Number of COM inside the polygon
        outside: float
            Number of COM outside the polygon
    """
    inside = []
    outside = []
    for day in range(len(data['area'])):
        print(day)
        for contour in range(data['contour_no'][day]):
            XY = list(zip((-180 + np.array(data["contour_pnts"][day][contour][1])*0.25), (90 - np.array(data["contour_pnts"][day][contour][0])*0.25)))
            #print(day,contour)
            if len(XY)>3:
                polygon = Polygon(XY) # create polygon
                point = Point( -180 + data['com'][day][contour][0][1]*0.25,90 - data['com'][day][contour][0][0]*0.25) # create point

                #Point touches or is inside the polygon
                if polygon.touches(point) or point.within(polygon):
                    inside.append(1)
            
                else:
                    #If it is slightly off the border it still counts
                    if polygon.exterior.distance(point)<0.001:
                        inside.append(1)
                    #Otherwise it is not inside
                    else:
                        #plt.plot(*polygon.exterior.xy)
                        #plt.scatter(-180 + data['com'][day][contour][0][1]*0.25,90 - data['com'][day][contour][0][0]*0.25)
                        #print(contour,day)
                        outside.append(1)
            #If there is less than 3 points, then it must be inside
            else:
                inside.append(1)
    return sum(inside),sum(outside)



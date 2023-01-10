import shapely
import shapely.ops
import numpy as np
import matplotlib.pyplot as plt
from elastica.external_forces import NoForces

def read_from_file(file_name):
    shapes=[]
    with open(file_name) as f:
        lines=f.readlines()
    for line in lines:
        myarray = np.fromstring(line, dtype=float, sep=',')
        myarray = np.reshape(myarray,(-1,2))
        shapes.append(shapely.Polygon(myarray))
    return shapes


class ObstacleForce(NoForces):
    def __init__(self,data):
        self.shapes=[]
        for d in data:
            d = np.reshape(d,(-1,2))
            self.shapes.append(shapely.Polygon(d))
        self.tree = shapely.STRtree(self.shapes)
        self.k = 100000.0


    def apply_forces(self, system, time: np.float64 = 0.0):

        positons = system.position_collection[0:2, :].T
        pts = [shapely.Point(p[0:2]) for p in positons]
        indice = self.tree.nearest(pts)

        for i in range(1,len(pts)-1):
            p1,p2 = shapely.ops.nearest_points(pts[i],self.shapes[indice[i]])
            distance = shapely.distance(p1,p2)
            direction = np.array([p1.x-p2.x,p1.y-p2.y,0.0])
            direction = (direction/np.linalg.norm(direction)).reshape(3,)
            system.external_forces[:,i]=self.k/distance/distance*direction



if __name__=='__main__':
    shapes = read_from_file('../data/map/refined_obstacle.txt')
    for poly in shapes:
        x,y = poly.exterior.xy
        plt.plot(x,y,'-g')

    tree = shapely.STRtree(shapes)

    pt = [shapely.Point([[0,1.2]]),shapely.Point([[-2,6]])]
    index = tree.nearest(pt)

    poly = shapes[index[0]]
    x,y = poly.exterior.xy
    plt.plot(x,y,'-r')

    print(shapely.distance(shapely.Point([-1,4]),poly))

    p1,p2 = shapely.ops.nearest_points(shapely.Point([-1,4]),poly)
    print(p2.x)


    plt.show()
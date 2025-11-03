import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
from Shapes import Point,Sphere
# Code for the class that we are building
class simulation_box:
    """ Class that represent the simulation box of the project"""
    def __init__(self,x_min,x_max,y_min,y_max,z_min,z_max):
        """ Initialize the simulation box with its size in each direction , TASK 0 """
        self.size_x_min=x_min
        self.size_x_max=x_max

        self.size_y_min=y_min
        self.size_y_max=y_max

        self.size_z_min=z_min
        self.size_z_max=z_max
        self.volume_of_the_box= (x_max - x_min)*(y_max - y_min)*(z_max - z_min)

        self.list_of_points=[]
        self.list_of_sphere=[]
        self.list_of_start_walker=[]
        self.set_of_valid_points= set() # Use of a set to avoid 

        self.number_of_points_in_sphere=0
        self.fig = None
        self.ax = None

    def setup_box(self):
        """ 
        Setup the 3D plot for the box, so that every sphere and point
        can be ploted on the same graph 
        """
        self.fig = plt.figure(figsize=(10,10))
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_xlim([self.size_x_min, self.size_x_max])
        self.ax.set_ylim([self.size_y_min, self.size_y_max])
        self.ax.set_zlim([self.size_z_min, self.size_z_max])
        self.ax.set_xlabel('X axis')
        self.ax.set_ylabel('Y axis')
        self.ax.set_zlabel('Z axis')
        self.ax.set_title('3D Sphere in Simulation Box')

    def random_point(self,n=1):
        """
        Generate n random points in the box, Function that correspond to TASK 1
        Input : n number of points to generate
        """
        for i in range(n):
            point=Point(np.random.uniform(self.size_x_min,self.size_x_max),
                        np.random.uniform(self.size_y_min,self.size_y_max),
                        np.random.uniform(self.size_z_min,self.size_z_max))
            self.list_of_points.append(point)

    #TODO ptet retirer cette fonction elle sert a rien 
    def show_points(self):
        """ Function that print all the points in the box, for debug purpose """
        print("Number of point in the box :", len(self.list_of_points))
        for p in self.list_of_points:
            print(p)

    def create_random_sphere(self):
        """
        Create the data for a random sphere in the box and add it the list of sphere of the class
        TASK 2
        """
        center= Point(np.random.uniform(self.size_x_min,self.size_x_max),
                      np.random.uniform(self.size_y_min,self.size_y_max),
                      np.random.uniform(self.size_z_min,self.size_z_max))
        radius= np.random.uniform(1,5)
        sphere=Sphere(center,radius)
        self.list_of_sphere.append(sphere)

    def create_specific_sphere(self,center,radius):
        """
        Create a sphere with wanted center and radius
        """
        sphere=Sphere(center,radius)
        self.list_of_sphere.append(sphere)
    def plot_sphere(self):
        """ Plot all the sphere in the list of sphere of the class """
        for sphere in self.list_of_sphere:
            sphere_x,sphere_y,sphere_z,radius= sphere.send_sphere_plot_data()
            self.ax.plot_surface(sphere_x,sphere_y,sphere_z, color='b', alpha=0.5)

    def plot_points(self):
        """ Plot all the points in the list of points of the class """
        for point in self.list_of_points:
            self.ax.scatter(point.x,point.y,point.z,color='r',s=50)
    def generate_sphere(self,n=1,specific=False):
        """ 
        Generate n random sphere in the box, Function that correspond to TASK 2 and TASK 6
        Input : n number of sphere to generate
        """
        for i in range(n):
            if specific:
                self.create_specific_sphere(center=Point(0,0,0),radius=1)
            else:
                self.create_random_sphere()
        print("Generating of spheres,number of spheres :",len(self.list_of_sphere))
    def test_volume_sphere(self):
        """ Function to test the volume calculation"""
        total_volume=0
        for sphere in self.list_of_sphere:
            coord_x,coord_y,coord_z,radius=sphere.send_sphere_center()
            volume=(4/3)*np.pi*radius**3
            total_volume+=volume
        print ("---Analytical Volume Test Results ---")
        print("Volume of the box is :",self.volume_of_the_box)
        print("The volume of the sphere are :",total_volume)
        return total_volume
    def check_point_spot(self):
        """
        Check if the points are inside the sphere or not, Function that correspond to TASK 4 and TASK 7 and TASK 10
        """
        for point in self.list_of_points:
            for sphere in self.list_of_sphere:
                coord_x,coord_y,coord_z,radius=sphere.send_sphere_center()
                if (point.x-coord_x)**2+(point.y - coord_y)**2+(point.z - coord_z)**2 <= radius**2:
                    self.number_of_points_in_sphere+=1
                    break
        print("--- Results ---")
        print("There is a total of ", len(self.list_of_points),
              "points, with ",self.number_of_points_in_sphere,"inside the spheres.")
        print("Fraction of point that are inside the spheres :",
              self.number_of_points_in_sphere/len(self.list_of_points))
        estimated_sphere_volume= self.volume_of_the_box *(self.number_of_points_in_sphere/len(self.list_of_points))
        print("estimated volume of sphere :", estimated_sphere_volume)

        return estimated_sphere_volume
    def estimate_pi(self):
        """ Function to estimate the value of pi using Monte Carlo Method, TASK 5"""
        number_of_point_in=0
        for point in self.list_of_points:
            if point.x**2 + point.y**2 + point.z**2 <= 1:
                number_of_point_in+=1

        pi_estimate = 6 * number_of_point_in / len(self.list_of_points)
        print("Estimated value of pi :", pi_estimate)
        print("Number of points inside the octant of sphere :", number_of_point_in)
        return pi_estimate


    def read_dna_coords(self):
        """ 
        Read the dna_coords.txt file and create spheres from the data, 
        Function that correspond to TASK 8  and TASK 9 
        """
        data=pd.read_csv("dna_coords.txt",sep="\s+",names=["atom","x","y","z","radius"])
        radius_data={"H": 1.20, "C" :1.70,"N": 1.55, "O": 1.50, "P": 1.80}
        data["radius"] = data["atom"].map(radius_data)
        self.data_radius=data["radius"]
        for atom in data.itertuples():
            point = Point(atom.x,atom.y,atom.z)
            sphere=Sphere(point,atom.radius)
            self.list_of_sphere.append(sphere)
        print("Number of spheres created from DNA data :", len(self.list_of_sphere))

    def generate_starting_walker(self,n=1):
        """
        Generate n random starting point for the random walker in the box
        Input : n number of starting point to generate
        """
        for i in range(n):
            point=Point(np.random.uniform(self.size_x_min,self.size_x_max),
                        np.random.uniform(self.size_y_min,self.size_y_max),
                        np.random.uniform(self.size_z_min,self.size_z_max))
            self.list_of_start_walker.append(point)
    def start_walker_simulation(self,steps=1000):
        """
        Function That start the walker simulation
        Input : steps -> maximum number of step a walker can do
        """
        directions = {
                    0: ("x", +1),
                    1: ("x", -1),
                    2: ("y", +1),
                    3: ("y", -1),
                    4: ("z", +1),
                    5: ("z", -1)
                }
        for walker in self.list_of_start_walker:
            for step in range(steps):
                number=np.random.randint(0,6)
                axis, delta = directions[number]
                setattr(walker, axis, getattr(walker, axis) + delta)
                bool_outside_box = (walker.x < self.size_x_min or walker.x >= self.size_x_max or
                                    walker.y < self.size_y_min or walker.y >= self.size_y_max or
                                    walker.z < self.size_z_min or walker.z >= self.size_z_max)
                if bool_outside_box:
                    break
                bool_inside_sphere = False
                for sphere in self.list_of_sphere:
                    coord_x,coord_y,coord_z,radius=sphere.send_sphere_center()
                    bool_inside_sphere =(walker.x-coord_x)**2+(walker.y - coord_y)**2+(walker.z - coord_z)**2 <= radius**2
                    if bool_inside_sphere:
                        break
                if bool_inside_sphere :
                    break
                self.set_of_valid_points.add((int(walker.x),int(walker.y),int(walker.z)))
        print("Number of valid  point for the walker :", len(self.set_of_valid_points))
        estimated_sphere_volume=self.volume_of_the_box - len(self.set_of_valid_points)
        print("Volume estimated of the sphere from the walker simulation :", self.volume_of_the_box - len(self.set_of_valid_points))
        return estimated_sphere_volume

def launch_simulation(box_size=(1,1,1,1,1,1),number_of_point=100,number_of_sphere=1,specific_sphere=False,dna_simulation=False):
    """
    Function to launch a box simulation, to avoid duplicating code in the jupyter notebook
    Input : box_size -> the values of the edges of the box 
            number_of_point -> number of random point to add to the simulation
            number_of_sphere -> number of random sphere to add to the simulation
            specific_sphere -> Bool value, specific condition to create a sphere for the estimation of pi 
            dna_simulation -> Bool value, specific condition to load or no the data form the "dna_coords.txt" file
    Output : box create in case we need to do other actions on it
    """
    box=simulation_box(box_size[0],box_size[1],box_size[2],box_size[3],box_size[4],box_size[5])
    box.setup_box()
    if dna_simulation:
        box.read_dna_coords()
        number_of_sphere=0  # No need to generate spheres if we read DNA data
    box.random_point(n=number_of_point)
    box.generate_sphere(n=number_of_sphere,specific=specific_sphere)
    box.plot_points()
    box.plot_sphere()
    box.check_point_spot()
    return box
import numpy as np 

# Code for the class Point
class Point:
    """ Class that represent a point in a 3D space, helpful for the simulation box"""
    def __init__(self,x,y,z):
        """ Initialize the point with x,y,z coordinates """
        self.x=x
        self.y=y
        self.z=z
    def __str__(self):
        """ String representation of the point """
        return f"Point({self.x}, {self.y}, {self.z})"
#Code for the class Sphere
class Sphere:
    """ Class that represent a sphere in a 3D space, helpful for the simulation box"""
    def __init__(self,center,radius):
        """ Initialize the sphere with center (Point) and radius """
        self.center=center
        self.radius=radius
    def __str__(self):
        """ String representation of the sphere """
        return f"Sphere(Center: {self.center}, Radius: {self.radius})"

    def send_sphere_plot_data(self):
        """ Create the data for the sphere, helpful for the ploting function """
        theta = np.linspace(0, 2 * np.pi, 100)
        phi = np.linspace(0, np.pi, 50)
        theta, phi = np.meshgrid(theta, phi) # Get a 2D grid
        sphere_x = self.radius * np.sin(phi) * np.cos(theta) + self.center.x
        sphere_y = self.radius * np.sin(phi) * np.sin(theta) + self.center.y
        sphere_z = self.radius * np.cos(phi) + self.center.z
        return sphere_x, sphere_y, sphere_z,self.radius
    def send_sphere_center(self):
        """ Send the center coordinates and radius of the sphere """
        return self.center.x, self.center.y, self.center.z, self.radius

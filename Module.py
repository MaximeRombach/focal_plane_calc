import numpy as np
np.set_printoptions(legacy='1.25')
from math import sqrt
from shapely.geometry import Polygon, MultiPolygon
from shapely.plotting import plot_polygon
from shapely.ops import unary_union
from Robot import Robot
import matplotlib.pyplot as plt
import CustomLegends as cl
import geopandas as gpd

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

class Module:
    """Class to uniquely define a module of the focal plane

    Contains the following attributes:
    - (int) nb_robots: number of robots in the module
    - (float) pitch: pitch of the module
    - (float) x0: x-coordinate of the center of the module in the global coordinate system
    - (float) y0: y-coordinate of the center of the module in the global coordinate system
    - (float) z0: z-coordinate of the center of the module in the global coordinate system
    - (int) module_id: unique identifier of the module
    """
    def __init__(self, nb_robots, pitch, **kwargs):
        
        self.pitch = pitch
        self.nb_robots = nb_robots # Number of robots in the module
        self.chamfer_size = kwargs.get('chamfer_size', 10.5)
        self.l_alpha = kwargs.get('l_alpha', 1.8)
        self.l_beta = kwargs.get('l_beta', 1.8)
        self.HR_l_alpha = kwargs.get('HR_l_alpha', 1.8)
        self.HR_l_beta = kwargs.get('HR_l_beta', 1.8)
        self.HR_fibers = kwargs.get('HR_fibers', [10, 12, 15, 17, 29, 34, 39, 50, 52, 59])
        # self.HR_fibers = kwargs.get('HR_fibers', [])

        self.x0 = kwargs.get('x0', 0)
        self.y0 = kwargs.get('y0', 0)
        self.z0 = kwargs.get('z0', 0)
        self.module_centroid = [self.x0, self.y0, self.z0]
        self.module_points_up = kwargs.get('module_points_up', True) # True if module oriented upward; False if module points down

        self.safety_margin = kwargs.get('safety_margin', 0.5)
        self.beta2fibre = kwargs.get('beta2fibre', 1)
        self.dist2wall = self.beta2fibre + self.safety_margin # The fiber is physically not at the edge of the beta arm which prevents it from reaching the wall + add a safety margin at which the positioner has to stay away from the wall
        self.LR_area = 0 # area covered by Low Resolution fibers
        self.HR_area = 0 # area covered by High Resolution fibers
        self.module_area = 0 # total area covered by all robots in a module
        self.module_id = kwargs.get('module_id', 1)
        
        self.__r_flat = kwargs.get('r_flat', 0)
        self.__theta_flat = kwargs.get('theta_flat', 0)
        self.__phi_flat = kwargs.get('phi_flat', 0)

        self.is_wall = kwargs.get('is_wall', True)
        self.remove_corners = kwargs.get('remove_corners', True)
        
        self.robots = []
        self.dataframe ={'module_id':[], 'robot_id':[], 'fiber_type':[], 'x0':[], 'y0':[], 'z0':[], 'color':[], 'geometry':[]}
        self.__workspace = None

    @property
    def module_side_length(self):
        """Computes the side length of the module.
        
        Arbitrary original side length set to 80 for 75 robots with 6.2 mm pitch (12 lines of robots).
        The side length is then adjusted according to the number of robots in the module.
        
        TEST: use only pitch*number of lines to scale better triangle size with pitch --> accomodate better change of pitch for WST for example
        
        """
        # return 80 + self.pitch * (self.nb_lines_of_robots() - 12)
        return self.pitch * (self.nb_lines_of_robots() + 1)
    
    @property
    def module_height(self):
        """ Computes the height of the equilateral triangular enveloppe of the module """

        return self.module_side_length * np.sqrt(3)/2
    
    @property
    def nbots(self):
        """ Computes the actual number of robots in the module. If for example the corners are accounted for"""
        if not self.remove_corners:
            return self.nb_robots + 3
        else:
            return self.nb_robots
    
    @property
    def module_boundaries(self):
        """Chamfers the corners of the original triangle by replacing sharp corners with a bevel."""
        t = self.module_side_length * np.sqrt(3) / 3
        # Adding the position of the module center as offset
        # This ways the module can be placed anywhere in the global coordinate system
        original_triangle = [[0, t, 0],
                             [t * np.cos(np.pi / 6), -t / 2, 0],
                             [-t * np.cos(np.pi / 6), -t / 2, 0]]
        
        if not self.module_points_up:
            """ If the module is not oriented upwards, flip the coordinates of the triangle """
            for i in range(len(original_triangle)):
                original_triangle[i] = self.flip_coordinates(original_triangle[i])

        original_triangle = np.array(original_triangle) + np.array(self.module_centroid)
        
        new_coords = []
        num_points = len(original_triangle)

        for i in range(num_points):
            p1 = np.array(original_triangle[i - 1])  # Previous point
            p2 = np.array(original_triangle[i])      # Current corner point
            p3 = np.array(original_triangle[(i + 1) % num_points])  # Next point

            # Compute direction vectors
            v1 = p1 - p2
            v2 = p3 - p2

            # Normalize and scale by chamfer size
            v1 = v1 / np.linalg.norm(v1) * self.chamfer_size
            v2 = v2 / np.linalg.norm(v2) * self.chamfer_size

            # New chamfered points
            new_coords.append(tuple(p2 + v1))
            new_coords.append(tuple(p2 + v2))

        return Polygon(new_coords)
        # return new_coords
    
    
    @property
    def fiber_type_assignment(self):
        """Assigns a fiber type to each robot in the module.
        
        Returns:
        - (list) fiber_types: list of fiber types assigned to the robots in the module
        
        """
        fiber_types = []
        for i in range(self.nbots):
            if i in self.HR_fibers:
                """ Assign High Resolution fibers to the robots in the module """
                fiber_types.append('HR')
            else:
                fiber_types.append('LR')
        return fiber_types
    
    @property
    def nb_of_HR_fibers(self):
        """Computes the number of High Resolution fibers in the module.
        
        Returns:
        - (int) nb_HR_fibers: number of High Resolution fibers in the module
        
        """
        return len(self.HR_fibers)
    
    @property
    def nb_of_LR_fibers(self):
        """Computes the number of Low Resolution fibers in the module.
        
        Returns:
        - (int) nb_LR_fibers: number of Low Resolution fibers in the module
        
        """
        return self.nbots - self.nb_of_HR_fibers
    
    @property
    def robots_layout(self) -> list:
        """Computes the layout of the robots in a module with respect to the number of robots lines.
        
        Returns:
        - (list) robots: list of objects Robot in the module
        
        """
        n = self.nb_lines_of_robots()
        fiber_types = self.fiber_type_assignment
        layout_center_x = -((n - 1) * self.pitch / 2)
        layout_center_y = -((n - 1) * self.pitch * np.sqrt(3) / 6)
        robot_index =  0
        "Starts from the bottom left corner of the module and builds each line of robots from left to right"
        for j in range(n):
            for i in range(n-j):
                if (j == n - 1 or (j == 0 and i == 0) or (j == 0 and i == n - 1)) and self.remove_corners:
                    """ Remove corner robots"""
                    continue
                else:
                    x = (i * self.pitch + 0.5 * self.pitch * j) + layout_center_x
                    y = (j * self.pitch * np.sqrt(3) / 2) + layout_center_y
                    z = self.z0

                    if not self.module_points_up:
                        x,y,z = self.flip_coordinates([x,y,z])

                    if robot_index in self.HR_fibers:
                        """Change the arm lengths for High Resolution fibers"""
                        new_l_alpha = self.HR_l_alpha
                        new_l_beta = self.HR_l_beta
                    else:
                        """Otherwise keep nominal"""
                        new_l_alpha = self.l_alpha
                        new_l_beta = self.l_beta

                    new_robot = Robot(robot_id = robot_index * self.module_id,
                                    module_id = self.module_id,
                                    l_alpha = new_l_alpha,
                                    l_beta = new_l_beta,
                                    fiber_type = fiber_types[robot_index],
                                    x0 = x + self.x0,
                                    y0 = y + self.y0,
                                    z0 = z + self.z0,)

                    self.dataframe['module_id'].append(new_robot.module_id)
                    self.dataframe['robot_id'].append(new_robot.robot_id)
                    self.dataframe['fiber_type'].append(new_robot.fiber_type)
                    self.dataframe['x0'].append(new_robot.x0)
                    self.dataframe['y0'].append(new_robot.y0)
                    self.dataframe['z0'].append(new_robot.z0)
                    self.dataframe['color'].append(new_robot.color)
                    self.dataframe['geometry'].append(new_robot.workspace)

                    self.robots.append(new_robot)
                robot_index += 1

        return self.robots
        
    @property
    def LR_robots(self):
        """
        Input: list of all the robots in the module
        Output: list of robots objects assigned to LR fibers
        """
        LR_list = []
        for robot in self.robots:

            if robot.fiber_type == 'LR':
                LR_list.append(robot)

        return LR_list
    
    @property
    def LR_coverage(self):
        """
        Input: list of all the LR robots in the module
        Output:
        - Polygon of coverage of the LR fibers in the module
        - Area covered by the LR fibers
        """

        LR_workspaces = [rob.workspace for rob in self.LR_robots]
        LR_union = unary_union(LR_workspaces)

        if self.is_wall:
            dist2wall = self.beta2fibre + self.safety_margin # The fiber is physically not at the edge of the beta arm which prevents it from reaching the wall + add a safety margin at which the positioner has to stay away from the wall
            lr_cov = LR_union.intersection(self.module_boundaries.buffer(-dist2wall))

        else:
            lr_cov = LR_union

        self.LR_area = lr_cov.area
                
        return lr_cov
    
    @property
    def HR_robots(self):
        """
        Input: list of all the robots in the module
        Output: list of robots objects assigned to HR fibers
        """
        HR_list = []
        for robot in self.robots:

            if robot.fiber_type == 'HR':
                HR_list.append(robot)
                
        return HR_list
    
    @property
    def HR_coverage(self):
        """
        Input: list of all the HR robots in the module
        Output: Polygon of coverage of the HR fibers in the module
        """
        HR_workspaces = [rob.workspace for rob in self.HR_robots]
        HR_union = unary_union(HR_workspaces)

        if self.is_wall:
            
            hr_cov = HR_union.intersection(self.module_boundaries.buffer(-self.dist2wall))

        else:
            hr_cov = HR_union

        self.HR_area = hr_cov.area
                
        return hr_cov
    
    @property
    def module_coverage(self):
        
        """
        Calculate the total module coverage by combining the LR and HR coverage polygons
        Input: LR_coverage, HR_coverage polygons 
        Output:
                - Polygon of total coverage of the module
                - Updates the module area attribute
        
        """
        
        mod_cov = unary_union([self.LR_coverage, self.HR_coverage])

        self.module_area = mod_cov.area

        return mod_cov
    
    def nb_lines_of_robots(self):
        """Computes the number of lines of robots in a module for a given amount of robots.
        Solves equation for n: n(n+1)/2 = nb_robots + 3 (see triangular number formula: https://en.wikipedia.org/wiki/Triangular_number)
        """
        discriminant = 1 + 4 * (2*(self.nb_robots + 3))
        # discriminant = 1 + 4 * (2*(self.nb_robots))
        if discriminant < 0 or discriminant == 0:
            """
            We know that any other value than a positive integer is not valid.
            """
            raise ValueError("Discriminant must be > 0. \n The number of robots is not valid. Examples of valid numbers 42, 52, 63, 75, 88, 102")
        else:
            n = (-1 + sqrt(discriminant)) / 2
            if n.is_integer():
                return int(n)
            else:
                """ We also know that only integer number of lines is valid."""
                raise ValueError("Number of lines must be an integer. \n The number of robots is not valid. Examples of valid numbers are 42, 52, 63, 75, 88, 102")

        return n
    
    def change_arms_lengths(self, new_l_alpha, new_l_beta):
            """Changes the lengths of the alpha and beta arms of the robots in the module.
            
            Args:
            - (float) l_alpha: new length of the alpha arm
            - (float) l_beta: new length of the beta arm
            
            """
            self.l_alpha = new_l_alpha
            self.l_beta = new_l_beta

    def flip_coordinates(self, coord):
        """ Flips coordinates 180° about z axis """
        R_theta = np.array([[np.cos(np.pi), -np.sin(np.pi), 0],
                             [np.sin(np.pi), np.cos(np.pi), 0],
                             [0, 0, 1]])
        coord = R_theta @ np.array(coord).T
        coord = coord.tolist()

        return coord
    
if __name__ == "__main__":

    mod = Module(63,
                6.2, 
                module_points_up = True,
                x0 = 0,
                y0 = 0,
                z0 = 0)

    robots = mod.robots_layout

    figure, ax = plt.subplots(figsize=(10, 10))
    """ Using geopandas to plot Polygon objects is way faster than looping plot_polygon """
    geo = gpd.GeoDataFrame(mod.dataframe)
    geo.plot(ax = ax, color=geo['color'], alpha=0.5, markersize=10, legend=True)
    plot_polygon(mod.module_boundaries, fill = False, add_points=False, color = 'black')
    plot_polygon(mod.module_boundaries.buffer(-mod.dist2wall), fill = False, add_points=False, color = 'red', linestyle = '--')
    for rob in robots:
        plt.text(rob.x0, rob.y0, str(rob.robot_id), fontsize=11, ha='center', va='center')
    plt.legend(handles = [cl.LR_handle(lab = f'LR fibers: {mod.nb_of_LR_fibers}'),
                          cl.HR_handle(lab = f'HR fibers: {mod.nb_of_HR_fibers}'),
                          cl.safety_margin_handle(lab = f'Dist2wall: {mod.dist2wall} mm')])
    plt.grid(visible=True)
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')
    plt.title(f'Raw workspaces of robots \n {len(robots)} robots per module')
    plt.show()
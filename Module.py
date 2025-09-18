import numpy as np
np.set_printoptions(legacy='1.25')
from math import sqrt
from shapely.geometry import Polygon, MultiPolygon
from shapely.plotting import plot_polygon
from shapely.ops import unary_union
from Robot import Robot
import matplotlib.pyplot as plt
plt.rc('axes', labelsize=13)    # fontsize of the x and y labels
plt.rc('figure', titlesize=16)  # fontsize of the figure title
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('xtick', labelsize=13)    # fontsize of the tick labels
plt.rc('ytick', labelsize=13)    # fontsize of the tick labels
plt.rc('legend', fontsize=13)    # legend fontsize
import CustomLegends as cl
import geopandas as gpd
import SavingResults as sr

import copy

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
        # self.HR_fibers = kwargs.get('HR_fibers', [10, 12, 15, 17, 29, 34, 39, 50, 52, 59])
        self.HR_fibers = kwargs.get('HR_fibers', [21, 25, 51])
        # self.HR_fibers = kwargs.get('HR_fibers', [12, 15, 29, 34, 50, 52])
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
        self._LR_coverage = None # Polygon of coverage of the LR fibers in the module
        self.HR_area = 0 # area covered by High Resolution fibers
        self._HR_coverage = None # Polygon of coverage of the HR fibers in the module
        self.module_coverage_area = 0 # total area covered by all robots in a module
        self.module_id = kwargs.get('module_id', 1)
        self.is_wall = kwargs.get('is_wall', True)
        self.remove_corners = kwargs.get('remove_corners', True)
        
        self.robots = []
        self.dataframe ={'module_id':[], 'robot_id':[], 'fiber_type':[], 'x0':[], 'y0':[], 'z0':[], 'color':[], 'geometry':[]}

    @property
    def module_side_length(self):
        """Computes the side length of the module.
        
        Arbitrary original side length set to 80 for 75 robots with 6.2 mm pitch (12 lines of robots).
        The side length is then adjusted according to the number of robots in the module.
        
        TEST: use only pitch*number of lines to scale better triangle size with pitch --> accomodate better change of pitch for WST for example
        
        """
        # return 80 + self.pitch * (self.nb_lines_of_robots() - 12)
        if self.nb_robots == 63:
            "Force value to 73.8 for 63 robots"
            return 73.8
        else:
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

    def module_boundaries_3D(self):
        """Computes the 3D boundaries of the module."""
        # Get the 2D polygon of the module boundaries
        polygon_2d = self.module_boundaries
        
        # Create a 3D polygon by adding z-coordinates
        points_3d = [(x, y, self.z0) for x, y in polygon_2d.exterior.coords]
        
        # Create a Poly3DCollection for the 3D polygon
        poly3d = Poly3DCollection([points_3d], alpha=0.5, edgecolor='k')
        
        return poly3d
    
    @property
    def module_boundaries_with_safety_margin(self):
        """Computes the limits inside the module where the robots can operats

        Input: 
        - (Polygon) module_boundaries: physical limits of the module walls
        - (float) dist2wall: beta2fiber + safety_margin

        Output:
        - (Polygon) module_boundaries_with_safety_margin: boundaries inside which robots workspace are limited
        
        """
        
        return self.module_boundaries.buffer(-self.dist2wall)
    
    @property
    def fiber_type_assignment(self):
        """Assigns a fiber type to each robot in the module.
        
        Input:
        - (list) predefined list of the HR fibers in the module
        
        Output:
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
        sqrt3 = np.sqrt(3)
        layout_center_x = -((n - 1) * self.pitch / 2)
        layout_center_y = -((n - 1) * self.pitch * sqrt3 / 6)
        robot_index =  0

        "Starts from the bottom left corner of the module and builds each line of robots from left to right"
        for j in range(n):
            for i in range(n-j):
                if (j == n - 1 or (j == 0 and i == 0) or (j == 0 and i == n - 1)) and self.remove_corners:
                    """ Remove corner robots"""
                    continue
                else:
                    x = (i * self.pitch + 0.5 * self.pitch * j) + layout_center_x
                    y = (j * self.pitch * sqrt3 / 2) + layout_center_y
                    z = self.z0
                    is_hr = False

                    if not self.module_points_up:
                        x,y,z = self.flip_coordinates([x,y,z])

                    if robot_index in self.HR_fibers:
                        is_hr =True
                        new_l_alpha = self.HR_l_alpha
                        new_l_beta = self.HR_l_beta
                    else:
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
        if self._LR_coverage is None:
            LR_workspaces = [rob.workspace for rob in self.LR_robots]
            LR_union = unary_union(LR_workspaces)

            if self.is_wall:
                lr_cov = LR_union.intersection(self.module_boundaries_with_safety_margin)
                # lr_cov = MultiPolygon(LR_workspaces)
            else:
                lr_cov = LR_union

            self.LR_area = lr_cov.area
            self._LR_coverage = lr_cov
                
        return self._LR_coverage
    
    @LR_coverage.setter
    def LR_coverage(self, coverage_polygon):
        """
        Setter for the LR_coverage property.
        This allows to set the LR_coverage polygon directly.
        
        Input:
        - (Polygon): coverage polygon of the LR fibers
        
        """
        self._LR_coverage = coverage_polygon
        self.LR_area = coverage_polygon.area
        self.dataframe['geometry'] = coverage_polygon
    
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
        if self._HR_coverage is None:
            HR_workspaces = [rob.workspace for rob in self.HR_robots]
            HR_union = unary_union(HR_workspaces)

            if self.is_wall:
                
                hr_cov = HR_union.intersection(self.module_boundaries_with_safety_margin)
                # hr_cov_raw = MultiPolygon(HR_workspaces).intersection(self.module_boundaries.buffer(-self.dist2wall))

            else:
                hr_cov = HR_union

            self.HR_area = hr_cov.area
            self._HR_coverage = hr_cov
                
        return self._HR_coverage
    
    @HR_coverage.setter
    def HR_coverage(self, coverage_polygon):
        """
        Setter for the HR_coverage property.
        This allows to set the HR_coverage polygon directly.
        
        Input:
        - (Polygon): coverage polygon of the HR fibers
        
        """
        self._HR_coverage = coverage_polygon
        self.HR_area = coverage_polygon.area
        self.dataframe['geometry'] = coverage_polygon

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

        self.module_coverage_area = mod_cov.area

        return mod_cov
    
    def raw_module_coverage(self, coverage_polygon):
        """
        Computes coverage polygon of the robots WITHOUT summing up all the similar reaches.
        Used for exporting DXFs of robots arrangements and visualization
        Input:

        - (Polygon): coverage polygon of corresponding HR or LR fibers
        - (Polygon): module boundaries reduced with dist2wall

        Output:
       
        - (Polygon): trimmed coverage polygon
        
        """
        
        return coverage_polygon.difference(self.module_boundaries_with_safety_margin)

    
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
    
    def plot_module(self, plot_rob_numbers = False):
        """Plots the module boundaries and the robots in the module."""
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)
        
        # Plot module boundaries
        plot_polygon(self.module_boundaries, ax=ax, add_points=False, fill=False, color='black', linestyle='--')
        
        # Plot LR coverage
        # if self.LR_coverage is not None:
        #     plot_polygon(self.LR_coverage, add_points=False, ax=ax, fill=True, color='C0', alpha=0.4)
        
        # Plot HR coverage
        # if self.HR_coverage is not None:
        #     plot_polygon(self.HR_coverage, add_points=False, ax=ax, fill=True, color='red', alpha=0.2)
        
        # Plot module boundaries with safety margin
        plot_polygon(self.module_boundaries_with_safety_margin, ax=ax, add_points=False, fill=False, color='red', linestyle='--')

        # Plot robots
        geo = gpd.GeoDataFrame(self.dataframe)
        geo.plot(ax = ax, color=geo['color'], alpha=0.5, markersize=10, legend=True)
        
        if plot_rob_numbers == True:
            """ Plot robot ID at the center of each robot """
            for rob in self.robots_layout:
                plt.text(rob.x0, rob.y0, str(rob.robot_id+1), fontsize=11, ha='center', va='center')


        if self.nb_of_HR_fibers != 0:
            plt.legend(handles = [cl.LR_handle(extra_lab = f'{self.nb_of_LR_fibers}'),
                                cl.HR_handle(extra_lab = f'{self.nb_of_HR_fibers}'),
                                cl.safety_margin_handle(lab = f'Dist2wall: {self.dist2wall} mm')])
        else:
            plt.legend(handles = [cl.LR_handle(extra_lab = f'{self.nb_of_LR_fibers} robots workspaces'),
                                cl.safety_margin_handle(lab = f'Dist2wall: {self.dist2wall} mm')])

        plt.xlabel('x [mm]')
        plt.ylabel('y [mm]')
        # plt.xlim([-500, 500])
        # plt.ylim([-500, 500])
        plt.title(cl.module_title(self.nb_robots, self.module_side_length, self.pitch, self.l_alpha, self.l_beta, self.HR_l_alpha, self.HR_l_beta))
        plt.grid()
    
if __name__ == "__main__":

    save = sr.SavingResults({"save_plots": True,
                            "save_txt": False},
                            project_name = 'test')
    mod = Module(63,
                6.2, 
                module_points_up = True,
                x0 = 100,
                y0 = 100,
                z0 = 0,
                HR_fibers = [21, 25, 39, 51])

    robots = mod.robots_layout
    print(robots[0].theta, robots[0].phi, robots[0].r, robots[0].r_flat)

    mod.plot_module(plot_rob_numbers=True)
    plt.scatter(robots[31].x0, robots[31].y0, color='blue')
    plt.show()
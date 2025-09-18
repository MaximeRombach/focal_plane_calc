import numpy as np
from shapely.geometry import Polygon
from shapely.plotting import plot_polygon
import matplotlib.pyplot as plt

class Robot:

    """Class to uniquely define a robot of the focal plane within a module
    
    Contains the following attributes:
    - (float) l_alpha: length of the first arm of the robot
    - (float) l_beta: length of the second arm of the robot
    - (string) fiber_type: type of fiber that the robot is carrying.
        - Default is None if the fibers are all the same in the project
        - Options are 'LR' for low resolution fibers and 'HR' for high resolution fibers
    - (int) robot_id: unique identifier of the robot
    - (int) module_id: unique identifier of the module
    - (float) x0: x-coordinate of the center of the robot workspace in the global coordinate system
    - (float) y0: y-coordinate of the center of the robot workspace in the global coordinate system
    - (float) z0: z-coordinate of the center of the robot workspace in the global coordinate system
    - (float) r_flat: radial position of the center of the robot in spherical coordinates
    - (float) theta_flat: azimuthal position of the center of the robot in spherical coordinates
    - (float) phi_flat: polar position of the center of the robot in spherical coordinates
    - (Polygon) workspace: shapely Polygon object representing the workspace of the robot

    """


    def __init__(self, l_alpha, l_beta, **kwargs):

        self.l_alpha = l_alpha
        self.l_beta = l_beta
        self._workspace_outer_radius = self.l_alpha + self.l_beta # [mm] outer radius of the workspace of the robot
        self._workspace_inner_radius = abs(self.l_alpha - self.l_beta)

        self.__robot_id = kwargs.get('robot_id', None)
        self.__module_id = kwargs.get('module_id', None)
        self.fiber_type = kwargs.get('fiber_type', None) # Type of fiber that the robot is carrying: Default is None if the fibers are all the same in the project
        self.__x0 = kwargs.get('x0', 0)
        self.__y0 = kwargs.get('y0', 0)
        self.__z0 = kwargs.get('z0', 0)
        self.r = np.sqrt(self.__x0**2 + self.__y0**2 + self.__z0**2) # [mm] radial distance of the center of the robot in 3D space
        self.phi = np.degrees(np.arctan2(self.__y0, self.__x0)) # [deg] azimuthal angle of the center of the robot in spherical coordinates
        self.theta = np.degrees(np.arccos(self.__z0 / self.r)) if self.r != 0 else 0 # [deg] polar angle of the center of the robot in spherical coordinates

    # Internal cache for expensive computation
        self._workspace = None


    @property
    def robot_id(self):
        return self.__robot_id
    
    @robot_id.setter
    def robot_id(self, robot_id):
        self.__robot_id = robot_id
    
    @property
    def module_id(self):
        return self.__module_id
    
    @module_id.setter
    def module_id(self, module_id):
        self.__module_id = module_id
    
    @property
    def x0(self):
        return self.__x0
    
    @x0.setter
    def x0(self, x0):
        self.__x0 = x0

    @property
    def y0(self):
        return self.__y0
    
    @y0.setter
    def y0(self, y0):
        self.__y0 = y0

    @property
    def z0(self):
        return self.__z0
    
    @z0.setter
    def z0(self, z0):
        self.__z0 = z0

    @property
    def r_flat(self):
        return np.sqrt(self.__x0**2 + self.__y0**2)

    @property
    def workspace(self):

        """Lazy-cached workspace property."""
        if self._workspace is not None:
            return self._workspace

        """Returns the workspace of the robot as a shapely Polygon object"""

        angle_outer = np.linspace(0, 2* np.pi, 20)
        angle_lower = np.linspace(0, 2* np.pi, 10)

        "Outer boundaries of positioner workspace"
        # CALCULATE THE those things only once in declaration of class!!!
        x_wk_ext = self.x0 + self._workspace_outer_radius * np.cos(angle_outer)
        y_wk_ext = self.y0 + self._workspace_outer_radius * np.sin(angle_outer)
        z_wk_ext = self.z0 * np.ones_like(x_wk_ext)
        wk_ext = [(x,y,z) for x,y,z in zip(x_wk_ext, y_wk_ext, z_wk_ext)]

        "Inner boundaries of positioner workspace; 0 if l_alpha = l_beta"
        if self.l_alpha == self.l_beta:
            wk_int = []
        else:
            x_wk_int = self.x0 + self._workspace_inner_radius * np.cos(angle_lower)
            y_wk_int = self.y0 + self._workspace_inner_radius * np.sin(angle_lower)
            z_wk_int = self.z0 * np.ones_like(x_wk_ext)
            wk_int = [(x,y,z) for x,y,z in zip(x_wk_int, y_wk_int, z_wk_int)]
            wk_int = wk_int[::-1] # Remove last point of inner boundary to avoid having a straight line in the center of the workspace

        workspace = Polygon(wk_ext, [wk_int])
        # workspace = wk_ext
        return workspace
    
    @property
    def color(self):
        if self.fiber_type == 'LR':
            return 'C0'
        elif self.fiber_type == 'HR':
            return 'red'
        else:
            return 'C0'
        
    def set_spherical_coords(self, r: float, theta: float, phi: float):
        self.r = r
        self.theta = theta
        self.phi = phi
        self.__x0 = r * np.sin(np.radians(theta)) * np.cos(np.radians(phi))
        self.__y0 = r * np.sin(np.radians(theta)) * np.sin(np.radians(phi))
        self.__z0 = r * np.cos(np.radians(theta))
        self.r_flat = np.sqrt(self.__x0**2 + self.__y0**2)

if __name__ == '__main__':
    robot = Robot(l_alpha = 1,
                  l_beta = 2,
                  robot_id=1,
                  module_id=1,
                  x0=0,
                  y0=0, 
                  z0=0)
    plot_polygon(robot.workspace) # Should print the workspace of the robot as a shapely Polygon object
    plt.show()
    print(robot.color) # Should print the color of the robot based on the fiber type
        

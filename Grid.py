import updown_tri as tri
import numpy as np
from shapely.geometry import Point
from shapely.ops import unary_union
from shapely.plotting import plot_polygon
import pandas as pd
from Focal_Suface import FocalSurf
from matplotlib import pyplot as plt
import CustomLegends as cl
import json
import geopandas as gpd

import logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-4s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

class Grid(FocalSurf):
    def __init__(self, project: str,
                 inner_gap: float, 
                 global_gap: float,
                 module_side_length: float,
                 **kwargs):

        self.project = project # project name    
        self.inner_gap = inner_gap
        self.global_gap = global_gap
        self.module_side_length = module_side_length
        self.centered_on_triangle = kwargs.get("centered_on_triangle", False)
        self.n = None

        self.flat_grid_dict = {'x': [], 'y': [], 'z': [], 'tri_points_up': []}
        self.fiducials = {'x': [], 'y': [], 'r': [], 'phi': [], 'z': [], 'geometry': []}

        super().__init__(project, **kwargs)

    @property
    def inter_triangle_side_length(self):

        """ Defines the triangular enveloppe of the "intermediate triangle"
        that locally brings closer 4 modules together
        """
        return 2 * self.single_triangle_side_length + 2*self.global_gap*np.cos(np.pi/6)
        # return 2*self.module_side_length + 2*self.inner_gap*np.cos(np.pi/6) + 2*self.global_gap*np.cos(np.pi/6)

    @property
    def single_triangle_side_length(self):
         
         " Defines the triangular enveloppe of a single module within the grid"
         
         return self.module_side_length + 2*self.inner_gap*np.cos(np.pi/6)

    def flat_grid(self):
        # Make global grid out of triangular grid method credited in updown_tri.py
          # Its logic is also explained in updown_tri.py
          # TODO: Use two diffrent n numbers for the two different grids: modules and ficudials OR find another optimization
          self.n = int(self.vigR/self.module_side_length) # first raw guess of number of lines of modules
          n = self.n
          center_coords = []
          a_max = n
          b_max = n
          c_max = n

          if self.centered_on_triangle:
               origin = 'triangle'
               valid = [0,1]
          else:
               origin = 'vertex'
               valid = [1,2]
          
          vigR_tresh = 150

          for a in np.arange(-a_max,a_max):
               for b in np.arange(-b_max,b_max):
                    for c in np.arange(-c_max,c_max):
                         sum_abc = a + b + c
                         # if valid == 1 or valid == 2: 
                         if valid.count(sum_abc): # check if sum abc corresponds to a valid triangle depending on the centering case
                              x,y = tri.tri_center(a,b,c,self.inter_triangle_side_length)
                              # Intermediate triangles placement
                              if Point(x,y).within(self.vignetting_disk): # allow centroid of inter modules to go out of vigR for further filling purpose
                                   self.flat_grid_dict['x'].append(x)
                                   self.flat_grid_dict['y'].append(y)
                                   # self.flat_grid_dict['geometry'].append(Point(x,y))
                                   intermediate_tri_points_up = tri.points_up(a,b,c, origin = origin) # check if INTERMEDIATE triangle is up or down; True: tri points up; False: tri points down
                                   self.flat_grid_dict['tri_points_up'].append(not intermediate_tri_points_up) # append inverse of points because center module of intermediate triangle is always pointing in the opposite directiob
                                   self.add_neighboring_modules(x,y,intermediate_tri_points_up) # add the three neighboring modules to the grid
                                   
                              #TODO: make smaller grid of fiducials to fill up more spaces
                              fiducial = tri.tri_corners(a,b,c,self.inter_triangle_side_length / 2)
                              # Fiducials placement
                              for fid in fiducial:
                                   if Point(fid[0],fid[1]).within(self.vignetting_disk.buffer(20)): # limit fiducials outside vigR to inter modules with center inside vigR
                                   
                                        self.fiducials['x'].append(np.round(fid[0],4))
                                        self.fiducials['y'].append(np.round(fid[1],4))
                                        self.fiducials['r'].append(np.round(np.sqrt(fid[0]**2 + fid[1]**2),4))
                                        self.fiducials['phi'].append(np.degrees(np.arctan2(fid[1],fid[0])))
                                        self.fiducials['z'].append(0)
                                        self.fiducials['geometry'].append(Point(fid))
                                   # self.fiducials['xyz'] = np.vstack((self.fiducials['xyz'], np.asarray(fiducial)))
          
          self.flat_grid_dict['z'] = np.zeros_like(self.flat_grid_dict['x'])
          self.flat_grid_dict = pd.DataFrame(self.flat_grid_dict)
          self.flat_grid_dict.drop_duplicates(subset=['x','y'],inplace=True)

          self.fiducials = pd.DataFrame(self.fiducials)
          self.fiducials.drop_duplicates(subset=['x','y'],inplace=True)

          return
    
    def add_neighboring_modules(self, x, y, points_up):
     
     x1 = []
     y1 =[]
     if not points_up:
          neighbours_angles = [30, 150, 270]
     else:
          neighbours_angles = [90, 210, 330]

     
     for angle in neighbours_angles:
          x_new = x + (2*self.module_side_length*np.sqrt(3)/6 + self.inner_gap) * np.cos(np.deg2rad(angle))
          y_new = y + (2*self.module_side_length*np.sqrt(3)/6 + self.inner_gap) * np.sin(np.deg2rad(angle))
          self.flat_grid_dict['x'].append(x_new)
          self.flat_grid_dict['y'].append(y_new)
          self.flat_grid_dict['tri_points_up'].append(points_up)
          x1.append(x_new)
          y1.append(y_new)

     return x1, y1


if __name__ == "__main__":
    # Example of how to use the Grid class
    

     project = "WST25"
     project_parameters = json.load(open('projects.json', 'r'))
     inner_gap = 0.5 # [mm] gap between two adjacent modules
     global_gap = 4 # [mm] gap between two adjacent modules
     from Module import Module
     mod = Module(nb_robots = 63, 
               pitch = 6.2,
               module_points_up = True)

     grid = Grid(project,
               inner_gap,
               global_gap,
               mod.module_side_length,
               **project_parameters[project])
     modules = []

     grid.flat_grid()
     print(grid.flat_grid_dict)

     LR_coverage_dict = {'geometry': [], 'color': [], 'label': []}
     HR_coverage_dict = {'geometry': [], 'color': [], 'label': []}
     boundaries = {'geometry': [], 'color': [], 'label': []}
     total_HR = 0
     total_LR = 0

     for mod_id,(x,y,z,points_up) in enumerate(zip(grid.flat_grid_dict['x'], grid.flat_grid_dict['y'], grid.flat_grid_dict['z'], grid.flat_grid_dict['tri_points_up'])):
          mod1 = Module(module_id = mod_id+1,
                    nb_robots = 63, 
                    pitch = 6.2,
                    module_points_up = points_up,
                    x0 = x,
                    y0 = y,
                    z0 = z)
          total_HR += mod1.nb_of_HR_fibers
          total_LR += mod1.nb_of_LR_fibers
          robots = mod1.robots_layout
          logging.info(f"Module {mod1.module_id}/{len(grid.flat_grid_dict['x'])}")
          modules.append(mod1)
          LR_coverage_dict['geometry'].append(mod1.LR_coverage)
          LR_coverage_dict['label'].append('LR coverage')
          LR_coverage_dict['color'].append('blue')

          HR_coverage_dict['geometry'].append(mod1.HR_coverage)
          HR_coverage_dict['label'].append('HR coverage')
          HR_coverage_dict['color'].append('red')

          boundaries['geometry'].append(mod1.module_boundaries)
          boundaries['label'].append('Module boundaries')
          boundaries['color'].append('green')

     figure, ax = plt.subplots(figsize=(10, 10))
     geo_HR = gpd.GeoDataFrame(HR_coverage_dict)
     geo_HR.plot(ax = ax, color=geo_HR['color'], alpha=0.4, legend=True)
     geo_LR = gpd.GeoDataFrame(LR_coverage_dict)
     geo_LR.plot(ax = ax, alpha=0.4, legend=True)
     geo_boundaries = gpd.GeoDataFrame(boundaries)
     geo_boundaries.plot(ax = ax, facecolor = 'None', edgecolor = boundaries['color'])
     plot_polygon(grid.vignetting_disk, ax = ax, fill = False, add_points=False, linestyle = '--', color = 'black')
     plt.legend(handles=[cl.HR_handle(), cl.LR_handle()], loc='upper right')
     plt.title(f'{grid.project} - {len(modules)} modules - {total_HR + total_LR} fibers')
     plt.xlabel('x [mm]')
     plt.ylabel('y [mm]')
     plt.grid(visible=True)
     plt.show()
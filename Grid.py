#%%
# -*- coding: utf-8 -*-
import updown_tri as tri
import numpy as np
from shapely.geometry import Point, MultiPolygon, MultiPoint, GeometryCollection
from shapely.ops import unary_union
from shapely.validation import make_valid
from shapely import concave_hull, is_valid
from shapely.plotting import plot_polygon
import pandas as pd
from Focal_Suface import FocalSurf
from matplotlib import pyplot as plt
from GFAs import GFA
import CustomLegends as cl
import json
import geopandas as gpd

import logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-4s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

class Grid():
    def __init__(self,
                 focal_surface: FocalSurf,
                 inner_gap: float, 
                 global_gap: float,
                 module_side_length: float,
                 **kwargs):

        self.surf = focal_surface
        self.project = self.surf.project_name # project name    
        self.inner_gap = inner_gap
        self.global_gap = global_gap
        self.module_side_length = module_side_length
        self.module_length = 600 # [mm] length of the module (default value, can be changed later)
        self.centered_on_triangle = kwargs.get("centered_on_triangle", False)
        self.GFA_polygon = kwargs.get("GFA_polygon", None)
        self.limiting_polygon = kwargs.get("limiting_polygon", None)
        self.trimming_angle = kwargs.get("trimming_angle", None) # [deg] angle to trim the grid in phi direction
        self.module_centroids_bounding_polygon = None
        self.n = None

        self.flat_grid_dict = {'x': [], 'y': [], 'z': [], 'tri_points_up': []}
        self.grid_3d_dict = {'x':[], 'y':[], 'z':[], 'r':[], 's':[], 'phi':[], 'theta':[], 'tri_spin':[], 'type':[], 'grid_pos':[]}
        self.fiducials = {'x': [], 'y': [], 'r': [], 'phi': [], 'z': [], 'geometry': []}

        self.R2Z, self.R2CRD, self.R2NORM, self.R2NUT, self.S2R = self.surf.transfer_functions()

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

    @property
    def fiducials_bounding_polygon(self):
        """ Returns the polygon that bounds the fiducials within the grid.
        """
        fid_lim_pol = None
        buffer = 0 # [mm] buffer to always include fiducials close to the edge
        # Check if limiting polygon already provided, otherwise use vignetting disk
        if self.limiting_polygon is not None:
            fid_lim_pol = self.limiting_polygon.difference(self.surf.donut_hole)
        else:
             fid_lim_pol = self.vignetting_disk
        # Then check if GFA polygon is provided, otherwise use vignetting disk
        if self.GFA_polygon is not None:
            fid_lim_pol = fid_lim_pol.buffer(buffer)
            fid_lim_pol = fid_lim_pol.difference(self.GFA_polygon)

        else:
            fid_lim_pol = fid_lim_pol.buffer(buffer)

         # buffer outwads to always include fiducials close to the edge

        return fid_lim_pol

    def flat_grid(self):
        # Make global grid out of triangular grid method credited in updown_tri.py
          # Its logic is also explained in updown_tri.py
          # TODO: Use two diffrent n numbers for the two different grids: modules and ficudials OR find another optimization
          self.n = int(self.surf.vigR/self.module_side_length) + 2 # first raw guess of number of lines of modules
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
          
          vigR_tresh = 10
          self.module_centroids_bounding_polygon = self.surf.vignetting_disk.buffer(vigR_tresh)
          self.module_centroids_bounding_polygon = self.module_centroids_bounding_polygon.difference(self.GFA_polygon)

          for a in np.arange(-a_max,a_max):
               for b in np.arange(-b_max,b_max):
                    for c in np.arange(-c_max,c_max):
                         sum_abc = a + b + c
                         # if valid == 1 or valid == 2: 
                         if valid.count(sum_abc): # check if sum abc corresponds to a valid triangle depending on the centering case
                              x,y = tri.tri_center(a,b,c,self.inter_triangle_side_length)
                              # Intermediate triangles placement
                              if Point(x,y).within(self.module_centroids_bounding_polygon): # allow centroid of inter modules to go out of vigR for further filling purpose
                                   self.flat_grid_dict['x'].append(x)
                                   self.flat_grid_dict['y'].append(y)
                                   # self.flat_grid_dict['geometry'].append(Point(x,y))
                                   intermediate_tri_points_up = int(tri.points_up(a,b,c, origin = origin)) # check if INTERMEDIATE triangle is up or down; True: tri points up; False: tri points down
                                   self.flat_grid_dict['tri_points_up'].append(int(not intermediate_tri_points_up)) # append inverse of points because center module of intermediate triangle is always pointing in the opposite directiob
                                   self.add_neighboring_modules(x,y,intermediate_tri_points_up) # add the three neighboring modules to the grid
                                   
                              #TODO: make smaller grid of fiducials to fill up more spaces
                              fiducial = tri.tri_corners(a,b,c,self.inter_triangle_side_length / 2)
                              # Fiducials placement
                              for fid in fiducial:
                                   if Point(fid[0],fid[1]).within(self.fiducials_bounding_polygon): # limit fiducials outside vigR to inter modules with center inside vigR
                                   
                                        self.fiducials['x'].append(np.round(fid[0],4))
                                        self.fiducials['y'].append(np.round(fid[1],4))
                                        self.fiducials['r'].append(np.round(np.sqrt(fid[0]**2 + fid[1]**2),4))
                                        self.fiducials['phi'].append(np.degrees(np.arctan2(fid[1],fid[0])))
                                        self.fiducials['z'].append(0)
                                        self.fiducials['geometry'].append(Point(fid))
                                   # self.fiducials['xyz'] = np.vstack((self.fiducials['xyz'], np.asarray(fiducial)))
          
          self.flat_grid_dict['z'] = np.zeros_like(self.flat_grid_dict['x'])
          self.flat_grid_dict = pd.DataFrame(self.flat_grid_dict)
          self.flat_grid_dict.drop_duplicates(subset=['x','y'],inplace=True) # remove duplicates to avoid problems with geometries
          self.flat_grid_dict = self.trim_grid(self.flat_grid_dict, trimming_angle=self.trimming_angle) if self.trimming_angle is not None else self.flat_grid_dict #trim grid if trimming angle provided

          self.fiducials = pd.DataFrame(self.fiducials)
          self.fiducials.drop_duplicates(subset=['x','y'],inplace=True)
          self.fiducials = self.trim_grid(self.fiducials, trimming_angle=self.trimming_angle) if self.trimming_angle is not None else self.fiducials

          return
    
    def grid_3d(self, x, y):
         """ Project flat grid to 3D to include curvature of the focal surface """
         grid_3d = {'x':[], 'y':[], 'z':[], 'r':[], 's':[], 'phi':[], 'theta':[], 'tri_spin':[], 'type':[], 'grid_pos':[]}
         grid_3d = pd.DataFrame.from_dict(grid_3d)  # Define the variable "grid_asph_pd"
         # Define transfer functions for aspherical surface i.e. Z position, theta angle and s position along the aspherical curve as function of radial position r


         grid_3d['s'] = np.sqrt(np.array(x)**2 + np.array(y)**2)
         grid_3d['phi']= np.rad2deg(np.arctan2(np.array(self.flat_grid_dict['y']), np.array(self.flat_grid_dict['x'])))
         r = self.S2R(grid_3d['s'])

         grid_3d['x'] = r * np.cos(np.deg2rad(grid_3d['phi']))
         grid_3d['dx_from_flat'] = grid_3d['x'] - np.array(self.flat_grid_dict['x'])
         grid_3d['y'] = r * np.sin(np.deg2rad(grid_3d['phi']))
         grid_3d['dy_from_flat'] = grid_3d['y'] - np.array(self.flat_grid_dict['y'])
         grid_3d['r'] = np.sqrt(grid_3d['x']**2 + grid_3d['y']**2)
         grid_3d['tri_points_up'] = np.asarray(self.flat_grid_dict['tri_points_up'])    
         grid_3d['z'] = self.R2Z(grid_3d['r'])
         grid_3d['theta'] = self.R2NUT(r)
         grid_3d['type'] = 'module' # add a column to specify the type of point (module or fiducial)
         grid_3d['grid_pos'] = 'front'
         grid_3d['geometry'] = [Point(x, y, z) for x, y, z in zip(grid_3d['x'], grid_3d['y'], grid_3d['z'])]
         grid_3d = grid_3d.round(3)
         
         return grid_3d
    
    def grid_3d_back(self, grid: pd):
         grid = grid.copy() # copy the grid to calculate the back grid position

         orientation_vectors = self.orientation_vector(np.deg2rad(grid['phi']), np.deg2rad(grid['theta'])) # get the orientation vectors of each module to project back grid in the right direction
         grid_xyz = np.vstack((grid['x'], grid['y'], grid['z'])).T # build numpy matrix for easier calculations
         grid_xyz_back = grid_xyz - self.module_length * orientation_vectors # calculate the back grid position projecting the front grid in the previously calculated orientation vectors

         grid['x'] = grid_xyz_back[:, 0]
         grid['y'] = grid_xyz_back[:, 1]
         grid['z'] = grid_xyz_back[:, 2]
         grid['r'] = np.sqrt(grid['x']**2 + grid['y']**2)
         grid['geometry'] = [Point(x, y, z) for x, y, z in zip(grid['x'], grid['y'], grid['z'])]
         grid['grid_pos'] = 'back'
         
         return grid
    
    def layout_concave_hull(self):

          """
          Caluclate the exterior polygon of the layout by taking the concave hull of the fiducials
          
          Input:
               - fiducials: [DataFrame] contains the x,y,z coordinates of the fiducials

          Output:
               - ch: [shapely Polygon] contains the exterior polygon of the layout  
          
          """
          
          self.fiducials.sort_values(by=['r'], inplace=True)
          fiducials = MultiPoint(self.fiducials['geometry'].to_list())
          if not is_valid(fiducials):
               fiducials = make_valid(fiducials)

          try :     
               ch = concave_hull(fiducials, ratio=0.2, allow_holes=True)
               if not is_valid(ch):
                    ch = make_valid(ch)
                    if type(ch) == GeometryCollection:
                         ch = ch.geoms[0]
               
               # if 'WST' in self.project:
               #      ch = ch.difference(self.donut_hole)
               return ch
          except: 
               return fiducials.convex_hull
          
          

    
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
          self.flat_grid_dict['tri_points_up'].append(int(points_up))
          x1.append(x_new)
          y1.append(y_new)

     return x1, y1
    
    
    def orientation_vector(self, phi: float, theta: float):

          """
          Input:

          - phi: [float] azimuthal angle in spherical coordinates
          - theta: [float] polar angle in spherical coordinates

          Output:

          - orientation_vector: [3x1 numpy array] contains the x,y,z coordinates of the orientation vector
          """

          x = - np.sin(theta) * np.cos(phi)
          y = - np.sin(theta) * np.sin(phi)
          z = np.cos(theta)

          return np.array([x,y,z]).T
    
    def trim_grid(self, grid : pd.DataFrame, trimming_angle: float = 0):
         #     index2drop = grid[(grid['phi'] > trimming_angle) & (grid['phi'] < 0)].index
         grid['phi'] =  np.rad2deg(np.arctan2(np.array(grid['y']), np.array(grid['x'])))
         trimmed_grid = grid[(0 <= grid['phi']) & (grid['phi'] <= trimming_angle)]
         
         return trimmed_grid


if __name__ == "__main__":
    # Example of how to use the Grid class
    

     PROJECT = "VLT_2030"
     project_parameters = json.load(open('projects.json', 'r'))
     INNER_GAP = 0.5 # [mm] gap between two adjacent modules
     GLOBAL_GAP = 4 # [mm] gap between two adjacent modules
     from Module import Module
     mod = Module(nb_robots = 63, 
               pitch = 6.2,
               module_points_up = True)
     
     nb_gfa = 6
     angle_offset = 30
     gfa_length = 60
     gfa_width = 60
     gfa = GFA(nb_gfa = nb_gfa,
          angle_offset = angle_offset,
          vigR = project_parameters[PROJECT]['vigD'] / 2,
          length = gfa_length,
          width = gfa_width)
     gdf_gfa = gfa.gdf_gfa
     polygon_gfa = MultiPolygon(list(gdf_gfa['geometry']))

     # limiting_polygon = surf.trimming_polygon(geometry='hex', trim_diff_to_vigR = 10)
     surf = FocalSurf(PROJECT, **project_parameters[PROJECT])
     limiting_polygon = surf.vignetting_disk

     grid = Grid(focal_surface = surf,
               inner_gap = INNER_GAP,
               global_gap = GLOBAL_GAP,
               module_side_length = mod.module_side_length,
               GFA_polygon = polygon_gfa,
               limiting_polygon = limiting_polygon,
               **project_parameters[PROJECT])

     modules = []

     grid.flat_grid()
     grid_3d = grid.grid_3d(grid.flat_grid_dict['x'], grid.flat_grid_dict['y'])
     grid_3d_back = grid.grid_3d_back(grid_3d)
     print(grid.fiducials)

#%%
     # project2 =  'VLT_2030'
     # grid2 = Grid(project2,
     #           inner_gap,
     #           global_gap,
     #           mod.module_side_length,
     #           GFA_polygon = polygon_gfa,
     #           **project_parameters[project2])
     
     # grid2.flat_grid()
     # grid2_3d = grid2.grid_3d(grid2.flat_grid_dict['x'], grid2.flat_grid_dict['y'])


     LR_coverage_dict = {'geometry': [], 'color': [], 'label': []}
     HR_coverage_dict = {'geometry': [], 'color': [], 'label': []}
     boundaries = {'geometry': [], 'color': [], 'label': []}
     total_HR = 0
     total_LR = 0
     # print(f'Number of modules: {len(grid.flat_grid_dict)}')
     print(f'Number of fiducials: {len(grid.fiducials["x"])}')

     figure, ax = plt.subplots(figsize=(10, 10))
     # geo_HR = gpd.GeoDataFrame(HR_coverage_dict)
     # geo_HR.plot(ax = ax, color=geo_HR['color'], alpha=0.4, legend=True)
     # geo_LR = gpd.GeoDataFrame(LR_coverage_dict)
     # geo_LR.plot(ax = ax, alpha=0.4, legend=True)
     # geo_boundaries = gpd.GeoDataFrame(boundaries)
     # geo_boundaries.plot(ax = ax, facecolor = 'None', edgecolor = boundaries['color'])
     plt.scatter(grid.flat_grid_dict['x'], grid.flat_grid_dict['y'], color='blue', alpha=0.4, label='Modules')
     plot_polygon(grid.surf.vignetting_disk, ax = ax, fill = False, add_points=False, linestyle = '--', color = 'black', label = 'Vignetting disk')
     plot_polygon(polygon_gfa, ax = ax, fill = True, add_points=False, alpha = 0.3, color = 'orange', label = 'GFAs')
     plot_polygon(grid.module_centroids_bounding_polygon, ax = ax, fill = False, add_points=False, linestyle = '--', color = 'green', label = 'Modules boundary')
     plot_polygon(grid.fiducials_bounding_polygon, ax = ax, fill = False, add_points=False, linestyle = '--', color = 'red', label = 'Fiducials boundary')
     plt.scatter(grid.fiducials['x'], grid.fiducials['y'], color='red', alpha=0.4, label='Fiducials')
     plt.legend(loc='upper right')
     plt.title(f'{grid.project} - {len(modules)} modules - {total_HR + total_LR} fibers')
     plt.xlabel('x [mm]')
     plt.ylabel('y [mm]')
     plt.grid(visible=True)
     ax.axis('equal')


     fig = plt.figure(figsize=(10, 10))
     ax = fig.add_subplot(projection='3d')

     ax.scatter(grid_3d['x'], grid_3d['y'], grid_3d['z'], c='blue', alpha=0.4, label=f'grid_front')
     ax.scatter(grid_3d_back['x'], grid_3d_back['y'], grid_3d_back['z'], c='red', alpha=0.4, label=f'grid_back')
     for x_start, y_start, z_start, x_end, y_end, z_end in zip(grid_3d['x'], grid_3d['y'], grid_3d['z'], grid_3d_back['x'], grid_3d_back['y'], grid_3d_back['z']):
          ax.plot([x_start, x_end], [y_start, y_end], [z_start, z_end], c='gray', alpha=0.2)
     ax.set_xlabel('X [mm]')
     ax.set_ylabel('Y [mm]')
     ax.set_zlabel('Z [mm]')
     ax.set_title(f'{grid.project} - 3D grid')
     ax.legend()
     ax.set_box_aspect((5,5,1))
     plt.show()
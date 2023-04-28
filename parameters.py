#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import logging
from shapely import affinity, MultiPolygon, MultiPoint, GeometryCollection
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from shapely.ops import unary_union
import time
from shapely.plotting import plot_polygon, plot_points
import os
from datetime import datetime
import alphashape # WARNING: check this git issue to solve compatibility problem between shapely 2.0 and alphashape: https://github.com/bellockk/alphashape/issues/36 

logging.basicConfig(level=logging.INFO)
tan30 = np.tan(np.deg2rad(30))

## Parameters ##

""" Focal plane parameters """

vigR = 613.27 # [mm] radius of the instrument
curve_radius = 11067 # [mm] curvature of the focal plane


# """ Intermediate frame parameters """

# intermediate_frame_thick =  1 # [mm] spacing between modules inside intermediate frame

# """ Global frame parameters """

# global_frame_thick = 3 # [mm] spacing between modules in global arrangement

""" Module parameters """ 

class Module:

     def __init__(self, nb_robots, width_increase = 0, chanfer_length = 7.5):

          """Robot parameters""" 

          self.nb_robots = nb_robots
          self.la = 1.8 # alpha arm length [mm] /!\ lb > la /!\
          self.lb = 1.8 # beta arm length [mm]
          self.pitch = 6.2 # pitch [mm]
          self.alpha = np.linspace(-180,180,180) # [deg] rotational range of alpha arm
          self.beta = np.linspace(-180,180,180) # [deg] rotational range of beta arm (holding the fiber)
          self.beta2fibre = 1 # [mm] distance fiber center to edge of beta arm

          """Module parameters"""

          self.chanfer_length = chanfer_length # [mm] Length of chanfer on module vertices
          self.width_increase = width_increase # [mm] Increase module side length w/o changing nb of pos
          
          self.x_first = 5.9 + self.width_increase/2 # [mm] Horizontal pos of first robot (bottom left)
          self.y_first = 3.41 + self.width_increase/2 * tan30 # [mm] Vertical pos of first robot (bottom left)
          
          self.x_inc = 3.1 # [mm] Horizontal increment at each row
          self.y_inc = 5.369 # [mm] Vertical increment at each row
          self.test_pitch = np.linalg.norm(np.array([self.x_inc,self.y_inc]))

          self.safety_distance = 0.4 # [mm] physical distance kept between shields and edge of beta arm
          self.offset_from_module_edges = self.safety_distance + self.beta2fibre
          self.start_offset_x = 6.2 # [mm]
          self.start_offset_y = 3.41 # [mm]

          self.is_wall = True # flag for protective shields or not on modules

          # 1 row addition from one case to another 
          if self.nb_robots == 52:

               self.module_width = 67.6 + self.width_increase # [mm] triangle side length
               self.nb_rows = 9 # number of rows of positioners
               self.key = 'n52'

          elif self.nb_robots == 63:
               
               self.module_width = 73.8 + self.width_increase# [mm] triangle side length
               self.nb_rows = 10 # number of rows of positioners
               self.key = 'n63'

          elif self.nb_robots == 75:

               self.module_width = 80 + self.width_increase # [mm] triangle side length
               self.nb_rows = 11 # number of rows of positioners
               self.key = 'n75'

          elif self.nb_robots == 88:
               self.module_width = 86.2 + self.width_increase # [mm] triangle side length
               self.nb_rows = 12 # number of rows of positioners
               self.key = 'n88'

          elif self.nb_robots == 102:
               
               self.module_width = 92.4 + self.width_increase # [mm] triangle side length
               self.nb_rows = 13 # number of rows of positioners
               self.key = 'n102'

          else:
               raise Exception('Error: only 52, 63, 75, 88, 102 robots per module supported')
          
     def equilateral_vertices(self, triangle_width, xoff = 0, yoff = 0):
          """
          Create the raw triangular base of a module w/ the correct width

          Output: sets of x,y coords of the triangle's vertices
          """
          x_vertices = np.ones(4)
          offset = [xoff, yoff]
          x_vertices[0] = offset[0]
          x_vertices[1] = x_vertices[0] + triangle_width
          x_vertices[2] = x_vertices[0] + triangle_width * np.cos(np.deg2rad(60))
          x_vertices[3] = x_vertices[0]


          y_vertices = np.ones(4)
          y_vertices[0] = offset[1]
          y_vertices[1] = y_vertices[0]
          y_vertices[2] = y_vertices[0] + triangle_width * np.sin(np.deg2rad(60))
          y_vertices[3] = y_vertices[0]

          return x_vertices,y_vertices

     def chanfered_base(self):
          """
          1) Creation of small triangles that will serve to cut the raw triangle for chanfering
          2) Cuts triangular base with small triangles to create the chanfers
          
          Output: sets of x,y coords of the chanfered module's vertices
          """
          
          x_vertices,y_vertices = self.equilateral_vertices(self.module_width)
          module_triangle = Polygon(to_polygon_format(x_vertices,y_vertices))

          x_chanfer,y_chanfer = self.equilateral_vertices(self.chanfer_length)

          # First chanfer created from the origin
          chanfer_triangle = Polygon(to_polygon_format(x_chanfer-0.0001,y_chanfer-0.0001)) # tiny diff for geometry to intersect in one point rather tha infinite points (there is a better solution)

          chanfers = [chanfer_triangle]
          angles = [120, 240]
          module_centroid = module_triangle.centroid
          
          for angle in angles:
               # Move the 2 remaining triangles to their corresponding location at the triangle's edges
               rot_chanfer_triangle = rotate_and_translate(chanfer_triangle, angle, 0, 0, origin = module_centroid)
               chanfers.append(rot_chanfer_triangle)
               chanfered_base = module_triangle.difference(rot_chanfer_triangle)
               

          multi_chanfers = MultiPolygon(chanfers)
          # Apply cut on triangular base
          chanfered_base = module_triangle.difference(multi_chanfers)
          
          return chanfered_base
     
     def remove_positioner(self, xx, yy, list_to_remove):
          """ Input:
               - xx, yy: grid of center of positioners
               - list_to_remove: list of positioners position to remove

               Output:
               - xx_sliced, yy_sliced: grid of center of positioners minus the ones removed
          """
          xx_sliced = np.delete(xx,list_to_remove)
          yy_sliced = np.delete(yy,list_to_remove)

          return xx_sliced, yy_sliced
     
     def create_module(self):
          
          module = self.chanfered_base() # Create chamfered module shape as shapely Polygon
          # coords_module_x, coords_module_y = module.exterior.coords.xy
          reference_centroid = module.centroid

          # Their norm corresponds to the pitch defined in parameters
          xx1 = np.ones(self.nb_rows + 1)
          yy1 = np.ones(self.nb_rows + 1)

          for idx in range(self.nb_rows): # Create the grid of positioner's center points
          
               # Place the first robot of each row
               start_offset_x = self.x_inc * idx + self.x_first
               start_offset_y = self.y_inc * idx + self.y_first

                    # Generate each row of robot: nb robots/row decreases as going upward in triangle
               xx_new = np.linspace(start_offset_x, (self.nb_rows - idx ) * self.pitch + start_offset_x, self.nb_rows - idx + 1)
               yy_new = start_offset_y * np.ones(len(xx_new))

               if idx == 0:
                    xx1 = xx_new
                    yy1 = yy_new
               else: 
                    xx1 = np.hstack((xx1, xx_new))
                    yy1 = np.hstack((yy1, yy_new))


               if idx == self.nb_rows - 1:
                    xx_new = np.array([self.x_inc * (idx + 1)])
                    yy_new = self.y_inc * (idx + 1) * np.ones(len(xx_new))
                    xx1 = np.hstack((xx1, xx_new))
                    yy1 = np.hstack((yy1, yy_new))

          list_to_remove = [0, self.nb_rows, -1] # remove positioners at the edges of the triangle
          # list_to_remove = [] # remove no positioner
          xx1, yy1 = self.remove_positioner(xx1, yy1, list_to_remove)
          nbots = len(xx1)


          triang_meshgrid = MultiPoint(to_polygon_format(xx1, yy1)) # convert the meshgrid into shapely standard for later manipulation

          """ 1)b) Define coverage for 1 modules """

          c1 = np.cos(np.deg2rad(self.alpha))
          s1 = np.sin(np.deg2rad(self.alpha))

          c2 = np.cos(np.deg2rad(self.beta))
          s2 = np.sin(np.deg2rad(self.beta))

          xa, ya, xab, yab = (self.lb-self.la)*c1, (self.lb-self.la)*s1, (self.lb+self.la)*c2, (self.la+self.lb)*s2

          wks_list = []

          for idx, (dx, dy) in enumerate(zip(xx1, yy1)):

               xa1, ya1, xab1, yab1 = xa + dx, ya + dy, xab + dx, yab + dy

               coords_int = to_polygon_format(xa1, ya1)
               coords_ext = to_polygon_format(xab1, yab1)

               interior = coords_int[::-1]
               poly_c1 = Polygon(coords_ext, [interior])
               wks_list.append(poly_c1)

          multi_wks = MultiPolygon(wks_list)
          total_positioners_workspace = unary_union(wks_list)
          module_w_beta_and_safety_dist = module.buffer(-self.offset_from_module_edges)

          if self.is_wall:
               effective_wks = module_w_beta_and_safety_dist.intersection(total_positioners_workspace)
          else:
               effective_wks = total_positioners_workspace


          module_collection = GeometryCollection([module, module_w_beta_and_safety_dist, effective_wks, triang_meshgrid])
          module_collection = affinity.translate(module_collection, xoff = -reference_centroid.x, yoff = -reference_centroid.y)
          multi_wks_list = affinity.translate(multi_wks, xoff = -reference_centroid.x, yoff = -reference_centroid.y)
          coverage_with_walls = round(effective_wks.area/module.area * 100,1)
          coverage_no_walls = round(total_positioners_workspace.area/module.area * 100,1)
          coverages = [coverage_with_walls, coverage_no_walls]

          logging.info(f'Module object created with width {self.module_width} & {nbots} robots')
          return module_collection, multi_wks_list, coverages

class IntermediateTriangle: 

     def __init__(self, module_collection, module_width, intermediate_frame_thick):

          self.module_collection = module_collection
          self.dist_inter = 2*module_width*np.sqrt(3)/6 + intermediate_frame_thick # distance between each neighbor from center module
          self.angles = np.array([-30, 90, 210])
          self.flip = [True,False,False,False]

          self.x_grid_inter = np.cos(np.deg2rad(self.angles))*self.dist_inter
          self.x_grid_inter = np.insert(self.x_grid_inter, 0, 0)
          self.y_grid_inter = np.sin(np.deg2rad(self.angles))*self.dist_inter
          self.y_grid_inter = np.insert(self.y_grid_inter, 0, 0)

     def create_intermediate_triangle(self):
          
          boundaries = []
          boundary_xx = []
          boundary_yy = []
          intermediate_collection = []
          inter_coverage = []
          covered_area = 0
          inter_boundaries_df = {'name':[],'geometry':[], 'color': []}
          inter_modules_df = {'geometry':[]}
          inter_coverage_df = {'geometry':[]}
          inter_robots_df = {'geometry':[]}
          

          for idx, (rotate, dx, dy) in enumerate(zip(self.flip, self.x_grid_inter, self.y_grid_inter)):

               if rotate:
                    angle = 180
               else:
                    angle = 0

               transformed_all = rotate_and_translate(self.module_collection, angle, dx, dy, origin = "centroid") # place module at the correct pos and orientation
               boundaries.append(transformed_all.geoms[0]) # log the exterior boundary points of the transformed module
               inter_modules_df['geometry'].append(transformed_all.geoms[0])
               inter_coverage.append(transformed_all.geoms[2])
               inter_coverage_df['geometry'].append(transformed_all.geoms[2])
               intermediate_collection.append(transformed_all)
               inter_robots_df['geometry'].append(transformed_all.geoms[3])
               
               xx,yy = transformed_all.geoms[0].exterior.coords.xy
               boundary_xx.append(xx.tolist())
               boundary_yy.append(yy.tolist())


          modules_polygon_intermediate = MultiPolygon(boundaries)
          # Convex hull of boundaries
          bounding_polygon_intermediate = modules_polygon_intermediate.convex_hull
          inter_boundaries_df['name'].append('Inter_convex_hull')
          inter_boundaries_df['geometry'].append(modules_polygon_intermediate.convex_hull)
          inter_boundaries_df['color'].append('green')
          # Concave hull of boundaries
          int_cent = np.array(bounding_polygon_intermediate.centroid.xy).reshape((1,2))
          boundary_xx = flatten_list(boundary_xx)
          boundary_yy = flatten_list(boundary_yy)
          pol_points = sort_points_for_polygon_format(boundary_xx,boundary_yy, int_cent)
          concave_hull = Polygon(pol_points)
          inter_boundaries_df['name'].append('Inter_concave_hull')
          inter_boundaries_df['geometry'].append(concave_hull)
          inter_boundaries_df['color'].append('orange')
          coverage_polygon_intermediate = MultiPolygon(inter_coverage) # save coverage as whole to speedup global calculation
          # plot_polygon(coverage_polygon_intermediate, add_points=False)
          covered_area_inter = coverage_polygon_intermediate.area
          intermediate_collection.append(bounding_polygon_intermediate)
          intermediate_collection = GeometryCollection(intermediate_collection)
          intermediate_collection_speed = GeometryCollection([bounding_polygon_intermediate, modules_polygon_intermediate, coverage_polygon_intermediate, concave_hull])
          available_intermediate_area = round(bounding_polygon_intermediate.area,1)
          intermediate_coverage = round(covered_area_inter/available_intermediate_area*100,1)
          inter_df = {'inter_boundaries': inter_boundaries_df, 'inter_modules': inter_modules_df, 'inter_coverage': inter_coverage_df, 'inter_robots': inter_robots_df, 'intermediate_coverage': intermediate_coverage}

          return intermediate_collection, intermediate_collection_speed, intermediate_coverage, inter_df
     
def to_polygon_format(x,y):
     """ Input:
          - x,y: 2 sets of coordinates for polygon creation

          Output:
          - coords: list of tuples for each polygon vertices (just for convenience) """
     coords = []
     for (i,j) in zip(x,y):
               coords.append((i,j))
               
     return coords

def rotate_and_translate(geom, angle, dx, dy, dz = None, origin = 'centroid'):

     rotated_geom = affinity.rotate(geom, angle, origin=origin)
     transformed_geom = affinity.translate(rotated_geom, dx, dy, dz)

     return transformed_geom


def save_figures_to_dir(save, suffix_name):

     if not save:
           return
     
     script_dir = os.path.dirname(__file__)
     results_dir = os.path.join(script_dir, 'Results/')

     now = datetime.now()
     today_filename = now.strftime("%Y-%m-%d--%H-%M-%S_") + suffix_name + ".png"


     if not os.path.isdir(results_dir):
          os.makedirs(results_dir)

     plt.savefig(results_dir + today_filename, bbox_inches = 'tight', format='png', dpi = 1200)

def make_vigR_polygon(pizza_angle = 360, r = vigR):
      
     n_vigR = 500
     vigR_lim_x = r * np.cos(np.deg2rad(np.linspace(0,pizza_angle,n_vigR)))
     vigR_lim_y = r * np.sin(np.deg2rad(np.linspace(0,pizza_angle,n_vigR)))
     if pizza_angle == 360:
          end_point = [vigR_lim_x[0], vigR_lim_y[0]]
     else:
          end_point = [0, 0]
     vigR_lim_x = np.insert(vigR_lim_x, 0, end_point[0])
     vigR_lim_y = np.insert(vigR_lim_y, 0, end_point[1])
     pizza = Polygon(to_polygon_format(vigR_lim_x, vigR_lim_y))

     return pizza

def plot_vigR_poly(pizza, label = None):
     plot_polygon(pizza, add_points = False, edgecolor = 'black', linestyle = '--', facecolor= 'None', label = label)

def plot_module(module_collection, label_coverage, label_robots, ignore_points):
     for jdx, geo in enumerate(module_collection.geoms):
          # plot_polygon(geometry[jdx], add_points=False, facecolor='None' , edgecolor='black')
          if (isinstance (module_collection.geoms[jdx], Polygon)) and jdx == 0:
               plot_polygon(module_collection.geoms[jdx], add_points=False, facecolor='None' , edgecolor='black')
          elif (isinstance (module_collection.geoms[jdx], Polygon)) and jdx == 1:
               plot_polygon(module_collection.geoms[jdx], add_points=False, facecolor='None', linestyle = '--')
          elif (isinstance (module_collection.geoms[jdx], Polygon)) and jdx == 2:
               plot_polygon(module_collection.geoms[jdx], add_points=False, alpha=0.2, edgecolor='black', label = label_coverage)
          elif (isinstance (module_collection.geoms[jdx], MultiPoint) and not ignore_points):
               plot_points(module_collection.geoms[jdx], marker='.', color='k', label = label_robots)

def plot_intermediate(intermediate_collection, nb_robots, ignore_points, intermediate_coverage=None, available_intermediate_area=None, draw_legend = False):
     for idx, mod_collection in enumerate(intermediate_collection.geoms):
          if idx == 0 and draw_legend:
               label_coverage = 'Coverage: {} %'.format(intermediate_coverage)
               # label_coverage = 'Coverage: {} %'.format(area_to_cover)
               label_robots = "{} robots".format(nb_robots)
          else: 
               label_coverage = None
               label_robots = None
          if (isinstance (mod_collection, Polygon)):
               plot_polygon(mod_collection, add_points=False, facecolor='None', linestyle = '-.', color = 'green')

               # plot_polygon(mod_collection, add_points=False, facecolor='None', linestyle = '-.', color = 'green', label = 'Available area: {} mm$^2$'.format(available_intermediate_area))
               continue
          plot_module(mod_collection, label_coverage, label_robots, ignore_points)

def plot_intermediate_speed(mod_collection, label_coverage):
          s11=time.time()
          if (isinstance (mod_collection.geoms[0], (Polygon, MultiPolygon))):
               plot_polygon(mod_collection.geoms[0], add_points=False, facecolor='None' , linestyle = '--')
          else:
               print(f'mod_collection.geoms[0] is {type(mod_collection.geoms[1])}')
          if (isinstance (mod_collection.geoms[1], MultiPolygon)):
               plot_polygon(mod_collection.geoms[1], add_points=False, facecolor='None' , edgecolor='black')
          else:
               print(f'mod_collection.geoms[1] is {type(mod_collection.geoms[2])}')
          if (isinstance (mod_collection.geoms[2], MultiPolygon)):
               plot_polygon(mod_collection.geoms[2], add_points=False, alpha=0.2, edgecolor='black', label = label_coverage)
          s12=time.time()
          # print(f"0: {s12-s11} s")
               
def final_title(nb_robots, total_modules, total_robots, inter_frame_thick, global_frame_thick, allow_small_out, out_allowance, disp_robots_info = True):


     if allow_small_out:
          small_out_info = f"Out allowance: {out_allowance * 100:0.1f} %"
     else:
          small_out_info = ''

     if disp_robots_info:
          robots_info = f"\n Total # modules: {total_modules} - Total # robots: {total_robots}"
          modules_info = f" {nb_robots} robots per module"
     else:
          robots_info = ''
          modules_info = ''
     if inter_frame_thick != global_frame_thick:
          figtitle = f"Semi frameless - {modules_info} \n Inner gap: {inter_frame_thick} mm - Global gap: {global_frame_thick} mm {robots_info} \n {small_out_info} \n"
     elif inter_frame_thick == global_frame_thick and global_frame_thick == 0:
          figtitle = f"Frameless - {modules_info} {robots_info} \n {small_out_info} \n"
     else:
          figtitle = f"Framed - {modules_info} \n Gap: {inter_frame_thick} mm {robots_info} \n {small_out_info} \n"

     return figtitle

def flatten_list(l:list):
     return [item for sublist in l for item in sublist]

def makeKey(key_int: int):
     return f'n{key_int}'

def sort_points_for_polygon_format(x: list, y: list, centroid):
     """
     Sort points (np array) counterclockwise for later polygon creation
     """
     points = []
     for (xi,yi) in zip(x, y):
          points.append([xi,yi])
     points = np.array(points)
     vec = points - centroid # get vector connecting centroid of points to them
     angles = np.arctan2(vec[:,1],vec[:,0]) # calculate angle of each vector
     d = np.hstack((vec, angles.reshape(len(angles),1))) # put everything in one array
     d = d[d[:, 2].argsort()] # sort the points by ascending order array
     return to_polygon_format(d[:,0] + centroid[0,0], d[:,1]+ centroid[0,1])

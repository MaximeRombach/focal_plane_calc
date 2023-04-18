#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from shapely import affinity, MultiPolygon, MultiPoint, GeometryCollection
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from shapely.ops import unary_union
import time
from shapely.plotting import plot_polygon, plot_points
import os
from datetime import datetime

## Parameters ##

""" Focal plane parameters """

vigR = 613.27 # [mm] radius of the instrument
curve_radius = 11067 # [mm] curvature of the focal plane

"""Positioner parameters""" 

la = 1.8 # alpha arm length [mm] /!\ lb > la /!\
lb = 1.8 # beta arm length [mm]
pitch = 6.2 # pitch [mm]
alpha = np.linspace(-180,180,180) # [deg] rotational range of alpha arm
beta = np.linspace(-180,180,180) # [deg] rotational range of beta arm (holding the fiber)
beta2fibre = 1 # [mm] distance fiber center to edge of beta arm

x_first = 5.9 # [mm] Horizontal pos of first robot (bottom left)
y_first = 3.41 # [mm] Vertical pos of first robot (bottom left)
x_inc = 3.1 # [mm] Horizontal increment at each row
y_inc = 5.369 # [mm] Vertical increment at each row
test_pitch = np.linalg.norm(np.array([x_inc,y_inc]))

""" Module parameters """ 

class Module:

     def __init__(self, nb_robots):
          
          self.nb_robots = nb_robots
          self.la = 1.8 # alpha arm length [mm] /!\ lb > la /!\
          self.lb = 1.8 # beta arm length [mm]
          self.pitch = 6.2 # pitch [mm]
          self.alpha = np.linspace(-180,180,180) # [deg] rotational range of alpha arm
          self.beta = np.linspace(-180,180,180) # [deg] rotational range of beta arm (holding the fiber)
          self.beta2fibre = 1 # [mm] distance fiber center to edge of beta arm

          self.x_first = 5.9 # [mm] Horizontal pos of first robot (bottom left)
          self.y_first = 3.41 # [mm] Vertical pos of first robot (bottom left)
          self.x_inc = 3.1 # [mm] Horizontal increment at each row
          self.y_inc = 5.369 # [mm] Vertical increment at each row
          self.test_pitch = np.linalg.norm(np.array([x_inc,y_inc]))

          self.safety_distance = 0.5 # [mm] physical distance kept between shields and edge of beta arm
          self.offset_from_module_edges = self.safety_distance + beta2fibre
          self.start_offset_x = 6.2 # [mm]
          self.start_offset_y = 3.41 # [mm]

          self.is_wall = True # flag for protective shields or not on modules

          if self.nb_robots == 75:

               self.module_width = 80 # [mm] triangle side length
               self.nb_rows = 11 # number of rows of positioners

          elif self.nb_robots == 63:
               
               self.module_width = 73.8# [mm] triangle side length
               self.nb_rows = 10 # number of rows of positioners

          elif self.nb_robots == 102:
               
               self.module_width = 92.4 # [mm] triangle side length
               self.nb_rows = 13 # number of rows of positioners

          else:
               raise Exception('Error: only 63, 75, 102 robots per module supported')
          
     def equilateral_vertices(self, triangle_width, xoff = 0, yoff = 0):
          
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

     def chanfered_base(self, chanfer_length = 7.5):
          
          x_vertices,y_vertices = self.equilateral_vertices(self.module_width)
          module_triangle = Polygon(to_polygon_format(x_vertices,y_vertices))

          x_chanfer,y_chanfer = self.equilateral_vertices(chanfer_length)


          chanfer_triangle = Polygon(to_polygon_format(x_chanfer-0.0001,y_chanfer-0.0001)) # tiny diff for geometry to intersect in one point rather tha infinite points (there is a better solution)

          chanfers = [chanfer_triangle]
          angles = [120, 240]
          module_centroid = module_triangle.centroid
          
          for angle in angles:
               rot_chanfer_triangle = self.rotate_and_translate(chanfer_triangle, angle, 0, 0, origin = module_centroid)
               chanfers.append(rot_chanfer_triangle)
               chanfered_base = module_triangle.difference(rot_chanfer_triangle)
               

          multi_chanfers = MultiPolygon(chanfers)

          plt.show()
          chanfered_base = module_triangle.difference(multi_chanfers)
          
          return chanfered_base
     
     def rotate_and_translate(self, geom, angle, dx, dy, dz = None, origin = 'centroid'):

          rotated_geom = affinity.rotate(geom, angle, origin=origin)
          transformed_geom = affinity.translate(rotated_geom, dx, dy, dz)

          return transformed_geom
     
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
     
     def to_polygon_format(self, x,y):
          """ Input:
               - x,y: 2 sets of coordinates for polygon creation

               Output:
               - coords: list of tuples for each polygon vertices (just for convenience) """
          coords = []
          for (i,j) in zip(x,y):
                    coords.append((i,j))
                    
          return coords
     
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
          nb_robots = len(xx1)

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

               coords_int = self.to_polygon_format(xa1, ya1)
               coords_ext = self.to_polygon_format(xab1, yab1)

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
          wks_list = affinity.translate(multi_wks, xoff = -reference_centroid.x, yoff = -reference_centroid.y)
          coverage_with_walls = round(effective_wks.area/module.area * 100,1)
          coverage_no_walls = round(total_positioners_workspace.area/module.area * 100,1)
          coverages = [coverage_with_walls, coverage_no_walls]

          return module_collection, wks_list, coverages


# Raw triangle
module_vertices_x = np.array([0,80,40,0])
module_vertices_y = np.array([0,0,69.3,0])

# Edge cut module
module_vertices_x = np.array([7.5, 72.5, 76.25, 43.75, 36.25, 3.75, 7.5]) # [mm]
module_vertices_y = np.array([0, 0, 6.5, 62.8, 62.8, 6.5, 0]) # [mm]



""" Intermediate frame parameters """

intermediate_frame_thick =  3 # [mm] spacing between modules inside intermediate frame

""" Global frame parameters """

global_frame_thick = 3 # [mm] spacing between modules in global arrangement


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

     plt.savefig(results_dir + today_filename, bbox_inches = 'tight')

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
               plot_polygon(mod_collection, add_points=False, facecolor='None', linestyle = '-.', color = 'green', label = 'Available area: {} mm$^2$'.format(available_intermediate_area))
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
               

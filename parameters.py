#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from shapely import affinity, MultiPolygon, MultiPoint
import matplotlib.pyplot as plt
from shapely.geometry import Polygon

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

class module_param:

     def __init__(self, nb_robots):
          
          self.nb_robots = nb_robots

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


# Raw triangle
module_vertices_x = np.array([0,80,40,0])
module_vertices_y = np.array([0,0,69.3,0])

# Edge cut module
module_vertices_x = np.array([7.5, 72.5, 76.25, 43.75, 36.25, 3.75, 7.5]) # [mm]
module_vertices_y = np.array([0, 0, 6.5, 62.8, 62.8, 6.5, 0]) # [mm]

safety_distance = 0.5 # [mm] physical distance kept between shields and edge of beta arm
offset_from_module_edges = safety_distance + beta2fibre
start_offset_x = 6.2 # [mm]
start_offset_y = 3.41 # [mm]

is_wall = True # flag for protective shields or not on modules

""" Intermediate frame parameters """

intermediate_frame_thick =  3 # [mm] spacing between modules inside intermediate frame

""" Global frame parameters """

global_frame_thick = 3 # [mm] spacing between modules in global arrangement

def remove_positioner(xx,yy, list_to_remove):
     """ Input:
          - xx, yy: grid of center of positioners
          - list_to_remove: list of positioners position to remove

          Output:
          - xx_sliced, yy_sliced: grid of center of positioners minus the ones removed
     """
     xx_sliced = np.delete(xx,list_to_remove)
     yy_sliced = np.delete(yy,list_to_remove)

     return xx_sliced, yy_sliced

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

def equilateral_vertices(module_width, offset = [0,0]):
      
     x_vertices = np.ones(4)
     x_vertices[0] = offset[0]
     x_vertices[1] = x_vertices[0] + module_width
     x_vertices[2] = x_vertices[0] + module_width * np.cos(np.deg2rad(60))
     x_vertices[3] = x_vertices[0]


     y_vertices = np.ones(4)
     y_vertices[0] = offset[1]
     y_vertices[1] = y_vertices[0]
     y_vertices[2] = y_vertices[0] + module_width * np.sin(np.deg2rad(60))
     y_vertices[3] = y_vertices[0]

     return x_vertices,y_vertices

def chanfered_base(module_width, chanfer_length = 7.5):

     x_vertices,y_vertices = equilateral_vertices(module_width)
     module_triangle = Polygon(to_polygon_format(x_vertices,y_vertices))
     x_vertices,y_vertices = equilateral_vertices(chanfer_length)
     chanfer_triangle = Polygon(to_polygon_format(x_vertices-0.0001,y_vertices-0.0001)) # tiny diff for geometry to intersect in one point rather tha infinite points (there is a better solution)
     chanfers = [chanfer_triangle]
     angles = [120, 240]
     module_centroid = module_triangle.centroid
     for angle in angles:
          rot_chanfer_triangle = rotate_and_translate(chanfer_triangle, angle, 0, 0, origin = module_centroid)
          chanfers.append(rot_chanfer_triangle)
          chanfered_base = module_triangle.difference(rot_chanfer_triangle)
          

     multi_chanfers = MultiPolygon(chanfers)
     chanfered_base = module_triangle.difference(multi_chanfers)
     
     return chanfered_base

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
               

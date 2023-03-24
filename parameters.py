#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from shapely import affinity, MultiPolygon
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.plotting import plot_polygon
import os
from datetime import datetime

## Parameters ##

""" Focal plane parameters """

vigR = 613.27
curve_radius = 11067

"""Positioner parameters""" 

la = 1.8 # alpha arm length [mm] /!\ lb > la /!\
lb = 1.8 # beta arm length [mm]
pitch = 6.2 # pitch [mm]
alpha = np.linspace(-180,180,180) # [deg] rotational range of alpha arm
beta = np.linspace(-180,180,180) # [deg] rotational range of beta arm (holding the fiber)
beta2fibre = 1 # [mm] distance fiber center to edge of beta arm


x_inc = 3.1 # [mm] Horizontal increment at each row
y_inc = 5.369 # [mm] Vertical increment at each row

""" Module parameters """ 

nb_robots = 102

if nb_robots == 75:

     module_width = 80 # [mm] triangle side length
     nb_rows = 11 # number of rows of positioners

elif nb_robots == 63:
      
     module_width = 73.8# [mm] triangle side length
     nb_rows = 10 # number of rows of positioners

elif nb_robots == 102:
      
     module_width = 92.4 # [mm] triangle side length
     nb_rows = 13 # number of rows of positioners

else:
      raise Exception('Invalid number of robots or number of robots not supported')


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
intermediate_frame_thick = 0 # [mm] spacing between modules inside intermediate frame

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

def rotate_and_translate(geom, angle, dx, dy, origin = 'centroid'):

     rotated_geom = affinity.rotate(geom, angle, origin=origin)
     transformed_geom = affinity.translate(rotated_geom, dx, dy)

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
     chanfer_triangle = Polygon(to_polygon_format(x_vertices-0.0001,y_vertices-0.0001))
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

     plt.savefig(results_dir + today_filename)

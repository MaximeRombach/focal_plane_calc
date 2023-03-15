#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from shapely import affinity
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


""" Module parameters """ 

module_width = 80 # [mm] triangle side length

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
intermediate_frame_thick = 3 # [mm] spacing between modules inside intermediate frame

def remove_positioner(xx,yy, list_to_remove):
     
     xx_sliced = np.delete(xx,list_to_remove)
     yy_sliced = np.delete(yy,list_to_remove)

     return xx_sliced, yy_sliced

def to_polygon_format(x,y):

    coords = []
    for (i,j) in zip(x,y):
            coords.append((i,j))

    return coords

def rotate_and_translate(geom, angle, dx, dy, origin = 'centroid'):

     rotated_geom = affinity.rotate(geom, angle, origin=origin)
     transformed_geom = affinity.translate(rotated_geom, dx, dy)

     return transformed_geom

def create_equilateral_module_base(module_width, offset = [0,0]):
      
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


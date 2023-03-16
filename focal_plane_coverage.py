#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import random
import array as array
import math
import json

from shapely.geometry import Polygon, Point
from shapely.ops import unary_union
from shapely.plotting import plot_polygon, plot_points
from shapely import affinity, MultiPoint, MultiPolygon, GeometryCollection, BufferCapStyle

import parameters as param

"""
This code is a tool for calculating the coverage performed for different focal plane arrangements.
It is split in 2 main steps:
1) Define the coverage of the positioners in 1 module unit
     a) Make the triangular meshgrid of the center position of each robot within a raft
     b) Derive their actual coverage 
2) Pattern module units together and calculate the total effective coverage

The effective coverage shall be calculated as the usable positioner area vs the actual reachable area of the positioners
"""
## Positioners numbering:
            # 0: bottom left
            # count a row then switch to left side of next one to continue
            # last positoner: top one

# with open(".vscode\launch.json","r") as file:
#      data = json.load(file)
# print(data)
"""1)a) Define triangular meshgrid of positioners' centeral """

module = param.chanfered_base(param.module_width, chanfer_length=7.5)
coords_module_x, coords_module_y = module.exterior.coords.xy
print(coords_module_x[0])

nb_rows = 11
x_inc = 3.1 # [mm] Horizontal increment at each row
y_inc = 5.369 # [mm] Vertical increment at each row
# Their norm corresponds to the pitch defined in parameters
xx = []
xx1 = np.ones(nb_rows + 1)
yy = []
yy1 = np.ones(nb_rows + 1)

for idx in range(nb_rows):
    
    start_offset_x = x_inc * idx 
    start_offset_y = y_inc * idx 

    xx_new = np.linspace(start_offset_x, (nb_rows - idx ) * param.pitch + start_offset_x, nb_rows - idx + 1)
    xx_new += 5.9
    yy_new = start_offset_y * np.ones(len(xx_new))
    yy_new += 3.41
    xx.append(xx_new.tolist())
    yy.append(yy_new.tolist())

    if idx == 0:
         xx1 = xx_new
         yy1 = yy_new
    else: 
         xx1 = np.hstack((xx1, xx_new))
         yy1 = np.hstack((yy1, yy_new))


    if idx == nb_rows - 1:
        xx_new = np.array([x_inc * (idx + 1)])
        xx_new += 5.9
        yy_new = y_inc * (idx + 1) * np.ones(len(xx_new))
        yy_new += 3.41
        xx.append(list(xx_new))
        yy.append(list(yy_new))
        xx1 = np.hstack((xx1, xx_new))
        yy1 = np.hstack((yy1, yy_new))

# list_to_remove = [0, 1,2,3,4,5,6,7,8,9,10, nb_rows,12, 22, -1]
list_to_remove = [0, nb_rows, -1]
# list_to_remove = []
xx1, yy1 = param.remove_positioner(xx1, yy1, list_to_remove)
nb_robots = len(xx1)

triang_meshgrid = MultiPoint(param.to_polygon_format(xx1, yy1)) # convert the meshgrid into shapely standard for later manipulation


"""1)b) Define coverage for 1 modules"""

c1 = np.cos(np.deg2rad(param.alpha))
s1 = np.sin(np.deg2rad(param.alpha))

c2 = np.cos(np.deg2rad(param.beta))
s2 = np.sin(np.deg2rad(param.beta))

xa, ya, xab, yab = (param.lb-param.la)*c1, (param.lb-param.la)*s1, (param.lb+param.la)*c2, (param.la+param.lb)*s2

wks_list = []

for idx, (dx, dy) in enumerate(zip(xx1, yy1)):

        xa1, ya1, xab1, yab1 = xa + dx, ya + dy, xab + dx, yab + dy

        coords_int = param.to_polygon_format(xa1, ya1)
        coords_ext = param.to_polygon_format(xab1, yab1)

        interior = coords_int[::-1]
        poly_c1 = Polygon(coords_ext, [interior])
        wks_list.append(poly_c1)

multi_wks = MultiPolygon(wks_list)
total_positioners_workspace = unary_union(wks_list)
module_w_beta_and_safety_dist = module.buffer(-param.offset_from_module_edges)

if param.is_wall:
     effective_wks = module_w_beta_and_safety_dist.intersection(total_positioners_workspace)
else:
     effective_wks = total_positioners_workspace

module_collection = GeometryCollection([module, module_w_beta_and_safety_dist, effective_wks, triang_meshgrid])
# print(module_collection.geoms)
coverage_with_walls = round(effective_wks.area/module.area * 100,1)
coverage_no_walls = round(total_positioners_workspace.area/module.area * 100,1)

"""2)a) Start meshing the grid for intermediate frame (4 triangles: 3 upwards + 1 downwards)"""

angles = np.array([-30, 90, 210])
flip = [True,False,False,False]
dist = 2*param.module_width*np.sqrt(3)/6 + param.intermediate_frame_thick

x_inter = np.cos(np.deg2rad(angles))*dist
# x_inter = np.insert(x_inter, 0, np.cos(np.deg2rad(30))*dist)
x_inter = np.insert(x_inter, 0, 0)
y_inter = np.sin(np.deg2rad(angles))*dist
# y_inter = np.insert(x_inter, 0, np.sin(np.deg2rad(30))*dist)
y_inter = np.insert(y_inter, 0, 0)

listing = []
intermediate_collection = []
covered_area = 0

for idx, (rotate, dx, dy) in enumerate(zip(flip, x_inter, y_inter)):

     if rotate:
          angle = 180
     else:
          angle = 0

     transformed_all = param.rotate_and_translate(module_collection, angle, dx, dy, origin = "centroid")
     listing.append(transformed_all.geoms[0])
     intermediate_collection.append(transformed_all)
     covered_area += transformed_all.geoms[2].area

bounding_polygon_intermediate_frame = unary_union(MultiPolygon(listing)).convex_hull
# to_origin = bounding_polygon_intermediate_frame.bounds[0:2]
# print(to_origin)
intermediate_collection = GeometryCollection(intermediate_collection)

# intermediate_collection = param.rotate_and_translate(intermediate_collection, 0, -to_origin[0], -to_origin[1])
# bounding_polygon_intermediate_frame = param.rotate_and_translate(bounding_polygon_intermediate_frame, 0, -to_origin[0], -to_origin[1])

available_intermediate_area = round(bounding_polygon_intermediate_frame.area,1)
intermediate_coverage = round(covered_area/available_intermediate_area*100,1)

"""2)b) Start meshing the grid for the whole thing"""

centroid_intermediate_frame = bounding_polygon_intermediate_frame.centroid


""" Plot plot time """ 

draw = True
is_timer = False
plot_time = 20 # [s] plotting time

fig = plt.figure(figsize=(8,8))
plt.title("Module coverage raw")
plot_polygon(module, facecolor='None', edgecolor='black', add_points=False)
plot_polygon(module_w_beta_and_safety_dist, facecolor='None', edgecolor='red', linestyle = '--', add_points=False
             , label = "Safety dist = {} mm".format(param.offset_from_module_edges))
plt.scatter(xx1,yy1, marker='.', color='k', label = "{} robots".format(nb_robots))
for idx, wks in enumerate(wks_list):
    if idx == 0:
         label_cov = "Coverage w/o walls: {} %".format(coverage_no_walls)
    else:
         label_cov = None
    plot_polygon(wks, add_points=False, alpha=0.2, facecolor='red', edgecolor='black', label = label_cov)
plt.xlabel('x position [mm]')
plt.ylabel('y position [mm]')
plt.legend(shadow = True)

plt.figure(figsize=(8,8))
plt.title("Module coverage with summed coverage & walls")
plot_polygon(module, facecolor='None', edgecolor='black', add_points=False)
plot_polygon(module_w_beta_and_safety_dist, facecolor='None', linestyle = '--', add_points=False
             , label = "Safety dist = {} mm".format(param.offset_from_module_edges))
plt.scatter(xx1,yy1, marker='.', color='k', label = "{} robots".format(nb_robots))
plot_polygon(effective_wks, add_points=False, alpha=0.2, edgecolor='black', label = "Coverage with walls: {} %".format(coverage_with_walls))

plt.xlabel('x position [mm]')
plt.ylabel('y position [mm]')
plt.legend(shadow = True)

plt.figure(figsize=(8,8))
plt.title("Intermediate frame coverage")
for idx, mod_collection in enumerate(intermediate_collection.geoms):
     if idx == 0:
          label_coverage = 'Coverage: {} %'.format(intermediate_coverage)
          label_robots = "{} robots".format(nb_robots)
     else: 
          label_coverage = None
          label_robots = None
     for jdx, geo in enumerate(mod_collection.geoms):
          # plot_polygon(geometry[jdx], add_points=False, facecolor='None' , edgecolor='black')
          if (isinstance (mod_collection.geoms[jdx], Polygon)) and jdx == 0:
               plot_polygon(mod_collection.geoms[jdx], add_points=False, facecolor='None' , edgecolor='black')
          elif (isinstance (mod_collection.geoms[jdx], Polygon)) and jdx == 1:
               plot_polygon(mod_collection.geoms[jdx], add_points=False, facecolor='None', linestyle = '--')
          elif (isinstance (mod_collection.geoms[jdx], Polygon)) and jdx == 2:
               plot_polygon(mod_collection.geoms[jdx], add_points=False, alpha=0.2, edgecolor='black', label = label_coverage)
          elif (isinstance (mod_collection.geoms[jdx], MultiPoint)):
               plot_points(mod_collection.geoms[jdx], marker='.', color='k', label = label_robots)

plot_polygon(bounding_polygon_intermediate_frame, add_points=False, facecolor='None', linestyle = '-.', color = 'green', label = 'Available area: {} mm$^2$'.format(available_intermediate_area))

plt.xlabel('x position [mm]')
plt.ylabel('y position [mm]')
plt.legend(shadow = True)

if draw and is_timer:

     plt.show(block=False)
     plt.pause(plot_time)
     plt.close('all')

elif draw and not is_timer:
     
     plt.show()
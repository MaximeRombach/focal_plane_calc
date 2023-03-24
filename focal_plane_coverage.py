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

"""1)a) Define triangular meshgrid of positioners' centeral """

module = param.chanfered_base(param.module_width, chanfer_length=7.5)
coords_module_x, coords_module_y = module.exterior.coords.xy
reference_centroid = module.centroid


# Their norm corresponds to the pitch defined in parameters
xx = []
xx1 = np.ones(param.nb_rows + 1)
yy = []
yy1 = np.ones(param.nb_rows + 1)

for idx in range(param.nb_rows):
    
    start_offset_x = param.x_inc * idx 
    start_offset_y = param.y_inc * idx 

    xx_new = np.linspace(start_offset_x, (param.nb_rows - idx ) * param.pitch + start_offset_x, param.nb_rows - idx + 1)
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


    if idx == param.nb_rows - 1:
        xx_new = np.array([param.x_inc * (idx + 1)])
        xx_new += 5.9
        yy_new = param.y_inc * (idx + 1) * np.ones(len(xx_new))
        yy_new += 3.41
        xx.append(list(xx_new))
        yy.append(list(yy_new))
        xx1 = np.hstack((xx1, xx_new))
        yy1 = np.hstack((yy1, yy_new))

list_to_remove = [0, param.nb_rows, -1] # remove positioners at the edges of the triangle
# list_to_remove = [] # remove no positioner
xx1, yy1 = param.remove_positioner(xx1, yy1, list_to_remove)
print(xx1)
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

""" 2)a) Start meshing the grid for intermediate frame (4 triangles: 3 upwards + 1 downwards) """

n_modules = 3

dist = 2*param.module_width*np.sqrt(3)/6 + param.intermediate_frame_thick
x_grid_inter = np.zeros(n_modules)
x_grid_inter[0:2]= np.linspace(reference_centroid.x, 2*reference_centroid.x, 2)
x_grid_inter[2]= np.linspace(reference_centroid.x , reference_centroid.x , 1)
y_grid_inter = np.zeros(n_modules)
y_grid_inter[0:2:2]+= np.sin(np.deg2rad(30))*dist
y_grid_inter[2]+= np.sin(np.deg2rad(30))*dist + 2*param.module_width*np.sqrt(3)/6
flip = np.zeros(n_modules)
flip[0:2:2] += 1

# n_row = 2
# dist = 2*param.module_width*np.sqrt(3)/6 + param.intermediate_frame_thick
# x_grid_inter = np.zeros(n_row)
# y_grid_inter = np.zeros(n_row)
# flip = np.zeros(n_row)

# for idx in range(n_row):

#      offset_x = reference_centroid.x + dist * np.cos(np.deg2rad(30)) * idx
#      offset_y = reference_centroid.y + dist * (np.sin(np.deg2rad(30)) + 1) * idx

#      new_row_x = np.linspace(offset_x, (n_row - idx ) * dist * np.cos(np.deg2rad(30)) + offset_x, n_row - idx + 1)
#      new_row_y = offset_y * np.ones(len(new_row_x))
#      new_flip = np.zeros(len(new_row_x))
     

#      if idx == 0:
#          x_grid_inter = new_row_x
#          y_grid_inter = new_row_y
#          new_flip[::2] += 1
#          flip = new_flip
#      else: 
#          x_grid_inter = np.hstack((x_grid_inter, new_row_x))
#          y_grid_inter = np.hstack((y_grid_inter, new_row_y))
#          new_flip[1::2] += 1
#          flip = np.hstack((flip, new_flip))

#      if idx == param.nb_rows - 1:
#         new_row_x = np.array([offset_x * (idx + 1)])
#         new_row_y = offset_y * (idx + 1) * np.ones(len(new_row_x))
#         new_flip = np.zeros(1)

#         x_grid_inter = np.hstack((x_grid_inter, new_row_x))
#         y_grid_inter = np.hstack((y_grid_inter, new_row_y))
#         flip = np.hstack((flip, new_flip))


# angles = np.array([-30, 90, 210])
# flip = [True,False,False,False]

# x_grid_inter = np.cos(np.deg2rad(angles))*dist
# # x_inter = np.insert(x_inter, 0, np.cos(np.deg2rad(30))*dist)
# # x_grid_inter = np.insert(x_grid_inter, 0, 0)
# y_grid_inter = np.sin(np.deg2rad(angles))*dist
# # y_inter = np.insert(x_inter, 0, np.sin(np.deg2rad(30))*dist)
# # y_grid_inter = np.insert(y_grid_inter, 0, 0)

listing = [module_collection.geoms[0]]
intermediate_collection = [module_collection]
covered_area = module_collection.geoms[2].area

# listing = []
# intermediate_collection = []
# covered_area = 0

for idx, (rotate, dx, dy) in enumerate(zip(flip, x_grid_inter, y_grid_inter)):

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

# # intermediate_collection = param.rotate_and_translate(intermediate_collection, 0, -to_origin[0], -to_origin[1])
# # bounding_polygon_intermediate_frame = param.rotate_and_translate(bounding_polygon_intermediate_frame, 0, -to_origin[0], -to_origin[1])

available_intermediate_area = round(bounding_polygon_intermediate_frame.area,1)
intermediate_coverage = round(covered_area/available_intermediate_area*100,1)
# intermediate_coverage = round(covered_area/area_to_cover*100,1)

"""2)b) Start meshing the grid for the whole thing"""

# centroid_intermediate_frame = bounding_polygon_intermediate_frame.centroid
## Pizza of module
# pizza_angle = 60
# vigR_lim_x = param.vigR * np.cos(np.deg2rad(np.linspace(0,pizza_angle,20)))
# vigR_lim_x = np.insert(vigR_lim_x, 0, 0)
# print(vigR_lim_x)
# vigR_lim_y = param.vigR * np.sin(np.deg2rad(np.linspace(0,pizza_angle,20)))
# vigR_lim_y = np.insert(vigR_lim_y, 0, 0)
# pizza = Polygon(param.to_polygon_format(vigR_lim_x, vigR_lim_y))
# area_to_cover = pizza.area


""" Plot plot time """ 

draw = True
is_timer = False
save_plots = False
plot_time = 20 # [s] plotting time

fig = plt.figure(figsize=(8,8))
figtitle = "Module_coverage_raw_" + str(param.nb_robots)
plt.title(figtitle)
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
param.save_figures_to_dir(save_plots, figtitle)

plt.figure(figsize=(8,8))
figtitle = "Module coverage with summed coverage & walls_" + str(param.nb_robots)
plt.title(figtitle)
plot_polygon(module, facecolor='None', edgecolor='black', add_points=False)
plot_polygon(module_w_beta_and_safety_dist, facecolor='None', linestyle = '--', add_points=False
             , label = "Safety dist = {} mm".format(param.offset_from_module_edges))
plt.scatter(xx1,yy1, marker='.', color='k', label = "{} robots".format(nb_robots))
plot_polygon(effective_wks, add_points=False, alpha=0.2, edgecolor='black', label = "Coverage with walls: {} %".format(coverage_with_walls))

plt.xlabel('x position [mm]')
plt.ylabel('y position [mm]')
plt.legend(shadow = True)
param.save_figures_to_dir(save_plots, figtitle)

plt.figure(figsize=(10,10))
# plt.figure()
figtitle = "Intermediate frame coverage" + str(param.nb_robots)
plt.title(figtitle)
for idx, mod_collection in enumerate(intermediate_collection.geoms):
     if idx == 0:
          label_coverage = 'Coverage: {} %'.format(intermediate_coverage)
          # label_coverage = 'Coverage: {} %'.format(area_to_cover)
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
     # plt.scatter(x_grid_inter,y_grid_inter, marker='o', color='r')

plot_polygon(bounding_polygon_intermediate_frame, add_points=False, facecolor='None', linestyle = '-.', color = 'green', label = 'Available area: {} mm$^2$'.format(available_intermediate_area))
plt.xlabel('x position [mm]')
plt.ylabel('y position [mm]')
plt.legend(shadow = True)
param.save_figures_to_dir(save_plots, figtitle)

if draw and is_timer:
     
     plt.show(block=False)
     plt.pause(plot_time)

     plt.close('all')

elif draw and not is_timer:
     
     plt.show()
     
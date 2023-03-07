#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import random
import array as array

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

"""1)a)"""
## Define triangular meshgrid ##

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

list_to_remove = [0, nb_rows, -1]
xx1, yy1 = param.remove_positioner(xx1, yy1, list_to_remove)

triang_meshgrid = MultiPoint(param.to_polygon_format(xx1, yy1)) # convert the meshgrid into shapely standard for later manipulation

"""1)b)"""
### Define area for 1 modules ###

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

coords_module = param.to_polygon_format(param.module_vertices_x, param.module_vertices_y)
module = Polygon(coords_module)

multi_wks = MultiPolygon(wks_list)
total_positioners_workspace = unary_union(wks_list)
module_w_beta_and_safety_dist = module.buffer(-param.offset_from_module_edges)
workspace_with_walls = module_w_beta_and_safety_dist.intersection(total_positioners_workspace)

module_collection = GeometryCollection([module, module_w_beta_and_safety_dist, workspace_with_walls, triang_meshgrid])
print(module_collection.geoms)
coverage_with_walls = round(workspace_with_walls.area/module.area * 100,1)
coverage_no_walls = round(total_positioners_workspace.area/module.area * 100,1)

"""2)a) Start meshing the grid for intermediate frame (4 triangles: 3 upwards + 1 downwards)"""

intermediate_frame_thick = 0
angles = np.array([-30, 90, 210])
flip = [True,False,False,False]
dist = 2*param.module_width*np.sqrt(3)/6 + intermediate_frame_thick

x_inter = np.cos(np.deg2rad(angles))*dist
x_inter = np.insert(x_inter, 0,0)
y_inter = np.sin(np.deg2rad(angles))*dist
y_inter = np.insert(y_inter, 0,0)
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
intermediate_collection = GeometryCollection(intermediate_collection)
available_intermediate_area = bounding_polygon_intermediate_frame.area
intermediate_coverage = round(covered_area/available_intermediate_area*100,1)

"""2)b) Start meshing the grid for the whole thing"""

### Plot plot time ###

draw = True
is_timer = True
plot_time = 20 # [s] plotting time

# fig = plt.figure(figsize=(8,8))
# plt.title("Module coverage raw")
# plot_polygon(module, facecolor='None', edgecolor='black', add_points=False)
# plot_polygon(module_w_beta_and_safety_dist, facecolor='None', edgecolor='red', linestyle = '--', add_points=False
#              , label = "Safety dist = {} mm".format(param.offset_from_module_edges))
# plt.scatter(xx1,yy1, marker='.', color='k')
# for idx, wks in enumerate(wks_list):
#     if idx == 0:
#          label_cov = "Coverage w/o walls: {} %".format(coverage_no_walls)
#     else:
#          label_cov = None
#     plot_polygon(wks, add_points=False, alpha=0.2, facecolor='red', edgecolor='black', label = label_cov)
# plt.xlabel('x position [mm]')
# plt.ylabel('y position [mm]')
# plt.legend(shadow = True)

# plt.figure(figsize=(8,8))
# plt.title("Module coverage with summed coverage & walls")
# plot_polygon(module, facecolor='None', edgecolor='black', add_points=False)
# plot_polygon(module_w_beta_and_safety_dist, facecolor='None', linestyle = '--', add_points=False
#              , label = "Safety dist = {} mm".format(param.offset_from_module_edges))
# plt.scatter(xx1,yy1, marker='.', color='k')
# plot_polygon(workspace_with_walls, add_points=False, alpha=0.2, edgecolor='black', label = "Coverage with walls: {} %".format(coverage_with_walls))

# plt.xlabel('x position [mm]')
# plt.ylabel('y position [mm]')
# plt.legend(shadow = True)

plt.figure(figsize=(8,8))
plt.title("Intermediate frame coverage")
for idx, mod_collection in enumerate(intermediate_collection.geoms):
     if idx == 0:
          label = 'Coverage: {} %'.format(intermediate_coverage)
     else: 
          label = None
     for jdx, geo in enumerate(mod_collection.geoms):
          # plot_polygon(geometry[jdx], add_points=False, facecolor='None' , edgecolor='black')
          if (isinstance (mod_collection.geoms[jdx], Polygon)) and jdx == 0:
               plot_polygon(mod_collection.geoms[jdx], add_points=False, facecolor='None' , edgecolor='black')
          elif (isinstance (mod_collection.geoms[jdx], Polygon)) and jdx == 1:
               plot_polygon(mod_collection.geoms[jdx], add_points=False, facecolor='None', linestyle = '--')
          elif (isinstance (mod_collection.geoms[jdx], Polygon)) and jdx == 2:
               plot_polygon(mod_collection.geoms[jdx], add_points=False, alpha=0.2, edgecolor='black', label = label)
          elif (isinstance (mod_collection.geoms[jdx], MultiPoint)):
               plot_points(mod_collection.geoms[jdx], marker='.', color='k')

plot_polygon(bounding_polygon_intermediate_frame, add_points=False, facecolor='None', linestyle = '-.', color = 'green', label = 'Available area')

plt.xlabel('x position [mm]')
plt.ylabel('y position [mm]')
plt.legend(shadow = True)

if draw and is_timer:

     plt.show(block=False)
     plt.pause(plot_time)
     plt.close('all')

elif draw and not is_timer:
     
     plt.show()
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import random
import array as array

from shapely.geometry import Polygon, Point
from shapely.ops import unary_union
from shapely.plotting import plot_polygon, plot_points
from shapely import affinity, MultiPoint, MultiPolygon, GeometryCollection

import parameters as param

# ## Parameters ##

# la = 1.8 # alpha arm length [mm] /!\ lb > la /!\
# lb = 1.8 # beta arm length [mm]
# p = 6.2 # pitch [mm]
# alpha = np.linspace(-180,180,180)
# beta = np.linspace(-180,180,180)
# ## Raw triangle
# # module_vertices_x = np.array([0,80,40,0])
# # module_vertices_y = np.array([0,0,69.3,0])

# ## Edge cut module
# module_vertices_x = np.array([7.5, 72.5, 76.25, 43.75, 36.25, 3.75, 7.5]) # [mm]
# module_vertices_y = np.array([0, 0, 6.5, 62.8, 62.8, 6.5, 0]) # [mm]
# module_width = 80 # [mm] triangle side length

# beta2fibre = 1 # [mm] distance fiber center to edge of beta arm
# safety_distance = 0.5 # [mm] self explicit for collision avoidance
# offset_from_module_edges = safety_distance + beta2fibre
# frame_thickness = 3 # [mm] thickness of frame between modules
# start_offset_x = 6.2 # [mm]
# start_offset_y = 3.41 # [mm]
# draw = True
# is_timer = True
# plot_time = 20 # [s] plotting time

## Positioners numbering:
            # 0: bottom left
            # count a row then switch to left side of next one to continue
            # last positoner: top one

def remove_positioner(xx,yy, list_to_remove):
     
     xx_sliced = np.delete(xx,list_to_remove)
     yy_sliced = np.delete(yy,list_to_remove)

     return xx_sliced, yy_sliced


def to_polygon_format(x,y):

    coords = []
    for (i,j) in zip(x,y):
            coords.append((i,j))

    return coords

def rotate_and_translate(geom, angle, dx, dy):

     rotated_geom = affinity.rotate(geom, angle, origin='centroid')
     transformed_geom = affinity.translate(rotated_geom, dx, dy)

     return transformed_geom

coords_module = to_polygon_format(param.module_vertices_x, param.module_vertices_y)
module = Polygon(coords_module)
## Define triangular meshgrid ##

nb_rows = 11
x_inc = 3.1 # [mm] Horizontal increment at each row
y_inc = 5.369 # [mm] Vertical increment at each row
# Their norm corresponds to the pitch defined above
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
xx1, yy1 = remove_positioner(xx1, yy1, list_to_remove)

triang_meshgrid = MultiPoint(to_polygon_format(xx1, yy1))


### Define area for 1 positioner ###


c1 = np.cos(np.deg2rad(param.alpha))
s1 = np.sin(np.deg2rad(param.alpha))

c2 = np.cos(np.deg2rad(param.beta))
s2 = np.sin(np.deg2rad(param.beta))

xa, ya, xab, yab = (param.lb-param.la)*c1, (param.lb-param.la)*s1, (param.lb+param.la)*c2, (param.la+param.lb)*s2

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
module_w_beta_and_safety_dist = module.buffer(-param.offset_from_module_edges)
workspace_with_walls = module_w_beta_and_safety_dist.intersection(total_positioners_workspace)

absolutely_all = MultiPolygon([module, module_w_beta_and_safety_dist, workspace_with_walls])

angle = 180
dist = param.module_width*np.sqrt(3)/6
dx = np.cos(np.deg2rad(30))*dist
dy = np.sin(np.deg2rad(30))*dist

DX = 2*dx 
DY = 2*dy

# DX = 2*dx + 40/69.3*3.1
# DY = 2*dy - 3.1

transformed_mesh = rotate_and_translate(triang_meshgrid, angle, DX, DY)

transformed_module = rotate_and_translate(module, angle, DX, DY)

transformed_wks = rotate_and_translate(multi_wks, angle, DX, DY)

transformed_all = rotate_and_translate(absolutely_all, angle, DX, DY)


coverage_with_walls = round(workspace_with_walls.area/module.area * 100,1)
coverage_no_walls = round(total_positioners_workspace.area/module.area * 100,1)

### Plot plot time ###

draw = True
is_timer = True
plot_time = 20 # [s] plotting time

fig = plt.figure(figsize=(8,8))
plt.title("Module with raw positioner coverage without walls")

# plot_polygon(transformed_module, add_points=False, facecolor='None', edgecolor='black')
# plot_points(transformed_mesh, marker='.', color='k')
# plot_polygon(transformed_wks, add_points = False, alpha=0.2, facecolor='red', edgecolor='black')
# plot_points(module.centroid)
plot_polygon(module, facecolor='None', edgecolor='black', add_points=False)
plot_polygon(module_w_beta_and_safety_dist, facecolor='None', edgecolor='red', linestyle = '--', add_points=False
             , label = "Safety dist = {} mm".format(param.offset_from_module_edges))
plt.scatter(xx1,yy1, marker='.', color='k')
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
plt.title("Module with summed coverage & walls")
plot_polygon(module, facecolor='None', edgecolor='black', add_points=False)
plot_polygon(module_w_beta_and_safety_dist, facecolor='None', linestyle = '--', add_points=False
             , label = "Safety dist = {} mm".format(param.offset_from_module_edges))
plt.scatter(xx1,yy1, marker='.', color='k')
plot_polygon(workspace_with_walls, add_points=False, alpha=0.2, edgecolor='black', label = "Coverage with walls: {} %".format(coverage_with_walls))
plot_polygon(transformed_all.geoms[0], add_points=False, alpha=0.2, edgecolor='black')
plot_points(transformed_mesh, marker='.', color='k')
plt.xlabel('x position [mm]')
plt.ylabel('y position [mm]')
plt.legend(shadow = True)

if draw and is_timer:

     plt.show(block=False)
     plt.pause(plot_time)
     plt.close('all')

elif draw and not is_timer:
     
     plt.show()
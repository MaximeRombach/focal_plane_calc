#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import random
import array as array

from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
from shapely.plotting import plot_polygon, plot_points

## Parameters ##
la = 1.8 # alpha arm length [mm]
lb = 1.8 # beta arm length [mm]
p = 6.2 # pitch [mm]
alpha = np.linspace(-180,180,180)
beta = np.linspace(-180,180,180)
module_vertices_x = np.array([0,80,40,0])
module_vertices_y = np.array([0,0,69.3,0])
start_offset_x = 6.2
start_offset_y = 3.41

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


coords_module = to_polygon_format(module_vertices_x, module_vertices_y)
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

    xx_new = np.linspace(start_offset_x, (nb_rows - idx ) * p + start_offset_x, nb_rows - idx + 1)
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


### Define area for 1 positioner ###

draw = True
c1 = np.cos(np.deg2rad(alpha))
s1 = np.sin(np.deg2rad(alpha))

c2 = np.cos(np.deg2rad(beta))
s2 = np.sin(np.deg2rad(beta))

xa, ya, xab, yab = la*c1, la*s1, (lb+la)*c2, (la+lb)*s2

wks_list = []

for idx, (dx, dy) in enumerate(zip(xx1, yy1)):

        xa1, ya1, xab1, yab1 = xa + dx, ya + dy, xab + dx, yab + dy

        coords_int = to_polygon_format(xa1, ya1)
        coords_ext = to_polygon_format(xab1, yab1)

        interior = coords_int[::-1]
        poly_c1 = Polygon(coords_ext, [interior])
        wks_list.append(poly_c1)

total_positioners_area = unary_union(wks_list)
positioners_area_with_walls = module.intersection(total_positioners_area)

coverage_with_walls = round(positioners_area_with_walls.area/module.area * 100,1)
coverage_no_walls = round(total_positioners_area.area/module.area * 100,1)

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot()
plt.title("Module with raw positioner coverage")
plot_polygon(module, facecolor='None', edgecolor='black', add_points=False)
plt.scatter(xx1,yy1, marker='.', color='k')
for idx, wks in enumerate(wks_list):
    if idx == 0:
         label_cov = "Coverage: {} %".format(coverage_with_walls)
    else:
         label_cov = None
    plot_polygon(wks, add_points=False, alpha=0.2, facecolor='red', edgecolor='black', label = label_cov)
plt.xlabel('x position [mm]')
plt.ylabel('y position [mm]')
plt.legend(shadow = True)

# plt.figure(figsize=(8,8))
# plt.title("Module with summed coverage & walls")
# plot_polygon(module, facecolor='None', edgecolor='black', add_points=False)
# plt.scatter(xx1,yy1, marker='.', color='k')
# plot_polygon(coverage_walls, add_points=False, alpha=0.2, edgecolor='black', label = "Coverage: {} %".format(coverage_with_walls))
# plt.legend()


print(coverage_no_walls)

if draw:
    plt.show()
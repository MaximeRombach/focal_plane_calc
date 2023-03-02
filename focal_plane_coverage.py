#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import random
import array as array

from shapely.geometry import Polygon
from shapely.plotting import plot_polygon, plot_points

## Parameters ##
la = 1.8 # alpha arm length [mm]
lb = 1.8 # beta arm length [mm]
p = 6.2 # pitch [mm]
alpha = np.linspace(-180,180,200)
beta = np.linspace(-180,180,200)
module_vertices_x = np.array([0,80,40,0])
module_vertices_y = np.array([0,0,69.3,0])
start_offset_x = 6.2
start_offset_y = 3.41

## Define triangular meshgrid ##

nb_rows = 11
x_inc = 3.1 # [mm] Horizontal increment at each row
y_inc = 5.369 # [mm] Vertical increment at each row
# Their norm corresponds to the pitch defined above
xx = []
yy = []

for idx in range(nb_rows):
    
    start_offset_x = x_inc * idx 
    start_offset_y = y_inc * idx 

    xx_new = np.linspace(start_offset_x, (nb_rows - idx ) * p + start_offset_x, nb_rows - idx + 1)
    xx_new += 5.9
    yy_new = start_offset_y * np.ones(len(xx_new))
    yy_new += 3.41
    xx.append(list(xx_new))
    yy.append(list(yy_new))
    plt.scatter(xx_new, yy_new, marker='o', color='k')

    if idx == nb_rows - 1:
        xx_new = np.array([x_inc * (idx + 1)])
        xx_new += 5.9
        yy_new = y_inc * (idx + 1) * np.ones(len(xx_new))
        yy_new += 3.41
        xx.append(list(xx_new))
        yy.append(list(yy_new))
        plt.scatter(xx_new, yy_new, marker='o', color='k')

xx = np.array(xx, dtype='object')
yy = np.array(yy, dtype='object')




c1 = np.cos(np.deg2rad(alpha))
s1 = np.sin(np.deg2rad(alpha))

c2 = np.cos(np.deg2rad(beta))
s2 = np.sin(np.deg2rad(beta))


### Define area for 1 positioner ###

xa, ya, xb, yb = la*c1, la*s1, (lb+la)*c2, (la+lb)*s2
coords_a = []
coords_b = []
for (x,y) in zip(xa,ya):
    coords_a.append((x,y))
for (x,y) in zip(xb,yb):
    coords_b.append((x,y))

poly_a = Polygon(coords_a)
poly_b = Polygon(coords_b)
workspace = poly_b.difference(poly_a)
x_ext,y_ext = workspace.exterior.coords.xy
x_wks_int, y_wks_int = workspace.interiors[0].coords.xy

# coord = Polygon(workspace.exterior.coords.xy+1)
plot_points(workspace.centroid)
plot_polygon(workspace)
plt.plot(module_vertices_x,module_vertices_y)
plt.grid()
plt.show()

# plt.plot(module_vertices_x,module_vertices_y)
# plt.show()
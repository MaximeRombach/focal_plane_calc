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
alpha = np.linspace(-180,180,180)
beta = np.linspace(-180,180,180)
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


### Define area for 1 positioner ###

draw = True
c1 = np.cos(np.deg2rad(alpha))
s1 = np.sin(np.deg2rad(alpha))

c2 = np.cos(np.deg2rad(beta))
s2 = np.sin(np.deg2rad(beta))

xa, ya, xab, yab = la*c1, la*s1, (lb+la)*c2, (la+lb)*s2


for idx, (x_shift, y_shift) in enumerate(zip(xx,yy)):
    for (x_shifted, y_shifted)in zip(x_shift, y_shift):
        xa1, ya1, xab1, yab1 = xa + x_shifted, ya + y_shifted, xab + x_shifted, yab + y_shifted
        coords_int = []
        coords_ext = []
        for (x,y) in zip(xa1,ya1):
            coords_int.append((x,y))
        for (x,y) in zip(xab1,yab1):
            coords_ext.append((x,y))

        interior = coords_int[::-1]
        poly_c1 = Polygon(coords_ext, [interior])
        plot_polygon(poly_c1)

# interior = coords_a2[::-1]
# poly_c2 = Polygon(coords_ab2, [interior])
# plot_polygon(poly_c2)

plt.plot(module_vertices_x,module_vertices_y, color='orange')
plt.grid()
if draw:
    plt.show()
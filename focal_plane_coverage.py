# %%
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import random
import array as array
import updown_tri as tri
import time
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
start_time = time.time()
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
module_collection = affinity.translate(module_collection, xoff = -reference_centroid.x, yoff = -reference_centroid.y)
coverage_with_walls = round(effective_wks.area/module.area * 100,1)
coverage_no_walls = round(total_positioners_workspace.area/module.area * 100,1)

""" 2)a) Start meshing the grid for intermediate frame (4 triangles: 3 upwards + 1 downwards) """

dist_inter = 2*param.module_width*np.sqrt(3)/6 + param.intermediate_frame_thick # distance between each neighbor from center module
angles = np.array([-30, 90, 210]) # angular positions of neighbors
flip = [True,False,False,False] # 

x_grid_inter = np.cos(np.deg2rad(angles))*dist_inter
x_grid_inter = np.insert(x_grid_inter, 0, 0)
y_grid_inter = np.sin(np.deg2rad(angles))*dist_inter
y_grid_inter = np.insert(y_grid_inter, 0, 0)

boundaries = []
intermediate_collection = []
covered_area = 0

for idx, (rotate, dx, dy) in enumerate(zip(flip, x_grid_inter, y_grid_inter)):

     if rotate:
          angle = 180
     else:
          angle = 0

     transformed_all = param.rotate_and_translate(module_collection, angle, dx, dy, origin = "centroid")
     boundaries.append(transformed_all.geoms[0])
     intermediate_collection.append(transformed_all)
     covered_area += transformed_all.geoms[2].area # add the net covered area of each module


bounding_polygon_intermediate_frame = unary_union(MultiPolygon(boundaries)).convex_hull
intermediate_collection.append(bounding_polygon_intermediate_frame)
intermediate_collection = GeometryCollection(intermediate_collection)

available_intermediate_area = round(bounding_polygon_intermediate_frame.area,1)
intermediate_coverage = round(covered_area/available_intermediate_area*100,1)
# intermediate_coverage = round(covered_area/area_to_cover*100,1)

# %% Global grid
"""2)b) Start meshing the grid for the whole thing"""

inter_frame_width = 2*param.module_width + 2*param.intermediate_frame_thick*np.cos(np.deg2rad(30))+2*param.global_frame_thick
dist_global = 2*inter_frame_width*np.sqrt(3)/6 + param.global_frame_thick
inter_centroid = inter_frame_width*np.sqrt(3)/3*np.array([np.cos(np.deg2rad(30)), np.sin(np.deg2rad(30))])

nb_max_modules = round(2*param.vigR / (inter_frame_width + param.global_frame_thick), 0)
print(nb_max_modules)

center_coords = []
start = time.time()
a_max = 6
b_max = 6
c_max = 6

x_grid = []
y_grid = []
flip_global = []
vigR_tresh = 50
for a in np.arange(-a_max,a_max):
     for b in np.arange(-b_max,b_max):
          for c in np.arange(-c_max,c_max):
               valid = a + b + c
               if valid == 1 or valid == 2:
                    x,y = tri.tri_center(a,b,c,inter_frame_width)
                    if np.sqrt(x**2 + y**2) < param.vigR - vigR_tresh: # check if module falls in vignetting radius
                         center_coords.append((a,b,c))
                         x_grid.append(x)
                         y_grid.append(y)
                         if tri.points_up(a,b,c):
                              flip_global.append(0)
                         else: 
                              flip_global.append(1)
                              

pizza = param.make_vigR_polygon()

grid = MultiPoint(param.to_polygon_format(x_grid, y_grid))
plt.show()

global_boundaries = []
global_collection = []
covered_area_global = 0

for idx, (rotate, dx, dy) in enumerate(zip(flip_global, x_grid, y_grid)):

     if rotate:
          angle = 180
     else:
          angle = 0

     transformed_all = param.rotate_and_translate(intermediate_collection, angle, dx, dy, origin = "centroid")
     global_collection.append(transformed_all)
     global_boundaries.append(transformed_all.geoms[4])
     # plot_polygon(transformed_all.geoms[0].geoms[0], add_points = False, facecolor = 'None', edgecolor = 'black')
     # plot_polygon(transformed_all.geoms[4], add_points = False, facecolor = 'None', linestyle = '--')
     # covered_area += transformed_all.geoms[2].area # add the net covered area of each module

global_bounding_polygon = unary_union(MultiPolygon(global_boundaries)).convex_hull
plot_polygon(global_bounding_polygon, add_points = False, facecolor = 'None', edgecolor = 'orange', linestyle = '--')
plt.scatter(x_grid,y_grid,color='red')

total_modules = len(x_grid)*len(x_grid_inter)
total_robots = total_modules*param.nb_robots

print(f"Total # modules: {total_modules} \n", f"Total # robots: {total_robots}")

# %% 

""" Plot plot time """ 

draw = True
is_timer = False
save_plots = False
plot_time = 20 # [s] plotting time
ignore_points = False

fig = plt.figure(figsize=(8,8))
figtitle = f"Module coverage raw - {param.nb_robots} robots per module"
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
figtitle = f"Module coverage with summed coverage & walls - {param.nb_robots} robots per module"
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
figtitle = f"Intermediate frame coverage - {param.nb_robots} per module"
plt.title(figtitle)
param.plot_intermediate(intermediate_collection, intermediate_coverage, available_intermediate_area, draw_legend = True)
plt.xlabel('x position [mm]')
plt.ylabel('y position [mm]')
plt.legend(shadow = True)
param.save_figures_to_dir(save_plots, figtitle)

plt.figure(figsize=(10,10))
if param.intermediate_frame_thick != param.global_frame_thick:
     figtitle = f"Semi frameless - {param.nb_robots} per module"
elif param.intermediate_frame_thick == param.global_frame_thick and param.global_frame_thick == 0:
     figtitle = f"Frameless - {param.nb_robots} per module"
else:
     figtitle = f"Framed - {param.nb_robots} per module"
plt.title(figtitle)
plot_polygon(pizza, add_points = False, facecolor = 'None', edgecolor = 'black', linestyle = '--')
for idx, inter_collection in enumerate(global_collection):
     param.plot_intermediate(inter_collection, ignore_points = True)
param.save_figures_to_dir(save_plots, figtitle)

if draw and is_timer:
     
     plt.show(block=False)
     plt.pause(plot_time)

     plt.close('all')

elif draw and not is_timer:
     
     plt.show()
end_time = time.time()

print(f'Elapsed time: {end_time-start_time} s')

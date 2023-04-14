# %%
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import geopandas as gpd
import array as array
import updown_tri as tri
import time
from tqdm import tqdm
from shapely.geometry import Polygon, Point
from shapely.ops import unary_union
from shapely.plotting import plot_polygon, plot_points
from shapely import affinity, MultiPoint, MultiPolygon, GeometryCollection, BufferCapStyle

import parameters as param
R = param.curve_radius
#%%

"""
This code is a tool for calculating the coverage performed for different focal plane arrangements.
It is split in 3 main steps:
1) Define the coverage of the positioners in 1 module unit
     a) Make the triangular meshgrid of the center position of each robot within a raft
     b) Derive their actual coverage 
2) Pattern intermediate modul units (4 of them) together
3) Derive the total grid coverage

The effective coverage shall be calculated as the usable positioner area vs the actual reachable area of the positioners
"""
## Positioners numbering:
            # 0: bottom left
            # count a row then switch to left side of next one to continue
            # last positoner: top one
start_time = time.time()
"""1)a) Define triangular meshgrid of positioners' centeral """

module = param.chanfered_base(param.module_width, chanfer_length=7.5) # Create chamfered module shape as shapely Polygon
# coords_module_x, coords_module_y = module.exterior.coords.xy
reference_centroid = module.centroid


# Their norm corresponds to the pitch defined in parameters
xx1 = np.ones(param.nb_rows + 1)
yy1 = np.ones(param.nb_rows + 1)

for idx in range(param.nb_rows): # Create the grid of positioner's center points
    
     # Place the first robot of each row
    start_offset_x = param.x_inc * idx + param.x_first
    start_offset_y = param.y_inc * idx + param.y_first

     # Generate each row of robot: nb robots/row decreases as going upward in triangle
    xx_new = np.linspace(start_offset_x, (param.nb_rows - idx ) * param.pitch + start_offset_x, param.nb_rows - idx + 1)
    yy_new = start_offset_y * np.ones(len(xx_new))

    if idx == 0:
         xx1 = xx_new
         yy1 = yy_new
    else: 
         xx1 = np.hstack((xx1, xx_new))
         yy1 = np.hstack((yy1, yy_new))


    if idx == param.nb_rows - 1:
        xx_new = np.array([param.x_inc * (idx + 1)])
        yy_new = param.y_inc * (idx + 1) * np.ones(len(xx_new))
        xx1 = np.hstack((xx1, xx_new))
        yy1 = np.hstack((yy1, yy_new))

list_to_remove = [0, param.nb_rows, -1] # remove positioners at the edges of the triangle
# list_to_remove = [] # remove no positioner
xx1, yy1 = param.remove_positioner(xx1, yy1, list_to_remove)
nb_robots = len(xx1)

triang_meshgrid = MultiPoint(param.to_polygon_format(xx1, yy1)) # convert the meshgrid into shapely standard for later manipulation
plot_points(triang_meshgrid)
# %%

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
inter_coverage = []
covered_area = 0

for idx, (rotate, dx, dy) in enumerate(zip(flip, x_grid_inter, y_grid_inter)):

     if rotate:
          angle = 180
     else:
          angle = 0

     transformed_all = param.rotate_and_translate(module_collection, angle, dx, dy, origin = "centroid")
     boundaries.append(transformed_all.geoms[0])
     inter_coverage.append(transformed_all.geoms[2])
     intermediate_collection.append(transformed_all)

modules_polygon_intermediate = MultiPolygon(boundaries)
bounding_polygon_intermediate = modules_polygon_intermediate.convex_hull
coverage_polygon_intermediate = unary_union(MultiPolygon(inter_coverage)) # save coverage as whole to speedup global calculation
plot_polygon(coverage_polygon_intermediate, add_points=False)
covered_area_inter = coverage_polygon_intermediate.area
intermediate_collection.append(bounding_polygon_intermediate)
intermediate_collection = GeometryCollection(intermediate_collection)
intermediate_collection_speed = GeometryCollection([bounding_polygon_intermediate, modules_polygon_intermediate, coverage_polygon_intermediate])
available_intermediate_area = round(bounding_polygon_intermediate.area,1)
intermediate_coverage = round(covered_area_inter/available_intermediate_area*100,1)
# intermediate_coverage = round(covered_area/area_to_cover*100,1)
# param.plot_intermediate_speed(intermediate_collection_speed)

"""2)b) Start meshing the grid for the whole thing"""

inter_frame_width = 2*param.module_width + 2*param.intermediate_frame_thick*np.cos(np.deg2rad(30))+2*param.global_frame_thick
dist_global = 2*inter_frame_width*np.sqrt(3)/6 + param.global_frame_thick
inter_centroid = inter_frame_width*np.sqrt(3)/3*np.array([np.cos(np.deg2rad(30)), np.sin(np.deg2rad(30))])

nb_max_modules = round(2*param.vigR / (inter_frame_width + param.global_frame_thick), 0)
print(nb_max_modules)

center_coords = []
a_max = 5
b_max = 5
c_max = 5

x_grid = []
y_grid = []
flip_global = []

for a in np.arange(-a_max,a_max):
     for b in np.arange(-b_max,b_max):
          for c in np.arange(-c_max,c_max):
               valid = a + b + c
               if valid == 1 or valid == 2:
                    x,y = tri.tri_center(a,b,c,inter_frame_width)
                    if np.sqrt(x**2 + y**2) < param.vigR - param.vigR_tresh: # check if module falls in vignetting radius
                         center_coords.append((a,b,c))
                         x_grid.append(x)
                         y_grid.append(y)
                         if tri.points_up(a,b,c): # flip the modules depending on there position on the grid
                              flip_global.append(0)
                         else: 
                              flip_global.append(1)                             

pizza = param.make_vigR_polygon()
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# grid = MultiPoint(param.to_polygon_format(x_grid, y_grid))

longitude = np.array(x_grid)/R
latitude = 2 * np.arctan(np.exp(np.array(y_grid)/R))-np.pi/2

x_sph = R * np.cos(latitude) * np.cos(longitude)
y_sph = R * np.cos(latitude) * np.sin(longitude)
z_sph = R * np.sin(latitude)

# plt.scatter(y_sph, z_sph, color='red')
# plt.figure()
# plt.scatter(x_sph, z_sph, color='red', label=f'{max(x_sph)-min(x_sph)}mm')
# plt.legend()
# ax.plot_trisurf(x_sph, y_sph, z_sph)
# ax.set_xlabel('X Label')
# ax.set_ylabel('Y Label')
# ax.set_zlabel('Z Label')
#%%
global_boundaries = []
global_collection = []
global_coverage = []
covered_area_global = 0
boundaries_df = {'geometry':[]}
modules_df = {'geometry':[]}
coverage_df ={'geometry':[]}
# Create module arrangement from the global grid
for idx, (rotate, dx, dy) in enumerate(zip(flip_global, x_grid, y_grid)):

     if rotate:
          angle = 180
     else:
          angle = 0

     transformed_all = param.rotate_and_translate(intermediate_collection_speed, angle, dx, dy, origin = "centroid")
     global_collection.append(transformed_all)
     global_boundaries.append(transformed_all.geoms[0])
     boundaries_df['geometry'].append(transformed_all.geoms[0])
     modules_df['geometry'].append(transformed_all.geoms[1])
     coverage_df['geometry'].append(transformed_all.geoms[2])
     covered_area += transformed_all.geoms[0].area # add the net covered area of each module


global_bounding_polygon = unary_union(MultiPolygon(global_boundaries)).convex_hull
instrumented_area = global_bounding_polygon.area
global_coverage = round(covered_area/instrumented_area*100,1)
gdf_bound = gpd.GeoDataFrame(boundaries_df)
gdf_modules = gpd.GeoDataFrame(modules_df)

gdf_coverage = gpd.GeoDataFrame(coverage_df)
gdf_coverage['label'] = f'Coverage: {global_coverage} %'

f, ax = plt.subplots()
gdf_modules.plot(ax=ax,facecolor='None')
gdf_bound.plot(ax=ax,facecolor='None', edgecolor='green')
gdf_coverage.plot(column='label',ax=ax, alpha=0.2, legend=True, label=gdf_coverage['label'])
ax.get_legend_handles_labels()
ax.legend()
plot_polygon(pizza, add_points=False, edgecolor='black', facecolor='None', linestyle='--')
plot_polygon(global_bounding_polygon, add_points=False, edgecolor='orange', facecolor='None', linestyle='--')

plt.xlabel('x position [mm]')
plt.ylabel('y position [mm]')
plt.show()


total_modules = len(x_grid)*len(x_grid_inter)
total_robots = total_modules*param.nb_robots

print(f"Total # modules: {total_modules} \n", f"Total # robots: {total_robots}")

# %% 

""" Plot plot time """ 

draw = True
is_timer = True
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
param.plot_intermediate(intermediate_collection, False, intermediate_coverage, available_intermediate_area, draw_legend = True)
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

for idx, inter_collection in enumerate(tqdm(global_collection)):
     if idx == 0:
          label_coverage = "Coverage: {} %".format(global_coverage)
     else: 
          label_coverage = None
     param.plot_intermediate_speed(inter_collection, label_coverage)
end_bisous = time.time()

plot_polygon(global_bounding_polygon, add_points = False, facecolor = 'None', edgecolor ='orange', linestyle = '--')
plt.xlabel('x position [mm]')
plt.ylabel('y position [mm]')
plt.legend(shadow = True)
plt.legend()
param.save_figures_to_dir(save_plots, figtitle)


end_time = time.time()

print(f'Elapsed time: {end_time-start_time:0.3f} s')

if draw and is_timer:
     
     plt.show(block=False)
     plt.pause(plot_time)

     plt.close('all')

elif draw and not is_timer:
     
     plt.show()


# %% For further information contact Maxime Rombach (EPFL - LASTRO) at maxime.rombach@epfl.ch
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import geopandas as gpd
import array as array
import updown_tri as tri
import time
from tqdm import tqdm
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.plotting import plot_polygon, plot_points, plot_line
from shapely import affinity, MultiPoint, MultiPolygon, GeometryCollection, BufferCapStyle

import parameters as param
R = param.curve_radius
vigR = param.vigR

"""
This code is a tool for calculating the coverage performed for different focal plane arrangements.
It is split in 3 main steps:
1) Define the coverage of the positioners in 1 module unit
     a) Make the triangular meshgrid of the center position of each robot within a raft
     b) Derive their actual coverage 
2) Pattern intermediate modul units (4 of them) together
3) Derive the total grid coverage

Main focal plane parameters that can be tuned in parameters.py:

- nb_robots : 63, 75, 102 (per module), all automatically associated to an adapted module size
- intermediate_frame_thick [mm]: gap between modules inside intermediate triangle
- global_frame_thick [mm]: gap between intermediate triangles in global frame (if inter = global, equivalent to framed cas)

The effective coverage is calculated as the usable positioner area vs total area of vigR (for comparison purposes
)
"""
#%% 1)a) Define triangular meshgrid of positioners' centeral

## Positioners numbering:
            # 0: bottom left
            # count a row then switch to left side of next one to continue
            # last positoner: top one


start_time = time.time()
nbots = [63, 75, 102]

draw = True
is_timer = False
save_plots = False
plot_time = 20 # [s] plotting time
ignore_points = False
keys = ['n63', 'n75', 'n102']

# Intermediate things

intermediate_df = {}

# Global things

nb_robots = 75
mod_param = param.Module(nb_robots)
module_width = mod_param.module_width
nb_rows = mod_param.nb_rows
module_collection, wks_list, coverages = mod_param.create_module()
module, module_w_beta_and_safety_dist, effective_wks, triang_meshgrid = module_collection.geoms
coverage_with_walls, coverage_no_walls = coverages


# %% 2)a) Meshing the grid for intermediate frame (4 triangles: 3 upwards + 1 downwards)

dist_inter = 2*module_width*np.sqrt(3)/6 + param.intermediate_frame_thick # distance between each neighbor from center module
angles = np.array([-30, 90, 210]) # angular positions of neighbors
flip = [True,False,False,False]

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
better_bound = unary_union(modules_polygon_intermediate).boundary
plot_line(better_bound)
coverage_polygon_intermediate = unary_union(MultiPolygon(inter_coverage)) # save coverage as whole to speedup global calculation
# plot_polygon(coverage_polygon_intermediate, add_points=False)
covered_area_inter = coverage_polygon_intermediate.area
intermediate_collection.append(bounding_polygon_intermediate)
intermediate_collection = GeometryCollection(intermediate_collection)
intermediate_collection_speed = GeometryCollection([bounding_polygon_intermediate, modules_polygon_intermediate, coverage_polygon_intermediate])
available_intermediate_area = round(bounding_polygon_intermediate.area,1)
intermediate_coverage = round(covered_area_inter/available_intermediate_area*100,1)
# plot_polygon(modules_polygon_intermediate.boundary)
# intermediate_coverage = round(covered_area/area_to_cover*100,1)
# param.plot_intermediate_speed(intermediate_collection_speed)

# %% 2)b) Start meshing the grid for the whole thing

inter_frame_width = 2*module_width + 2*param.intermediate_frame_thick*np.cos(np.deg2rad(30))+2*param.global_frame_thick
rho = inter_frame_width * np.sqrt(3)/6
dist_global = 2*rho + param.global_frame_thick
inter_centroid = inter_frame_width*np.sqrt(3)/3*np.array([np.cos(np.deg2rad(30)), np.sin(np.deg2rad(30))])

nb_max_modules = round(2*vigR / (inter_frame_width + param.global_frame_thick), 0)

center_coords = []
n = 6
a_max = n
b_max = n
c_max = n

x_grid = []
y_grid = []
flip_global = []
vigR_tresh = 20
# Make global grid out of triangular grid method credited in updown_tri.py
for a in np.arange(-a_max,a_max):
     for b in np.arange(-b_max,b_max):
          for c in np.arange(-c_max,c_max):
               valid = a + b + c
               if valid == 1 or valid == 2:
                    x,y = tri.tri_center(a,b,c,inter_frame_width)
                    if np.sqrt(x**2 + y**2) < vigR + vigR_tresh: # allow centroid of inter modules to go out of vigR for further filling purpose
                         center_coords.append((a,b,c))
                         x_grid.append(x)
                         y_grid.append(y)
                         if tri.points_up(a,b,c): # flip the modules depending on there position on the grid
                              flip_global.append(0)
                         else: 
                              flip_global.append(1)                             

pizza = param.make_vigR_polygon()

grid = MultiPoint(param.to_polygon_format(x_grid, y_grid))
plot_polygon(pizza, add_points = False, facecolor = 'None', edgecolor = 'black')
plot_points(grid)
longitude = np.array(x_grid)/R
latitude = 2 * np.arctan(np.exp(np.array(y_grid)/R))-np.pi/2

x_sph = R * np.cos(latitude) * np.cos(longitude)
y_sph = R * np.cos(latitude) * np.sin(longitude)
z_sph = R * np.sin(latitude)

x_grid = y_sph
y_grid = z_sph
grid_sph = MultiPoint(param.to_polygon_format(x_grid, y_grid))
plot_points(grid_sph, color='red')
min_pizz = param.make_vigR_polygon(r = vigR-rho)
max_pizz = param.make_vigR_polygon(r = vigR+rho)
plot_polygon(min_pizz, add_points =False, edgecolor = 'red', linestyle='--', facecolor='None')
plot_polygon(max_pizz, add_points =False, edgecolor = 'blue', linestyle='--', facecolor='None')

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(x_sph, y_sph, z_sph , label=f'{max(x_sph)-min(x_sph)}mm')
ax.set_box_aspect((1, 1, 1))
plt.legend()
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
#%% 2)c) Place the intermediate triangles accordingly on the grid

fill_empty = True
covered_area_global = 0 
total_modules = 0
boundaries_df = {'geometry':[], 'color': []}
modules_df = {'geometry':[]}
coverage_df ={'geometry':[]}
# Create module arrangement from the global grid
for idx, (rotate, dx, dy) in enumerate(zip(flip_global, x_grid, y_grid)):

     if rotate:
          angle = 180
     else:
          angle = 0
     # Apply rotation around the inter module centroid and translation to its corresponding position on the grid
     transformed_all = param.rotate_and_translate(intermediate_collection_speed, angle, dx, dy, origin = "centroid")

     new_boundary = transformed_all.geoms[0]
     color_boundary = 'green'
     sticks_out = new_boundary.overlaps(pizza)
     if sticks_out and not fill_empty:
          color_boundary = 'red'
          # continue
     new_modules = transformed_all.geoms[1]
     new_coverage = transformed_all.geoms[2]
     
     # Check if intermediate module goes out from vigR
     if fill_empty and sticks_out:
     # if new_boundary.overlaps(pizza) and fill_empty:
          # Keep the indiv triangle that are within vigR     
               new_modules = MultiPolygon([coucou for coucou in new_modules.geoms if not coucou.overlaps(pizza) and norm(coucou.centroid.xy)<vigR-10])
               new_coverage = MultiPolygon([coucou for coucou in new_coverage.geoms if not coucou.overlaps(pizza) and norm(coucou.centroid.xy)<vigR-10])
               new_boundary = new_modules.convex_hull
               color_boundary = 'blue'
     
     boundaries_df['geometry'].append(new_boundary)
     boundaries_df['color'].append(color_boundary)
     modules_df['geometry'].append(new_modules)
     coverage_df['geometry'].append(new_coverage)
     # print(len(list(new_modules.geoms)))
     total_modules += len(list(new_modules.geoms))
     covered_area += new_coverage.area # add the net covered area of each module
global_bounding_polygon = MultiPolygon(boundaries_df['geometry']).convex_hull
instrumented_area = global_bounding_polygon.area
# global_coverage = round(covered_area/instrumented_area*100,1)
global_coverage = round(covered_area/pizza.area*100,1)
gdf_bound = gpd.GeoDataFrame(boundaries_df)
gdf_modules = gpd.GeoDataFrame(modules_df)
gdf_coverage = gpd.GeoDataFrame(coverage_df)
gdf_coverage['label'] = f'Coverage vigR: {global_coverage} %'

total_robots = total_modules*nb_robots

print(f"Total # modules: {total_modules} \n", f"Total # robots: {total_robots}")
# %% Plot plot time 

fig = plt.figure(figsize=(8,8))
figtitle = f"Module coverage raw - {nb_robots} robots per module"
plt.title(figtitle)
plot_polygon(module, facecolor='None', edgecolor='black', add_points=False)
plot_polygon(module_w_beta_and_safety_dist, facecolor='None', edgecolor='red', linestyle = '--', add_points=False
          , label = "Safety dist = {} mm".format(mod_param.offset_from_module_edges))
plot_points(triang_meshgrid, marker='.', color='k', label = "{} robots".format(nb_robots))

for idx, wks in enumerate(wks_list.geoms):
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
figtitle = f"Module coverage with summed coverage & walls \n {nb_robots} robots per module"
plt.title(figtitle)
plot_polygon(module, facecolor='None', edgecolor='black', add_points=False)
plot_polygon(module_w_beta_and_safety_dist, facecolor='None', linestyle = '--', add_points=False
          , label = "Safety dist = {} mm".format(mod_param.offset_from_module_edges))
plot_points(triang_meshgrid, marker='.', color='k', label = "{} robots".format(nb_robots))
plot_polygon(effective_wks, add_points=False, alpha=0.2, edgecolor='black', label = "Coverage with walls: {} %".format(coverage_with_walls))

plt.xlabel('x position [mm]')
plt.ylabel('y position [mm]')
plt.legend(shadow = True)
param.save_figures_to_dir(save_plots, figtitle)

plt.figure(figsize=(10,10))
figtitle = f"Intermediate frame - {nb_robots} robots per module \n Inner gap: {param.intermediate_frame_thick} mm \n Total # modules: 4 - Total # robots: {nb_robots*4}"
plt.title(figtitle)
param.plot_intermediate(intermediate_collection, nb_robots, False, intermediate_coverage, available_intermediate_area, draw_legend = True)
plt.xlabel('x position [mm]')
plt.ylabel('y position [mm]')
plt.legend(shadow = True)
param.save_figures_to_dir(save_plots, figtitle)

if param.global_frame_thick > 0:
     frame=pizza.difference(GeometryCollection(list(boundaries_df['geometry'])))
     plt.figure(figsize=(10,10))   
     figtitle = f'Frame to manufacture - {nb_robots} per module'
     plt.title(figtitle)
     plot_polygon(global_bounding_polygon, add_points = False, facecolor = 'None', edgecolor ='orange', linestyle = '--')
     plot_polygon(frame, add_points=False, facecolor='red', alpha = 0.2, edgecolor = 'black', label=f'Wall thickness = {param.global_frame_thick} mm')
     plt.xlabel('x position [mm]')
     plt.ylabel('y position [mm]')
     plt.legend(shadow = True)
     plt.legend()
     param.save_figures_to_dir(save_plots, figtitle)


if param.intermediate_frame_thick != param.global_frame_thick:
     figtitle = f"Semi frameless - {nb_robots} robots per module \n Inner gap: {param.intermediate_frame_thick} mm - Global gap: {param.global_frame_thick} mm \n Total # modules: {total_modules} - Total # robots: {total_robots}"
elif param.intermediate_frame_thick == param.global_frame_thick and param.global_frame_thick == 0:
     figtitle = f"Frameless - {nb_robots} robots per module"
else:
     figtitle = f"Framed - {nb_robots} robots per module \n Gap: {param.intermediate_frame_thick} mm \n Total # modules: {total_modules} - Total # robots: {total_robots}"
f, ax= plt.subplots(figsize=(10, 10), sharex = True, sharey=True)
f.suptitle(figtitle)
gdf_modules.plot(ax=ax,facecolor='None')
gdf_bound.plot(ax=ax,facecolor='None', edgecolor=gdf_bound['color'])
gdf_coverage.plot(column='label',ax=ax, alpha=0.2, legend=True, legend_kwds={'loc': 'upper right'})

plot_polygon(pizza, ax=ax, add_points=False, edgecolor='black', facecolor='None', linestyle='--')
plot_polygon(global_bounding_polygon, ax=ax, add_points=False, edgecolor='orange', facecolor='None', linestyle='--',label='Instrumented area')

param.save_figures_to_dir(save_plots, figtitle)

if param.intermediate_frame_thick != param.global_frame_thick:
     figtitle = f"Semi frameless - {nb_robots} robots per module \n Inner gap: {param.intermediate_frame_thick} mm - Global gap: {param.global_frame_thick} mm \n Total # modules: {total_modules} - Total # robots: {total_robots}"
elif param.intermediate_frame_thick == param.global_frame_thick and param.global_frame_thick == 0:
     figtitle = f"Frameless - {nb_robots} robots per module"
else:
     figtitle = f"Framed - {nb_robots} robots per module \n Gap: {param.intermediate_frame_thick} mm \n Total # modules: {total_modules} - Total # robots: {total_robots}"
f, axes= plt.subplots(nrows=2,ncols=2, figsize=(10, 10), sharex = True, sharey=True)
f.suptitle(figtitle)
[ax.grid() for ax in axes.flatten()]
gdf_modules.plot(ax=axes[0,0],facecolor='None')
gdf_bound.plot(ax=axes[0,0],facecolor='None', edgecolor=gdf_bound['color'])
gdf_coverage.plot(column='label',ax=axes[0,0], alpha=0.2, legend=True, legend_kwds={'loc': 'upper right'})

plot_polygon(pizza, ax=axes[0,0], add_points=False, edgecolor='black', facecolor='None', linestyle='--')
plot_polygon(global_bounding_polygon, ax=axes[0,0], add_points=False, edgecolor='orange', facecolor='None', linestyle='--',label='Instrumented area')

[ax.set_xlabel('x position [mm]') for ax in axes.flatten()]
[ax.set_ylabel('y position [mm]') for ax in axes.flatten()]
end_time = time.time()

print(f'Elapsed time: {end_time-start_time:0.3f} s')

if draw and is_timer:
     
     plt.show(block=False)
     plt.pause(plot_time)

     plt.close('all')

elif draw and not is_timer:
     
     plt.show()


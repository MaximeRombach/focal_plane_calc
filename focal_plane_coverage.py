# %% For further information contact Maxime Rombach (EPFL - LASTRO) at maxime.rombach@epfl.ch
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import array as array
import updown_tri as tri
import time
from tqdm import tqdm
from shapely.geometry import Polygon
from shapely.ops import unary_union, polygonize
from shapely.plotting import plot_polygon, plot_points, plot_line
from shapely import affinity, MultiPoint, MultiPolygon, GeometryCollection, BufferCapStyle
import logging
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning) 
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

The effective coverage is calculated as the usable positioner area vs total area of vigR (for comparison purposes)
"""
#%% 1)a) Define triangular meshgrid of positioners' centeral

## Positioners numbering:
            # 0: bottom left
            # count a row then switch to left side of next one to continue
            # last positoner: top one


start_time = time.time()

""" Global variables """
nbots = [63, 75, 88, 102]
# nbots = [63]
width_increase = 1 # [mm] How much we want to increase the base length of a module
chanfer_length = 8.2 # [mm] Size of chanfers of module vertices (base value: 7.5)
centered_on_triangle = False

""" Intermediate frame parameters """

intermediate_frame_thick =  1 # [mm] spacing between modules inside intermediate frame

""" Global frame parameters """

global_frame_thick = 3 # [mm] spacing between modules in global arrangement

""" Drawing parameters """
draw = True
is_timer = False
save_plots = False
plot_time = 20 # [s] plotting time
ignore_points = False

""" Data storage """
keys = []
global_dict = {'n52': {}, 'n63': {}, 'n75': {}, 'n88': {}, 'n102': {},
               'Overall_results': {'nb_robots_list':[], 'total_robots_list' : [],
                                    'total_modules_list':[], 'coverages_list':[]}}


for nb_robots in nbots:
     
     mod_param = param.Module(nb_robots, width_increase, chanfer_length)
     key = mod_param.key
     keys.append(key)
     module_width = mod_param.module_width
     nb_rows = mod_param.nb_rows
     module_collection, wks_list, coverages = mod_param.create_module()
     module, module_w_beta_and_safety_dist, effective_wks, triang_meshgrid = module_collection.geoms
     coverage_with_walls, coverage_no_walls = coverages

     # plot_polygon(module, add_points= False, facecolor='None', edgecolor='black')
     # plot_polygon(effective_wks, add_points= False, alpha = 0.2, edgecolor='black', label=f'Coverage = {coverage_with_walls} %')
     # plot_points(triang_meshgrid, color='black')
     # plt.legend()
     # plt.show()

     # %% 2)a) Meshing the grid for intermediate frame (4 triangles: 3 upwards + 1 downwards)

     inter_param = param.IntermediateTriangle(module_collection, mod_param.module_width, intermediate_frame_thick)
     intermediate_collection, intermediate_collection_speed, intermediate_coverage, inter_df = inter_param.create_intermediate_triangle()

     # %% 2)b) Start meshing the grid for the whole thing

     inter_frame_width = 2*module_width + 2*intermediate_frame_thick*np.cos(np.deg2rad(30))+2*global_frame_thick
     rho = inter_frame_width * np.sqrt(3)/6
     dist_global = 2*rho + global_frame_thick
     inter_centroid = inter_frame_width*np.sqrt(3)/3*np.array([np.cos(np.deg2rad(30)), np.sin(np.deg2rad(30))])

     nb_max_modules = round(2*vigR / (inter_frame_width + global_frame_thick), 0)
     
     n = 6
     center_coords = []
     # if centered_on_triangle:
     #      n+= 2  
     a_max = n
     b_max = n
     c_max = n

     x_grid = []
     y_grid = []
     flip_global = []
     if centered_on_triangle:
          origin = 'triangle'
          valid = [0,1]
     else:
          origin = 'vertex'
          valid = [1,2]
     vigR_tresh = 150
     # Make global grid out of triangular grid method credited in updown_tri.py
     # It's logic is also explained in updown_tri.py
     for a in np.arange(-a_max,a_max):
          for b in np.arange(-b_max,b_max):
               for c in np.arange(-c_max,c_max):
                    sum_abc = a + b + c
                    # if valid == 1 or valid == 2: 
                    if valid.count(sum_abc): # x,y "valid" if sum(a,b,c) == 1 or 2 for origin == 'vertice'
                         x,y = tri.tri_center(a,b,c,inter_frame_width) 
                         if np.sqrt(x**2 + y**2) < vigR + vigR_tresh: # allow centroid of inter modules to go out of vigR for further filling purpose
                              center_coords.append((a,b,c))
                              x_grid.append(x)
                              y_grid.append(y)
                              if tri.points_up(a,b,c, origin = origin): # flip the modules depending on there position on the grid
                                   flip_global.append(0)
                              else: 
                                   flip_global.append(1)                             

     pizza = param.make_vigR_polygon()

     # if centered_on_triangle:
     #      x_grid -= inter_centroid[0]
     #      y_grid -= inter_centroid[1]
     grid = MultiPoint(param.to_polygon_format(x_grid, y_grid))

     param.plot_vigR_poly(pizza)
     plot_points(grid, label=f"{nb_robots}")
     plt.legend()
     # plt.show()
     longitude = np.array(x_grid)/R
     latitude = 2 * np.arctan(np.exp(np.array(y_grid)/R))-np.pi/2

     x_sph = R * np.cos(latitude) * np.cos(longitude)
     y_sph = R * np.cos(latitude) * np.sin(longitude)
     z_sph = R * np.sin(latitude)

     x_grid = y_sph
     y_grid = z_sph
     grid_sph = MultiPoint(param.to_polygon_format(x_grid, y_grid))
     # plot_points(grid_sph, color='red')
     # min_pizz = param.make_vigR_polygon(r = vigR-rho)
     # max_pizz = param.make_vigR_polygon(r = vigR+rho)
     # plot_polygon(min_pizz, add_points =False, edgecolor = 'red', linestyle='--', facecolor='None')
     # plot_polygon(max_pizz, add_points =False, edgecolor = 'blue', linestyle='--', facecolor='None')

     # fig = plt.figure()
     # ax = fig.add_subplot(projection='3d')
     # ax.scatter(x_sph, y_sph, z_sph , label=f'{max(x_sph)-min(x_sph)}mm')
     # ax.set_box_aspect((1, 1, 1))
     # plt.legend()
     # ax.set_xlabel('X Label')
     # ax.set_ylabel('Y Label')
     # ax.set_zlabel('Z Label')
     #%% 2)c) Place the intermediate triangles accordingly on the grid

     fill_empty = True # Fill empty spaces by individual modules
     allow_small_out = True # allow covered area of module to stick out of vigR (i.e. useless covered area because does not receive light)
     out_allowance = 0.07 # percentage of the covered area of a module authorized to stick out of vigR
     covered_area = 0 
     total_modules = 0
     boundaries_df = {'geometry':[], 'color': []}
     modules_df = {'geometry':[]}
     coverage_df ={'geometry':[]}
     # Create module arrangement from the global grid
     logging.info(f'Arranging focal plane for {nb_robots} robot case')
     for idx, (rotate, dx, dy) in enumerate(zip(flip_global, x_grid, y_grid)):

          if rotate:
               angle = 180
          else:
               angle = 0
          # Apply rotation around the inter module centroid and translation to its corresponding position on the grid
          transformed_all = param.rotate_and_translate(intermediate_collection_speed, angle, dx, dy, origin = "centroid")

          new_boundary = transformed_all.geoms[3]
          new_modules = transformed_all.geoms[1]
          new_coverage = transformed_all.geoms[2]
          color_boundary = 'green'

          sticks_out = new_boundary.overlaps(pizza) # a portion of the int triangle sticks out
          int_centroid_out = np.sqrt(dx**2 + dy**2) > vigR
          if not fill_empty and sticks_out:
               color_boundary = 'red'
               continue
          
          # Check if intermediate module goes out from vigR and keep the modules that are either right inside and, if allowed,
          # the ones that only stick out of a certain amount
          elif fill_empty and sticks_out:

               temp_mod = []
               temp_cov = []
               xx = []
               yy = []
               
               for mod, cov in zip(new_modules.geoms, new_coverage.geoms):
                    # plot_polygon(pizza,add_points=False,facecolor='None')
                    # plot_polygon(mod)
                    cent = np.array(cov.centroid.xy)
                    cov_out = cov.difference(pizza)

                    centroid_out = np.sqrt(cent[0]**2 + cent[1]**2) > vigR
                    if not cov.overlaps(pizza) and not centroid_out:
                    # If the coverage area of a single module does not stick out of vigR and is within vigR, we keep it
                         temp_mod.append(mod)
                         temp_cov.append(cov)
                         color_boundary = 'blue'

                         xx.append(mod.exterior.coords.xy[0].tolist())
                         yy.append(mod.exterior.coords.xy[1].tolist())
                         

                    elif allow_small_out and cov_out.area/cov.area < out_allowance and not centroid_out:
                         # If the coverage area of a module sticks out by less than the authorized amount, we keep it
                         remaining_cov = cov.intersection(pizza)
                         temp_mod.append(mod)
                         temp_cov.append(remaining_cov)
                         color_boundary = 'blue'

               if not temp_mod: 
                    #check if empty list; skip case if True (corrects huge bug)
                    continue
               if len(temp_mod) == 4:
                    color_boundary = 'green'
                        
               new_modules = MultiPolygon(temp_mod)
               new_coverage = MultiPolygon(temp_cov)
               convex_hull_modules = new_modules.convex_hull
               temp_cent = np.array(convex_hull_modules.centroid.xy).reshape((1,2))
               xx = param.flatten_list(xx)
               yy = param.flatten_list(yy)
               if not xx or not yy:
                    new_boundary = convex_hull_modules
               else:
                    new_boundary = Polygon(param.sort_points_for_polygon_format(xx, yy, temp_cent))

          elif int_centroid_out and not sticks_out:
               # If int_triangle is not in vigR AND does not overlap with it, no point keeping it
               continue  

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

     global_dict[key]['boundaries_df'] = boundaries_df
     global_dict[key]['global_bounding_polygon'] = global_bounding_polygon
     global_dict[key]['modules_df'] = modules_df
     global_dict[key]['coverage_df'] = coverage_df
     global_dict[key]['coverage_df']['label'] = f'Coverage vigR: {global_coverage} %'
     global_dict[key]['instrumented_area'] = instrumented_area
     global_dict[key]['global_coverage'] = global_coverage
     global_dict[key]['nb_robots'] = nb_robots
     global_dict[key]['total_modules'] = total_modules
     global_dict[key]['total_robots'] = total_robots

     global_dict['Overall_results']['nb_robots_list'].append(nb_robots)
     global_dict['Overall_results']['total_robots_list'].append(total_robots)
     global_dict['Overall_results']['total_modules_list'].append(total_modules) 
     global_dict['Overall_results']['coverages_list'].append(global_coverage)
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
figtitle = f"Module coverage with summed coverage + walls \n {nb_robots} robots per module"
filename = f"Module cov w walls - {nb_robots} robots per mod"
plt.title(figtitle)
plot_polygon(module, facecolor='None', edgecolor='black', add_points=False)
plot_polygon(module_w_beta_and_safety_dist, facecolor='None', linestyle = '--', add_points=False
          , label = "Safety dist = {} mm".format(mod_param.offset_from_module_edges))
plot_points(triang_meshgrid, marker='.', color='k', label = "{} robots".format(nb_robots))
plot_polygon(effective_wks, add_points=False, alpha=0.2, edgecolor='black', label = "Coverage with walls: {} %".format(coverage_with_walls))

plt.xlabel('x position [mm]')
plt.ylabel('y position [mm]')
plt.legend(shadow = True)
param.save_figures_to_dir(save_plots, filename)

plt.figure(figsize=(10,10))
figtitle = f"Intermediate frame - {nb_robots} robots per module \n Inner gap: {intermediate_frame_thick} mm \n Total # modules: 4 - Total # robots: {nb_robots*4}"
filename = f"Intermediate plot - {nb_robots} robots per mod"
plt.title(figtitle)
param.plot_intermediate(intermediate_collection, nb_robots, False, intermediate_coverage, draw_legend = True)
gdf_inter_bound = gpd.GeoDataFrame(inter_df)
plt.xlabel('x position [mm]')
plt.ylabel('y position [mm]')
plt.legend(shadow = True)
param.save_figures_to_dir(save_plots, filename)

if global_frame_thick > 0:
     frame=pizza.buffer(20).difference(GeometryCollection(list(global_dict['n63']['boundaries_df']['geometry'])))
     plt.figure(figsize=(10,10))   
     figtitle = f'Frame to manufacture - {nb_robots} per module'
     filename = figtitle
     plt.title(figtitle)
     # plot_polygon(global_bounding_polygon, add_points = False, facecolor = 'None', edgecolor ='orange', linestyle = '--')
     plot_polygon(frame, add_points=False, facecolor='red', alpha = 0.2, edgecolor = 'black', label=f'Wall thickness = {global_frame_thick} mm')
     param.plot_vigR_poly(pizza, label = 'vigR')
     plt.xlabel('x position [mm]')
     plt.ylabel('y position [mm]')
     plt.legend(shadow = True)
     plt.legend()
     param.save_figures_to_dir(save_plots, filename)

figtitle = param.final_title(nb_robots, total_modules, total_robots, intermediate_frame_thick, global_frame_thick, allow_small_out, out_allowance)
filename = f"Coverage global - {nb_robots} rob - Inner {intermediate_frame_thick} mm - Global {global_frame_thick} mm"
f, ax= plt.subplots(figsize=(10, 10), sharex = True, sharey=True)
f.suptitle(figtitle)
gdf_modules.plot(ax=ax,facecolor='None')
gdf_bound.plot(ax=ax,facecolor='None', edgecolor=gdf_bound['color'])
gdf_coverage.plot(column='label',ax=ax, alpha=0.2, legend=True, legend_kwds={'loc': 'upper right'})
ax.scatter(0,0,s=7,color='red')
plot_polygon(pizza, ax=ax, add_points=False, edgecolor='black', facecolor='None', linestyle='--')
plot_polygon(global_bounding_polygon, ax=ax, add_points=False, edgecolor='orange', facecolor='None', linestyle='--',label='Instrumented area')
param.save_figures_to_dir(save_plots, filename)

figtitle = param.final_title(nb_robots, total_modules, total_robots, intermediate_frame_thick, global_frame_thick, allow_small_out, out_allowance, disp_robots_info=False)
filename = f"Summary of coverages for inner {intermediate_frame_thick} and global {global_frame_thick}"
f, axes= plt.subplots(nrows=2,ncols=2, figsize=(12, 12), sharex = True, sharey=True)
f.suptitle(figtitle)
axes = axes.flatten()
# ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
# ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
# ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
# ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
# ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)
# axes = [ax1, ax2, ax3, ax4, ax5]
for idx, (k,ax) in tqdm(enumerate(zip(keys, axes))):

     gdf_bound = gpd.GeoDataFrame(global_dict[k]['boundaries_df'])
     gdf_modules = gpd.GeoDataFrame(global_dict[k]['modules_df'])
     gdf_coverage = gpd.GeoDataFrame(global_dict[k]['coverage_df'])
     
     gdf_modules.plot(ax=ax,facecolor='None')
     gdf_bound.plot(ax=ax,facecolor='None', edgecolor=gdf_bound['color'])
     gdf_coverage.plot(column='label',ax=ax, alpha=0.2, legend=True, legend_kwds={'loc': 'upper right'})

     plot_polygon(pizza, ax=ax, add_points=False, edgecolor='black', facecolor='None', linestyle='--')
     plot_polygon(global_dict[k]['global_bounding_polygon'], ax=ax, add_points=False, edgecolor='orange', facecolor='None', linestyle='--',label='Instrumented area')

     ax.scatter(0,0,s=7,color='red')
     ax.set_title(f"{global_dict[k]['nb_robots']} robots / module \n # modules: {global_dict[k]['total_modules']} - # robots: {global_dict[k]['total_robots']}")
     ax.set_xlabel('x position [mm]')
     ax.set_ylabel('y position [mm]')
param.save_figures_to_dir(save_plots, filename)

fig,ax = plt.subplots(figsize = (8,8))
figtitle = f"Coverages and # robots for {keys} cases"
filename = figtitle
fig.suptitle(figtitle)
ax.plot(nbots, global_dict['Overall_results']['coverages_list'],
        color="red", marker="o")
ax.set_xlabel("# robots/module")
ax.set_ylabel("Coverage [% of vigR area]",
              color="red")

ax2=ax.twinx()
# make a plot with different y-axis using second axis object
ax2.plot(nbots, global_dict['Overall_results']['total_robots_list'],
         color="blue",marker="o")
ax2.set_ylabel("Total # robots",color="blue")
ax.grid()
param.save_figures_to_dir(save_plots, filename)

end_time = time.time()

print(f'Elapsed time: {end_time-start_time:0.3f} s')

if draw and is_timer:
     
     plt.show(block=False)
     plt.pause(plot_time)

     plt.close('all')

elif draw and not is_timer:
     
     plt.show()


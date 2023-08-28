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
from shapely.geometry import Polygon, Point
from shapely.plotting import plot_polygon, plot_points, plot_line
from shapely import MultiPoint, MultiPolygon, GeometryCollection
from shapely.ops import unary_union
import logging
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning) 
import parameters as param
from datetime import datetime

from astropy.coordinates import cartesian_to_spherical


logging.basicConfig(level=logging.INFO)

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

project_surface = 'MUST'
surf = param.FocalSurf(project = project_surface)
R = abs(surf.R)
vigR = surf.vigR
BFS = surf.BFS

""" Global variables """
nbots = [63, 75, 88, 102]
# nbots = [75]
out_allowances = np.arange(0, 0.9, 0.05)
# out_allowances = [0]
width_increase = 0 # [mm] How much we want to increase the base length of a module
chanfer_length = 10 # [mm] Size of chanfers of module vertices (base value: 7.5)
centered_on_triangle = False # move the center of the grid (red dot) on the centroid on a triangle instead of the edge
full_framed = False # flag to check wether we are in semi frameless or in full framed case (turns True if inter_gap = global_gap = 0 mm)

""" Intermediate frame parameters """

intermediate_frame_thick =  1 # [mm] spacing between modules inside intermediate frame

""" Global frame parameters """

global_frame_thick = 3 # [mm] spacing between modules in global arrangement

""" Protective shields on module """

is_wall = True # [bool] puts walls or not on modules (changes greatly coverage)

if intermediate_frame_thick == global_frame_thick and intermediate_frame_thick != 0:
     full_framed = True

""" Drawing parameters """
draw = True
is_timer = False

plot_time = 20 # [s] plotting time
ignore_robots_positions = False

save_plots = False # Save final plots drew with *draw* flag
save_all_plots = False
save_frame_as_dxf = False # Save the outline of the frame for Solidworks integration
save_csv = False # Save position of robots (flat for now, TBI: follow focal surface while staying flat in modules)
save_txt = False # Save positions of modules along focal surface
saving_df = {"save_plots": save_plots, "save_dxf": save_frame_as_dxf, "save_csv": save_csv, "save_txt": save_txt}
saving = param.SavingResults(saving_df)

"""GFA stuff"""
gfa_tune = 1
nb_gfa = 8
gfa = param.GFA(length = 33.3*gfa_tune, width = 61*gfa_tune, nb_gfa = nb_gfa, vigR=vigR, saving_df=saving_df)
gdf_gfa = gfa.gdf_gfa
polygon_gfa = MultiPolygon(list(gdf_gfa['geometry']))

logging.info(f'Loaded parameters: \n - Surface:  {project_surface} \n - Inter gap: {intermediate_frame_thick} mm & Global gap: {global_frame_thick} mm \n - {nbots} robots/module')
""" Data storage """
keys = []

global_dict = {'n42': {'local_total_robots_list' : [], 'local_total_modules_list':[], 'local_coverages_list':[], 'useless_robots_list': [], 'useful_robots_list': [], 'efficiency_list': [] },
               'n52': {'local_total_robots_list' : [], 'local_total_modules_list':[], 'local_coverages_list':[], 'useless_robots_list': [], 'useful_robots_list': [], 'efficiency_list': [] },
               'n63': {'local_total_robots_list' : [], 'local_total_modules_list':[], 'local_coverages_list':[], 'useless_robots_list': [], 'useful_robots_list': [], 'efficiency_list': [] },
               'n75': {'local_total_robots_list' : [], 'local_total_modules_list':[], 'local_coverages_list':[], 'useless_robots_list': [], 'useful_robots_list': [], 'efficiency_list': [] },
               'n88': {'local_total_robots_list' : [], 'local_total_modules_list':[], 'local_coverages_list':[], 'useless_robots_list': [], 'useful_robots_list': [], 'efficiency_list': [] },
               'n102': {'local_total_robots_list' : [], 'local_total_modules_list':[], 'local_coverages_list':[], 'useless_robots_list': [], 'useful_robots_list': [], 'efficiency_list': []},
               'Overall_results': {'nb_robots_list':[], 'total_robots_list' : [],
                                    'total_modules_list':[], 'coverages_list':[]}}
to_dxf_dict = {}


for nb_robots in nbots: # iterate over number of robots/module cases

     
     mod_param = param.Module(nb_robots, saving_df, is_wall, width_increase, chanfer_length)
     key = mod_param.key
     keys.append(key)
     global_dict[key]['nbots/module'] = nb_robots
     module_width = mod_param.module_width
     nb_rows = mod_param.nb_rows
     module_collection, wks_list, coverages = mod_param.module_collection, mod_param.multi_wks_list, mod_param.coverages
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

     # x_inter = inter_df['inter_modules']['centroids'][:,0]
     # print(np.reshape(x_inter, (len(x_inter), 1)))
     # y_inter = inter_df['inter_modules']['centroids'][:,1]
     # z_inter = -np.sqrt(R**2 - (vigR)**2)*np.ones(len(x_inter))
     # inter_grid = np.vstack((x_inter.T, y_inter.T, z_inter.T))
     # print(inter_grid)

     # proj_front = grid.project_grid_on_sphere()

     # %% 2)b) Start meshing the grid for the whole thing

     # inter_frame_width = 2*module_width + 2*intermediate_frame_thick*np.cos(np.deg2rad(30))+2*global_frame_thick
     inter_frame_width = 2*module_width + 2*intermediate_frame_thick*np.cos(np.deg2rad(30)) + 2*global_frame_thick*np.cos(np.deg2rad(30))
     rho = inter_frame_width * np.sqrt(3)/6
     dist_global = 2*rho
     inter_centroid = inter_frame_width*np.sqrt(3)/3*np.array([np.cos(np.deg2rad(30)), np.sin(np.deg2rad(30))])

     grid = param.Grid(module_width, intermediate_frame_thick, global_frame_thick, vigR, R)
     grid_df = grid.grid_df
     grid.plot_2D_grid()

     flip_global = grid_df['flip_global']

     pizza = surf.make_vigR_polygon()
     pizza_with_GFA = pizza.difference(polygon_gfa)
     # plot_polygon(pizza_with_GFA)
     # plt.show()

     #%% 2)c) Place the intermediate triangles accordingly on the grid

     for outage in out_allowances: #iterate over how much we allow a module coverage to be out of the vignetting radius

          fill_empty = True # Fill empty spaces by individual modules
          allow_small_out = True # allow covered area of module to stick out of vigR (i.e. useless covered area because does not receive light)
          out_allowance = outage # percentage of the covered area of a module authorized to stick out of vigR
          covered_area = 0 # total covered area of vigR
          total_modules = 0 # total number of modules in vigR
          useless_robots = 0 # robots that are not in the vigR
          boundaries_df = {'geometry':[], 'color': []}
          modules_df = {'geometry':[], 'color': []}
          coverage_df = {'geometry':[]}
          robots_df = {'geometry':[]}
          final_grid = {'inter': {'x': [], 'y': [], 'z': [], 'xyz': [], 'geometry': []}, 
                    'indiv': {'x': [], 'y': [], 'z': [], 'xyz': [], 'geometry': [],
                              'robots': {'x': [], 'y': [], 'z': [], 'xyz': [], 'geometry': []}
                              }
                         }
          # Create module arrangement from the global grid
          logging.info(f'Arranging focal plane for {nb_robots} robot case')
          for idx, (rotate, dx, dy, dz) in enumerate(zip(flip_global, grid.x_grid, grid.y_grid, grid.z_grid)):

               if rotate:
                    angle = 180
               else:
                    angle = 0
               # Apply rotation around the inter module centroid and translation to its corresponding position on the grid
               transformed_all = param.rotate_and_translate(intermediate_collection_speed, angle, dx, dy, origin = "centroid")

               new_boundary = transformed_all.geoms[3]
               new_modules = transformed_all.geoms[1]
               new_coverage = transformed_all.geoms[2]
               new_robots = transformed_all.geoms[4]
               # print(new_robots.geoms)
               new_centroid = [dx, dy]
               out = []

               color_boundary = 'green'
               color_module = 'black'

               # Flags 
               sticks_out = new_boundary.overlaps(pizza) # a portion of the int triangle sticks out
               # int_centroid_out = np.sqrt(dx**2 + dy**2) > vigR # centroid of intermediate triangle out of vigR
               int_centroid_out = not Point(dx,dy).within(pizza)

               # Is on GFA?
               int_overlaps_GFA = new_boundary.overlaps(polygon_gfa)
               
               # closest_gfa_index = gfa.closest_gfa(new_centroid) # determine closest gfa to the current grid position
               # closest_gfa = gdf_gfa.at[closest_gfa_index, "geometry"] # extract closest gfa polygon for overlapping testing

               # int_tri_overlaps_gfa = new_boundary.overlaps(closest_gfa)
               # plot_polygon(new_modules)
               # param.plot_vigR_poly(pizza)
               # plt.show()

               if not fill_empty and sticks_out:
                    color_boundary = 'red'
                    continue

               elif int_centroid_out and not sticks_out:
                    # If int_triangle is not in vigR AND does not overlap with it, no point keeping it
                    continue
               
               # Check if intermediate module goes out from vigR and keep the modules that are either right inside and, if allowed,
               # the ones that only stick out of a certain amount
               elif fill_empty and sticks_out or int_overlaps_GFA:

                    temp_mod = []
                    temp_cov = []
                    temp_rob = []
                    xx = []
                    yy = []
                    
                    for mod, cov, robs in zip(new_modules.geoms, new_coverage.geoms, new_robots.geoms):
                         # plot_polygon(pizza,add_points=False,facecolor='None')
                         # plot_polygon(mod)
                         cent = np.array(cov.centroid.xy)
                         cov_out = cov.difference(pizza)

                         centroid_out = not Point(cent[0], cent[1]).within(pizza)
                         mod_overlaps_GFA = mod.overlaps(polygon_gfa)
                         if not cov.overlaps(pizza) and not centroid_out and not mod_overlaps_GFA:
                         # If the coverage area of a single module does not stick out of vigR and is within vigR, we keep it
                              temp_mod.append(mod)
                              temp_cov.append(cov)
                              temp_rob.append(robs)
                              color_boundary = 'blue'

                              # Log coordinates of boundaries of individual modules to create concave hull later
                              xx.append(mod.exterior.coords.xy[0].tolist())
                              yy.append(mod.exterior.coords.xy[1].tolist())
                              
                         elif allow_small_out and cov_out.area/cov.area < out_allowance and not mod_overlaps_GFA:
                              # If the coverage area of a module sticks out by less than the authorized amount, we keep it
                              plot_polygon(cov)
                              remaining_cov = cov.intersection(pizza)

                              # Sometimes coverage polygon is split in 2 polygons by pizza at the level of the arc de cercle given by the workspaces
                              # That gives a MultiPolygon which is better to split in separate polygons for further manipulations
                              if isinstance(remaining_cov, MultiPolygon):
                                   remaining_temp = []
                                   for poly in remaining_cov.geoms:
                                        remaining_temp.append(poly)
                                   remaining_cov = remaining_temp

                              remaining_robs = robs     
                              robs_in_vigR = robs.intersection(pizza)
                              useless_robots += len(list(robs.geoms)) - len(list(robs_in_vigR.geoms)) # counts the number of robots that are not in vigR
                              # pizz = surf.make_vigR_polygon()
                              # surf.plot_vigR_poly(pizz)
                              # plot_points(remaining_robs)
                              # plot_polygon(mod)
                              # plt.show()
                              temp_mod.append(mod)
                              temp_cov.append(remaining_cov)
                              temp_rob.append(remaining_robs)

                              # Log coordinates of boundaries of individual modules to create concave hull later
                              xx.append(mod.exterior.coords.xy[0].tolist())
                              yy.append(mod.exterior.coords.xy[1].tolist())
                              color_boundary = 'blue'

                    # Check if empty list; skip case if True (corrects huge bug)
                    if not temp_mod: 
                         
                         continue
                    if len(temp_mod) == 4:
                         color_boundary = 'green'
                         
                    new_modules = MultiPolygon(temp_mod)
                    new_coverage = MultiPolygon(temp_cov)
                    # new_robots = unary_union(temp_rob)
                    new_robots = GeometryCollection(temp_rob)
                    convex_hull_modules = new_modules.convex_hull

                    # Create concave hull in case of individual module
                    temp_cent = np.array(convex_hull_modules.centroid.xy).reshape((1,2))
                    xx = param.flatten_list(xx)
                    yy = param.flatten_list(yy)
                    if not xx or not yy:
                         new_boundary = new_modules
                    else:
                         # new_boundary = Polygon(param.sort_points_for_polygon_format(xx, yy, temp_cent))
                         new_boundary = new_modules

               if full_framed: # boundaries are the individual modules frontiers for the full framed case
                    new_boundary = new_modules

               # Store individual xy location for each module
               for new_mod, new_rob in zip(new_modules.geoms, new_robots.geoms) :
                    
                    final_grid['indiv']['x'].append(new_mod.centroid.x)
                    final_grid['indiv']['y'].append(new_mod.centroid.y)
                    final_grid['indiv']['z'].append(dz)
                    final_grid['indiv']['xyz'].append([new_mod.centroid.x, new_mod.centroid.y, dz])
                    final_grid['indiv']['geometry'].append(new_mod.centroid)
                    final_grid['indiv']['robots']['geometry'].append(new_rob)

               # Store individual xy location for each intermediate frame
               # final_grid['intermediate']['x'].append(dx)
               # final_grid['intermediate']['y'].append(dy)
               # final_grid['intermediate']['z'].append(dz)
               # final_grid['intermediate']['xyz'].append([dx, dy, dz])
               # final_grid['intermediate']['geometry'].append(Point(dx,dy))

               final_grid['inter']['x'].append(new_boundary.centroid.x)
               final_grid['inter']['y'].append(new_boundary.centroid.y)
               final_grid['inter']['z'].append(dz)
               final_grid['inter']['xyz'].append([dx, dy, dz])
               final_grid['inter']['geometry'].append(Point(new_boundary.centroid.x,new_boundary.centroid.y))

               boundaries_df['geometry'].append(new_boundary)
               boundaries_df['color'].append(color_boundary)
               modules_df['geometry'].append(new_modules)
               modules_df['color'].append(color_module)
               coverage_df['geometry'].append(new_coverage)

               robots_df['geometry'].append(new_robots)

               total_modules += len(list(new_modules.geoms))
               covered_area += new_coverage.area # add the net covered area of each module

          global_coverage = round(covered_area/pizza_with_GFA.area*100,1)
          gdf_bound = gpd.GeoDataFrame(boundaries_df)
          gdf_modules = gpd.GeoDataFrame(modules_df)
          gdf_coverage = gpd.GeoDataFrame(coverage_df)
          gdf_coverage['label'] = f'Coverage vigR: {global_coverage} %'
          robots_df['geometry'] = [unary_union(robots_df['geometry'])]
          gdf_robots = gpd.GeoDataFrame(robots_df)
          gdf_robots['markersize'] = 0.05
          total_robots = total_modules*nb_robots

          gdf_final_grid_int = gpd.GeoDataFrame(final_grid['inter'])

          global_dict[key]['boundaries_df'] = boundaries_df

          global_dict[key]['modules_df'] = modules_df
          global_dict[key]['coverage_df'] = coverage_df
          global_dict[key]['coverage_df']['label'] = f'Coverage vigR: {global_coverage} %'

          global_dict[key]['global_coverage'] = global_coverage
          global_dict[key]['nb_robots'] = nb_robots
          global_dict[key]['total_modules'] = total_modules
          global_dict[key]['total_robots'] = total_robots

          global_dict['Overall_results']['nb_robots_list'].append(nb_robots)
          global_dict['Overall_results']['total_robots_list'].append(total_robots)
          global_dict['Overall_results']['total_modules_list'].append(total_modules) 
          global_dict['Overall_results']['coverages_list'].append(global_coverage)

          global_dict[key]['local_coverages_list'].append(global_coverage)
          global_dict[key]['local_total_robots_list'].append(total_robots)
          global_dict[key]['local_total_modules_list'].append(total_modules)
          global_dict[key]['useless_robots_list'].append(useless_robots)
          global_dict[key]['useful_robots_list'].append(total_robots - useless_robots)
          global_dict[key]['efficiency_list'].append((total_robots - useless_robots)/total_robots)
          print(f"Out allowance: {outage} \n", f"Total # modules: {total_modules} \n", f"Total # robots: {total_robots} \n", f"Covergae: {global_coverage}")

#%% 2)d) Project grid on focal surface (BFS as a first try, aspherical will come later)

projection = {'front': {'x': [], 'y': [], 'z': [], 'xyz': [], 'theta': [], 'phi': []},
              'back': {'x': [], 'y': [], 'z': [], 'xyz': [], 'theta': [], 'phi': []}}

trim_angle = 60 # [deg], put 360° for full grid


if trim_angle == 360:
     grid_points = np.asarray(final_grid['indiv']['xyz'])
     grid_points_inter = np.asarray(final_grid['inter']['xyz'])
else:
     pizza_slice = surf.make_vigR_polygon(pizza_angle = trim_angle)
     grid_points = []
     grid_points_inter = []
     for idx, point in enumerate(final_grid['indiv']['geometry']):
          if point.within(pizza_slice):
               grid_points.append(final_grid['indiv']['xyz'][idx])
     grid_points = np.asarray(grid_points)
     for idx, point in enumerate(final_grid['inter']['geometry']):
          if point.within(pizza_slice):
               grid_points_inter.append(final_grid['inter']['xyz'][idx])
     grid_points_inter = np.asarray(grid_points_inter)

front_proj = grid.project_grid_on_sphere(grid_points, BFS, 'front_proj')
back_proj = grid.project_grid_on_sphere(grid_points, BFS + mod_param.module_length, 'back_proj')
proj = np.vstack((front_proj, back_proj))
r, lat, lon = cartesian_to_spherical(front_proj[:,0], front_proj[:,1], front_proj[:,2],)
projection['front']['x'] = front_proj[:,0]
projection['front']['y'] = front_proj[:,1]
projection['front']['z'] = front_proj[:,2]
projection['front']['xyz'] = front_proj
projection['front']['theta'] = np.rad2deg(lat.value)
projection['front']['phi'] = np.rad2deg(lon.value)
cols = ['x', 'y', 'z', 'theta', 'phi']
df = pd.DataFrame(projection['front'], columns = cols)

proj[:,2] = proj[:,2] + BFS
saving.save_grid_to_txt(proj, f'grid_indiv_{nb_robots}')
saving.save_grid_to_txt(front_proj, f'front_grid_indiv_{nb_robots}')
saving.save_grid_to_txt(back_proj, f'back_grid_indiv_{nb_robots}')
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
if save_all_plots:
     saving.save_figures_to_dir(figtitle)

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
if save_all_plots:
     saving.save_figures_to_dir(filename)

plt.figure(figsize=(10,10))
figtitle = f"Intermediate frame - {nb_robots} robots per module \n Inner gap: {intermediate_frame_thick} mm \n Total # modules: 4 - Total # robots: {nb_robots*4}"
filename = f"Intermediate plot - {nb_robots} robots per mod"
plt.title(figtitle)
param.plot_intermediate(intermediate_collection, nb_robots, False, intermediate_coverage, draw_legend = True)
gdf_inter_bound = gpd.GeoDataFrame(inter_df)
plt.xlabel('x position [mm]')
plt.ylabel('y position [mm]')
plt.legend(shadow = True)
if save_all_plots:
     saving.save_figures_to_dir(filename)

if global_frame_thick > 0:
     
     key_frame = keys[-1]
     # key_frame ='n75'
     robots = global_dict[key_frame]['nb_robots']
     modules = global_dict[key_frame]['total_modules']
     extra_material_for_frame = 50 # [mm] amount of material added on each side of vigR to make the frame structure
     figtitle = f'Frame to manufacture - {robots} robots per module - {modules} modules \n Vignetting diam: {2*vigR} mm - Extra: {extra_material_for_frame} mm\n Total diam: {2*(vigR + extra_material_for_frame)} mm'
     filename = f'Frame to manufacture - {robots} robots per module - {modules} modules'

     f, ax = plt.subplots(figsize=(10, 10))
     f.suptitle(figtitle)

     
     frame=pizza.buffer(extra_material_for_frame).difference(GeometryCollection(list(global_dict[key_frame]['boundaries_df']['geometry']))) 
     frame_ishish = unary_union(global_dict[key_frame]['boundaries_df']['geometry'])
     gdf_gfa.plot(ax=ax,facecolor = 'None', edgecolor=gdf_gfa['color'], linestyle='--', legend = True, label = 'GFA')
     # gdf_final_grid_int.plot(ax=ax)
     
     plot_polygon(frame, ax=ax, add_points=False, facecolor='red', alpha = 0.2, edgecolor = 'black', label=f'Wall thickness = {global_frame_thick} mm')
     surf.plot_vigR_poly(pizza, ax=ax, label = f'vigD = {2*vigR} mm')
     ax.set_xlabel('x position [mm]')
     ax.set_ylabel('y position [mm]')
     ax.legend(shadow = True)
     saving.save_figures_to_dir(filename)

     figtitle = f'Frame to manufacture - {robots} robots per module - {modules} modules'
     filename = figtitle

     f, ax = plt.subplots(figsize=(10, 10))     
     plot_polygon(frame, ax=ax, add_points=False, facecolor='red', alpha = 0.2, edgecolor = 'black')
     ax.spines['top'].set_visible(False)
     ax.spines['right'].set_visible(False)
     ax.spines['left'].set_visible(False)
     ax.spines['bottom'].set_visible(False)
     ax.axes.get_xaxis().set_visible(False)
     ax.axes.get_yaxis().set_visible(False)
     to_dxf_dict['frame'] = frame_ishish
     saving.save_dxf_to_dir(to_dxf_dict, f'frame_{robots}_robots_{modules}_modules')     
     

figtitle = param.final_title(nb_robots, total_modules, total_robots, intermediate_frame_thick, global_frame_thick, allow_small_out, out_allowance)
filename = f"Coverage global - {nb_robots} rob - Inner {intermediate_frame_thick} mm - Global {global_frame_thick} mm"
f, ax= plt.subplots(figsize=(10, 10), sharex = True, sharey=True)
f.suptitle(figtitle)
gdf_modules.plot(ax=ax,facecolor='None',edgecolor=gdf_modules['color'])
gdf_bound.plot(ax=ax,facecolor='None', edgecolor=gdf_bound['color'])
gdf_coverage.plot(column='label',ax=ax, alpha=0.2, legend=True, legend_kwds={'loc': 'upper right'})
ax.set_xlabel('x position [mm]')
ax.set_ylabel('y position [mm]')

if not ignore_robots_positions:
     gdf_robots.plot(ax = ax, markersize = 0.05)

if save_csv:
     indiv_pos_df = {'x [mm]': [], 'y [mm]' :[], 'geometry' : []}
     for idx, point in enumerate(robots_df['geometry'][0].geoms):
          indiv_pos_df['x [mm]'].append(point.x)
          indiv_pos_df['y [mm]'].append(point.y)
          indiv_pos_df['geometry'].append(point)
     now = datetime.now()
     info_case = f"{nb_robots}_robots-per-module_{total_robots}_robots_{intermediate_frame_thick}_inner_gap_{global_frame_thick}_global_gap"
     csv_filename = now.strftime("%Y-%m-%d-%H-%M-%S_") + info_case + ".csv"
     gpd.GeoDataFrame(indiv_pos_df).to_csv(saving.results_dir_path() + csv_filename, index_label = 'robot_number', sep = ";", decimal = ".", header='BONJOUR')
     # TBI: robot pos accounting for focal plane curvature
     logging.info(f'Robots positions saved to .csv file')

ax.scatter(0,0,s=7,color='red')
plot_polygon(pizza, ax=ax, add_points=False, edgecolor='black', facecolor='None', linestyle='--')

gdf_gfa.plot(ax=ax,facecolor = 'None', edgecolor=gdf_gfa['color'], linestyle='--')

saving.save_figures_to_dir(filename)

if len(nbots)>1: # Useless to do multiple plots for only one case
     figtitle = param.final_title(nb_robots, total_modules, total_robots, intermediate_frame_thick, global_frame_thick, allow_small_out, out_allowance, disp_robots_info=False)
     filename = f"Summary of coverages for inner {intermediate_frame_thick} and global {global_frame_thick}"
     f, axes= plt.subplots(nrows=2,ncols=2, figsize=(12, 12), sharex = True, sharey=True)
     f.suptitle(figtitle)
     axes = axes.flatten()

     for idx, (k,ax) in tqdm(enumerate(zip(keys, axes))):

          gdf_bound = gpd.GeoDataFrame(global_dict[k]['boundaries_df'])
          gdf_modules = gpd.GeoDataFrame(global_dict[k]['modules_df'])
          gdf_coverage = gpd.GeoDataFrame(global_dict[k]['coverage_df'])
          
          gdf_modules.plot(ax=ax,facecolor='None')
          gdf_bound.plot(ax=ax,facecolor='None', edgecolor=gdf_bound['color'])
          gdf_coverage.plot(column='label',ax=ax, alpha=0.2, legend=True, legend_kwds={'loc': 'upper right'})
          gdf_gfa.plot(ax=ax,facecolor = 'None', edgecolor=gdf_gfa['color'], linestyle='--')

          plot_polygon(pizza, ax=ax, add_points=False, edgecolor='black', facecolor='None', linestyle='--')
          # plot_polygon(global_dict[k]['global_bounding_polygon'], ax=ax, add_points=False, edgecolor='orange', facecolor='None', linestyle='--',label='Instrumented area')

          ax.scatter(0,0,s=7,color='red')
          ax.set_title(f"{global_dict[k]['nb_robots']} robots / module \n # modules: {global_dict[k]['total_modules']} - # robots: {global_dict[k]['total_robots']}")
          ax.set_xlabel('x position [mm]')
          ax.set_ylabel('y position [mm]')
     saving.save_figures_to_dir(filename)

     # fig,ax = plt.subplots(figsize = (8,8))
     # figtitle = f"Coverages and # robots for {keys} cases"
     # filename = figtitle
     # fig.suptitle(figtitle)
     # ax.plot(nbots, global_dict['Overall_results']['coverages_list'],
     #      color="red", marker="o")
     # ax.set_xlabel("# robots/module")
     # ax.set_ylabel("Coverage [% of vigR area]",
     #           color="red")

     # ax2=ax.twinx()
     # # make a plot with different y-axis using second axis object
     # ax2.plot(nbots, global_dict['Overall_results']['total_robots_list'],
     #      color="blue",marker="o")
     # ax2.set_ylabel("Total # robots",color="blue")
     # ax.grid()
     # saving.save_figures_to_dir(filename)


fig,ax = plt.subplots(figsize = (8,8))
gdf_robots_indiv = gpd.GeoDataFrame(geometry=final_grid['indiv']['robots']['geometry'])
gdf_robots_indiv.plot(ax=ax,markersize=0.1)

if len(out_allowances) > 1:
     fig = plt.figure(figsize=(8,8))
     filename = "Coverage_VS_out_allowance"
     plt.title(f'Focal plane coverage VS out allowance of modules \n Inner gap {intermediate_frame_thick} mm - Global gap {global_frame_thick} mm')
     for key in keys:
          plt.plot(np.array(out_allowances)*100, global_dict[key]['local_coverages_list'], '.-', label = f"{global_dict[key]['nbots/module']} robots/module")
     plt.xlabel('Out allowance [%]')
     plt.ylabel('Coverage [%]')
     plt.grid()
     plt.legend(shadow = True)
     saving.save_figures_to_dir(filename)

     fig = plt.figure(figsize=(8,8))
     filename = "Robots_VS_out_allowance"
     plt.title(f'Total robots VS out allowance of modules \n Inner gap {intermediate_frame_thick} mm - Global gap {global_frame_thick} mm')
     for key in keys:
          plt.plot(np.array(out_allowances)*100, global_dict[key]['local_total_robots_list'], '.-', label = f"{global_dict[key]['nbots/module']} robots/module")
     plt.xlabel('Out allowance [%]')
     plt.ylabel('Total robots [-]')
     plt.grid()
     plt.legend(shadow = True)
     saving.save_figures_to_dir(filename)

     fig = plt.figure(figsize=(8,8))
     filename = "Robots_VS_Coverage"
     plt.title(f'Focal plane coverage VS total robots \n Inner gap {intermediate_frame_thick} mm - Global gap {global_frame_thick} mm')
     for key in keys:
          plt.plot(global_dict[key]['local_total_robots_list'], global_dict[key]['local_coverages_list'], '.-', label = f"{global_dict[key]['nbots/module']} robots/module")
     plt.xlabel('Total robots [-]')
     plt.ylabel('Coverage [%]')
     plt.grid()
     plt.legend(shadow = True)
     saving.save_figures_to_dir(filename)

     fig = plt.figure(figsize=(8,8))
     filename = "Useless_VS_total_robots"
     plt.title(f'Useless robots VS total robots \n Inner gap {intermediate_frame_thick} mm - Global gap {global_frame_thick} mm')
     for key in keys:
          plt.plot(global_dict[key]['useless_robots_list'], global_dict[key]['local_total_robots_list'], '.-', label = f"{global_dict[key]['nbots/module']} robots/module")
     plt.xlabel('Useless robots [-]')
     plt.ylabel('Total robots [-]')
     plt.grid()
     plt.legend(shadow = True)
     saving.save_figures_to_dir(filename)

     fig = plt.figure(figsize=(8,8))
     filename = "Useful_robots_VS_total_robots"
     plt.title(f'Useful robots VS total robots \n Inner gap {intermediate_frame_thick} mm - Global gap {global_frame_thick} mm')
     for key in keys:
          plt.plot(global_dict[key]['efficiency_list'], global_dict[key]['local_coverages_list'], '.-', label = f"{global_dict[key]['nbots/module']} robots/module")
     plt.xlabel('Useful robots / Total robots [-]')
     plt.ylabel('Coverage [%]')
     plt.grid()
     plt.legend(shadow = True)
     saving.save_figures_to_dir(filename)
     

end_time = time.time()

print(f'Elapsed time: {end_time-start_time:0.3f} s')

if draw and is_timer:
     
     plt.show(block=False)
     plt.pause(plot_time)

     plt.close('all')

elif draw and not is_timer:
     
     plt.show()
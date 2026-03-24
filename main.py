#%%
# For further information contact Maxime Rombach (EPFL - LASTRO) at maxime.rombach@epfl.ch
# This code is a tool for calculating the focal plane layouts of future highly multiplexed MOS telescope instruments.
# It is inspired from the work of Joseph H. Silber (LBNL): https://github.com/joesilber/raft-design 
#!/usr/bin/env python
# -*- coding: utf-8 -*-
from Module import Module
from Focal_Suface import FocalSurf
from Grid import Grid
from GFAs import GFA
from SavingResults import SavingResults
import CustomLegends as cl
import time
from progress.bar import Bar

from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import geopandas as gpd
import pandas as pd
import numpy as np

import json

from shapely.plotting import plot_polygon
from shapely.ops import unary_union
from shapely.geometry import Point, MultiPolygon

import logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-4s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

timesstamp0 = time.time()

""" Available projects: MUST, Spec-S5, WST25, WST27, VLT_2030"""
PROJECT = "MUST"

""" Saving results """
save = SavingResults({"save_plots": True,
                      "save_txt": False,
                      "save_csv": False,
                      "save_dxf": False,
                      "project_name": PROJECT},
                      project_name = PROJECT)

""" Load focal surface data """
project_parameters = json.load(open('projects.json', 'r'))
trimming_angle = None # [deg] angle to trim the grid in phi direction
surf = FocalSurf(PROJECT, trimming_angle = trimming_angle, **project_parameters[PROJECT])
vigR = surf.vigD / 2
vignetting_area = surf.vignetting_disk.area
limit_pol = True
if limit_pol:
    limiting_polygon = surf.trimming_polygon(geometry='circle', trim_diff_to_vigR = 15)
else:
    limiting_polygon = surf.vignetting_disk

""" Module parameters"""

nb_robots_per_module = 63 # number of robots per module
pitch = 6.2 # [mm] distance between two adjacent robots
# HR_fibers = [10, 12, 15, 17, 29, 34, 39, 50, 52, 59] # position of the HR fibers in the module
# HR_fibers = [12, 15, 29, 34, 50, 52]
# HR_fibers = [21, 25, 39, 51]
HR_fibers = []
is_HR = len(HR_fibers) != 0

HR_l_beta = 3.6 # [mm] length of the HR fibers in beta direction
HR_l_alpha = 3.6 # [mm] length of the HR fibers in alpha direction
arms_length_tol = 0 # [mm] max error on the arms lengths of the robots, used to compute the worst case coverage of the module
is_wall = False

mod0 = Module(nb_robots = 63, 
                pitch = 6.2,
                HR_fibers = HR_fibers,
                HR_l_beta = HR_l_beta,
                HR_l_alpha = HR_l_alpha,
                is_wall = is_wall,
                arms_length_tol = arms_length_tol)
robots0 = mod0.robots_layout

INNER_GAP = 0.3 # [mm] gap between two modules within an intermediate triangle
GLOBAL_GAP = 4.4 # [mm] gap between two intermediate triangles; if inner = global, all modules are eqaully spaced
OUT_ALLOWANCE = 0 # fraction of the module coverage that is allowed to stick out of the vignetting disk

""" GFA parameters """

nb_gfa = 6
angle_offset = 0
gfa_length = 100
gfa_width = 100
#TODO: fix warning appearing twice --> FocalSurf called twice at begining of main AND within Grid class
#FIX: input surf as parameter to Grid class
gfa = GFA(nb_gfa = nb_gfa,
          angle_offset = angle_offset,
          vigR = surf.vigR,
          length = gfa_length,
          width = gfa_width)
gdf_gfa = gfa.make_GFA_array()
gfa_polygon = MultiPolygon(list(gdf_gfa['geometry'])) # GFA polygon

""" Grid parameters """

grid = Grid(focal_surface = surf,
            inner_gap = INNER_GAP,
            global_gap = GLOBAL_GAP,
            module_side_length = mod0.module_side_length,
            GFA_polygon = gfa_polygon,
            limiting_polygon = limiting_polygon,
            trimming_angle = trimming_angle,
            **project_parameters[PROJECT])

grid.flat_grid() # create the grid of modules
grid_3d = grid.grid_3d(grid.flat_grid_dict['x'], grid.flat_grid_dict['y'])
concave_hull = grid.layout_concave_hull() # concave hull of the module layout
concave_hull_area = concave_hull.area # area of the concave hull


modules = [] # list of modules
robots_workspaces = {'Module ID': [], 'Robot ID': [], 'x':[], 'y':[], 'z':[], 'geometry': []}
# robots_workspaces = {'geometry': []}
# robots_workspaces = {'Module ID': [], 'Robot ID': [], 'x':[], 'y':[]}

total_HR_fibers = 0
total_HR_area = 0
total_LR_fibers = 0
total_LR_area = 0
times = []
geom_times = []
assign_times = []
for_times = []
index2drop = []


with Bar('Aranging focal plane modules', max = len(grid_3d['x'])) as bar:
    # for mod_id,(x,y,z,points_up) in enumerate(zip(grid.flat_grid_dict['x'], grid.flat_grid_dict['y'], grid.flat_grid_dict['z'], grid.flat_grid_dict['tri_points_up'])):
    for mod_id,(x,y,z,points_up) in enumerate(zip(grid_3d['x'], grid_3d['y'], grid_3d['z'], grid_3d['tri_points_up'])):

        time1 = time.time()

        """ Create module object"""
        
        mod = Module(module_id = mod_id+1,
                        nb_robots = 63, 
                        pitch = 6.2,
                        module_points_up = points_up,
                        x0 = x,
                        y0 = y,
                        z0 = z,
                        HR_fibers = HR_fibers,
                        HR_l_beta = HR_l_beta,
                        HR_l_alpha = HR_l_alpha,
                        is_wall = is_wall,
                        arms_length_tol = arms_length_tol)

        robots = mod.robots_layout
        time2 = time.time()
        times.append(time2 - time1)
        """ Flags """
        geom_stamp = time.time()
        # if np.linalg.norm([x,y]) > surf.vigR-100:
        sticks_out = mod.module_boundaries.overlaps(limiting_polygon) # check if the module coverage sticks out of the vignetting disk
        area_sticking_out = mod.module_coverage.difference(limiting_polygon) # area of the module coverage that sticks out of the vignetting disk
        centroid_out = not Point(x,y).within(limiting_polygon)
        mod_overlaps_gfa = mod.module_boundaries.intersects(gfa_polygon) # check if the module coverage overlaps with the GFA polygon
        geom_times.append(time.time() - geom_stamp)

        assign_stamp = time.time()
        if (sticks_out and area_sticking_out.area/mod.module_coverage_area > OUT_ALLOWANCE) or mod_overlaps_gfa:                
            index2drop.append(mod_id)
            continue
        elif centroid_out and not sticks_out:
            """ If the centroid of the module coverage is outside the vignetting disk AND does NOT
            overlaps with vignetting --> no contribution to coverage, drop the module """
            index2drop.append(mod_id)
            continue
        elif mod_overlaps_gfa:
            """ GFAs are hard limits, no overlap allowed """
            index2drop.append(mod_id)
            continue
        else:
            if sticks_out:
                """ Update HR/LR coverages if part of the module coverage sticks out of the vignetting disk """
                LR_coverage_new = mod.LR_coverage.intersection(limiting_polygon)
                mod.LR_coverage = LR_coverage_new
                HR_coverage_new = mod.HR_coverage.intersection(limiting_polygon)
                mod.HR_coverage = HR_coverage_new
            modules.append(mod)

            total_HR_area += mod.HR_coverage.area
            total_LR_area += mod.LR_coverage.area

            total_HR_fibers += mod.nb_of_HR_fibers
            total_LR_fibers += mod.nb_of_LR_fibers


            robots_workspaces['geometry'].append(MultiPolygon(mod.dataframe['geometry']))
            # robots_workspaces['geometry'].append(mod.module_boundaries)
            # robots_workspaces['color'].append('yellow')
            robots_workspaces['Module ID'].extend(mod.dataframe["module_id"])
            robots_workspaces['Robot ID'].extend(mod.dataframe['robot_id'])
            robots_workspaces['x'].extend(mod.dataframe['x0'])
            robots_workspaces['y'].extend(mod.dataframe['y0'])
            robots_workspaces['z'].extend(mod.dataframe['z0'])

        assign_times.append(time.time() - assign_stamp)

        for_times.append(time.time() - time1)
        bar.next()

#%% 

# robots_workspaces = pd.DataFrame(robots_workspaces)
# # Renumber the 'Module ID' column in robots_workspaces to be consecutive starting from 1
# unique_ids = robots_workspaces['Module ID'].unique()
# id_map = {old_id: new_id for new_id, old_id in enumerate(sorted(unique_ids), start=1)}
# robots_workspaces['Module ID'] = robots_workspaces['Module ID'].map(id_map)
# # Renumber 'Robot ID' within each module from 1 to 63
# robots_workspaces['Robot ID'] = robots_workspaces.groupby('Module ID').cumcount() + 1

#%%
grid_3d.drop(index2drop, inplace = True) # drop the modules that do not contribute to the coverage
grid_3d = grid.trim_grid(grid_3d, trimming_angle = 360) # trim the grid to remove modules with phi < 0
grid_3d_back = grid.grid_3d_back(grid_3d)

save.save_grid_to_txt2(grid_3d, filename = f"Grid_{nb_robots_per_module}_rob__Inner_{INNER_GAP}_mm__Global_{GLOBAL_GAP}_mm", columns = ['x', 'y', 'z', 'tri_points_up'])
save.save_grid_to_txt2(grid_3d_back, filename = f"Grid_back_{nb_robots_per_module}_rob__Inner_{INNER_GAP}_mm__Global_{GLOBAL_GAP}_mm", columns = ['x', 'y', 'z', 'tri_points_up'])
save.save_grid_to_csv(grid_3d, filename = f"Grid_{nb_robots_per_module}_rob__Inner_{INNER_GAP}_mm__Global_{GLOBAL_GAP}_mm", results_string = f"Grid with {len(grid_3d)} modules, {len(modules)} contributing to the coverage")
save.save_grid_to_csv(robots_workspaces, filename = f"Robots_positions_{nb_robots_per_module}_rob__Inner_{INNER_GAP}_mm__Global_{GLOBAL_GAP}_mm")

LR_coverage_dict = {'geometry': [mod.LR_coverage for mod in modules]}
HR_coverage_dict = {'geometry': [mod.HR_coverage for mod in modules]}
boundaries = {'geometry': [mod.module_boundaries for mod in modules]}

HR_vignetting_coverage = 100 * total_HR_area/vignetting_area
HR_layout_coverage = 100 * total_HR_area/concave_hull_area

LR_vignetting_coverage = 100 * total_LR_area/vignetting_area
LR_layout_coverage = 100 * total_LR_area/concave_hull_area

#%%

mod0.plot_module(plot_rob_numbers = True)
filename = f"Module_{nb_robots_per_module}_rob_pitch_{mod.pitch}_mm"
save.save_figures_to_dir(filename)

figure, ax = plt.subplots(figsize=(10, 10))
if len(HR_fibers) != 0:
    geo_HR = gpd.GeoDataFrame(HR_coverage_dict)
    geo_HR.plot(ax = ax, color = 'red', alpha=0.4, legend=True)
geo_LR = gpd.GeoDataFrame(LR_coverage_dict)
geo_LR.plot(ax = ax, alpha=0.4, legend=True, color = 'C0')
geo_boundaries = gpd.GeoDataFrame(boundaries)
geo_boundaries.plot(ax = ax, facecolor = 'None', edgecolor = 'green')
gdf_gfa = gfa.make_GFA_array()
gdf_gfa.plot(ax = ax, facecolor = 'None', edgecolor = gdf_gfa['color'], linestyle = '--')
geo_fiducials = gpd.GeoDataFrame(grid.fiducials)
geo_fiducials.plot(ax = ax, facecolor = 'red', edgecolor = 'w', markersize = 10, label='Fiducials')
gdf_robots = gpd.GeoDataFrame({}, geometry = [Point(xy) for xy in zip(robots_workspaces['x'], robots_workspaces['y'])])
gdf_robots.plot(ax = ax, facecolor = 'blue', edgecolor = 'None', alpha =0.5, markersize = 1)
# for i, mod in enumerate(modules):
#     plt.text(mod.x0, mod.y0, f"{i+1}", fontsize=8, ha='center', va='center', color='white')
plot_polygon(grid.surf.vignetting_disk, ax = ax, fill = False, add_points=False, linestyle = '--', color = 'black')
plot_polygon(grid.layout_concave_hull(), ax = ax, fill = False, add_points=False, color = 'orange')

# plot_polygon(grid.fiducials_bounding_polygon, ax = ax, fill = False, add_points=False, color = 'purple')
# geo_grid = gpd.GeoDataFrame(grid_3d)
# geo_grid.plot(ax = ax, facecolor= 'None', markersize = 16, edgecolor = 'orange')
nb_fiducials = len(grid.fiducials['x'])
if len(HR_fibers) != 0:
    plt.legend(handles=[cl.HR_handle(extra_lab = f'HR vig : {HR_vignetting_coverage: .1f} %, layout : {HR_layout_coverage: .1f} %'),
                        cl.LR_handle(f'LR vig : {LR_vignetting_coverage: .1f} %, layout : {LR_layout_coverage: .1f} %'),
                        cl.GFA_handle(lab = f'GFAs: {nb_gfa}; {gfa_length}x{gfa_width} mm'), 
                        cl.fiducials_handle(lab = f'Fiducials: {nb_fiducials}')], loc='upper right')
else:
    plt.legend(handles=[cl.LR_handle(f'Vignetting coverage: {LR_vignetting_coverage: .1f} % \nLayout coverage: {LR_layout_coverage: .1f} %'),
                        cl.GFA_handle(lab = f'GFAs: {nb_gfa}'), 
                        cl.fiducials_handle(lab = f'Fiducials: {nb_fiducials}')], loc='upper right')
plt.title(cl.final_layout_title(PROJECT, surf.vigD, 
                                nb_robots_per_module, 
                                len(modules), 
                                nb_robots_per_module*len(modules), 
                                INNER_GAP, GLOBAL_GAP, 
                                OUT_ALLOWANCE,
                                total_HR_fibers, total_LR_fibers))
plt.xlabel('x [mm]')

if "WST" in PROJECT:
    def custom_formatter(x, pos):
        return f"{surf.mm2deg(x):.1f}"  # Converts y axis from mm to arcmin for WST case
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(custom_formatter))
    ax.set_ylabel('y [deg]')
    default_ticks = figure.gca().get_yticks()
    desired_ticks = [-vigR, -2/3*vigR, -vigR/3,0, vigR/3, vigR*2/3, vigR, -surf.arcmin2mm(surf.donut_diam/2), surf.arcmin2mm(surf.donut_diam/2)]
    ax.yaxis.set_ticks(desired_ticks)
else:
        ax.set_ylabel('y [mm]')

save.save_dxf_to_dir(geometries = robots_workspaces['geometry'], suffix_name = f"Robots_boundaries_{nb_robots_per_module}_rob__Inner_{INNER_GAP}_mm__Global_{GLOBAL_GAP}_mm")

plt.grid(visible=True)
print(f"Mean time for module: {sum(times)/len(times)} s")
print(f"Mean time for geom check: {sum(geom_times)/len(geom_times)} s")
# print(f"Mean time for assigning: {sum(assign_times)/len(assign_times)} s")
print(f"Mean time for 1 for loop: {sum(for_times)/len(for_times)} s")

filename = f"Coverage_global_{nb_robots_per_module}_rob__Inner_{INNER_GAP}_mm__Global_{GLOBAL_GAP}_mm"
save.save_figures_to_dir(filename, dpi = 800)

timesstamp1 = time.time()
print(f"Total time: {timesstamp1 - timesstamp0} s")
plt.show()
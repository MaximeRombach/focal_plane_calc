#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
# Or the path to the focal_plane_calc folder on your machine
# I know it is not the best way to do it but it works for now, To Be Fixed Later
sys.path.insert(0,'C:/Users/rombach/Documents/Astrobots/Inosuisse/focal_plane_calc')

import numpy as np
import matplotlib.pyplot as plt
plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
plt.rc('figure', titlesize=17)  # fontsize of the figure title
plt.rc('axes', titlesize=17)     # fontsize of the axes title
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsi
import math
from astropy.table import Table, vstack
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import parameters as param
import pandas as pd
from datetime import datetime

## Parameters ##

project = 'MUST'
surf = param.FocalSurf(project)

saving_df = {"save_plots": False, "save_txt": False}
save = param.SavingResults(saving_df, project_name=project)

## Functions ##

circles_diam = np.array([7.159, 14.318, 18.941, 25.813, 28.637, 31.206, 35.796, 37.883, 39.860, 43.547, 50.114, 51.625, 61.168]) # Diameter of the circles passing trough all pos centers [mm] --> from CAD model of module chassis

radial_pos = circles_diam / 2 # Radial pos used on BFS [mm]
radial_pos_diff = np.diff(radial_pos) # Radial pos differences [mm]
R2BFS = lambda r: (surf.BFS - np.sqrt(surf.BFS**2 - r**2))
r_plot = np.linspace(0,radial_pos[-1] * 1.1 ,500) # Define radius vector for focal plane curve
z_focus = R2BFS(radial_pos) # Calculate focal plane curve from csv data
z_diff = np.diff(z_focus) # Calculate differences between circle centers
z_plot = R2BFS(r_plot) # Calculate focal plane curve from csv data

results = {'Circles #': np.arange(1,len(circles_diam)+1),
            'Circles diameter [mm]': circles_diam,
            'Radial Position [mm]': radial_pos,
           r'Nominal focus position [$\mu$m]': z_focus * 10**3, 
           r"Focus Difference [$\mu$m]": np.concat((np.zeros(1), z_diff * 10**3))}
pd_results = pd.DataFrame(results)

save.save_grid_to_csv(pd_results, 'Module_focuses_table')

plt.figure(figsize=(10,6))
plt.plot(r_plot,z_plot * 10**3, label=f'BFS radius = {surf.BFS} mm')
plt.scatter(radial_pos,z_focus * 10**3, label='Fiber radial from module axis', color='red')
plt.xlabel('Radial Position [mm]')
plt.ylabel(r'Nominal focus position [$\mu$m]')
plt.title('Nominal fiber focus positions on BFS within module')
plt.legend()
plt.grid()
save.save_figures_to_dir('Nominal fiber focus positions on BFS within module', save_eps=False, dpi=400)

plt.figure(figsize=(10,6))
plt.scatter(range(len(z_diff)),z_diff * 10**3, color='red')
plt.xlabel('Neighbouring fibers')
plt.ylabel(r'Height difference [$\mu$m]')
plt.title('Focus differences between neighboring circles')
# plt.legend()
plt.grid()
save.save_figures_to_dir('Focus differences between neighboring circles', save_eps=False, dpi=400)
plt.show()


# For further information contact Maxime Rombach (EPFL - LASTRO) at maxime.rombach@epfl.ch
# This code is a tool for calculating the focal plane layouts of future highly multiplexed MOS telescope instruments.
# It is inspired from the work of Joe Silber (LBNL): https://github.com/joesilber/raft-design 
#!/usr/bin/env python
# -*- coding: utf-8 -*-
from Module import Module
from Focal_Suface import FocalSurf
from Grid import Grid

from matplotlib import pyplot as plt
import geopandas as gpd
import pandas as pd

import json

from shapely.plotting import plot_polygon
from shapely.ops import unary_union
from shapely.plotting import plot_polygon

import logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-4s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

""" Available projects: MUST, Spec-S5, WST25, WST27 """
project = "WST25"

""" Load focal surface data """
project_parameters = json.load(open('projects.json', 'r'))
surf = FocalSurf(project, **project_parameters[project])

""" Module parameters"""

nb_robots_per_module = 63 # number of robots per module
pitch = 6.2 # [mm] distance between two adjacent robots

mod = Module(nb_robots = nb_robots_per_module, 
             pitch = pitch)

""" Grid parameters """

inner_gap = 0.5 # [mm] gap between two modules within an intermediate triangle
global_gap = 4 # [mm] gap between two intermediate triangles; if inner = global, all modules are eqaully spaced

grid = Grid(project,
            inner_gap,
            global_gap,
            mod.module_side_length)

robots = mod.robots_layout
data = pd.DataFrame(mod.dataframe)
print(data)
geo = gpd.GeoDataFrame(data)
# geo.plot(color=geo['color'], alpha = 0.4, legend = True)
plot_polygon(mod.LR_coverage, fill = True, add_points=False, label = f'LR fibers: {len(mod.LR_robots)}')
plot_polygon(mod.HR_coverage, fill = True, add_points=False, color = 'red', label = f'HR fibers: {len(mod.HR_robots)}')
plt.scatter(mod.dataframe['x0'], mod.dataframe['y0'])
plot_polygon(mod.module_boundaries, fill = False, add_points=False, color = 'black')
plt.text(mod.dataframe['x0'][0], mod.dataframe['y0'][0], '1', fontsize=12, ha = 'center', va = 'center')
plt.scatter(mod.x0, mod.y0, color = 'red', s = 50, marker = 'o')
plt.grid(visible=True)
plt.legend()
plt.show()
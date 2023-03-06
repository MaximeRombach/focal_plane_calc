#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
## Parameters ##

la = 1.8 # alpha arm length [mm] /!\ lb > la /!\
lb = 1.8 # beta arm length [mm]
pitch = 6.2 # pitch [mm]
alpha = np.linspace(-180,180,180)
beta = np.linspace(-180,180,180)
## Raw triangle
# module_vertices_x = np.array([0,80,40,0])
# module_vertices_y = np.array([0,0,69.3,0])

## Edge cut module
module_vertices_x = np.array([7.5, 72.5, 76.25, 43.75, 36.25, 3.75, 7.5]) # [mm]
module_vertices_y = np.array([0, 0, 6.5, 62.8, 62.8, 6.5, 0]) # [mm]
module_width = 80 # [mm] triangle side length

beta2fibre = 1 # [mm] distance fiber center to edge of beta arm
safety_distance = 0.5 # [mm] self explicit for collision avoidance
offset_from_module_edges = safety_distance + beta2fibre
frame_thickness = 3 # [mm] thickness of frame between modules
start_offset_x = 6.2 # [mm]
start_offset_y = 3.41 # [mm]
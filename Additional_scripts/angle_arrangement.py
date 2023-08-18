#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
import math

## Parameters ##

Rc = 11067 # Curve radius [mm]
vigR = 613.2713
d_plane = 1200 # Focal plane diameter [mm]
module_side = 80 # Module width [mm]
triangle_height = math.cos(math.pi/6)*module_side # weight of triangle [mm]
print(triangle_height)
nb_modules = 8 #number of modules of one side of the axis
tolerance_envelop_width = 0.1 # [mm] tolerance in positioning around nominal focal surface 

t = Table.read('MM1536-cfg1-20210910.csv', comment='#')
Z_data = -t['Z']
CRD_data = t['CRD']
R_data = t['R']

x = np.linspace(-triangle_height/2,vigR,1000)
# begin=np.array([0])
x_disc = np.arange(-triangle_height/2,nb_modules*triangle_height,triangle_height)
# x_disc = np.concatenate((begin,x_disc))

BFS = lambda x : -np.sqrt(Rc**2-np.square(x)) + Rc
angle = np.tan(np.diff(BFS(x_disc))/np.diff(x_disc))
print(np.degrees(angle))
plt.plot(x,BFS(x),'g--',label='BFS')
plt.plot(x_disc,BFS(x_disc),'ro-')
plt.grid()
plt.show()
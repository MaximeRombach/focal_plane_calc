#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
# Or the path to the focal_plane_calc folder on your machine
# I know it is not the best way to do it but it works for now, To Be Fixed Later
sys.path.insert(0,'C:/Users/rombach/Documents/Astrobots/Inosuisse/focal_plane_calc')
from shapely.geometry import Polygon
from shapely.plotting import plot_polygon
from shapely.affinity import translate
from parameters import to_polygon_format
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

fiber1Diam = 0.1
fiber2Diam = 0.1
fiber1R = fiber1Diam/2
fiber2R = fiber2Diam/2  

# Define the polygons
def make_circle(r, x_center = 0, y_center = 0):

    n = 100
    x = x_center + r * np.cos(np.linspace(0, 2 * np.pi, n))
    y = y_center + r * np.sin(np.linspace(0, 2 * np.pi, n))

    return Polygon(to_polygon_format(x, y))

def dbtok(db):
    return 10**(db/10)
def ktodb(k):
    return 10*np.log10(k)
    

fiber1 = make_circle(fiber1Diam/2)
fiber2 = make_circle(fiber2Diam/2)
offsets = np.linspace(0, fiber1R+fiber2R, 100)
throughputs = []
losses = []

for t in offsets:
    fiber3 = translate(fiber2, xoff = 0, yoff = t)
    throughput = fiber1.intersection(fiber3).area / fiber1.area
    loss = 1-throughput

    throughputs.append(throughput)
    losses.append(loss)
    
fig, ax = plt.subplots()
ax.plot(offsets, np.asarray(losses))
ax.set_xlabel('Offset [mm]')
ax.set_ylabel('Losses [-]')
secay = ax.secondary_yaxis('right', functions=(dbtok, ktodb))
secay.set_ylabel('Losses [dB]')
ax.grid()
plt.show()

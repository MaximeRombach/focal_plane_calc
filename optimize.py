#%% 1) Define asph surface

import pandas as pd
from scipy.interpolate import interp1d
import parameters as param
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import logging

# Choose project
# Available: MUST, Megamapper
project = 'MUST'
surf = param.FocalSurf(project = project)
logging.basicConfig(level=logging.INFO)
logging.info('Project loaded: %s', project)

draw = True

if surf.asph_formula: # Check if focal surface defined directly by ashperic coefficients (analytical solution)
    R2Z = surf.asph_R2Z()

else: # If no coefficients load data from csv file and interpolate to get focal plane curve

    filename = "MM1536-cfg1-20210910.csv" # optics data from Zemax
    comment_character = "#"  # The character that indicates a commented line
    # Read CSV file and ignore commented lines
    optics_data = pd.read_csv(filename, comment=comment_character)

    # Set headers using the first non-comment line
    with open(filename, 'r') as file:
        for line in file:
            if not line.startswith(comment_character):
                headers = line.strip().split(',')
                break

    optics_data.columns = headers

    Z = optics_data['Z']
    CRD = optics_data['CRD']
    R = optics_data['R']

    R2Z = interp1d(R,Z,kind='cubic', fill_value = "extrapolate") #leave 'cubic' interpolation for normal vectors calculations
    R2CRD = interp1d(R,CRD,kind='cubic')


r = np.linspace(0,surf.vigR,500) # Define radius vector for focal plane curve
z = R2Z(r) # Calculate focal plane curve from csv data

#%% 2) Define BFS

BFS = surf.calc_BFS(r,z)
z_BFS = np.sqrt(BFS**2 - r**2) - BFS # Define z coordinate of BFS

plt.figure()
plt.title('Focal plane aspherical curve and Best Fit Sphere')
plt.plot(r,z, label = 'Aspherical curve')
plt.xlabel('r [mm]')
plt.ylabel('z [mm]')
plt.plot(r,z_BFS, '--', label = 'BFS')
plt.grid()
plt.legend()

plt.figure()
plt.title('Z error between aspherical curve and BFS')
plt.plot(r, z-z_BFS, label = 'Difference between asph and BFS')
plt.xlabel('r [mm]')
plt.ylabel('Error [mm]')
plt.grid()
plt.legend()

if draw:
    plt.show()

#%% 3) Load module pos and angles from focal_plane_coverage.py

#%% 4) Optimize module z pos and angles with loss functions
# 4.1) Define loss functions
#%% 4.2) Optimize

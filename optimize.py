#%% 1) Define asph surface

import pandas as pd
from scipy.interpolate import interp1d
import parameters as param
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import logging
import circle_fit as cf

# Choose project
# Available: MUST, Megamapper
project = 'MUST'
surf = param.FocalSurf(project = project)
logging.basicConfig(level=logging.INFO)
logging.info('Project loaded: %s', project)

save_txt = False
saving_df = {'save_txt': save_txt}
saving = param.SavingResults(saving_df)

draw = True
r = np.linspace(0,surf.vigR,500) # Define radius vector for focal plane curve
if surf.asph_formula:
    z_analytical = surf.analytical_R2Z(r)
    BFS_analytical = surf.calc_BFS(r,z_analytical)

# elif not surf.asph_formula and project == 'Spec-s5':
#     filename = "./Data_focal_planes/2023_10_16_Spec-s5.txt" # optics data from Zemax
#     Z,CRD,R = surf.read_focal_plane_data(filename) # Read data from csv file

filename = f"./Data_focal_planes/{project}.txt" # optics data from Zemax
optics_data = surf.read_focal_plane_data() # Read data from csv file

R2Z, R2CRD, R2NORM, R2S = surf.transfer_functions(optics_data)


z = R2Z(r) # Calculate focal plane curve from csv data


#%% Plot 3D focal plane from curve

# n = 100
# plt.rcParams["figure.figsize"] = [7.00, 3.50]
# plt.rcParams["figure.autolayout"] = True
# fig = plt.figure(figsize=(12,4))
# # ax1 = fig.add_subplot(121)
# ax2 = plt.subplot(111, projection='3d')

# x = np.linspace(0,surf.vigR,n)
# y = R2Z(x)
# t = np.linspace(0, np.pi*2, n)

# curve = np.array([x,y,np.zeros_like(x)]).reshape(3,n)
# saving.save_grid_to_txt(curve.T, 'curve')

# xn = np.outer(x, np.cos(t))
# yn = np.outer(x, np.sin(t))
# zn = np.zeros_like(xn)

# for i in range(len(x)):
#     zn[i:i+1,:] = np.full_like(zn[0,:], y[i])

# # ax1.plot(x, y)
# ax2.plot_surface(xn, yn, zn)
# ax2.set_zlim3d(-20, 0)    

#%% 2) Define BFS

BFS = surf.calc_BFS(r,z)

z_BFS = np.sqrt(BFS**2 - r**2) - BFS # Define z coordinate of BFS
tol = 50e-3
z_tol_up = z + tol*np.ones(len(z))
z_tol_down = z - tol*np.ones(len(z))

plt.figure()
plt.title('Focal plane aspherical curve and Best Fit Sphere')
plt.plot(r,z, label = 'Aspherical curve')
plt.plot(r,z_analytical, label = 'Aspherical analytical curve')
plt.xlabel('r [mm]')
plt.ylabel('z [mm]')
plt.plot(r,z_BFS, '--', label = f'BFS = {BFS :.1f} mm')
plt.plot(r,z_tol_up, '--', color = 'g', label = fr'Tolerance = $\pm$ {tol} mm')
plt.plot(r,z_tol_down, '--' , color = 'g',)
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

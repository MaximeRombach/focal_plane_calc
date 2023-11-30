#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
# Or the path to the focal_plane_calc folder on your machine
# I know it is not the best way to do it but it works for now, To Be Fixed Later
sys.path.insert(0,'C:/Users/rombach/Documents/Astrobots/Inosuisse/focal_plane_calc')

import numpy as np
import matplotlib.pyplot as plt
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

Rc = surf.curvature_R # Curve radius [mm]
vigR = surf.vigR
d_plane = 2*vigR # Focal plane diameter [mm]
nb_robots = 63
saving_df = {"save_plots": False, "save_dxf": False, "save_csv": False, "save_txt": False }
display_BFS = False
save = param.SavingResults(saving_df, project_name=project)
mod_param = param.Module(nb_robots, saving_df)
module_width = mod_param.module_width * np.sqrt(3)/2 #= triangle HEIGHT and NOT side
nb_modules = 9 #number of modules of one side of the axis
tolerance_envelop_width = surf.focus_tolerance_width # [mm] tolerance in positioning around nominal focal surface
start_offset = module_width/2

t = Table.read(f'Data_focal_planes/{project}.csv', comment='#', delimiter=';')
Z_data = -t['Z']
CRD_data = t['CRD']
R_data = t['R']
r = np.linspace(0,vigR,500)

BFS = surf.calc_BFS(R_data, Z_data)
R2Z_BFS = interp1d(r,-np.sqrt(BFS**2-np.square(r)) + BFS,kind='cubic', fill_value = "extrapolate")

R2Z = interp1d(R_data,Z_data,kind='cubic', fill_value = "extrapolate") #leave 'cubic' interpolation for normal vectors calculations
R2CRD = interp1d(R_data,CRD_data,kind='cubic')

z = R2Z(r) # Calculate focal plane curve from csv data

saving = param.SavingResults(saving_df)
to_csv_dict = {}

def calc_modules_pos(module_width, nb_modules, R2Z, on_BFS = False):

    r_modules = np.arange(start_offset,nb_modules*module_width+start_offset,module_width)
    if on_BFS: # Calculate module pos on Best Fit Sphere
        z_modules = R2Z_BFS(r_modules)
    else: # Calcualte module pos on interpolation of values from Berkeley optical team (csv file)
        z_modules = R2Z(r_modules) 

    r_diff = np.diff(r_modules)
    z_diff = np.diff(z_modules)

    print("x difference btw modules [mm]: ", np.around(r_diff,2))
    print("y difference btw modules [mm]: ", np.around(z_modules,2))

    to_csv_dict["x difference btw modules [mm]"] = np.around(r_diff,2)
    to_csv_dict["y difference btw modules [mm]"] = np.around(z_diff,2)

    #x_modules = np.concatenate((np.flip(-x_modules,0), x_modules), axis = 0)
    #y_modules = np.concatenate((np.flip(y_modules,0), y_modules), axis = 0)

    return r_modules, z_modules

def get_normals(r_modules, R2Z, h=1e-5):
    """Calculates the normal vectors to the curve
    Their coordinates are stored in "normal" element-wise ex: normal[0]=[nx0,ny0]"""
    normal = []
    for idx in range(len(r_modules)):
        x0, y0, xa, ya = r_modules[idx], R2Z(r_modules[idx]), r_modules[idx]+h, R2Z(r_modules[idx]+h)
        dx = h
        dy = ya - y0
        norm = math.hypot(dx, dy)
        dx /= norm
        dy /= norm
        normal.append([-dy, dx])

    return np.array(normal)


def get_normals_angles(normal, print_data = True):
    """Calculates the angle between the normal vectors and the vertical"""
    angles = []

    for idx in range(len(normal)):
        theta = np.arccos(np.dot(normal[idx],np.array([0,1])))
        angles.append(theta)
    
    angles = np.array(angles)
    angles_deg = np.around(np.degrees(angles),3)
    to_csv_dict["angle to normal [deg]"] = 0
    if print_data:
        print("angle to normal [deg]: ", angles_deg)

    

    return angles
  

def rotate_modules_tangent_2_surface(x_modules, y_modules, angles):
    new_left_module = []
    new_right_module = []

    for px,py,theta in zip(x_modules,y_modules,angles):
        p = np.array([px,py,1]).reshape(3,1)
        left_module = p - np.array([module_width/2,0,0]).reshape(3,1)
        right_module = p + np.array([module_width/2,0,0]).reshape(3,1)

        c,s = np.cos(theta), np.sin(theta)
        R_theta = np.array([[c,-s,0],[s,c,0],[0,0,1]]) #
        Rt_minus= np.array([[1,0,-px],[0,1,-py],[0,0,1]])
        Rt_plus= np.array([[1,0,px],[0,1,py],[0,0,1]])
        new_left_module.append(Rt_plus @ R_theta @ Rt_minus @ left_module)
        new_right_module.append(Rt_plus @ R_theta @ Rt_minus @ right_module)

    return np.array(new_left_module), np.array(new_right_module)

def translate_module_in_envelop(dist, module_number, draw = True):

    if not draw:
        return

    new_left_module_translated = []
    new_right_module_translated = []
    for (theta,module_l,module_r) in zip(angles, new_left_module, new_right_module):
        dx, dy= np.sin(theta)*dist, np.cos(theta)*dist
        Rt= np.array([[1,0,dx],[0,1,dy],[0,0,1]])
        new_left_module_translated.append(Rt @ module_l)
        new_right_module_translated.append(Rt @ module_r)
    
    plt.plot((new_left_module_translated[module_number][0],new_right_module_translated[module_number][0]),
                        (new_left_module_translated[module_number][1],new_right_module_translated[module_number][1]),'k',
                        label="Module translated ({} $\mu$m)".format(int(dist*1000)))

def tolerance_envelop(r,R2Z):

    normal_vectors = get_normals(r,R2Z)
    normal_angles = get_normals_angles(normal_vectors, print_data = False)
    dx, dy= np.sin(normal_angles)*tolerance_envelop_width/2, np.cos(normal_angles)*tolerance_envelop_width/2
    
    r_envelop_plus = r + dx
    z_envelop_plus = z + dy

    r_envelop_minus = r - dx
    z_envelop_minus = z - dy

    return r_envelop_plus, z_envelop_plus, r_envelop_minus, z_envelop_minus

def draw_modules_oriented(plot_axis=0, draw=True, draw_one = False, module_number = None):
    if not draw:
        return
    if not draw_one:
        for i in range(len(new_left_module)): # draw modules width
            if i == 0:
                plt.plot((new_left_module[i][0],new_right_module[i][0]),
                        (new_left_module[i][1],new_right_module[i][1]), color='b', label="Tangent modules", linewidth=1.25)
            else:
                plt.plot((new_left_module[i][0],new_right_module[i][0]),
                        (new_left_module[i][1],new_right_module[i][1]))
    else:
        plt.plot((new_left_module[module_number][0],new_right_module[module_number][0]),
                        (new_left_module[module_number][1],new_right_module[module_number][1]), color='b', label="Tangent module", linewidth=1.25)

def draw_normals(x_modules, y_modules, normal, length=10, draw=True):
    if not draw:
        return
    else:
        for idx in range(len(x_modules)):
            x0, y0 = x_modules[idx], y_modules[idx]
            nx, ny = normal[idx][0]*length, normal[idx][1]*length
            plt.plot((x0, x0+nx), (y0, y0+ny))
            # plt.quiver(x0, y0,  y0-ny, x0+nx)
            # axs[0].plot((x0, x0), (y0, y0+1*length),'r--')def draw_normals(x_modules, y_modules, normal, length=10, draw=True):
    if not draw:
        return
    else:
        for idx in range(len(x_modules)):
            x0, y0 = x_modules[idx], y_modules[idx]
            nx, ny = normal[idx][0]*length, normal[idx][1]*length
            plt.plot((x0, x0+nx), (y0, y0+ny))
            # plt.quiver(x0, y0,  y0-ny, x0+nx)
            # axs[0].plot((x0, x0), (y0, y0+1*length),'r--')
def draw_modules_horizontal(plot_axis=0, draw=True, draw_one = False, module_number = None):
    if not draw:
        return
    
    if not draw_one:
        # Draw all modules
        for i,data in enumerate(x_modules): # draw modules width
            x_min = x_modules[i]-module_width/2
            x_max = x_modules[i]+module_width/2

            if i == 0: # Draw all the modules
                plt.plot((x_min,x_max),(y_modules[i], y_modules[i]), color='r', label="Horizontal modules", linewidth=1.25)
            else:
                plt.plot((x_min,x_max),(y_modules[i], y_modules[i]), color='r')
    else: # Draw only module of interest
        x_min, x_max = x_modules[module_number]-module_width/2, x_modules[module_number]+module_width/2

        plt.plot((x_min,x_max),(y_modules[module_number], y_modules[module_number]), color='r', label="Horizontal modules", linewidth=1.25)

def draw_BFS(BFS,vigR, display_BFS=False, full_curve = False):
    """Draws the Best Fit Sphere of the focal surface for comparison with the interpolated one
    full_curve criteria draws both sides of the curve"""
    if not display_BFS:
        return
    if full_curve: 
        x = np.linspace(-vigR,vigR,100)
    else:
        x = np.linspace(0,vigR,100)

    y = R2Z_BFS(x)

    plt.plot(x,y,label=f'BFS = {int(BFS)} mm',color='orange',linestyle='-.')
    plt.legend()

def draw_R2CRD(r,R_data,CRD_data,R2CRD, draw=True):
    if not draw:
        return
    fig = plt.figure(num='R2CRD')
    plt.plot(r,R2CRD(r),label='interpolation')
    plt.scatter(R_data,CRD_data,marker='.',color='red',label='data')
    plt.legend()
    plt.title('Chief Ray Deviation in terms of radial pos on focal surface')
    plt.xlabel('R [mm]')
    plt.ylabel('CRD [deg]')
    plt.grid()

def draw_envelop(plot_axis = 0,draw = True, draw_legend=False):
    
    envelop_looks = 'r--'
    if not draw:
        return
    if draw_legend:
        plt.plot(r_envelop_minus, z_envelop_minus, envelop_looks, linewidth = 0.75,
                            label='Tolerance envelop = {} $\mu$m'.format(int(tolerance_envelop_width*10**3)))
    else:
        plt.plot(r_envelop_minus, z_envelop_minus, envelop_looks, linewidth = 0.75)

    plt.plot(r_envelop_plus, z_envelop_plus, envelop_looks, linewidth = 0.75)

def zoom_in_1module(module_number, xbound, ybound, draw=True):

    if not draw:
        return
    
    figtitle = 'Project: $\\bf{'+f'{project}''}$'+'\nZoom in module #{} - {} robots'.format(module_number+1, nb_robots)
    fig = plt.figure('Zoom in module', figsize=(15,5))
    x_center, y_center = x_modules[module_number], y_modules[module_number]
    xlim_min, xlim_max, ylim_min, ylim_max = x_center - xbound, x_center + xbound,  y_center - ybound, y_center + ybound

    plt.plot(r,z,'--g',label="Focal surface")
    draw_modules_horizontal(plot_axis=1, draw=True, draw_one=True, module_number=module_number)
    draw_envelop(plot_axis=1, draw_legend=True)
    draw_modules_oriented(plot_axis=1, draw_one=True, module_number=module_number)
    # draw_normals(x_modules,y_modules,normal,length=1)
    translate_module_in_envelop(dist=0.025, module_number=module_number)
    draw_BFS(BFS,vigR, display_BFS=display_BFS, full_curve=False)

    plt.xlim(xlim_min, xlim_max)
    plt.ylim(ylim_min, ylim_max)
    plt.title(figtitle)
    plt.legend(shadow=True)
    plt.xlabel('R [mm]')
    plt.ylabel('Z [mm]')
    plt.grid()

def makemoduleline(x, start_left, end_right):
    a = (end_right[1]-start_left[1])/(end_right[0]-start_left[0])
    b = start_left[1] - a*start_left[0]

    return a*x + b
    

def normIntersects(normals, module_number, module_left, module_right):
    r_curve = np.array([])
    a = normals[module_number][1]/normals[module_number][0]
    x_iterable = np.linspace(module_left[module_number][0], module_right[module_number][0], 100)
    y = makemoduleline(x_iterable, module_left[module_number], module_right[module_number])
    b_iterable = y - a*x_iterable

    for b,x in zip(b_iterable, x_iterable):

        def fun(r1):
            return a*r1 + b - R2Z(r1)
    
        initial_guess = x

        r_solved = fsolve(fun, initial_guess)

        r_curve = np.append(r_curve, r_solved, axis=0)
    
    z_curve = R2Z(r_curve)

    return x_iterable, y, r_curve, z_curve

def dist2curve(r_module,z_module,r_curve,z_curve):
    """Calculate the distance between a point and a curve"""
    dr = r_curve-r_module.reshape(1,len(r_curve))
    dz = z_curve-z_module.reshape(1,len(z_curve))
    return np.sqrt(dr**2+dz**2)


x_modules, y_modules = calc_modules_pos(module_width, nb_modules, R2Z, on_BFS=False)
normals = get_normals(x_modules, R2Z)
angles = get_normals_angles(normals)
new_left_module, new_right_module = rotate_modules_tangent_2_surface(x_modules, y_modules, angles)
r_envelop_plus, z_envelop_plus, r_envelop_minus, z_envelop_minus = tolerance_envelop(r,R2Z)

module_number = 3
r_module, z_module, r_curve, z_curve = normIntersects(normals, module_number, new_left_module, new_right_module)
d = dist2curve(r_module,z_module,r_curve,z_curve)


## Plotting time ##

# fig = plt.figure(figsize=(12,8))
fig = plt.figure('Full focal surface',figsize=(15,5))
# fig, ax = plt.subplots()

draw = True
plot_time = 20 #seconds
is_timer = False
save_fig = True

draw_normals(x_modules,y_modules,normals,length=vigR)
# plt.plot(x_radius, y_radius, '-.',label="BFS")
draw_modules_horizontal()
draw_modules_oriented(draw = True)
plt.plot(r,z,'--g',label="Focal surface")
draw_BFS(BFS,vigR, display_BFS=display_BFS)
plt.scatter(r_curve,z_curve, marker="o", color="C1")
# plt.plot((r1,new_left_module[3][0]),(z1, new_left_module[3][1]),'--ok')
# z_intersect = normIntersects(normal, 3, new_left_module)
# print(z_intersect)
# plt.plot(z, np.interp(z_intersect, r, z), marker="o", ls="", ms=4, color="C1")
plt.legend(shadow=True)
plt.title('Project: $\\bf{'+f'{project}''}$'+f'\nFocal surface and module arrangement with normal vectors - {nb_robots} robots/module')
plt.xlabel('R [mm]')
plt.ylabel('Z [mm]')
plt.grid()
save.save_figures_to_dir('Focal surface and module arrangement with normal vectors')

module_number = 3
if nb_robots == 102:
    xbound=45
    ybound=1.2
elif nb_robots == 63:
    xbound=35
    ybound=0.75
elif nb_robots == 75:
    xbound=40
    ybound=0.9
zoom_in_1module(module_number = module_number, xbound=xbound, ybound=ybound, draw=True)
save.save_figures_to_dir(f'Zoomed_in_module_{module_number+1}')

plt.figure('dis2curve')
plt.scatter(r_module, d*10**3)
plt.plot((r_module[0], r_module[-1]), (tolerance_envelop_width/2*10**3, tolerance_envelop_width/2*10**3), 'r--')
plt.grid()
plt.xlabel('Position on fiber tips plane [mm]')
plt.ylabel(r'Normal distance to focal surface [$\mu$m]')
plt.title(f'Distance between module {module_number+1} and focal surface')

draw_R2CRD(r,R_data,CRD_data,R2CRD, draw=False)




if is_timer:
    plt.show(block=False)
    plt.pause(plot_time)
    plt.close('all')
else:
    plt.show()
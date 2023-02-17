#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.table import Table, vstack
from scipy.interpolate import interp1d

## Parameters ##

Rc = 11067 # Curve radius [mm]
vigR = 613.2713
d_plane = 1200 # Focal plane diameter [mm]
module_width = 80 # Module width [mm]
nb_modules = 8 #number of modules of one side of the axis
tolerance_envelop_width = 0.1 # [mm] tolerance in positioning around nominal focal surface 

t = Table.read('MM1536-cfg1-20210910.csv', comment='#')
Z_data = -t['Z']
CRD_data = t['CRD']
R_data = t['R']

R2Z = interp1d(R_data,Z_data,kind='cubic', fill_value = "extrapolate") #leave 'cubic' interpolation for normal vectors calculations
R2CRD = interp1d(R_data,CRD_data,kind='cubic')
r = np.linspace(0,vigR,1000)
z = R2Z(r) # Calculate focal plane curve from csv data

def calc_modules_pos(Rc, module_width, nb_modules, R2Z, BFS = False):

    r_modules = np.array(range(0,nb_modules*module_width,module_width))
    if BFS: # Calculate module pos on Best Fit Sphere
        z_modules = -np.sqrt(Rc**2-np.square(r_modules)) + Rc
    else: # Calcualte module pos on interpolation of values from Berkeley optical team (csv file)
        z_modules = R2Z(r_modules) 

    r_diff = np.diff(r_modules)
    z_diff = np.diff(z_modules)

    print("x difference btw modules [mm]: ", np.around(r_diff,2))
    print("y difference btw modules [mm]: ", np.around(z_diff,2))

    #x_modules = np.concatenate((np.flip(-x_modules,0), x_modules), axis = 0)
    #y_modules = np.concatenate((np.flip(y_modules,0), y_modules), axis = 0)

    return r_modules, z_modules

def get_normals(r_modules, R2Z, h=0.001):
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
    if print_data:
        print("angle to normal [deg]: ",   np.around(np.degrees(angles),3))

    return angles
  
def draw_normals(x_modules, y_modules, normal, length=10, plot_axis=0, draw=True):
    if not draw:
        return
    else:
        for idx in range(len(x_modules)):
            x0, y0 = x_modules[idx], y_modules[idx]
            nx, ny = normal[idx][0]*length, normal[idx][1]*length
            axs[plot_axis].plot((x0, x0+nx), (y0, y0+ny))
            # axs[0].plot((x0, x0), (y0, y0+1*length),'r--')

def rotate_modules_2_surface(x_modules, y_modules, angles):
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

def draw_modules(plot_axis=0, draw=True, draw_one = False, module_number = None):
    if not draw:
        return
    if not draw_one:
        for i in range(len(new_left_module)): # draw modules width
            if i == 0:
                axs[plot_axis].plot((new_left_module[i][0],new_right_module[i][0]),
                        (new_left_module[i][1],new_right_module[i][1]), color='b', label="tangent modules", linewidth=1.25)
            else:
                axs[plot_axis].plot((new_left_module[i][0],new_right_module[i][0]),
                        (new_left_module[i][1],new_right_module[i][1]))
    else:
        axs[plot_axis].plot((new_left_module[module_number][0],new_right_module[module_number][0]),
                        (new_left_module[module_number][1],new_right_module[module_number][1]), color='b', label="tangent module", linewidth=1.25)

def tolerance_envelop(r,R2Z):

    normal_vectors = get_normals(r,R2Z)
    normal_angles = get_normals_angles(normal_vectors, print_data = False)
    dx, dy= np.sin(normal_angles)*tolerance_envelop_width/2, np.cos(normal_angles)*tolerance_envelop_width/2
    
    r_envelop_plus = r + dx
    z_envelop_plus = z + dy

    r_envelop_minus = r - dx
    z_envelop_minus = z - dy

    return r_envelop_plus, z_envelop_plus, r_envelop_minus, z_envelop_minus

def draw_modules_width(plot_axis=0, draw=True, draw_one = False, module_number = None):
    if not draw:
        return
    
    if not draw_one:
        # Draw all modules
        for i,data in enumerate(x_modules): # draw modules width
            x_min = x_modules[i]-module_width/2
            x_max = x_modules[i]+module_width/2

            
            if i == 0: # Draw all the modules
                axs[plot_axis].plot((x_min,x_max),(y_modules[i], y_modules[i]), color='r', label="fiber tips line", linewidth=1.25)
            else:
                axs[plot_axis].plot((x_min,x_max),(y_modules[i], y_modules[i]), color='r')
    else: # Draw only module of interest
        x_min, x_max = x_modules[module_number]-module_width/2, x_modules[module_number]+module_width/2

        axs[plot_axis].plot((x_min,x_max),(y_modules[module_number], y_modules[module_number]), color='r')

def draw_BFS(Rc,vigR, draw=False, full_curve = False):
    """Draws the Best Fit Sphere of the focal surface for comparison with the interpolated one
    full_curve criteria draws both sides of the curve"""
    if not draw:
        return
    if full_curve: 
        x = np.linspace(-vigR,vigR,100)
    else:
        x = np.linspace(0,vigR,100)

    y = -np.sqrt(Rc**2-np.square(x)) + Rc

    axs[0].plot(x,y,label='BFS',color='orange',linestyle='-.')

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
        axs[plot_axis].plot(r_envelop_minus, z_envelop_minus, envelop_looks, linewidth = 0.75,
                            label='Tolerance envelop = {} $\mu$m'.format(tolerance_envelop_width*10**3))
    else:
        axs[plot_axis].plot(r_envelop_minus, z_envelop_minus, envelop_looks, linewidth = 0.75)

    axs[plot_axis].plot(r_envelop_plus, z_envelop_plus, envelop_looks, linewidth = 0.75)

def zoom_in_1module(module_number, xbound, ybound, plot_axis=1, draw=True):

    if not draw:
        return
    
    x_center, y_center = x_modules[module_number], y_modules[module_number]
    xlim_min, xlim_max, ylim_min, ylim_max = x_center - xbound, x_center + xbound,  y_center - ybound, y_center + ybound

    axs[plot_axis].plot(r,z,'--g',label="focal surface")
    draw_modules_width(plot_axis=1, draw=False, draw_one=True, module_number=module_number)
    draw_envelop(plot_axis=1, draw_legend=True)
    draw_modules(plot_axis=1, draw_one=True, module_number=module_number)

    axs[plot_axis].set_xlim(xlim_min, xlim_max)
    axs[plot_axis].set_ylim(ylim_min, ylim_max)
    axs[plot_axis].legend()
    axs[plot_axis].set_xlabel('R [mm]')
    axs[plot_axis].set_ylabel('Z [mm]')
    axs[plot_axis].grid()


x_modules, y_modules = calc_modules_pos(Rc, module_width, nb_modules, R2Z)
normal = get_normals(x_modules, R2Z)
angles = get_normals_angles(normal)
new_left_module, new_right_module = rotate_modules_2_surface(x_modules, y_modules, angles)
r_envelop_plus, z_envelop_plus, r_envelop_minus, z_envelop_minus = tolerance_envelop(r,R2Z)
# print(new_left_module)

## Plotting time ##

fig, axs = plt.subplots(2,1,figsize=(12,8))
# fig, ax = plt.subplots()

plot_time = 20 #seconds
is_timer = False

draw_normals(x_modules,y_modules,normal,length=10)
# plt.plot(x_radius, y_radius, '-.',label="BFS")
draw_modules_width()
draw_modules(draw = True)
axs[0].plot(r,z,'--g',label="focal surface")
draw_BFS(Rc,vigR, draw=False)
draw_envelop()
zoom_in_1module(module_number = 3, xbound=45, ybound=1, draw=True)
axs[0].legend(shadow=True)
axs[0].set_title('Focal surface and module arrangement')
axs[0].set_xlabel('R [mm]')
axs[0].set_ylabel('Z [mm]')
axs[0].grid()
draw_R2CRD(r,R_data,CRD_data,R2CRD, draw=False)


if is_timer:
    plt.show(block=False)
    plt.pause(plot_time)
    plt.close('all')
else:
    plt.show()
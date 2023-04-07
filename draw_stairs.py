#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.table import Table, vstack
from scipy.interpolate import interp1d
import parameters as param

## Parameters ##

Rc = 11067 # Curve radius [mm]
vigR = 613.2713
d_plane = 1200 # Focal plane diameter [mm]
module_width = param.module_width * np.sqrt(3)/2 #= triangle HEIGHT and NOT side
# module_width = 74*np.sqrt(3)/2 # Module width [mm] = triangle HEIGHT and NOT side
nb_modules = 9 #number of modules of one side of the axis
tolerance_envelop_width = 0.1 # [mm] tolerance in positioning around nominal focal surface
start_offset = 0

t = Table.read('MM1536-cfg1-20210910.csv', comment='#')
Z_data = -t['Z']
CRD_data = t['CRD']
R_data = t['R']

R2Z = interp1d(R_data,Z_data,kind='cubic', fill_value = "extrapolate") #leave 'cubic' interpolation for normal vectors calculations

R2CRD = interp1d(R_data,CRD_data,kind='cubic')
r = np.linspace(0,vigR,1000)
z = R2Z(r) # Calculate focal plane curve from csv data

def calc_modules_pos(Rc, module_width, nb_modules, R2Z, BFS = False):

    r_modules = np.arange(start_offset,nb_modules*module_width+start_offset,module_width)
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
  
def draw_normals(x_modules, y_modules, normal, length=10, draw=True):
    if not draw:
        return
    else:
        for idx in range(len(x_modules)):
            x0, y0 = x_modules[idx], y_modules[idx]
            nx, ny = normal[idx][0]*length, normal[idx][1]*length
            plt.plot((x0, x0+nx), (y0, y0+ny))
            # plt.quiver(x0, y0,  y0-ny, x0+nx)
            # axs[0].plot((x0, x0), (y0, y0+1*length),'r--')

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


def draw_modules(plot_axis=0, draw=True, draw_one = False, module_number = None):
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
                plt.plot((x_min,x_max),(y_modules[i], y_modules[i]), color='r', label="Horizontal modules", linewidth=1.25)
            else:
                plt.plot((x_min,x_max),(y_modules[i], y_modules[i]), color='r')
    else: # Draw only module of interest
        x_min, x_max = x_modules[module_number]-module_width/2, x_modules[module_number]+module_width/2

        plt.plot((x_min,x_max),(y_modules[module_number], y_modules[module_number]), color='r', label="Horizontal modules", linewidth=1.25)

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

    plt.plot(x,y,label='BFS',color='orange',linestyle='-.')

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

def zoom_in_1module(module_number, xbound, ybound, save_fig, draw=True):

    if not draw:
        return
    
    figtitle = 'Zoom in module #{} - {} robots'.format(module_number+1, param.nb_robots)
    fig = plt.figure(figtitle,figsize=(15,5))
    x_center, y_center = x_modules[module_number], y_modules[module_number]
    xlim_min, xlim_max, ylim_min, ylim_max = x_center - xbound, x_center + xbound,  y_center - ybound, y_center + ybound

    plt.plot(r,z,'--g',label="Focal surface")
    draw_modules_width(plot_axis=1, draw=True, draw_one=True, module_number=module_number)
    draw_envelop(plot_axis=1, draw_legend=True)
    draw_modules(plot_axis=1, draw_one=True, module_number=module_number)
    translate_module_in_envelop(dist=0.025, module_number=module_number)

    plt.xlim(xlim_min, xlim_max)
    plt.ylim(ylim_min, ylim_max)
    plt.title(figtitle)
    plt.legend(shadow=True)
    plt.xlabel('R [mm]')
    plt.ylabel('Z [mm]')
    plt.grid()
    param.save_figures_to_dir(save_fig, figtitle)


x_modules, y_modules = calc_modules_pos(Rc, module_width, nb_modules, R2Z)
normal = get_normals(x_modules, R2Z)
angles = get_normals_angles(normal)
new_left_module, new_right_module = rotate_modules_tangent_2_surface(x_modules, y_modules, angles)
r_envelop_plus, z_envelop_plus, r_envelop_minus, z_envelop_minus = tolerance_envelop(r,R2Z)
# print(new_left_module)

## Plotting time ##

# fig = plt.figure(figsize=(12,8))
fig = plt.figure('Full focal surface',figsize=(15,5))
# fig, ax = plt.subplots()

draw = True
plot_time = 20 #seconds
is_timer = False
save_fig = True

draw_normals(x_modules,y_modules,normal,length=10)
# plt.plot(x_radius, y_radius, '-.',label="BFS")
draw_modules_width()
draw_modules(draw = True)
plt.plot(r,z,'--g',label="Focal surface")
draw_BFS(Rc,vigR, draw=False)
draw_envelop()
plt.legend(shadow=True)
plt.title('Focal surface and module arrangement with normal vectors')
plt.xlabel('R [mm]')
plt.ylabel('Z [mm]')
plt.grid()
param.save_figures_to_dir(save_fig, 'Focal surface and module arrangement with normal vectors')

zoom_in_1module(module_number = 3, xbound=45, ybound=1, save_fig = save_fig, draw=True)
draw_R2CRD(r,R_data,CRD_data,R2CRD, draw=False)


if is_timer:
    plt.show(block=False)
    plt.pause(plot_time)
    plt.close('all')
else:
    plt.show()
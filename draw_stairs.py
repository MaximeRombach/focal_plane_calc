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

t = Table.read('MM1536-cfg1-20210910.csv', comment='#')
t1 = -t['Z']
CRD = t['CRD']

R2Z = interp1d(t['R'],t1,kind='cubic') #leave 'cubic' interpolation for normal vectors calculations
r = np.linspace(0,vigR,1000)
z = R2Z(r) # Calculate focal plane curve from csv data

def get_normals(r_modules, R2Z, h=0.001):
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

def get_normals_angles(normal):
    print(np.linalg.norm(normal[0]))



def draw_normals(x_modules, y_modules, normal, length=10, draw=True):
    if not draw:
        return
    else:
        for idx in range(len(x_modules)):
            x0, y0 = x_modules[idx], y_modules[idx]
            nx, ny = normal[idx][0]*length, normal[idx][1]*length
            ax.plot((x0, x0+nx), (y0, y0+ny))
            # ax.plot((x0, x0), (y0, y0+1*length),'r--')


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

def draw_BFS(Rc,vigR, draw=False, full_curve = False):

    if not draw:
        return
    if full_curve: 
        x = np.linspace(-vigR,vigR,100)
    else:
        x = np.linspace(0,vigR,100)

    y = -np.sqrt(Rc**2-np.square(x)) + Rc

    ax.plot(x,y,label='BFS',color='orange',linestyle='-.')

x_modules, y_modules = calc_modules_pos(Rc, module_width, nb_modules, R2Z)
normal = get_normals(x_modules, R2Z)
get_normals_angles(normal)

## Plotting time ##

fig, ax = plt.subplots(1,1,figsize=(15,1.5))

plot_time = 20 #seconds
is_timer = False
if is_timer:
    timer = fig.canvas.new_timer(interval = plot_time*1000)
    timer.add_callback(plt.close)

draw_normals(x_modules,y_modules,normal)
# plt.plot(x_radius, y_radius, '-.',label="BFS")
for i in range(len(x_modules)): # draw modules width
    x_min = x_modules[i]-module_width/2
    x_max = x_modules[i]+module_width/2
    if i == 0:
        ax.hlines(y=y_modules[i], xmin=x_min, xmax=x_max, color='r', label="fiber tips line", linewidth=1.25)
    else:
        ax.hlines(y=y_modules[i], xmin=x_min, xmax=x_max, color='r')
plt.plot(r,z,'--g',label="focal surface")
draw_BFS(Rc,vigR, draw=False)
plt.legend()
plt.title('Focal surface and module arrangement')
plt.xlabel('R [mm]')
plt.ylabel('Z [mm]')
plt.grid()

if is_timer:
    timer.start()

plt.show()
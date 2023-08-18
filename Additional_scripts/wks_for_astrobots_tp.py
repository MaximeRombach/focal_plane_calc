import numpy as np
import parameters as param
from shapely.geometry import Polygon, MultiPolygon
from shapely.plotting import plot_polygon
import matplotlib.pyplot as plt


la = 1.8 # alpha arm length [mm] /!/ lb > la /!/
lb = 3 # beta arm length [mm]
pitch = 6.2 # pitch [mm]
alpha = np.linspace(-180,180,180) # [deg] rotational range of alpha arm
beta = np.linspace(-180,180,180) # [deg] rotational range of beta arm (holding the fiber)

size = la + lb
angles_hex = np.linspace(60,360,6)
x_hex = np.zeros(7)
y_hex = np.zeros(7)
x_hex[1:] = size * np.cos(np.deg2rad(angles_hex))
y_hex[1:] = size * np.sin(np.deg2rad(angles_hex))

c1 = np.cos(np.deg2rad(alpha))
s1 = np.sin(np.deg2rad(alpha))

c2 = np.cos(np.deg2rad(beta))
s2 = np.sin(np.deg2rad(beta))

xa, ya, xab, yab = (lb-la)*c1, (lb-la)*s1, (lb+la)*c2, (la+lb)*s2

wks_list = []

for idx, (dx, dy) in enumerate(zip(x_hex, y_hex)):

        xa1, ya1, xab1, yab1 = xa + dx, ya + dy, xab + dx, yab + dy

        coords_int = param.to_polygon_format(xa1, ya1)
        coords_ext = param.to_polygon_format(xab1, yab1)

        interior = coords_int[::-1]
        poly_c1 = Polygon(coords_ext, [interior])
        wks_list.append(poly_c1)

multi_wks = MultiPolygon(wks_list)

if la == lb:
      option = "l_alpha = l_beta"
elif la < lb:
      option = "l_alpha < l_beta"
elif la > lb:
    option = "l_alpha > l_beta"

fig = plt.figure(figsize=(8,8))
plt.title("Workspace for 7 positioners " + option)
plt.scatter(x_hex,y_hex, marker='o', color='k')

for idx, wks in enumerate(wks_list):
    plot_polygon(wks, add_points=False, alpha=0.2, facecolor='red', edgecolor='black')

plt.grid()
plt.savefig('C:/Users/rombach/Dropbox/Apps/Overleaf/TP_Astrobots/figures/02_SCARA_theory/solution_q13.png', format='png')
plt.show()

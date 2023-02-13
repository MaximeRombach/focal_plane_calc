import numpy as np
import matplotlib.pyplot as plt
import math

## This code is for legacy. When printing the normal to a curve, it may not
## appear normal to thesaid curve. That might be due to the plot window
## dimensions that don't show it clearly
## Try playing with the following plot's dimensions to see the issue
## And also check your maths first

def get_normals(length=10, indexing = True):
    
    for idx in range(len(x)-1):
        x0, y0, xa, ya = x[idx], y[idx], x[idx+1], y[idx+1]
        dx, dy = xa-x0, ya-y0
        norm = math.hypot(dx, dy) * 1/length
        dx /= norm
        dy /= norm
        
        ax.plot((x0, x0-dy), (y0, y0+dx))    # plot the normals

fig, ax = plt.subplots()
plt.rcParams["figure.figsize"] = [8, 8]
x=np.linspace(-1,1,100)
y=x**2
plt.plot(x,y)
get_normals(length=-.3)
plt.show()
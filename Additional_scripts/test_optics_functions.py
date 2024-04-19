import numpy as np
from numpy.polynomial import Polynomial
import matplotlib.pyplot as plt

blur2loss = Polynomial([0, -0.000141553, 0.000373672, -1.76888E-06, -8.23219E-08, 8.72644E-10])
tilt2loss = Polynomial.fit(x=[0.0, 0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975],
                                    y=1-np.array([1.0, 0.99995, 0.99975, 0.99945, 0.999, 0.99845, 0.99775, 0.9969, 0.99595, 0.9949, 0.99365, 0.99225, 0.9907, 0.98895, 0.98705, 0.985, 0.98275, 0.9803, 0.9776, 0.97465, 0.9715]),
                                    deg=6,
                                    ),  # input units deg, I am fitting here since the polynomial coefficients provided by Excel chart in Fiber Tilt tab are poor
f_number = 3.6
defocus2blur = lambda dz_mm: (dz_mm*1000) / 2 / f_number / 3
dz = np.linspace(-0.5, 0.5, 1000)
d_deg = np.linspace(-1, 1, 1000)
blur = defocus2blur(dz)
loss_defocus = blur2loss(blur)
plt.figure()
plt.title('Defocus2loss')
plt.plot(dz, loss_defocus)
plt.xlabel('Defocus (mm)')
plt.ylabel('Losses (%)')
plt.grid()

plt.figure()
plt.title('Tilt2loss')
plt.plot(d_deg, tilt2loss[0](d_deg))
plt.xlabel('dAngle (°)')
plt.ylabel('Losses (%)')
plt.grid()
plt.show()
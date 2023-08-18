# # Series of tests for focal_plane coverage
# import parameters as par
# import numpy as np

# def test_pitch(pitch, test_pitch):
#     tol = 1e-3
#     if pitch - test_pitch <= tol:
#         print('Pitch test passed')
#         print(test_pitch)
#     else:
#         print('Pitch test failed: review x_inc and y_inc definition')

# def main():
#     test_pitch(par.pitch, par.test_pitch)

# if __name__ == '__main__':
#     main()

import numpy as np
import matplotlib.pyplot as plt

# Define the center of the sphere
center = np.array([0, 0, 0])

# Define the radius of the sphere
radius = 110

# Define the points to project onto the spherical surface
points = np.array([[1, 2, 0], [4, 5, 0], [7, 8, 0]])

# Normalize the points to have unit length
norm_points = points / np.linalg.norm(points, axis=1)[:, np.newaxis]

# Compute the projected points on the spherical surface
projected_points = center + radius * norm_points

print(projected_points)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(points[:,0], points[:,1], points[:,2] , label="original", color='red')
ax.scatter(projected_points[:,0], projected_points[:,1], projected_points[:,2] , label=f'projected', color='blue')
ax.set_box_aspect((1, 1, 1))
plt.legend()
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()
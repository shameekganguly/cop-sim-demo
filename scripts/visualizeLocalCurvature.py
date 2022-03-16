# visualizeLocalCurvature.py

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

'''
CAUTION: This code visualizes a model that we were considering for local curvature
for a spatial surface. But we have discontinued that model since it produces a
non-C2 continuous local surface. Instead we are now using a torus approximation.
'''

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    # ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


RA = 3
RB = 1
def R(theta):
	return (RA-RB)*np.sin(theta) + RA

thetas = np.linspace(0, np.pi, 101)
betas = np.linspace(-np.pi/2, np.pi/2, 201)

num_points = thetas.size * betas.size
xarr = np.zeros(num_points)
yarr = np.zeros(num_points)
zarr = np.zeros(num_points)

for i in range(thetas.size):
	for j in range(betas.size):
		theta = thetas[i]
		beta = betas[j]
		r = R(theta)
		plane_pt = np.array([r*np.sin(beta),0,r*(1.0 - np.cos(beta))])
		ct = np.cos(theta)
		st = np.sin(theta)
		space_pt = np.dot(np.array([[ct, -st, 0],[st, ct, 0],[0, 0, 1.0]]), plane_pt)
		# space_pt = plane_pt
		pt_ind = i*betas.size + j
		xarr[pt_ind] = space_pt[0]
		yarr[pt_ind] = space_pt[1]
		zarr[pt_ind] = space_pt[2]

# ax = plt.axes(projection='3d')
# ax.plot(xarr, yarr, zarr, '.')
# set_axes_equal(ax)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_trisurf(xarr, yarr, zarr, cmap=cm.jet, linewidth=0)
fig.colorbar(surf)

ax.xaxis.set_major_locator(MaxNLocator(5))
ax.yaxis.set_major_locator(MaxNLocator(6))
ax.zaxis.set_major_locator(MaxNLocator(5))
set_axes_equal(ax)
fig.tight_layout()

plt.show()
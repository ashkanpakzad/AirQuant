import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.interpolate import interpn

# Set up grid and array of values
x1 = np.arange(10)
x2 = np.arange(10)
arr = x1 + x2[:, np.newaxis] # broadcast to get 10x10 gradient

# Set up grid for plotting
X, Y = np.meshgrid(x1, x2)

# Plot the values as a surface plot to depict
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, arr, rstride=1, cstride=1, cmap=cm.jet,
                       linewidth=0, alpha=0.8)
fig.colorbar(surf, shrink=0.5, aspect=5)

## interpolation

interp_x = 3.5           # Only one value on the x1-axis
interp_y = np.arange(10) # A range of values on the x2-axis

# Note the following two lines that are used to set up the
# interpolation points as a 10x2 array!
interp_mesh = np.array(np.meshgrid(interp_x, interp_y))
interp_points = np.rollaxis(interp_mesh, 0, 3)
interp_points.reshape((10, 2))

# Perform the interpolation
interp_arr = interpn((x1, x2), arr, interp_points)

# Plot the result
ax.scatter(interp_x * np.ones(interp_y.shape), interp_y, interp_arr, s=20,
           c='k', depthshade=False)
plt.xlabel('x1')
plt.ylabel('x2')

plt.show()


### 3D example

x1 = np.arange(10)
x2 = np.arange(10)
x3 = np.arange(10)

# Set up grid for plotting
X, Y, Z = np.meshgrid(x1, x2, x3)

v = X*2+Y-Z

# Plot the values as a plot to depict
fig = plt.figure()
plt.imshow(v[1,:,:])

## suppose you want interpolation of plane normal to vector n = [1; 1; 1], that point origin lies in ..

origin = np.array(([4.25,4.5,5.0])) # point in plane
normal = np.array(([1, 1, 1]))
d = np.dot(-1*origin, normal)


perp_slice_size = 10 # even
precision = 0.25 # all values must be round
bounds = (perp_slice_size*precision)/2

lowerx1 = origin[0] - bounds
upperx1 = origin[0] + bounds

lowerx2 = origin[1] - bounds

upperx2 = origin[1] + bounds

interp_x1 = np.arange(lowerx1,upperx1,precision)
interp_x2 = np.arange(lowerx2,upperx2,precision)

XX, YY = np.meshgrid(interp_x1, interp_x2)
ZZ = (-normal[0]*XX - normal[1]*YY - d)/normal[2];

interp_mesh = np.array((XX,YY,ZZ))
interp_points = np.rollaxis(interp_mesh, 0, 3)

# Perform the interpolation
interp_v = interpn((x1, x2, x3), v, interp_points)

## final view
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X,Y,Z)
ax.scatter(XX,YY,ZZ)
# plotNonPenetrationMomentLineContactSlip.py

import numpy as np
from matplotlib import pyplot as plt

'''
parameters:
'''
# d = distance from last COP point
d_min = -0.1
d_max = 0.1

mu = 0.1 # friction coefficient
fz = 0.79907 # normal force evaluated at last COP point with current rhs, speed and other input params
fx_dist_to_com = 0.1 # this is the moment arm of the x axis force to the y moment. NEGATIVE if fx is positive.
vx = 1.3e-5 # slip speed in x direction. remains constant irrespective
vy = 2.47e-5 # slip speed in y direction at the last COP point.
omega = 0.0619413  # slip rotation speed in z direction
my0 = -0.00153 # y moment evaluated at last COP point

'''
plotting:
'''
# array of distance values
a = np.linspace(d_min, d_max, 1000)
b = fz*mu*vx/np.sqrt(vx**2 + (vy - omega*a)**2)
c = -b*fx_dist_to_com + fx_dist_to_com*fz*mu*vx/np.sqrt(vx**2 + (vy)**2)
d = c + a*fz + my0

plt.plot(a,d)
plt.grid()
plt.show()
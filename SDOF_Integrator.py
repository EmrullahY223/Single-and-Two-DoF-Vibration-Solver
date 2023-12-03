"""
This file is created as a newer version of the older Numerical Integration Algorithm, for a clean start
"""

import numpy as np
from numpy.linalg import inv
from matplotlib import pyplot as plt
#Variables
m = 2.0
k = 2.0
c = 1.0     # critical damping = 2 * SQRT(m*k)

F0 = 1.0
delta_t = 0.001
omega = 1.0
time = np.arange(0.0,60.0,delta_t)

#Initial Conditions

y = np.array([0,0])     #[velocity, displacement]

A = np.array([[m,0],[0,1]])
B = np.array([[c,k],[-1,0]])
F = np.array([0.0,0.0])

Y = []
force = []

#Time-Stepping Solution
for t in time:
    
    if t <= 15:
        F[0] = F0 * np.cos(omega*t)
    else:
        F[0] = 0.0
    
    y = y + delta_t*inv(A).dot(F - B.dot(y))
    Y.append(y[1])
    force.append(F[0])

    KE = 0.5 * m * y[0]**2
    PE = 0.5 * k * y[1]**2

    if t % 1 <= 0.01:

        print(f'Total Energy: {KE+PE}')

#Plot the Results

t = [i for i in time]

plt.plot(t,Y)
plt.plot(t,force)
plt.grid(True)
plt.legend(['Displacement','Force'], loc='lower right')

plt.show()

print(f'Critical Damping  : {np.sqrt(m*k)*2}')
print(f'Natural Frequency : {np.sqrt(k/m)}')





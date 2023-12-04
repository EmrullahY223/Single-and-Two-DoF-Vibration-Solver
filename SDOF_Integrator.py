"""
This file is created as a newer version of the older Numerical Integration Algorithm, for a clean start

Using State-Space Form

"""

import numpy as np
from numpy.linalg import inv
from matplotlib import pyplot as plt

def F(t):
    F = np.array([0.0,0.0])

    if t <= 15:
        F[0] = F0 * np.cos(omega*t)
    else:
        F[0] = 0.0
    
    return F

def G(y,t):
    return A_inv.dot(F(t) - B.dot(y))

def RK4_step(y,t,delta_t):

    K1 = G(y,t)
    K2 = G(y+(K1*delta_t/2), t+delta_t/2)
    K3 = G(y+(K2*delta_t/2), t+delta_t/2)
    K4 = G(y+(K3*delta_t), t+delta_t)

    return delta_t/6*(K1 + 2*K2 + 2*K3 + K4)

#Variables
m = 2.0
k = 2.0
c = 0.0     # critical damping = 2 * SQRT(m*k)

F0 = 0.0
delta_t = .1
omega = 1.0
time = np.arange(0.0,60.0,delta_t)

#Initial Conditions

y = np.array([0,1])     #[velocity, displacement]

A = np.array([[m,0],[0,1]])
B = np.array([[c,k],[-1,0]])
A_inv = inv(A)

Y = []
force = []

#Time-Stepping Solution
for t in time:
    
    y = y + RK4_step(y,t,delta_t)
    Y.append(y[1])
    force.append(F(t)[0])

    KE = 0.5 * m * y[0]**2
    PE = 0.5 * k * y[1]**2

    if t % 1 <= 0.01:

        print(f'Total Energy: {KE+PE}')

#Plot the Results

plt.plot(time,Y)
plt.plot(time,force)
plt.grid(True)
plt.legend(['Displacement','Force'], loc='lower right')

plt.show()

print(f'Critical Damping  : {np.sqrt(m*k)*2}')
print(f'Natural Frequency : {np.sqrt(k/m)}')





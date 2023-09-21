import numpy as np 
from numpy import diff
from numpy.linalg import inv
from matplotlib import pyplot as plt
from math import cos, sin, tan, pi
import sympy as smp

# defining state vector
def G(y,t): 
    theta_d, phi_d, theta, phi = y[0], y[1], y[2], y[3]

    theta_dd = phi_d**2 * cos(theta) * sin(theta) - g/l * sin(theta)
    phi_dd = -2.0 * theta_d * phi_d / tan(theta)

    return np.array([theta_dd, phi_dd, theta_d, phi_d])

# Setting up Runge-Kutta fourth order
def RK4(y, t, dt):
    k1 = G(y,t)
    k2 = G(y+0.5*k1*dt, t+0.5*dt)
    k3 = G(y+0.5*k2*dt, t+0.5*dt)
    k4 = G(y+k3*dt, t+dt)

    return dt * (k1 + 2*k2 + 2*k3 + k4)/6

# variables
m = 1.0
l = 1.0
g = 9.81

delta_t = 0.01 # time step
time = np.arange(0.0, 4.0, delta_t) # start time, duration, time step

# initial state
y = np.array([0, 2, pi/2 , 0])   # [velocity(theta, phi), displacement(theta, phi)]

Y1 = []
Y2 = []

# time-stepping solution
for t in time:
    y = y + RK4(y, t, delta_t) 

    Y1.append(y[2])
    Y2.append(y[3])

# plot the result
plt.figure(figsize=(7, 4))
plt.xlabel('Time in seconds')

plt.ylabel('Angle in rad')

plt.plot(time, Y1, color = 'black') # plots theta, angle with the normal
plt.plot(time, Y2, linestyle = 'dashed', color = 'black') # plots phi, azimuthal angle
plt.grid(True)
plt.legend(['θ', 'ϕ'], loc='lower right')
plt.show()

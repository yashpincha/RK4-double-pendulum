import numpy as np 
from numpy.linalg import inv
from matplotlib import pyplot as plt
from math import cos, sin, pi

def G(y,t): 
    a1d, a2d = y[0], y[1]
    a1, a2 = y[2], y[3]

    m11, m12 = (m1+m2)*l1, m2*l2*cos(a1-a2)
    m21, m22 = l1*cos(a1-a2), l2
    m = np.array([[m11, m12],[m21, m22]])

    f1 = -m2*l2*a2d*a2d*sin(a1-a2) - (m1+m2)*g*sin(a1)
    f2 = l1*a1d*a1d*sin(a1-a2) - g*sin(a2)
    f = np.array([f1, f2])
    accel = inv(m).dot(f)

    return np.array([accel[0], accel[1], a1d, a2d])

def RK4_step(y, t, dt):
    k1 = G(y,t)
    k2 = G(y+0.5*k1*dt, t+0.5*dt)
    k3 = G(y+0.5*k2*dt, t+0.5*dt)
    k4 = G(y+k3*dt, t+dt)

    return dt * (k1 + 2*k2 + 2*k3 + k4) /6

# variables
m1, m2 = 4.0, 1.0
l1, l2 = 2.0, 1.0
g = 9.81

delta_t = 0.01 # increment
time = np.arange(0, 30.0, delta_t)

# initial state
y = np.array([0,0,0.8,0])   # [velocity, displacement]

Y1 = []
Y2 = []

# time-stepping solution
for t in time:
    y = y + RK4_step(y, t, delta_t) 

    Y1.append(y[2])
    Y2.append(y[3])

# plot the result
plt.figure(figsize=(14, 6))
plt.plot(time,Y1, color = 'black')
plt.plot(time,Y2, linestyle='dashed', color = 'black')
ax = plt.gca()  
ax.set_ylim([min(Y2) - 0.5, max(Y2)+0.5])

plt.xlabel('Time in seconds', fontsize=15)
plt.ylabel('Angle in radians', fontsize=15)

plt.grid(True)
plt.legend(['$\\theta_1$', '$\\theta_2$'], loc='lower right', prop={'size': 15})
plt.show()

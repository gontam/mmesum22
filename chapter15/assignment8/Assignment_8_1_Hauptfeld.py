import numpy as np
import matplotlib.pyplot as plt
from math import ceil

def RK4_step(t,y,dt,f):
    k1 = dt * f(t, y)
    k2 = dt * f(t + dt/2, y + k1/2)
    k3 = dt * f(t + dt/2, y + k2/2)
    k4 = dt * f(t + dt, y + k3)
    # t(+1), y(+1)
    return t+dt, y + k1/6 + k2/3 + k3/3 + k4/6

def RK4integrator(y0,f,N,t0,dt):
    # Starting values
    y = y0
    t = t0
    # Arrays
    ta = np.zeros(N+1)
    ya = np.zeros(N+1)
    ta[0] = t
    ya[0] = y
    # Run through N iterations of RK4
    for i in range(1, N+1):
        t, y = RK4_step(t, y, dt, f)
        ta[i] = t
        ya[i] = y
    # Return the fruits of our labour
    return ta, ya

def test_integration():
    # Function
    def derive_a(t, y):
        return ((t - y) / 2)
    def derive_b(t, y):
        return -2*y
    # Parameters
    t_from = 0
    t_to = 10
    dt = .1
    # Go through all the steps
    t, y = RK4integrator(1, derive_a, ceil(t_to/dt), t_from, dt)
    plt.plot(t, y)
    plt.show()

if __name__ == "__main__":
    test_integration()

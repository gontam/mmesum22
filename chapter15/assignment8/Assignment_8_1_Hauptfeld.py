import numpy as np
import matplotlib.pyplot as plt
from math import ceil


# a) General IVP solving using RK4

def RK4_step(t,y,dt,f):
    d = len(y)

    k1 = dt * np.array(f(t,        [y[i] for i in range(d)]))
    k2 = dt * np.array(f(t + dt/2, [y[i] + k1[i]/2 for i in range(d)]))
    k3 = dt * np.array(f(t + dt/2, [y[i] + k2[i]/2 for i in range(d)]))
    k4 = dt * np.array(f(t + dt,   [y[i] + k3[i] for i in range(d)]))

    new_y = [y[i] + k1[i]/6 + k2[i]/3 + k3[i]/3 + k4[i]/6 for i in range(d)]

    # t(+1), y(+1)
    return t+dt, new_y

def RK4integrator(y0,f,N,t0,dt):
    # Starting values
    y = y0
    t = t0
    # Arrays
    ta = [t]
    ya = [y]
    # Run through N iterations of RK4
    for i in range(1, N+1):
        t, y = RK4_step(t, y, dt, f)
        ta.append(t)
        ya.append(y)
    # Return the fruits of our labour
    return ta, ya

# b) Specific IVP
# y'' = -4y-y'
# y(0) = 1; y'(0) = -1
# 1st order ODE system: (y1') = (y2)
#                       (y2')   (-4y1 - y2)
# 
# Initial Conditions:   (y1(0)) = (1)
#                       (y2(0))   (-1)

def test_integration():
    # Function
    #def derive_a(t, y):
    #    return ((t - y) / 2)
    #def derive_b(t, y):
    #    return -2*y
    def derive(t, y):
        return [
            y[1],
            -4*y[0] - y[1]
        ]
    # Parameters
    t_from = 0
    t_to = 10
    dt = .1
    # Go through all the steps
    t, y = RK4integrator([1, -1], derive, ceil(t_to/dt), t_from, dt)
    print(y)
    plt.plot(t, y)
    plt.show()

if __name__ == "__main__":
    test_integration()

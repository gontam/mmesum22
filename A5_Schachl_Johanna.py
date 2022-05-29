# Assignment 5 - Modeling of infectious diseases
# Author: Johanna Schachl
# Date: 15.05.2022
# Last change:

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population of germany in 2009, N.
N = 81900000

# Initial number of infected and recovered individuals, I0 and R0.
I0 = 1
R0 = 0

# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0

# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta = 0.3
beta_Int = 0.28
gamma = 1. / 7

# A grid of time points (in days)
t = np.linspace(0, 270, 270)


# The SIR model differential equations.
def deriv(y, t, N, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I
    dRdt = gamma * I
    return dSdt, dIdt, dRdt


# Initial conditions vectors
y0 = S0, I0, R0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma))
S, I, R = ret.T

# Integrate the SIR equations with intervention over the time grid, t.
ret_Int = odeint(deriv, y0, t, args=(N, beta_Int, gamma))
S_Int, I_Int, R_Int = ret_Int.T

# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S / 1000000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I / 1000000, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R / 1000000, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
ax.set_title('SIR Model: A/H1N1 pandemic, Germany 2009')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number (in Million)')
ax.set_ylim(0, 90)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()

# Plot the data with intervention on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S_Int / 1000000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I_Int / 1000000, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R_Int / 1000000, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
ax.set_title('SIR Model: A/H1N1 pandemic, Germany 2009 with intervention')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number (in Million)')
ax.set_ylim(0, 90)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()


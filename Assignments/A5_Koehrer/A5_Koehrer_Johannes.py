import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population in Europe in 2000
N = 726352883
# Initial number of infected and recovered individuals, I0 and R0.
# The values assume that the infection started with one patient and that
# the infection never occurred before in Europe and has therefore no recovered patients
I0, R0 = 1, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
V0 = 0
# Contact rate, beta (in 1/days).
# with beta = 0.5, around 50% of the population will be infected on the peak
beta = 0.3
# Mean recovery rate gamma between 7 and 14 days (in 1/days)
gamma = 1./10.5
# Vaccination rate
delta = 0.001
# A grid of time points (in days)
t = np.linspace(0, 200, 200)


# The SIR model differential equations.
def deriv(y, t, N, beta, gamma, delta):
    S, I, R, V = y
    dSdt = -beta * S * I / N - delta * N
    dIdt = beta * S * I / N - gamma * I
    dRdt = gamma * I + delta * N
    dVdt = delta * N
    return dSdt, dIdt, dRdt, dVdt


# Initial conditions vector
y0 = S0, I0, R0, V0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma, delta))
S, I, R, V = ret.T

# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S / 1000000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I / 1000000, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R / 1000000, 'g', alpha=0.5, lw=2, label='Recovered with immunity + Vaccinated')
ax.plot(t, V / 1000000, 'y', alpha=0.5, lw=2, label='Vaccinated')
ax.set_title('SIR Model: Diphtheria')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number (10^6)')
ax.set_ylim(0, 800)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(visible=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()

#Assignment 5 - Modeling of infectious diseases - HEIN Raphael
# Detail description of the simulation -> see A5_Hein_Raphael.pdf
# Github: https://github.com/gontam/mmesum22/blob/Upload_Hein/Assignment%205/main.py


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

N = 8000000 # Total population of Austria
I0, R0 = 4, 474 # infected and recovered individuals, I0 and R0.
S0 = N - I0 - R0 # S0 susceptible to infection initially.
beta, gamma = 10, 1/30 # Contact rate = beta / recovery = gamma -> needs to be adjusted when applying medical interventions
t = np.linspace(0,360,360) # Grid for plot diagram

# --------------------------------- SIR model differential equations ---------------------------------
def deriv(y, t, N, beta, gamma):
    S, I, R = y
    Sdt = -beta * S * I / N
    Idt = (beta * S * I / N - gamma * I) * 0.05 # 5 % adults fall ill after an infection
    Rdt = gamma * I
    return Sdt, Idt, Rdt

y0 = S0, I0, R0

# --------------------------------- Integrate the SIR equations ---------------------------------
ret = odeint(deriv, y0, t, args=(N, beta, gamma))
S, I, R = ret.T

# --------------------------------- Visuall output / plot diagram ---------------------------------
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S / 10000, 'b', alpha=1, lw=2.5, label='Susceptible') # per 10.000 people for visual reasons
ax.plot(t, I / 10000 , 'r', alpha=1, lw=2.5, label='Infected')   # per 10.000 people for visual reasons, needs to be adjusted e.g. medical interventions
ax.plot(t, R / 10000, 'g', alpha=1, lw=2.5, label='Recovered')   # per 10.000 people for visual reasons
ax.set_xlabel('Time /days')
ax.set_ylabel('Number per 10.0000')
ax.set_ylim(0, 850)                                            # Limiting X axis to 8.500.000 people
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(visible=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(1)
for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)
plt.show()

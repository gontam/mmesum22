import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
'''
------------------- RESEARCH -------------------------
What is R0? How is it different to Re? 
Answer: 
    R0 is dependent on disease
    Re is dependent on population (effective)

Disease: Malaria
Pathogen: Plasmodium Vivax (Parasite)

Case Study Population: 
Malaria is an endemic disease to sub-saharan and tropical regions, Ghana is no exception.

Accra, Ghana -> Prevalence of Malaria ranges from 5% to 31% in the nation. For Accra, a value of 25% is used

Considering the incubation time of the parasite 8 Days
Typical recovery time from Malaria is 14 Days

A time frame of 100 Days is selected, to allow for stabilisation, with following compartments

S - Susceptible                     Every member of population
I - Infected                        Malaria (8 Days)
R - Recovered w. Immunity           After Malaria (14 Days)

---Outbreak Factors---

Humidity, Malaria-Vector prevalence (Fem. Anophales Mosquitoe), standing water, population density
   
Malaria is a treatable disease if diagnosed early enough, but left alone is deadly especially in children.
Accessibility to medication, vaccination, safe shelter with working plumbing and nets are all measures to 
control an outbreak and protect people. With time, some indivduals develop partial immunity which may be lost.

---Modelling---

To Model an outbreak with the above information;

A duration of 100 Days and an initial infection rate beta of 0.25 are used.
Once more than 20% of the population has Malaria, there is a higher prevalence of the parasite in human bodies
and hence in the vector animal, the infection rate beta increases to 0.4
To model acquired immunity, 10% of individuals recover with immunity, whilst 90% are again susceptible

To model the described interventions in a different plot, once 8% of the population is infected, 
the contact rate is reduced to 0.15, by the assumed effectiveness of government programs.    

---Analysis---

The produced plots show that during an outbreak, a large population is infected and may die!
Only few individuals recover with immunity, other individuals are susceptible again or may 
soon get it again without any intervention. With the intervention, part of of the population is 
still infected with malaria, but it depicts a situation where majority of the population remains
healthy or recovers with immunity over the course of the 'outbreak' and the affected communities
are provided with medications over the next months. 

Ghana has been able to reduce the prevalance of Malaria to 5% in Accra in 2021.

These models do not make considerations for other parameters and are based on informed estimates.
'''
# The SIR model differential equations for Outbreak.
def deriv1(y, t, N, beta, gamma):
    S, I, R = y
    if (I > 200):
        beta = 0.4
    else:
        beta = 0.25
    dSdt = -beta * S * I / N + (9 * gamma * I)
    dIdt = beta * S * I / N - gamma * I - (9 * gamma * I)
    dRdt = (gamma * I) 
    return dSdt, dIdt, dRdt

# The SIR model differential equations for Intervention.
def deriv2(y, t, N, beta, gamma):
    S, I, R = y
    if(I > 80):
        beta = 0.15
    dSdt = -beta * S * I / N + (9 * gamma * I)
    dIdt = beta * S * I / N - gamma * I- (9 * gamma * I)
    dRdt = (gamma * I) 
    return dSdt, dIdt, dRdt

# Total population, N.
N = 1000
I0, R0 = 1, 0
S0 = N - I0 - R0

# A grid of time points (in days)
t = np.linspace(0, 100, 100)
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 0.25, 1/14 * 0.1

# Initial conditions vector
y0 = S0, I0, R0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv1, y0, t, args=(N, beta, gamma))
S, I, R = ret.T

# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S / 1000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I / 1000, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R / 1000, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
ax.set_title('SIR Model: Malaria Outbreak')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number (1000s)')
ax.set_ylim(0, 1.2)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()


# Initial conditions 
y0 = S0, I0, R0
beta, gamma = 0.25, 1/14 * 0.1

ret = odeint(deriv2, y0, t, args=(N, beta, gamma))
S, I, R = ret.T

# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S / 1000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I / 1000, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R / 1000, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
ax.set_title('SIR Model: Prevention of Outbreak')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number (1000s)')
ax.set_ylim(0, 1.2)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()

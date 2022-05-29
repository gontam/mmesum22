#Assignment 5

#Research: choose and gather information about an infection or a virus
#    - gather information about a population where an outbreak has occurred or could happen
#            – identify a reliable source for health information for that region
#            - choose a time frame (in days) for the model and/or the forecast
#            - depending on the selected infectious disease, what compartments are there in the population?
#    - define the factors that influence an outbreak either by using credible sources or by providing them by yourself with justification
#            - primarily: susceptibility, infection rate, immunity
#    - define the interventions that could have been/were undertaken to control the outbreak
#        - how did you define the intervention factors in the equations? what are the values based on?
#    - determine which variables are available and can be used to model this event
#    - What is R0? How is it different to Re? Note the answer down.

# The chosen disease is Covid-19
# I gathered information on how this is modelled, but did not choose a specific population, but rather played with the
# capacities of the model
# and the ranges of numbers in which the parameters alpha, beta, gamma and epsilon usually occur
# As this is covid, the SIER (and including Deaths) model was implemented based on the SIR model of Prof. Gonta

#as a source this website was used, where the specific steps are explained
# https://towardsdatascience.com/forecasting-the-covid-19-trend-using-the-seir-model-90979abb9e64
# but with Python instead of Java

# as an intervention a vacination was chosen
# the difference between R0 and Re is that Re ist the effective reproduction number, meaning how many people in the
# population can be infected by an individual at any specific time -> how many people are immuninized?


#example
# small village with 2000 people, after first infections people are not allowed to leave the village

#Tried to change the formulas of https://mdpi-res.com/d_attachment/mathematics/mathematics-09-00636/article_deploy/mathematics-09-00636.pdf?version=1615969616
# accordnigly (without quarantine, newborns, normal deaths, etc. -> therefore alternative version of the shown equations in the paper)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population, N.
N = 2000
# Initial number of infected and recovered individuals, I0 and R0. Deaths D0 and Exposed E0
I0, R0 = 10, 0
D0 = 0
E0 = 30
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0 - D0 - E0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 1, 0.25 #beta typical values: 0.05 to 1.0; gamma: 1/9 to 1/4
eps = 1/5 #average 1/5.2 -> rounded to 1/5
alpha = 0.0001 #Fatality rate: between 0.000005 and 0.0001
# A grid of time points (in days)
t = np.linspace(0, 100)


# The SEIR model differential equations.
def deriv(y, t, N, beta, gamma, eps, alpha):
    S, I, R, E, D = y
    dSdt = -(beta * S * I) / N
    dEdt = (beta * S *I) /N - eps * E
    dIdt = eps * E - (gamma + alpha) * I
    dRdt = gamma * I
    dDdt = alpha * I
    return dSdt, dIdt, dRdt, dEdt, dDdt


# Initial conditions vector
y0 = S0, I0, R0, E0, D0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma, eps, alpha))
S, I, R, E, D= ret.T

# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S / 100, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I / 100, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(t, E / 100, 'y', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, R / 100, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
ax.plot(t, D / 100, 'm', alpha=0.5, lw=2, label='Dead')
ax.set_title('SEIR Model: Covid in a small village with lockdown')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number (100s)')
ax.set_ylim(0, 21)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(visible=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()

# with intervention of vaccination

# Total population, N.
N = 2000
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 10, 0
D0 = 0
E0 = 30
V0 = 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0 - D0 - E0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 1/N, 0.25 #beta typical values: 0.05 to 1.0; gamma: 1/9 to 1/4
eps = 1/5 # Time until infected; average 1/5.2 -> rounded to 1/5
alpha = 0.0001 #Fatality rate: between 0.000005 and 0.0001
vac = 0.2 #Vaccination rate
sigma = 0.8#efficiency of vaccination
# A grid of time points (in days)
t = np.linspace(0, 100)
#V=10

#change alpha to vac

# The SEIR model differential equations.
def deriv(y, t, N, beta, gamma, eps, alpha, vac, sigma):
    S, I, R, E, D, V = y
    dSdt = -(beta * S * I) - vac * S
    dEdt = (beta * S *I) - eps * E + beta * sigma * V * I
    dIdt = eps * E - (gamma + alpha) * I
    dRdt = gamma * I
    dDdt = alpha * I
    dVdt = vac * S - sigma * beta * V * I
    return dSdt, dIdt, dRdt, dEdt, dDdt, dVdt


# Initial conditions vector
y0 = S0, I0, R0, E0, D0, V0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma, eps, alpha, vac, sigma))
S, I, R, E, D, V= ret.T

# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S / 100, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I / 100, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(t, E / 100, 'y', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, R / 100, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
ax.plot(t, D / 100, 'm', alpha=0.5, lw=2, label='Dead')
ax.plot(t, V / 100, 'c', alpha=0.5, lw=2, label='Vaccinated')
ax.set_title('SEIR Model: Covid in a small village with lockdown (with vaccinations)')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number (100s)')
ax.set_ylim(0, 21)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(visible=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()
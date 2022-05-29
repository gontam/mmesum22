'''
---------------------- 1: Reasearch ----------------------------
Info-------Source: https://www.who.int/health-topics/coronavirus#tab=tab_1
Disease: Coronavirus disease (COVID-19)
Pathogen: SARS-CoV-2 virus

Infectious disease cause by the virus SARS_CoV-2

 - gather information about a population where an outbreak has occurred or could happen
https://www.worldometers.info/coronavirus/country/austria/

 Worldwide
 Place choosen: Vienna, Austria




Time frame: Beginning of Pandemic Duration: 150 Days



Incubation Time: 5 days
Typical recovery time from mild symptoms is ca. 14 days

Compartments:

S - Susceptible         All members of population
I - Infected            Covid-19:  5 Days
R - Recovered           After Covid-19: 14 Days


Outbreak Factors:
    define the factors that influence an outbreak either by using credible sources or by providing them by yourself
    with justification:


        susceptibility: The more people susceptible, the more people can get it and spread it
        infection rate: If the infection rate is higher, more people are prone to an infection (number of infection/number of people at risk of infection)
        immunity: The higher the immunity (natural/vaccine), the less infections will occur as successful transmission and possible host are less uncommon

        population density: The more people live in one place, the more people meet each other every day and therefore a higher chance of infection occurs. Also, in densely populateted and poor regions
        the average access to clean supplies (example Water) and medical care is low and people live more cramped together-->even higher risk of disease outbrake



    - define the interventions that could have been/were undertaken to control the outbreak
        - how did you define the intervention factors in the equations? what are the values based on?
        Intervention choosen is vaccination, the factors are based on Austrian Covid-19 statistics and https://towardsdatascience.com/forecasting-the-covid-19-trend-using-the-seir-model-90979abb9e64

    - determine which variables are available and can be used to model this event
        Population Count
        Exposed count
        Infectious count
        Recovered with immunity count
        Death count
        Vaccination rate and effectivity


What is R0? How is it different to Re? Source https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7751056/

    R0: is the basic reproduction number, also known as basic reproduction ratio or rate which is an epidemiological metric used to measure the transmissibility of infectious agents
    Re: is the number of people in a population who can be infected by an individual at any specific time

To have a more precise model and as a lot of data about Covid-19 is found about Austria, the SEIR model was choosen.
The whole program is based on the example of Ms. Gontam and extended to fit the SEIR model and homework task
Sources to help build the programm with the SEIR model:
Explanation and formulas: https://towardsdatascience.com/forecasting-the-covid-19-trend-using-the-seir-model-90979abb9e64
Formulas with vac:        https://mdpi-res.com/d_attachment/mathematics/mathematics-09-00636/article_deploy/mathematics-09-00636.pdf?version=1615969616
                     The formulas are changed as quarantined, new born and trough other causes deceased people are not considered in this code, therefore a more basic model


---------------------- 3: Analysis ----------------------------
The choosen place for the model was Vienna, as the IO Value and E0 Value are not known, an estimated value was choosen
People were not allowed to leave the city in this scenario, therefore people from outside Vienna, people leaving Vienna and normal deaths, new borns and etc. are not taken into account to keep the scenario simple

The vaccination rate simulation is biased due to the fact the numbers are lower in the beginning than in real life, therefore just a showcase of effectiviness of the intervention if it would have happend in the beginning

To showcase the vaccination in an example as it happened in early 2021 I used the data from May 2021 to showcase the beginning of the vaccination in a 3rd plot

Comparison:

With a vaccination, the population susceptible gets lower accordingly to the number of people vaccinated and infacted and therefore much faster than in the scenario without the vaccination
Also the general number of infected, exposed, dead and recovered people is lower in the vaccination scenario
Therefore, the vaccination was successful, even tho with a better functioning vaccine it would be even better
'''

#---------------------- 2: Program -----------------------------


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


#------- No vaccination---------

# Total population, N.
N = 1900000

# Initial number of infected and recovered individuals, I0 and R0. Deaths D0 and Exposed E0 (infected put not contagious yet)
I0, R0 = 50, 0
D0 = 0
E0 =100
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0 - E0 - D0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 0.9, 0.25          #beta typical ranges: 0.05 to 1.0; gamma: from 1/9 to 1/4
eps = 1/5.2                         #average 1/5.2
alpha = 0.0001                      #fatality rate: between 0.000005 and 0.0001
# A grid of time points (in days)
t = np.linspace(0, 150)


# The SEIR model differential equations.
def deriv(y, t, N, beta, gamma, eps, alpha):
    S, I, R, E, D = y
    dSdt = -(beta * S * I) / N
    dIdt = eps * E - (alpha + gamma) * I
    dRdt = gamma * I
    dEdt = (beta * S * I) / N - eps * E
    dDdt = alpha * I
    return dSdt, dIdt, dRdt, dEdt, dDdt


# Initial conditions vector
y0 = S0, I0, R0, E0, D0
# Integrate the SEIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma, eps, alpha))
S, I, R, E, D= ret.T

# Plot the data on five separate curves
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S / 1000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I / 1000, 'g', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R / 1000, 'r', alpha=0.5, lw=2, label='Recovered with immunity')
ax.plot(t, E / 1000, 'm', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, D / 1000, 'aqua', alpha=0.5, lw=2, label='Dead')
ax.set_title('SEIR Model: Vienna with Lockdown')
ax.set_xlabel('Time / Days')
ax.set_ylabel('Population (1000s)')
ax.set_ylim(0, 2200)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(visible=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()

# -----------------With intervention of vaccination----------------------

# Total population, N.
N = 1900000

# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 50, 0
D0 = 0
E0 = 150
V0 = 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0 - E0 - D0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 0.9/N, 0.25               #beta typical ranges: 0.05 to 1.0; gamma: from 1/9 to 1/4
eps = 1/5.2             #Time until infected; average 1/5.2
alpha = 0.0001          #fatality rate: between 0.000005 and 0.0001
vac = 0.15               #vaccination rate
sigma = 0.8             #efficiency of vaccination-->choosen 80% as average
# A grid of time points (in days)
t = np.linspace(0, 150)
#V=10

#change alpha to vac

# The SEIR model differential equations.
def deriv(y, t, N, beta, gamma, eps, alpha, vac, sigma):
    S, I, R, E, D, V = y

    dEdt = (beta * S *I) - eps * E + beta * sigma * V * I
    dIdt = eps * E - (alpha + gamma) * I
    dRdt = gamma * I
    dSdt = -(beta * S * I) - vac * S
    dDdt = alpha * I
    dVdt = vac * S - sigma * beta * V * I
    return dSdt, dIdt, dRdt, dEdt, dDdt, dVdt


# Initial conditions vector
y0 = S0, I0, R0, E0, D0, V0
# Integrate the SEIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma, eps, alpha, vac, sigma))
S, I, R, E, D, V = ret.T

# Plot the data on six separate curves
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S / 1000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I / 1000, 'g', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R / 1000, 'r', alpha=0.5, lw=2, label='Recovered with immunity')
ax.plot(t, E / 1000, 'm', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, D / 1000, 'aqua', alpha=0.5, lw=2, label='Dead')
ax.plot(t, V / 1000, 'gold', alpha=0.5, lw=2, label='Vaccinated')
ax.set_title('SEIR Model: Vienna with lockdown and vaccination)')
ax.set_xlabel('Time / Days')
ax.set_ylabel('Population (1000s)')
ax.set_ylim(0, 2200)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(visible=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()

#--------With vaccination from the middle of the pandemic to be closer to reality----------



# Total population, N.
N = 1900000

# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 15600, 550000
D0 = 1000
E0 = 50000
V0 = 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0 - E0 - D0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 0.9/N, 0.25               #beta typical ranges: 0.05 to 1.0; gamma: from 1/9 to 1/4
eps = 1/5.2             #Time until infected; average 1/5.2
alpha = 0.0001          #fatality rate: between 0.000005 and 0.0001
vac = 0.15               #vaccination rate
sigma = 0.8             #efficiency of vaccination-->choosen 80% as average
# A grid of time points (in days)
t = np.linspace(0, 150)
#V=10



# The SEIR model differential equations.
def deriv(y, t, N, beta, gamma, eps, alpha, vac, sigma):
    S, I, R, E, D, V = y

    dEdt = (beta * S *I) - eps * E + beta * sigma * V * I
    dIdt = eps * E - (alpha + gamma) * I
    dRdt = gamma * I
    dSdt = -(beta * S * I) - vac * S
    dDdt = alpha * I
    dVdt = vac * S - sigma * beta * V * I
    return dSdt, dIdt, dRdt, dEdt, dDdt, dVdt


# Initial conditions vector
y0 = S0, I0, R0, E0, D0, V0
# Integrate the SEIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma, eps, alpha, vac, sigma))
S, I, R, E, D, V= ret.T

# Plot the data on six separate curves
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S / 1000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I / 1000, 'g', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R / 1000, 'r', alpha=0.5, lw=2, label='Recovered with immunity')
ax.plot(t, E / 1000, 'm', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, D / 1000, 'aqua', alpha=0.5, lw=2, label='Dead')
ax.plot(t, V / 1000, 'gold', alpha=0.5, lw=2, label='Vaccinated')
ax.set_title('SEIR Model: Vienna with lockdown and vaccination in middle of pandemic')
ax.set_xlabel('Time / Days')
ax.set_ylabel('Population (1000s)')
ax.set_ylim(0, 2200)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(visible=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()


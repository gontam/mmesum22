import math

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# following papers/links were used to answer the research questions:
# https://www.who.int/news-room/fact-sheets/detail/ebola-virus-disease
# https://www.nhs.uk/conditions/ebola/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8039563/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7751056/
#
# The model is mostly based on the paper written by M. Juga et al. The data for the model is based on outbreaks in
# North Kivu und South Kivu within May 2019 and June 2020 which leads to approximately 395 days. The model uses the
# SIRDH model, meaning there are 6 compartments within the population: Susceptible; Infectious; Recovered; Deceased;
# Hospitalized (The paper by M. Juga et al. suggests a 6 compartments system using as well a pathogens in the
# environment compartment. However, for this simple model, the SIRD model is enough.) There are multiple factors that
# influence an outbreak, some of them are listed below (taken from the M. Juga et al. paper, for the full table
# please check Table 3 of the paper):
# Recruitment rate: 200
# Effective human to human contact rate	0.67
# Effective human to pathogen contact rate	0.64
# Disease related death of the infected	0.6
# Disease related death of the hospitalized	0.2
#
# Possible interventions include hospitalisation, quarantine (meaning decrease in effective human to human contact rate)
# The equations are taken from the M. Juga et al. paper. Hospitalisation influences the following two equations:
# dH/dt=σ2I−ν2H,
# dD/dt=σ3I+γ2H−ρD,dDdt=σ3I+γ2H−ρD,
# These equations can be seen in a visually better form in chapter 2 of said paper.
#
# R0 is the basic reproduction number, it is an estimation of new infections. Re is the effective reproduction
# number, meaning it is number of humans who are able to be infected, therefore it changes as the population gets
# immunized e.g. due to past infections
#
# The 2 scenarios are based on the research numbers for South and North Kivu and are used for Austria.
# Some of the values are adapted to better reflect our health system (meaning less deaths and higher hospitalisation)
#
#

def plot(data, h_included):
    S = data[:, 0]
    I = data[:, 1]
    R = data[:, 2]
    D = data[:, 3]
    plt.figure()
    plt.plot(t, S, "m", label="S(t): Susceptibles")
    plt.plot(t, I, 'r', label="I(t): Infected")
    plt.plot(t, R, 'g', label="R(t): Recovered")
    plt.plot(t, D, 'k', label="D(t): Deceased")
    if h_included:
        H = data[:, 4]
        plt.plot(t, H, 'b', label="R(t): Hospitalized")
    plt.xlim([0, 180])
    plt.ylim([0, 9.5e6])
    plt.ylabel("Population [millions]")
    plt.xlabel("Time [days]")
    plt.title('Ebola model in Austria within 1 year')
    plt.legend()
    plt.grid(True)
    plt.show()


# function to calculate the parameters dSdt, dIdt and dRdt, dHdt and dDdt for the ODE
def SIRDHfunction(y,t):
    S, I, R, D, H = y
    dSdt = pi-(lamda+u)*S
    dIdt = lamda*S-v1*I
    dRdt = o1*I + y1*H
    dHdt = o2*I - v2*H
    dDdt = o3*I + y2*H
    return dSdt,dIdt,dRdt, dDdt, dHdt

# Adapted function without hospitalisation
def SIRDfunction(y,t):
    S, I, R, D = y
    dSdt = pi-(lamda+u)*S
    dIdt = lamda*S-v1alt*I
    dRdt = o1*I
    dDdt = o3*I
    return dSdt,dIdt,dRdt, dDdt

# defined variables
N = 8917000  # Population in Austria
I0 = 80
R0 = 0
D0 = 0
H0 = 3
S0 = N - I0 - R0 - H0 - D0  # susceptible individuals in Austria
y0 = [S0, I0, R0, D0, H0]
y02 = [S0, I0, R0, D0]
# y1 = Rate of recovery of the hospitalised
# y2 = disease related death of the hospitalized
y1 = 0.7
y2 = 0.1
# o1 = recovery rate
# o2 hospitalised rate
# 03 dying rate
# v1 = o1 + o2 + o3 + u
o1 = 0.33
o2 = 0.149
o3 = 0.5
# u = natural death rate
u = 0.000296
v1 = o1 + o2 + o3 + u
v1alt = o1 + o3 + u
v2 = u + y1 +y2
p = 0.009
pi = 200
lamda = 0.67
t = np.linspace(0, 180, 1000)
y = odeint(SIRDHfunction, y0, t)
y2 = odeint(SIRDfunction, y02, t)
plot(y, True)
plot(y2, False)



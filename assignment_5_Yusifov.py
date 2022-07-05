#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Made by:
#   Ukhnalyov Andrey
#   Yusifov Tamerlan

# Improting libraries
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

'''
------------------------------------RESEARCH------------------------------------
Parameters R0 and Re:
    R0 depends on disease
    Re depends on population


Disease: Tuberculosis
Location: Austria


GENERAL INFORMATION ABOUT DESEASE:
(according to WHO https://www.who.int/health-topics/tuberculosis)

TB is caused by bacteria (Mycobacterium tuberculosis) and it most often affects
the lungs. TB is spread through the air when people with lung TB cough, sneeze
or spit. A person needs to inhale only a few germs to become infected.

Every year, 10 million people fall ill with tuberculosis (TB). Despite being
a preventable and curable disease, 1.5 million people die from TB each year –
making it the world’s top infectious killer.

TB is the leading cause of death of people with HIV and also a major
contributor to antimicrobial resistance. Most of the people who fall ill
with TB live in low- and middle-income countries, but TB is present all
over the world. About half of all people with TB can be found in 8 countries:
Bangladesh, China, India, Indonesia, Nigeria, Pakistan, Philippines and
South Africa. About one-quarter of the world’s population is estimated to be
infected by TB bacteria. Only 5-15% of these people will fall ill with active
TB disease. The rest have TB infection but are not ill and cannot transmit the
disease. Both TB infection and disease are curable using antibiotics.

Symptoms:
- Prolonged cough
- Chest pain
- Weakness or fatigue
- Weight loss
- Fever
- Night sweats


INFORMATION ABOUT POPULATION:

Source: WHO for 2017-2018

Total population: 8 772 865

Number of cases: 570

Infection rate (beta): 15/365
People with active TB can infect 5–15 other people through close contact
over the course of a year.

Recovery (gama): from 6 Month = c.a. 1/180

Time period = 365 days


INTERVENTIONS:
1. Masks -> reduce contact rate (beta)
2. Isolation of a patient -> reduce contact rate (beta)
3. Application of antibiotics -> reduce recovery time and contact rate
  (beta and gamma)

According to the paper:
"Surgical Face Masks Worn by Patients with Multidrug-Resistant Tuberculosis"
PMID: 22323300, the risk of transmission (beta) dicreases by 56%
if patients ware masks


ANALYSIS:
Based on the graphs provided, we can clearly see the following trend:
the rate of tuberculosis drastically deacreases when masks are on.
Although, it does not fully prevent from getting tuberculosis, however it does
minimaise chanses to almost zero. Such slow rates of contamination allow to
completely prevent the epidemy by applying
some additional intervention methods.

The values used for calculation are based on:
- WHO https://www.who.int/health-topics/tuberculosis
- "Surgical Face Masks Worn by Patients with Multidrug-Resistant Tuberculosis"
  PMID: 22323300

'''


# The SIR model differential equations for Outbreak
def deriv(y, t, N, beta, gamma):
    S, I, R = y

    dSdt = (-beta*S*I / N) + gamma*I
    dIdt = (beta*S*I / N) - gamma*I
    dRdt = gamma*I

    return dSdt, dIdt, dRdt


# Plotting function
def make_graph(Susceptible, Infected, Recovered, time):
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
    ax.plot(time, Susceptible / 10000,
            'b', alpha=0.5, lw=2, label='Susceptible')
    ax.plot(time, Infected / 10000,
            'r', alpha=0.5, lw=2, label='Infected')
    ax.plot(time, Recovered / 10000,
            'g', alpha=0.5, lw=2, label='Recovered with immunity')
    ax.set_title('SIR Model: Tuberculosis Outbreak')
    ax.set_xlabel('Time /days')
    ax.set_ylabel('Number (10000s)')
    ax.set_ylim(0, 880)
    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which='major', c='w', lw=2, ls='-')
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)
    plt.show()


if __name__ == '__main__':
    # Initial conditions
    N, I0, R0 = 8772865, 570, 0
    S0 = N - I0 - R0

    # Initial conditions vector
    y0 = S0, I0, R0

    # Time vector
    t = np.linspace(0, 365, 365)

    # Parameters without interventions
    beta, gamma = 15/365, 1/180

    # Integrate the SIR equations over t
    ret = odeint(deriv, y0, t, args=(N, beta, gamma))
    S, I, R = ret.T

    # Plot the data on three separate curves for S(t), I(t) and R(t)
    make_graph(S, I, R, t)

    # Beta with masks intervention
    beta_new = beta*(1-0.56)

    # Integrate the SIR equations over t with beta_new
    ret = odeint(deriv, y0, t, args=(N, beta_new, gamma))
    S, I, R = ret.T

    # Plot the data on three separate curves for S(t), I(t) and R(t)
    make_graph(S, I, R, t)

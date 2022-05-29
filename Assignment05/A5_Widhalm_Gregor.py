# author: Gregor Widhalm

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


# The SIR model differential equations.
def deriv(y, t, N, beta, gamma, v, mu):
    S, I, R = y
    # v ... vaccination rate
    # beta ... transmission rate
    # gamma ... recovery rate
    # mu ... birth/mortality rate
    dSdt = -beta * S * I / N - v*S - mu*S + mu*N
    dIdt = beta * S * I / N - gamma * I - mu*I
    dRdt = gamma * I + v*S - mu*R
    return dSdt, dIdt, dRdt

def plotSIR(S, I, R, title, t):
    # Plot the data on three separate curves for S(t), I(t) and R(t)
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
    ax.plot(t, S / 1000000, 'b', alpha=0.5, lw=2, label='Susceptible')
    ax.plot(t, I / 1000000, 'r', alpha=0.5, lw=2, label='Infected')
    ax.plot(t, R / 1000000, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
    ax.set_title(title)
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Number [Millions]')
    # ax.set_ylim(0, 1.2)
    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(which='major', c='w', lw=2, ls='-')
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)
    plt.show()

if __name__ == '__main__':

    # COIVD-19 situation starting at 01.02.2022.


    t = np.linspace(0, 100, 100)

    # actual numbers on February 1st 2022.
    # information gained from https://covid19-dashboard.ages.at/ download section
    # Total population, N. in 8,932,664
    N = 8932664
    # Initial number of infected and recovered individuals, I0 and R0. in 1000: I0 380,851; R0 1,552,182
    I0 = 380851
    R0 = 1552182
    # Everyone else, S0, is susceptible to infection initially.
    S0 = N - I0 - R0

    print('initial condition: ')
    print(' total population N: ', N, '\n', 'susceptible S0: ', S0, '\n', 'recovered R0 : ', R0, '\n', 'infected I0: ', I0, '\n')


    # Contact rate, beta; Italy as reference; according to Cooper et al.
    # Cooper, Ian, Argha Mondal, und Chris G. Antonopoulos. „A SIR Model Assumption for the Spread of COVID-19 in Different Communities“. Chaos, Solitons & Fractals 139 (Oktober 2020): 110057. https://doi.org/10.1016/j.chaos.2020.110057.
    beta = 0.18

    # recovery duration assumed to be in mean 17 days, according to Espinosa et al.
    # Espinosa, Pablo, Paulina Quirola‐Amores, und Enrique Teran. „Application of a Susceptible, Infectious, and/or Recovered (SIR) Model to the COVID-19 Pandemic in Ecuador“. Frontiers in Applied Mathematics and Statistics 6 (30. November 2020): 571544. https://doi.org/10.3389/fams.2020.571544.
    # Kaul, D. An overview of coronaviruses including the SARS-2 coronavirus—molecular biology, epidemiology and clinical implications. Curr Med Res Pract (2020). 10(2):54–64. doi:10.1016/j.cmrp.2020.04.001.
    # Rothan, HA, and Byrareddy, SN. The epidemiology and pathogenesis of coronavirus disease (COVID-19) outbreak. J Autoimmun (2020). 109, 102433. doi:10.1016/j.jaut.2020.102433.
    # mean recovery rate, gamma, (in 1/days)
    gamma = 1 / 17

    # vaccination rate v
    # 3rd vaccination on 01.02.2022: 23.586 (according to BMI, Gesundheitsministerium: https://info.gesundheitsministerium.at/opendata#COVID19_vaccination_doses_timeline)
    # vaccination rate on 01.02.2022:
    v = 23586/N
    print('vaccination rate v: ', v)

    # birth/death ration µ
    # 1/µ set to 76 years
    # according to Bjørnstad et al.
    # Bjørnstad, Ottar N., Katriona Shea, Martin Krzywinski, und Naomi Altman. „The SEIRS Model for Infectious Disease Dynamics“. Nature Methods 17, Nr. 6 (Juni 2020): 557–58. https://doi.org/10.1038/s41592-020-0856-2.
    mu = 1/(76*365)
    print('birth/mortality rate mu: ', mu)

    # Initial conditions vector
    y0 = S0, I0, R0
    # Integrate the SIR equations over the time grid, t.
    ret = odeint(deriv, y0, t, args=(N, beta, gamma, v, mu))
    S, I, R = ret.T

    np.max(I)

    plotSIR(S, I, R, 'SIR Model: COVID-19; Starting on 02nd of March 2022 w/o intervention', t)


    # vaccination boost --> vaccination rate increase to 1% --> assumed by Law et al.
    # Law, Kian Boon, Kalaiarasu M. Peariasamy, Hishamshah Mohd Ibrahim, und Noor Hisham Abdullah. „Modelling Infectious Diseases with Herd Immunity in a Randomly Mixed Population“. Scientific Reports 11, Nr. 1 (Dezember 2021): 20574. https://doi.org/10.1038/s41598-021-00013-2.
    v_intervention = 0.01
    ret_intervention = odeint(deriv, y0, t, args=(N, beta, gamma, v_intervention, mu))
    S_int, I_int, R_int = ret_intervention.T
    plotSIR(S_int, I_int, R_int, 'SIR Model: COVID-19; Starting on 02nd of March 2022 w incr. vacc. rate', t)


    print('maximum I w/o intervention: ', int(np.max(I)))
    print('maximum I with intervention: ', int(np.max(I_int)))
    print('delta maximum I: ', int(np.max(I))-int(np.max(I_int)))


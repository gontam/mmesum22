# Author: Oliver Jovanović
# Monte Carlo PI Implementation
# Based on: https://blogs.sas.com/content/iml/2016/03/14/monte-carlo-estimates-of-pi.html#:~:text=To%20compute%20Monte%20Carlo%20estimates,quarter%20circle%20of%20unit%20radius.
# Ansatz: average value method

# Libraries:
from scipy import integrate
import numpy as np
import math
import random

# Set seed for reproducibility
random.seed(1, 0.1)


def mc_pi(N):
    # Generate N times random numbers between 0 and 1.
    for x in range(N):
        a = random.uniform(0, 1)
        y = math.sqrt(1 - a**2)
        approx_pi = 4 * np.mean(y)
    return approx_pi


def mc_pi_stat(N, M):
    mc_pi(N)


#def mc_pi_plt():


# Results:
N = 10**5
print(mc_pi(N))

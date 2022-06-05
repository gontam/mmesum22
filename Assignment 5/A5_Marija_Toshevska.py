import numpy as np
import matplotlib.pyplot as plt


N = 59554657
I0 = 36.880
R0 = 24.882
#S0 - other people.
S0 = N - I0 - R0
#transmission rate
beta = 0.18
#recovery duration 14 days
gamma = 1 / 14
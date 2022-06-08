import numpy as np
import math as math
import matplotlib.pyplot as plt


# Calculates pi based on Monte-Carlo Algorithm using sample size N
def mc_pi(n):
    count = 0
    for i in range(0, n):
        x = np.random.uniform(0, 1)
        y = np.random.uniform(0, 1)
        if math.sqrt((x * x) + (y * y)) <= 1.0:
            count += 1
    # returning pi approximation
    return (count / n) * 4


# Calculates mean and standard deviation for a number of m approximations of pi
def mc_pi_stat(n, m):
    approx_pi = []
    for i in range(0, m):
        approx_pi.append(mc_pi(n))
    # returning mean (first value) and standard deviation (second value)
    return np.mean(approx_pi), np.std(approx_pi, ddof=1)


# Illustrates standard deviation and mean
def mc_pi_plt():
    mean = []
    std = []

    # calculating for 10^2 to 10^7
    for i in range(2, 7):
        n = 10 ** i
        u, o = mc_pi_stat(n, 10)
        mean.append(u)
        std.append(o)
    x = np.array([2, 3, 4, 5, 6])
    pi = []
    for i in range(2, 7):
        pi.append(np.pi)

    # Creating first plot showing our mean and the errorbar
    plt.figure(0)
    plt.plot(x, mean, 'b*')
    plt.errorbar(x, pi, yerr=std, ecolor="b", color="#FFA500")
    plt.xlabel("log(N)")
    plt.ylabel("μ")

    # We ignore the polyfit warning here because we do not need the third value that it would return
    k1, d1 = np.polyfit(x, np.log10(std), 1)
    # As requested, printing alpha
    print(f'Alpha: {k1}')

    # Creating the second plot for the standard deviation
    plt.figure(1)
    plt.plot(x, np.log10(std), 'x')
    plt.plot(x, k1 * x + d1, '-')
    plt.xlabel("log(N)")
    plt.ylabel("log(σ)")

    # Showing the plots
    plt.show()


mc_pi_plt()

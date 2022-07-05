#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Made by:
#   Ukhnalyov Andrey
#   Yusifov Tamerlan

# Improting libraries
import numpy as np
import math
import matplotlib.pyplot as plt


def mc_pi(n: int) -> float:
    '''
    Calculates an approximation of pi using the Monte-Carlo Algorithm
    using sample size n
    '''
    count = 0
    for i in range(n):
        x = np.random.uniform(0, 1)
        y = np.random.uniform(0, 1)
        if math.sqrt((x**2) + (y**2)) <= 1:
            count += 1
    # Return a pi approximation
    return 4*count/n


def mc_pi_stat(n: int, m: int) -> tuple[float, float]:
    '''
    Calculates mean and standard deviation for a number
    of m approximations of pi
    '''
    # print("Inside of stat function")
    pi_values = []
    for i in range(m):
        pi_values.append(mc_pi(n))
    # Returns a tuple (mean, standard deviation)
    # print("Out of stat function")
    return np.mean(pi_values), np.std(pi_values, ddof=1)


def mc_pi_plt(N: int, Base=10) -> None:
    """
    Make plots for log(N)-μ and log(N)-log(σ)
    """

    # print("Starting the plt proc")

    mean_array = []
    std_array = []

    for i in range(1, N+1):
        # print(f"Making the mean and std arrays {i}")
        n = Base**i
        mean, std = mc_pi_stat(n, Base)
        mean_array.append(mean)
        std_array.append(std)

    # print(mean_array)
    # print(std_array)
    pi = []
    for i in range(N):
        pi.append(np.pi)

    x = np.array([i for i in range(1, N+1)])

    # Make a plot log(N)-μ
    plt.figure(0)
    plt.plot(x, mean_array, 'b*')
    plt.errorbar(x, pi, yerr=std_array, ecolor="b", color="#FFA500")
    plt.xlabel("log(N)")
    plt.ylabel("μ")

    # Find coeffitients k*x + d
    k, d = np.polyfit(x, [math.log(value, Base) for value in std_array], 1)
    # Print alpha
    print(f'Alpha: {k}')

    # Make a plot log(N)-log(σ)
    plt.figure(1)
    plt.plot(x, [math.log(value, Base) for value in std_array], 'x')
    plt.plot(x, k*x + d, '-')
    plt.xlabel("log(N)")
    plt.ylabel("log(σ)")

    # Showing the plots
    plt.show()


if __name__ == '__main__':
    size = int(input("Enter the size of a sample: "))
    base = input("Enter the base: ")

    print(size, base)

    if base == '':
        mc_pi_plt(size)
    elif base != '':
        mc_pi_plt(size, Base=int(base))

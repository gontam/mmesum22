#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Made by:
#   Ukhnalyov Andrey
#   Yusifov Tamerlan

# Improting libraries
import numpy as np
import random as rnd
import matplotlib.pyplot as plt
from time import time


# Sorting function that takes a list of floats and returns a list of floats
# Applied recursively
def myquicksort(list_to_sort: list[float]) -> list[float]:
    if list_to_sort == []:
        # If list is empty, stop the execution
        return list_to_sort
    else:
        # If list is not empty, take the element 0
        list_copy = list_to_sort.copy()
        pivot = list_copy[0]
        list_smaller = []
        list_higher = []
        # Sort all other elements in two lists:
        # smaller or higher then the element 0
        for i in range(len(list_copy)):
            if list_copy[i] < pivot:
                list_smaller.append(list_copy[i])
            elif list_copy[i] > pivot:
                list_higher.append(list_copy[i])
        # Make concatination of lists
        return myquicksort(list_smaller) + [pivot] + myquicksort(list_higher)


# Based on the bioinf_sort.py
def check_my_sort() -> None:
    # general initialization
    N_pts = 5  # number of different lentghs for which sort time is measured
    Nt = 10  # number of repetitions to measure time
    results = np.zeros([2, N_pts, 2])

    # range of lengths for built in method
    Emin = 5  # log10 of minimal sample size
    Emax = 6  # log10 of maximal sample size
    pts = np.linspace(Emin, Emax, N_pts)  # contains exponents for list lengths
    for i in range(N_pts):
        N = int(10 ** (pts[i]))  # lentgh of list
        num_list = np.zeros(N)
        print(N)
        # inititalization for timing of list with length N
        for j in range(N):
            num_list[j] = rnd.uniform(0, 100)
        # timing for built in method
        t1 = time()
        for t in range(Nt):
            sorted_num_list = np.sort(num_list, kind='mergesort')
        t2 = time()
        # write time for each list length in results
        results[0, i, 0] = N
        results[0, i, 1] = (t2 - t1) / Nt

    # range of lengths for built in method
    Emin = 2  # log10 of minimal sample size
    Emax = 3  # log10 of maximal sample size
    pts = np.linspace(Emin, Emax, N_pts)  # contains exponents for list lengths
    for i in range(N_pts):
        N = int(10 ** (pts[i]))  # lentgh of list
        num_list = np.zeros(N)
        print(N)
        # inititalization for timing of list with length N
        for j in range(N):
            num_list[j] = rnd.uniform(0, 100)
        # timing for myquicksort
        t1 = time()
        for t in range(Nt):
            sorted_num_list = my_sort(num_list)
        t2 = time()
        # write time for each list length in results
        results[1, i, 0] = N
        results[1, i, 1] = (t2 - t1) / Nt

    print(results[0, :, :])
    print(results[1, :, :])
    x1 = np.log10(results[0, :, 0])
    y1 = np.log10(results[0, :, 1])
    x2 = np.log10(results[1, :, 0])
    y2 = np.log10(results[1, :, 1])
    k1, d1 = np.polyfit(x1, y1, 1)
    k2, d2 = np.polyfit(x2, y2, 1)

    plt.figure(0)
    plt.plot(x1, y1, 'x')
    plt.plot(x1, k1 * x1 + d1, '-')
    plt.ylabel('log(t)')
    plt.xlabel('log(N)')
    print('Slope and intercept of built-in in log-log plot are ' + str(k1)
          + ' and ' + str(d1))

    plt.figure(1)
    plt.plot(x2, y2, 'x')
    plt.plot(x2, k2 * x2 + d2, '-')
    print('Slope and intercept of myquicksort in log-log plot are ' + str(k2)
          + ' and ' + str(d2))
    plt.show()


if __name__ == '__main__':
    check_my_sort()

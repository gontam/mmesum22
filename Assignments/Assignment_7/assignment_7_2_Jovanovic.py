# Author: Oliver Jovanović
# Quicksort
# Libraries for "check_my_sort()".
import numpy as np
import random as rnd
import matplotlib.pyplot as plt
from time import time

# Quicksort Algorithm
def qs(low, high, lst):
    # If length of list is 1, end algorithm.
    if len(lst) == 1:
        return lst
    # If the smallest value is smaller than the highest, sort.
    if low < high:
        p = partition(low, high, lst)
        qs(low, p - 1, lst)
        qs(p + 1, high, lst)
    # return the sorted list.
    return lst

# Partition (Pointer) of Quicksort Algorithm
def partition(low, high, lst):
    # Begin with Pointer at value 0
    pointer = low
    pivot = lst[high]
    # Iterate through whole list:
    for j in range(low, high):
        # if current number in list is smaller than Pivot value
        if lst[j] <= pivot:
            # Switch the numbers
            lst[pointer], lst[j] = lst[j], lst[pointer]
            # Go further with pointer
            pointer = pointer + 1
    lst[pointer], lst[high] = lst[high], lst[pointer]
    return pointer


def myquicksort(lst):
    qs(0, len(lst) - 1, lst)
    print(lst)


# By Dr. Ramberger and a bit grammatically corrected and programmatically adapted by Oliver Jovanović, B.A., B.Sc.:
def check_my_sort():
    # general initialization
    N_pts = 5  # number of different lengths for which sort time is measured
    Nt = 10  # number of repetitions to measure time
    results = np.zeros([2, N_pts, 2])

    # range of lengths for build in method
    Emin = 5  # log10 of minimal sample size
    Emax = 6  # log10 of maximal sample size
    pts = np.linspace(Emin, Emax, N_pts)  # contains exponents for list lengths
    for i in range(N_pts):
        N = int(10 ** (pts[i]))  # length of list
        num_list = np.zeros(N)
        print(N)
        # initialization for timing of list with length N
        for j in range(N):
            num_list[j] = rnd.uniform(0, 100)
        # timing for build in method
        t1 = time()
        for t in range(Nt):
            sorted_num_list = np.sort(num_list, kind='mergesort')
        t2 = time()
        # write time for each list length in results
        results[0, i, 0] = N
        results[0, i, 1] = (t2 - t1) / Nt

    # range of lengths for build in method
    Emin = 2  # log10 of minimal sample size
    Emax = 3  # log10 of maximal sample size
    pts = np.linspace(Emin, Emax, N_pts)  # contains exponents for list lengths
    for i in range(N_pts):
        N = int(10 ** (pts[i]))  # length of list
        num_list = np.zeros(N)
        print(N)
        # initialization for timing of list with length N
        for j in range(N):
            num_list[j] = rnd.uniform(0, 100)
        # timing for myquicksort
        t1 = time()
        for t in range(Nt):
            sorted_num_list = myquicksort(num_list)
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
    plt.plot(x1, k1 * x1 + d1, '-')
    plt.plot(x2, k2 * x2 + d2, '-')
    plt.ylabel('log(t)')
    plt.xlabel('log(N)')
    plt.title('myquicksort Algorithm vs. np.sort() Algorithm')
    plt.legend(['myquicksort Algorithm', 'np.sort() Algorithm'])
    print('Slope and intercept of built-in in log-log plot are ' + str(k1)
          + ' and ' + str(d1))
    print('Slope and intercept of myquicksort in log-log plot are ' + str(k2)
          + ' and ' + str(d2))
    plt.show()

# Execute the check_my_sort() script.
check_my_sort()

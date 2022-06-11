# MME2- BI - Assignment 7.2 - Quicksort - Leonhard Hauptfeld
import numpy as np
from time import time
import matplotlib.pyplot as plt
import random as rnd

# Loosely based off pseudocode from https://www.geeksforgeeks.org/quick-sort/
def myquicksort(unsorted_list, low=None, high=None):
    # Copy the list we're working with
    sorted_list = np.copy(unsorted_list)
    # If low and high aren't given, assume them as the first and last element
    if low is None:
        low = 0
    if high is None:
        high = len(sorted_list)-1
    # Internal Quicksort Helpers
    def list_swap(list, i, j):
        temp = list[i]
        list[i] = list[j]
        list[j] = temp
    def list_part(list, low, high):
        pivot = list[high]
        i = low - 1
        for j in range(low,high):
            if list[j] < pivot:
                i += 1
                list_swap(list, i, j)
        list_swap(list, i+1, high)
        return i+1
    # Perform quicksort
    if low < high:
        partIndex = list_part(sorted_list, low, high)
        myquicksort(sorted_list, low, partIndex - 1)
        myquicksort(sorted_list, partIndex + 1, high)
    return sorted_list

# Compare scaling of different sort algos
# Taken from bioinf_sort.py
def check_my_sort():
    #general initialization
    N_pts = 5 # number of different length for which sort time is measured
    Nt = 10 #number of repetitions to measure time
    results = np.zeros([2, N_pts, 2])
    
    #range of lengths for built in method
    Emin = 5  # log10 of minimal sample size
    Emax = 6  # log10 of maximal sample size
    pts = np.linspace(Emin, Emax, N_pts) #contains exponents for list lengths
    for i in range(N_pts):
        N = int(10 ** (pts[i])) # length of list
        num_list = np.zeros(N)
        print(N)
        #inititalization for timing of list with length N
        for j in range(N):
            num_list[j] = rnd.uniform(0, 100)
        #timing for built in method
        t1 = time()
        for t in range(Nt):
            sorted_num_list = np.sort(num_list, kind='mergesort')
        t2 = time()
        #write time for each list length in results
        results[0, i, 0] = N
        results[0, i, 1] = (t2 - t1) / Nt
        
    #range of lengths for built in method
    Emin = 2  # log10 of minimal sample size
    Emax = 3  # log10 of maximal sample size
    pts = np.linspace(Emin, Emax, N_pts) #contains exponents for list lengths
    for i in range(N_pts):
        N = int(10 ** (pts[i])) # length of list
        num_list = np.zeros(N)
        print(N)
        #inititalization for timing of list with length N
        for j in range(N):
            num_list[j] = rnd.uniform(0, 100)
        #timing for myquicksort
        t1 = time()
        for t in range(Nt):
            sorted_num_list = myquicksort(num_list)
        t2 = time()
        #write time for each list length in results
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
                
if __name__ == "__main__":
    check_my_sort()
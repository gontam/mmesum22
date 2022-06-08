import numpy as np
import matplotlib.pyplot as plt
from time import time
import random as rnd


# checks if the list between these parameters is sorted correctly
def check_sort(list, start, pivot):
    for x in range(start, pivot - 1, 1):
        if list[x + 1] < list[x]:
            return True
    return False


def myquicksort(list):
    sort_algorythm(list, 0, len(list) - 1)


# Sorting algorythm based on quick sort, recursively calling itself if not finished with sorting
def sort_algorythm(list, start, pivot):
    while check_sort(list, start, pivot):
        old_pivot = pivot
        pivot = pivot - 1
        for x in range(start, pivot, 1):
            if x < pivot:
                while list[x] > list[pivot]:
                    if pivot - 1 == x:
                        list[[x, pivot]] = list[[pivot, x]]
                        pivot = pivot - 1
                    else:
                        list[[pivot, pivot - 1]] = list[[pivot - 1, pivot]]
                        list[[x, pivot]] = list[[pivot, x]]
                        pivot = pivot - 1
        if check_sort(list, start, pivot):
            sort_algorythm(list, start, pivot)
        if check_sort(list, pivot, old_pivot):
            sort_algorythm(list, pivot, old_pivot)
    return pivot


# based on bioinf_sort.py from moodle
def check_my_sort():
    # general initialization
    N_pts = 5  # number of different lentghs for which sort time is measured
    Nt = 10  # number of repetitions to measure time
    results = np.zeros([2, N_pts, 2])

    # range of lengths for built in method
    Emin = 5  # log10 of minimal sample size
    Emax = 6  # log10 of maximal sample size
    pts = np.linspace(Emin, Emax, N_pts)  # contains exponents for list lengths

    for i in range(N_pts):
        N = int(10 ** (pts[i])) #lentgh of list
        num_list = np.zeros(N)
        # inititalization for timing of list with length N
        for j in range(N):
            num_list[j] = rnd.uniform(1, 100)
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
        # inititalization for timing of list with length N
        for j in range(N):
            num_list[j] = rnd.uniform(0, 100)
        # timing for myquicksort
        t1 = time()
        for t in range(Nt):
            sorted_num_list = myquicksort(num_list)
        t2 = time()
        # write time for each list length in results
        results[1, i, 0] = N

        # For the first value, we have the problem that t2 - t1 results in 0 which is why we add a very small number
        # to not run into problems later when calculating the log
        if ((t2 - t1) / Nt)==0:
            t2 += 0.000001
        results[1, i, 1] = (t2 - t1) / Nt

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
    plt.title('Performance np.sort')
    print(f'Performance np.sort :\n  Slope: {str(k1)}\n  Intercept: {str(d1)}')

    plt.figure(1)
    plt.plot(x2, y2, 'x')
    plt.plot(x2, k2 * x2 + d2, '-')
    plt.ylabel('log(t)')
    plt.xlabel('log(N)')
    plt.title('Performance custom quicksort')
    print(f'Performance custom quicksort :\n  Slope: {str(k2)}\n  Intercept: {str(d2)}')
    plt.show()


check_my_sort()

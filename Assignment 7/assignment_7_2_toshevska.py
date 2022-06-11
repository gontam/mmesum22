import numpy as np
from time import time
import random
import timeit
import matplotlib.pyplot as plt

#reference: https://stackoverflow.com/questions/20175380/quick-sort-python-recursion

def partition(lst, start, end):
    pos = start

    for i in range(start, end):
        if lst[i] < lst[end]:
            lst[i],lst[pos] = lst[pos],lst[i]
            pos += 1

    lst[pos],lst[end] = lst[end],lst[pos]
    return pos

def myquicksort(lst, start, end):
    if start < end:
        pos = partition(lst, start, end)
        myquicksort(lst, start, pos - 1)
        myquicksort(lst, pos + 1, end)
        return lst

def func_2(inp):
    return np.sort(inp)[10]

def check_my_sort():
    N_pts = 4 # number of lengths for which sort time is measured
    Emin = 5  # log10 of minimal sample size
    Emax = 6  # log10 of maximal sample size
    Nt = 8
    pts = np.linspace(Emin, Emax, N_pts)
    for i in range(N_pts):
        N = int(10 ** (pts[i]))
        print(N)
        num_list = np.random.rand(N);
    print(num_list)
    sorted_num_list = myquicksort(num_list, 0, len(num_list)-1)
    print(sorted_num_list)
    #we get result that there is no error with the comparison
    np_sort = np.sort(num_list)
    print(sum(sorted_num_list-np_sort))
    #check that it works
    #array = [38, 44, 78, 87, 19, 21, 76, 58, 88, 73, 97, 16, 21]
    #myquicksort(array, 0, len(array) - 1)
    #print(array)
    plt.figure()
    plt.xlabel("log(t)")
    plt.ylabel("log(N)")
    plt.title('Quick Sort vs. np.sort()')
    plt.loglog(sorted_num_list, color="purple")
    plt.loglog(np_sort, color="blue")
    plt.show()
    #the results say that both algorithms have the same complexity cause the curves are overlapping

check_my_sort()

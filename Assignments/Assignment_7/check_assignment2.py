# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 11:53:59 2020
check assignment programs
@author: Ramberger
"""

import numpy as np


def test1(lu_solve):
    """
    Check if the lu_solve algorithm really solves a linear system
    should print out the sum of residuals, i.e. a small number
    """
    Emin=2 # log10 of minimal sample size
    Emax=2 # log10 of maximal sample size
    N_pts=1
    pts=np.linspace(Emin,Emax,N_pts)
    for i in range(N_pts):
        N=int(10 ** pts[i])
        A=np.random.randint(1,10,[N,N]).astype('float')
        b=np.random.randint(1,10,N).astype('float')
        A_orig=A.copy()
        x=lu_solve(A,b)
        print(sum(np.matmul(A,x)[:]-b))
        print(sum(np.matmul(A_orig,x)[:]-b))

def test2(myquicksort):
    """
    Check if the myquicksort algorithm really sorts a list of floats
    should print out the sum of differences between the sorted list using
    my quicksort and a built in sort
    """
    N_pts=1
    Emin = 4 # log10 of minimal sample size
    Emax = 4  # log10 of maximal sample size
    pts = np.linspace(Emin, Emax, N_pts)
    for i in range(N_pts):
        N = int(10 ** (pts[i]))
        print(N)
        num_list = np.random.rand(N);
    print(sum(myquicksort(num_list)-np.sort(num_list)))
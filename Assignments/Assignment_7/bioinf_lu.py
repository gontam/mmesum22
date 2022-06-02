# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 18:58:46 2020
Module contains the function lu() that carries out LU factorization of a square matrix.
@author: Ramberger
"""

import numpy as np

def lu(A,piv=True,Tol=10**(-7)):
    """Performs LU decomposition of A, overwrites original matrix with L-E+U

    Parameters
    ----------
    A : numpy.array
        contains the matrix to be decomposed in LU.
    piv : Boolean, optional
        allows to turn off row pivoting. The default is True.
    Tol : float, optional
        tolerance to determine if matrix is (nearly) degenerate.
        The default is 10**(-7).

    Returns
    -------
    A : numpy.array
        Contains the matrix L-E below the diagonal, and the matrix U on
        and above the diagonal. The diagonal elements of L are all 1 by
        construction and don't need to be stored. Overwrites the original
        matrix.
    permu : numpy.array (vector)
        Stores the permutation vector permu[:]: The i-th row of the pivotised LU
        decomposition ( (LU)[i,:] ) corresponds to the permu[i]-th row of the
        original matrix A ( A[permu[i],:] ).
        In other words the permutation matrix is 1 at all positions
        P[i,permu[i]] and 0 else. Then LU=PA.

    """
    
    Nr,Nc = A.shape # Number of rows and columns
    if Nr != Nc:
        print('Matrix is not square')
        return None
    N = Nr
    permu = np.array(range(N)) # Initialize permutation vector p
    for j in range(N-1):
        p_val = abs(A[j,j])
        p_ind=j
        for i in range(j+1,N):
            if abs(A[i,j]) > p_val:
                p_val = abs(A[i,j])
                p_ind=i
        if p_val < Tol:
            print('Matrix nearly degenerate (tolereance = '+str(Tol)+').')
            return None
        if piv:
            #Row pivoting, i.e. put row with largest element in current
            #column on top
            permu[j], permu[p_ind] = permu[p_ind], permu[j]
            A[[j,p_ind],:]=A[[p_ind,j],:]
        
        for i in range(j+1,N):
            A[i,j] = A[i,j] / A[j,j] #A[i,j] now stores L[i,j]
            for k in range(j+1,N):
                A[i,k] = A[i,k] - A[i,j]*A[j,k]
            # A[i,j+1:N] = A[i,j+1:N] - A[i,j]*A[j,j+1:N] #faster assignment
    return A,permu
# LU Linear equation solver:
# @Author: Jovanović
# Import necessary libraries:
import numpy as np
import sys

# Give input for matrix:
print('Please give a number for columns:')
col = input()
col = int(col)
print('Please give a number for rows:')
row = input()
row = int(row)

# Randomized Matrix:
N = np.random.rand(row, col)

# Checking if rows and columns are equal:
if row == col:
    # If Columns and Rows are equal, do the LU decomposition.
    print("It is a quadratic matrix!")
    print("The current quadratic Matrix looks like:", N)
    Lsub(L, b)

else:
    # If Columns and Rows are unequal, exit the script.
    print("It is not a quadratic matrix. The Program will be shut down.")
    sys.exit(0)

# Lower triangular system:
def Lsub(L, b):
    L = N
    b = np.ones(row, dtype="float")

# Upper triangular system:
def Usub(U, y):
    pass

def lu(A, piv=True, Tol=10 ** (-7)):
    Nr, Nc = A.shape  # Number of rows and columns
    if Nr != Nc:
        print('Matrix is not square')
        return None
    N = Nr
    permu = np.array(range(N))  # Initialize permutation vector p
    for j in range(N - 1):
        p_val = abs(A[j, j])
        p_ind = j
        for i in range(j + 1, N):
            if abs(A[i, j]) > p_val:
                p_val = abs(A[i, j])
                p_ind = i
        if p_val < Tol:
            print('Matrix nearly degenerate (tolereance = ' + str(Tol) + ').')
            return None
        if piv:
            # Row pivoting, i.e. put row with largest element in current
            # column on top
            permu[j], permu[p_ind] = permu[p_ind], permu[j]
            A[[j, p_ind], :] = A[[p_ind, j], :]

        for i in range(j + 1, N):
            A[i, j] = A[i, j] / A[j, j]  # A[i,j] now stores L[i,j]
            for k in range(j + 1, N):
                A[i, k] = A[i, k] - A[i, j] * A[j, k]
            # A[i,j+1:N] = A[i,j+1:N] - A[i,j]*A[j,j+1:N] #faster assignment
    return A, permu
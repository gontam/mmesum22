# LU Linear equation solver:
# @Author: Jovanović
# Based on: https://www.youtube.com/watch?v=UlWcofkUDDU&ab_channel=Mathispower4u
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
tmp = 0.0

# Compute U by given L.
def compute(A, L):
    k = 1
    l = 0
    m = 0
    if A[row - k, l] == A[0, l]:
        pass
    else:
        if A[row - k, l] < 0:
            if A[0, l] > A[row - k, l]:
                tmp = A[row - k, l] / A[0, l]
                L[row - k, l] = tmp
                for y in range(l, row - 1):
                    A[row - k, l] = A[row - k, l] + tmp * A[0, l]
                    print(A)
                    print(L)
            elif A[0, l] < A[row - k, l]:
                tmp = A[0, l] / A[row - k, l]
                L[row - k, l] = tmp
                for y in range(l, row - 1):
                    A[row - k, l] = A[row - k, l] + tmp * A[0, l]
                    print(A)
                    print(L)

        elif A[row - k, l] > 0:
            if A[0, l] > A[row - k, l]:
                tmp = A[row - k, l] / A[0, l]
                L[row - k, l] = tmp
                for y in range(l, row - 1):
                    A[row - k, l] = A[row - k, l] - tmp * A[0, l]
                    print(A)
                    print(L)
            elif A[0, l] < A[row - k, l]:
                tmp = A[0, l] / A[row - k, l]
                L[row - k, l] = tmp
                for y in range(l, row - 1):
                    A[row - k, l] = A[row - k, l] - tmp * A[0, l]
                    print(A)
                    print(L)



# Randomized Matrix:
A = np.random.rand(row, col)
# Checking if rows and columns are equal:
if row == col:
    # If Columns and Rows are equal, do the LU decomposition.
    print("The quadratic Matrix A looks like:\n", A)
    L = np.eye(row)
    print("Matrix L looks like:\n", L)
    compute(A, L)
else:
    # If Columns and Rows are unequal, exit the script.
    print("It is not a quadratic matrix. The Program will be shut down.")
    sys.exit(0)



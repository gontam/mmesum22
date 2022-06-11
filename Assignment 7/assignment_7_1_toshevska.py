import random
import numpy as np
# reference: https://numpy.org/doc/stable/reference/generated/numpy.linalg.solve.html
#we must have an equal number for rows and columns

print("Enter number of linear equations/variables:")
n = int(input())

# #creating randomized matrix A
A = np.random.randint(-5, 5, size=(n, n))

#define matrix B
B = np.array([4, 5, 6])
print(A)

#linalg.solve is the function of NumPy to solve a system of linear scalar equations
x = np.linalg.solve(A, B)
print("Solution:\n",x)

#check if solution is correct
print(np.allclose(np.dot(A, x), B))






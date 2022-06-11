import bioinf_lu as bioinf
import numpy as np

# LU solver
def lu_solve(A,b,piv=True):
    A, permu = bioinf.lu(A, piv)
    y = Lsub(L, b)
    return Usub(U, y)

# Forward substitution
def Lsub(L,b):
    n = L.shape[0]
    y = np.zeros_like(b, dtype=np.double);
    y[0] = b[0] / L[0, 0]
    for i in range(1, n):
        y[i] = (b[i] - np.dot(L[i,:i], y[:i])) / L[i,i]
    return y

# Back substitution
def Usub(U,y):
    n = U.shape[0]
    x = np.zeros_like(y, dtype=np.double);

    x[-1] = y[-1] / U[-1, -1]
    for i in range(n-2, -1, -1):
        x[i] = (y[i] - np.dot(U[i,i:], x[i:])) / U[i,i]
        
    return x
import numpy as np

M = int(input("M: "))
N = int(input("N: "))

if M < 100 or N < 100:
    a = np.random.randint(10, size=(N, M))
    print(a)
else:
    print("You have selected a too high M or N!")
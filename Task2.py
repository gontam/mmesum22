from random import randrange
import numpy as np

array = np.ndarray()
for i in range(3):
    for j in range(3):
        random_number = randrange(0, 100, 1)
        array[i, j] = random_number

print(array)

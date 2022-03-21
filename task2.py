import random

import numpy as np

random_matrix = [[random.random() for e in range(2)] for e in range(3)]
for i in range(3):
    for j in range(2):
        print(random_matrix[i][j])
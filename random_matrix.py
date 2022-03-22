#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import modules
import numpy as np

# lambda funktion to generate a 3x3 matrix filled with random integers
random_3x3_matrix = lambda: np.random.randint(100, size=(3, 3))

def main():
    # save a generated matrix in a variable and print it
    generated_matrix = random_3x3_matrix()
    print(generated_matrix)


if __name__ == '__main__':
    main()


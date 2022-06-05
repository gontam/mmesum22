# Bioinformatics assignment 6
# authors: Theodor Abart, Gregor Widhalm

import numpy as np
from random import randint
import time
import matplotlib.pyplot as plt


def create_glider_gun_grid():
    # creates glider gun grid
    # according to https://www.unioviedo.es/compnum/labs/PYTHON/conway.html

    glider_gun = \
        [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
         [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]

    glider_gun_grid = np.zeros((50, 70))
    glider_gun_grid[1:10, 1:37] = glider_gun

    return glider_gun_grid


def create_grid(size_grid, n_alive):
    s = [size_grid, size_grid]
    grid = np.zeros(s)
    grid = create_start(grid, size_grid, n_alive)
    return grid


def create_start(grid, size_grid, n_alive):
    for i in range(n_alive):
        grid[randint(0, size_grid - 1), randint(0, size_grid - 1)] = 1
    return grid


def update_grid(grid):
    new_grid = grid.copy()

    # loop through every cell and count all neighbours in order to apply the rules
    for i in range(new_grid.shape[0] - 1):
        for j in range(new_grid.shape[0] - 1):
            this_sum = int((grid[i, (j - 1)] + grid[i, (j + 1)] +
                            grid[(i - 1), j] + grid[(i + 1), j] +
                            grid[(i - 1), (j - 1)] + grid[
                                (i - 1), (j + 1)] +
                            grid[(i + 1), (j - 1)] + grid[
                                (i + 1), (j + 1)]))
            if grid[i, j] == 1:
                if (this_sum < 2) or (this_sum > 3):
                    new_grid[i, j] = 0
            else:
                if this_sum == 3:
                    new_grid[i, j] = 1

    return new_grid


def show_grid(grid):
    # print(grid)
    plt.ion()
    fig, ax = plt.subplots()
    ax.matshow(grid)
    plt.show()


def main():
    print("Hello")
    input_selection = 0
    while input_selection != 1 and input_selection != 2:
        input_selection = (
            int(input("Do you prefer \"Gosper glider gun\" (guaranteed 100 generations) (1) or a random-based initialization (2)?: ")))

    # glider gun
    if input_selection == 1:
        n_gen = 100
        grid = create_glider_gun_grid()
    # random-based initialisation: 20% living cells
    elif input_selection == 2:
        size_grid = int(input("Enter Field Dimension: "))
        n_gen = int(input("Enter Generations to simulate: "))
        n_alive = round((size_grid * size_grid) * 0.20)
        grid = create_grid(size_grid, n_alive)

    show_grid(grid)
    for i in range(n_gen):
        grid = update_grid(grid)
        show_grid(grid)
        time.sleep(0.2)


if __name__ == '__main__':
    main()

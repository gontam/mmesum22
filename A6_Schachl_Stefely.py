#Assignment 6: Conway's Game of Life

#Authors: Johanna Schachl, Lisa Stefely

#Rules:
#    A living cell dies if it has less than two living neighboring cells.
#    A living cell with two or three living neighbors continues to live.
#    A living cell with more than three living neighbor cells dies in the next time step.
#    A dead cell is revived if it has exactly three living neighbor cells.

#Import libraries
import pygame
import numpy as np


#Define Colors of cells near death, alive, background and grid
near_death = (200, 200, 225)
alive = (144, 26, 52)
background = (10, 10, 10)
grid = (30, 30, 60)

#Define Dimensions of pygame window
dimx = 120
dimy = 90
cellsize = 8

#Define funtion for  one cycle, returns updated pattern
def cycle(surface_cycle, cells_cycle, size_cycle):
    new_pattern = np.zeros((cells_cycle.shape[0], cells_cycle.shape[1]))

    for r, c in np.ndindex(cells_cycle.shape):
        num_alive = np.sum(cells_cycle[r - 1:r + 2, c - 1:c + 2]) - cells_cycle[r, c]

        #Implementation of rules
        if cells_cycle[r, c] == 1 and num_alive < 2 or num_alive > 3:
            color = near_death
        elif (cells_cycle[r, c] == 1 and 2 <= num_alive <= 3) or (cells_cycle[r, c] == 0 and num_alive == 3):
            new_pattern[r, c] = 1
            color = alive

        color = color if cells_cycle[r, c] == 1 else background
        pygame.draw.rect(surface_cycle, color, (c * size_cycle, r * size_cycle, size_cycle - 1, size_cycle - 1))

    return new_pattern

#Define initialization function with given dimensions and pattern matrix
def init(dimx, dimy):
    cells = np.zeros((dimy, dimx))
    pattern = np.array([[0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,1,1,0,0,0,0,0,0],
                        [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                        [0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0],
                        [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,1,0,0,1,1,1,0,0,0,0,0,0,0,1,1,0,0,0],
                        [1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
                        [1,1,0,0,0,1,0,0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,1,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0],
                        [0,0,0,0,0,1,0,1,1,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                        [1,1,0,0,0,1,0,1,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0],
                        [1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0]]);
    pos = (3,3)
    cells[pos[0]:pos[0]+pattern.shape[0], pos[1]:pos[1]+pattern.shape[1]] = pattern
    return cells

def main(generations):

    #Initialize pygame
    pygame.init()
    surface = pygame.display.set_mode((dimx * cellsize, dimy * cellsize))
    pygame.display.set_caption("Conway's Game of Life")

    cells = init(dimx, dimy) #Call of function init for returning initial pattern

    #Checks if amount of generations is reached
    #If yes, quits the game
    for ii in range (1,generations):
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                return

        surface.fill(grid)

        #Cycle funciton is called -> pattern is updated
        cells = cycle(surface, cells, cellsize)
        pygame.display.update()

if __name__ == "__main__":
    print('Enter number of generations:')
    generations = input()
    print('Caution! The window will open in the background')
    main(int(generations))
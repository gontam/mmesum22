from filecmp import DEFAULT_IGNORES
import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation

DEFAULT_WIDTH = 50
DEFAULT_HEIGHT = 50
DEFAULT_GENERATIONS = 100

class GameOfLife:
    def __init__(self, size):
        self.field = np.zeros(size, dtype=int)
        self.generation = 1
        self.on(0,0)
    def progress(self):
        if self.at(0,0) == 0:
            self.on(0,0)
        else:
            self.off(0,0)
        self.generation += 1
    def get_field(self):
        return self.field
    def at(self, x, y):
        return self.field[x, y]
    def on(self, x, y):
        self.field[x, y] = 1
    def off(self, x, y):
        self.field[x, y] = 0
    def neighbors(self, x, y):
        n = self.field[x-1:x+2, y-1:y+2]
        return n.sum() - n[1,1]


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description="Run conway's game of life")
    parser.add_argument('--width', type=int, default=DEFAULT_WIDTH, help='Width of field')
    parser.add_argument('--height', type=int, default=DEFAULT_HEIGHT, help='Height of field')
    parser.add_argument('--generations', type=int, default=DEFAULT_GENERATIONS, help='Generations to simulate')
    args = parser.parse_args()
    # Game
    game = GameOfLife((args.width, args.height))
    # Display
    def update(frame):
        game.progress()
        #print(game.get_field())
        im.set_data(game.get_field())
        return [im]
    fig, ax = plt.subplots()
    im = ax.imshow(game.get_field(), cmap=cm.binary, interpolation='nearest')
    ani = FuncAnimation(fig, update, frames=np.linspace(0, 2*np.pi, 128), blit=True)
    plt.show()

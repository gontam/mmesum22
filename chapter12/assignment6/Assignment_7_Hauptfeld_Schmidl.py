import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import itertools
import sys

DEFAULT_WIDTH = 50
DEFAULT_HEIGHT = 50
DEFAULT_INTERVAL = 500
DEFAULT_PATTERN = 'pentomino'

PATTERNS = {
    'pentomino': [[ 0, 1, 1 ],
                  [ 1, 1, 0 ],
                  [ 0, 1, 0 ]],
    'tetris':    [[ 0, 1, 0 ],
                  [ 1, 1, 1 ]]
}

class GameOfLife:
    def __init__(self, size):
        self.field = np.zeros(size, dtype=int)
        self.generation = 1
    def progress(self):
        prev_field = np.copy(self.field)
        for y,x in itertools.product(range(self.field.shape[0]), range(self.field.shape[1])):
            cell_state = self.at(x, y, field=prev_field)
            alive_neighbors = self.neighbors(x, y, field=prev_field).sum() - cell_state
            if cell_state and alive_neighbors not in (2,3):
                self.off(x, y)
            elif not cell_state and alive_neighbors == 3:
                self.on(x, y)
        self.generation += 1
    def get_field(self):
        return self.field
    def insert_pattern(self, x, y, pattern, field=None):
        if field is None: field = self.field
        for py,px in itertools.product(range(pattern.shape[0]), range(pattern.shape[1])):
            field[y + py, x + px] = pattern[py, px]
    def at(self, x, y, field=None):
        if field is None: field = self.field
        return field[y, x]
    def on(self, x, y, field=None):
        if field is None: field = self.field
        field[y, x] = 1
    def off(self, x, y, field=None):
        if field is None: field = self.field
        field[y, x] = 0
    def neighbors(self, x, y, field=None):
        if field is None: field = self.field
        return field[
            max(y-1, 0):min(y+2, field.shape[0]),
            max(x-1, 0):min(x+2, field.shape[1])
        ]


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description="Run Conway's game of life")
    parser.add_argument('--width', type=int, default=DEFAULT_WIDTH, help='Width of field')
    parser.add_argument('--height', type=int, default=DEFAULT_HEIGHT, help='Height of field')
    parser.add_argument('--generations', type=int, required=True, help='Generations to simulate')
    parser.add_argument('--pattern', type=str, default=DEFAULT_PATTERN, choices=PATTERNS.keys(), help='Starting pattern. One of: ' + ', '.join(PATTERNS.keys()))
    parser.add_argument('--interval', type=str, default=DEFAULT_INTERVAL, help='At what interval a step is processed, in ms (animation speed, lower = faster)')
    args = parser.parse_args()
    # Convert starting pattern
    pattern = np.array(PATTERNS[args.pattern])
    # Game
    game = GameOfLife((args.width, args.height))
    game.insert_pattern(round(args.width / 2 - pattern.shape[1] / 2), round(args.height / 2 - pattern.shape[0] / 2), pattern)
    # Display
    def update(frame):
        game.progress()
        txt.set_text(f"Generation {game.generation}")
        im.set_data(game.get_field())
        return [im, txt]
    fig, ax = plt.subplots()
    fig.canvas.manager.set_window_title("Game of Life")
    im = ax.imshow(game.get_field(), cmap=cm.binary, interpolation='nearest')
    txt = ax.text(10, 10, f"Generation {game.generation}", bbox=dict(fill=False, edgecolor='red', linewidth=2))
    ani = FuncAnimation(fig, update, frames=np.linspace(0, 2*np.pi, 128), blit=True, interval=args.interval)
    plt.show()

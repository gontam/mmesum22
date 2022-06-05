import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import numpy as np
import time

# This function calculates the changes on the board for the next generation, this includes births as well as deaths simultaneously.
def calculateBirthAndDeath(index, column):
    # The criteria for birth and death is on purpose changed by 0.5 due to us including the stone itself when counting neighbors.
    # This way we can increase performance by ignoring fields without stones & neighbors
    if boardNeighbors.loc[index][column] == 3 or boardNeighbors.loc[index][column] == 3.5:
        board.loc[index][column] = 1
    elif boardNeighbors.loc[index][column] >= 4 or boardNeighbors.loc[index][column] <= 1.5:
        board.loc[index][column] = 0

# A board with a reduced amount of fields (for optimization purposes) is handed to this function.
# This function then iterates over the relevant fields and uses the checkNeighboringStones() function to calculate the neighbours.
# The values are saved into the boardNeighbors, which is later on used to calculates births and deaths.
def calculateNeighborsBoard(livingCells):
    for index in range(livingCells.index.min()-1, livingCells.index.max()+2, 1):
        for column in range(livingCells.columns.min()-1, livingCells.columns.max()+2, 1):
            if index < 0 or index > 199:
                continue
            if column < 0 or column > 199:
                continue
            boardNeighbors.loc[index][column] = checkNeighboringStones(index, column)

# This function calculates the amount of living neighbours for each field.
def checkNeighboringStones(index, column):
    numberOfLivingNeighbors = 0
    for x_check in range(-1, 2, 1):
        for y_check in range(-1, 2, 1):
            if index+x_check < 0 or index+x_check > 199:
                continue
            if column+y_check < 0 or column+y_check > 199:
                continue
            if board.loc[index+x_check][column+y_check] == 1:
                if x_check == 0 and y_check == 0:
                    numberOfLivingNeighbors += 0.5
                else:
                    numberOfLivingNeighbors += 1
    return numberOfLivingNeighbors

# This function calculates a subboard which contains alive cells and then calls the calculateNeighborsBoard() with this subboard.
# The boardNeighbors is then used to calculate the the births and deaths by calling the calculateBirthAndDeath() function.
def gameIteration():
    livingCells = board.where(board == 1).dropna(how='all').dropna(how='all',axis=1)
    calculateNeighborsBoard(livingCells)
    #print(livingCells.to_string())
    #print(boardNeighbors.where(boardNeighbors > 0).dropna(how='all').dropna(how='all',axis=1).to_string())

    for index in boardNeighbors.where(boardNeighbors > 0).dropna(how='all').dropna(how='all',axis=1).index:
        for column in boardNeighbors.where(boardNeighbors > 0).dropna(how='all').dropna(how='all',axis=1).columns:
            #print(f'{index}|{column} -> {boardNeighbors.loc[index][column]}')
            calculateBirthAndDeath(index, column)

    #print(board.where(board > 0).dropna(how='all').dropna(how='all', axis=1).to_string())

# Here the initial board formation is set.
# This means placing "alive" cells (1s) on the otherwise "dead" (0s) board
def setInitalFormation(n):
    for i in range(0,n,1):
        selectedFormation = np.random.randint(1,7)
        selectedPosition = np.random.randint(10,190, size=2)
        #Square
        if selectedFormation == 1:
            board.loc[selectedPosition[0]][selectedPosition[1]] = 1
            board.loc[selectedPosition[0]+1][selectedPosition[1]+1] = 1
            board.loc[selectedPosition[0]+0][selectedPosition[1]+1] = 1
            board.loc[selectedPosition[0]+1][selectedPosition[1]+0] = 1
        #Glider
        elif selectedFormation == 2:
            board.loc[selectedPosition[0]][selectedPosition[1]] = 1
            board.loc[selectedPosition[0] + 1][selectedPosition[1] + 1] = 1
            board.loc[selectedPosition[0] - 1][selectedPosition[1]] = 1
            board.loc[selectedPosition[0] - 1][selectedPosition[1] + 1] = 1
            board.loc[selectedPosition[0] - 1][selectedPosition[1] + 2] = 1
        #F-Pentomino
        elif selectedFormation == 3:
            board.loc[selectedPosition[0]][selectedPosition[1]] = 1
            board.loc[selectedPosition[0] + 1][selectedPosition[1]] = 1
            board.loc[selectedPosition[0]][selectedPosition[1] - 1] = 1
            board.loc[selectedPosition[0] - 1][selectedPosition[1]] = 1
            board.loc[selectedPosition[0] - 1][selectedPosition[1] + 1] = 1
        #Boat
        elif selectedFormation == 4:
            board.loc[selectedPosition[0]][selectedPosition[1]] = 1
            board.loc[selectedPosition[0]][selectedPosition[1] + 1] = 1
            board.loc[selectedPosition[0] + 1][selectedPosition[1]] = 1
            board.loc[selectedPosition[0] + 2][selectedPosition[1] + 1] = 1
            board.loc[selectedPosition[0] + 1][selectedPosition[1] + 2] = 1
        # Acorn
        elif selectedFormation == 5:
            board.loc[selectedPosition[0]][selectedPosition[1]] = 1
            board.loc[selectedPosition[0]][selectedPosition[1] + 1] = 1
            board.loc[selectedPosition[0] - 2][selectedPosition[1] + 1] = 1
            board.loc[selectedPosition[0] - 1][selectedPosition[1] + 3] = 1
            board.loc[selectedPosition[0]][selectedPosition[1] + 4] = 1
            board.loc[selectedPosition[0]][selectedPosition[1] + 5] = 1
            board.loc[selectedPosition[0]][selectedPosition[1] + 6] = 1
        elif selectedFormation == 6:
            board.loc[selectedPosition[0]][selectedPosition[1]] = 1
            board.loc[selectedPosition[0]][selectedPosition[1] + 1] = 1
            board.loc[selectedPosition[0]][selectedPosition[1] + 2] = 1
            board.loc[selectedPosition[0] - 1][selectedPosition[1] + 1] = 1


# Tried to work with FuncAnimation, however, we did not succeed and therefore decided to create a plot for each
# iterration
def update(frame):
    gameIteration()
    im.set_data(board)
    return [im]

# This function plots the pandas Dataframe
def plotting():
    dg = board[0]
    fig, ax = plt.subplots()
    ax.set_xlim(0, len(dg))
    ax.set_ylim(0, len(dg))
    plt.gca().invert_yaxis()
    global im
    im = ax.imshow(board, cmap=cm.binary, interpolation='nearest')
    gameIteration()
    #ani = FuncAnimation(fig, update, frames=np.linspace(0, 2*np.pi, 128), blit=True)
    plt.show()


boardSizeX = 200
boardSizeY = 200
board = pd.DataFrame(index=range(boardSizeY), columns=range(boardSizeX), dtype=int)
boardNeighbors = pd.DataFrame(index=range(boardSizeY), columns=range(boardSizeX), dtype=int)
for col in board.columns:
    board[col].values[:] = 0
    boardNeighbors[col].values[:] = 0
numberOfLoops = input('Please input the number of generations')
numberOfInitialGenerations = input('Please input the number of initial objects that you want to place on the board')
setInitalFormation(int(numberOfInitialGenerations))
for loopIteration in range(0,int(numberOfLoops)+1,1):
    plotting()
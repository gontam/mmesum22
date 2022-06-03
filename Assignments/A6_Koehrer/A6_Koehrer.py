def prompt_user():
    # Prompt the user to provide the number of generations to be simulated (n > 100)
    return int(input("The simulation becomes stable after generation 424 with the used initial generation\n"
                     "Pleaser enter how many generations should be simulated: "))


def output_grid(grid_, row_, column_, i_):
    # Print the grid on the console
    print(f"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
          f"+++++++++++++++++++++++++++++\nGeneration {i_}")
    for i in range(row_):
        for j in range(column_):
            if grid_[i][j] == 0:
                print(f".  ", end="")
            else:
                print(f"#  ", end="")
            if j % (column_-1) == 0 and j != 0:
                print("")


def initial_generation(grid_):
    # Define the initial generation with which the simulation starts
    # Initial data points for a tiny squid - becomes stable with generation 424 by a grid size of 40 x 40
    row_index = [15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 18, 18, 18,
                    18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20, 20, 21, 21,
                    21, 22, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 24]
    column_index = [18, 19, 20, 21, 22, 15, 16, 18, 19, 20, 21, 22, 24, 25, 16, 18, 20, 22, 24, 16, 18,
                       19, 20, 21, 22, 24, 16, 18, 19, 20, 21, 22, 24, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                       18, 20, 22, 15, 16, 18, 20, 22, 24, 25, 16, 18, 20, 22, 24, 16, 17, 18, 20, 22, 23, 24]

    for i in range(len(row_index)):
        grid_[row_index[i]][column_index[i]] = 1


def simulate_generation(grid_, row_, column_):
    # Simulate the next generation
    grid_temp = [[0 for j in range(row)] for i in range(column)]
    row_index = [-1, 0,  +1, +1, +1,  0, -1, -1]
    column_index = [+1, +1, +1,  0, -1, -1, -1,  0]

    count = 0
    for i in range(row_):
        for j in range(column_):
            for k in range(len(row_index)):
                # Check if neighbour is out of bound
                if i + row_index[k] < 0 or i + row_index[k] > row_ - 1 \
                        or j + column_index[k] < 0 or j + column_index[k] > column_ - 1:
                    continue
                elif grid_[i+row_index[k]][j+column_index[k]] == 1:
                    count += 1
            # Check if grid point has enough living neighbours to become alive
            if grid_[i][j] == 0:
                if count == 3:
                    grid_temp[i][j] = 1
                else:
                    grid_temp[i][j] = 0
            # Otherwise, check if grid point has enough living neighbours to stay alive
            else:
                if count == 2 or count == 3:
                    grid_temp[i][j] = 1
                else:
                    grid_temp[i][j] = 0
            count = 0
    return grid_temp


''' MAIN PROGRAM '''

row = column = 39
# Generate and fill the grid with 0's
grid = [[0 for j in range(row)] for i in range(column)]

# Prompt the user to provide the number of generations to be simulated
numb_gen = prompt_user()

# Define the initial generation with which the simulation starts
initial_generation(grid)

# Print the grid with generation 0 on the console
output_grid(grid, row, column, 0)

for i in range(numb_gen):
    # Simulate the next generation
    grid = simulate_generation(grid, row, column)
    # Print the grid on the console
    output_grid(grid, row, column, i+1)

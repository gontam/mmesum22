import random
import pandas as pd
import fileinput


def kmeans(data, n, seeds, iteration):
    data = data.astype('float')
    # calculating Manhattan distance
    data['seed1_distance'] = abs(data[0] - seeds.loc[0, "pos_x"]) + abs(data[1] - seeds.loc[0, "pos_y"])
    data['seed2_distance'] = abs(data[0] - seeds.loc[1, "pos_x"]) + abs(data[1] - seeds.loc[1, "pos_y"])
    data['seed3_distance'] = abs(data[0] - seeds.loc[2, "pos_x"]) + abs(data[1] - seeds.loc[2, "pos_y"])

    # checks which Manhatten distance is the smallest and assigns according cluster
    for index, row in data.iterrows():
        if (data.loc[index, 'seed1_distance'] < data.loc[index, 'seed2_distance'] and data.loc[
            index, 'seed1_distance'] < data.loc[index, 'seed3_distance']):
            data.loc[index, 'assigned_seed'] = 1
        elif (data.loc[index, 'seed2_distance'] < data.loc[index, 'seed1_distance'] and data.loc[
            index, 'seed2_distance'] < data.loc[index, 'seed3_distance']):
            data.loc[index, 'assigned_seed'] = 2
        elif (data.loc[index, 'seed3_distance'] < data.loc[index, 'seed1_distance'] and data.loc[
            index, 'seed3_distance'] < data.loc[index, 'seed2_distance']):
            data.loc[index, 'assigned_seed'] = 3

    # copy old seeds to check later
    old_seeds = seeds.copy()
    for x in range(3):
        count = 0
        new_x = 0
        new_y = 0
        # calculates new cluster seeds
        for index, row in data.loc[data['assigned_seed'] == x + 1].iterrows():
            count += 1
            new_x += data.loc[index, 0]
            new_y += data.loc[index, 1]
        seeds.loc[x, 'pos_x'] = new_x / count
        seeds.loc[x, 'pos_y'] = new_y / count

    # if old seeds and new seeds are equal, k means clustering is done
    if seeds.equals(old_seeds):
        with open('output.csv', 'a') as file:
            file.write(f'{iteration};;\n')
            file.write(f'{rows};{columns};\n')
        data['assigned_seed'] -= 1
        data['assigned_seed'] = data['assigned_seed'].astype('int32')
        data[['assigned_seed', 0, 1]].to_csv('output.csv', mode='a', index=False, header=False)
        with fileinput.FileInput('output.csv', inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace(',', ';'), end='')
        print('Finished, check output.csv')
        return
    else:
        kmeans(data, n, seeds, iteration + 1)
        return


# get input data
raw_input_data = pd.read_csv('input.csv', sep=';', header=None)

# split input data in data & information like number of clusters
number_of_clusters = raw_input_data.iloc[0, 0]
rows = raw_input_data.iloc[1, 0]
columns = raw_input_data.iloc[1, 1]
input_data = raw_input_data.iloc[2:, :]

# file is not in international format therefore dot based floating point needs to be converted
input_data = input_data.apply(lambda x: x.str.replace(',', '.'))

# get boundaries of original data to generate first seeds
col_x_max = input_data[0].max()
col_x_min = input_data[0].min()
col_y_max = input_data[1].max()
col_y_min = input_data[1].min()

# write first line of the field (and also delete possible old ones by using 'w' mode)
file = open('output.csv', 'w')
file.write(f'{number_of_clusters};;\n')
file.close()

# create seed dataframe to hold the 3 seeds of the clusters, the additional empty column is in order to create the
# needed output for the output.csv which is not a valid csv file
original_seeds = pd.DataFrame(columns=['pos_x', 'pos_y', 'empty'])
for x in range(3):
    original_seeds.loc[x, "pos_x"] = random.uniform(float(col_x_min), float(col_x_max))
    original_seeds.loc[x, "pos_y"] = random.uniform(float(col_y_min), float(col_y_max))

# writes seed data into the file
original_seeds.to_csv('output.csv', mode='a', index=False, header=False)
kmeans(input_data, number_of_clusters, original_seeds, 1)

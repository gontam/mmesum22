# MME2021 - Bioinformatics
# Assignment 3 - K-means
# Date: 08.04.2022
# Latest change: 01.05.2022
# Authors: Johannes Köhrer, Lisa Stefely, Johanna Schachl

import numpy as np
import random as rd


def data_preparation():
    # Read data from 'input.csv'
    input_file2 = open('input.csv')
    data = input_file2.read()
    # Replaces for string modification
    while ',' in data or '\n' in data:
        data = data.replace(',', '.')
        data = data.replace('\n', ';')
    # Transform string to list
    data = data.split(';')
    # Get number of clusters from list
    cluster_num = [int(s) for s in data[0] if s.isdigit()]
    # Remove all non-data points from list
    data = data[4:]
    return data, cluster_num


def generate_matrix(data):
    data_arr = np.ones((18, 2))
    count = 0
    for i in range(0, 18):
        for j in range(0, 2):
            data_arr[i][j] = data[count]
            count = count + 1
    return data_arr


def calculate_new_center(center, cluster):
    if len(cluster) > 0:
        x_sum = 0
        y_sum = 0
        center_new = [0, 0]
        for index in range(0, len(cluster)):
            x_sum += cluster[index][0]
            y_sum += cluster[index][1]
        center_new[0] = x_sum / len(cluster)
        center_new[1] = y_sum / len(cluster)
        return center_new
    else:
        return center


# --- MAIN PROGRAM---

# Get data input and prepare it for clustering
data_input, cluster_number = data_preparation()

# Transform list to 2d matrix
data_array = generate_matrix(data_input)

# Generate one empty matrix for each cluster (0 row, 2 columns)
cluster_1 = np.empty((0, 2))
cluster_2 = np.empty((0, 2))
cluster_3 = np.empty((0, 2))

# Generate center for clusters with random x-,y-coordinates between 0 and 1
center_1 = [rd.random(), rd.random()]
center_2 = [rd.random(), rd.random()]
center_3 = [rd.random(), rd.random()]

# Generate variable to check if position of centers change between two loops
center_1_old = [0, 0]
center_2_old = [0, 0]
center_3_old = [0, 0]

loop_count = 1

while True:
    print(f"Loop run {loop_count} ...")
    loop_count += 1
    for i in range(0, len(data_array)):
        # Calculate distance between cluster center and data point
        dist_center_1 = abs(center_1[0] - data_array[i][0]) + abs(center_1[1] - data_array[i][1])
        dist_center_2 = abs(center_2[0] - data_array[i][0]) + abs(center_2[1] - data_array[i][1])
        dist_center_3 = abs(center_3[0] - data_array[i][0]) + abs(center_3[1] - data_array[i][1])

        # Add data point to cluster with smallest distance
        if (dist_center_1 < dist_center_2) and (dist_center_1 < dist_center_3):
            cluster_1 = np.vstack((cluster_1, data_array[i]))
        elif (dist_center_2 < dist_center_1) and (dist_center_2 < dist_center_3):
            cluster_2 = np.vstack((cluster_2, data_array[i]))
        else:
            cluster_3 = np.vstack((cluster_3, data_array[i]))

    # Calculate new x-, y-coordinates for cluster center 1 to 3
    center_1 = calculate_new_center(center_1, cluster_1)
    center_2 = calculate_new_center(center_2, cluster_2)
    center_3 = calculate_new_center(center_3, cluster_3)

    # Check if the coordinates of the centers are no longer changing
    if (center_1_old == center_1) and (center_2_old == center_2) and (center_3_old == center_3):
        print("\nFinished.\nThe centres are no longer changing and will be saved in 'output.csv'")
        break

    # Save coordinates from centers for comparison
    center_1_old = center_1
    center_2_old = center_2
    center_3_old = center_3

    # Define the cluster arrays new, so they are empty
    cluster_1 = np.empty((0, 2))
    cluster_2 = np.empty((0, 2))
    cluster_3 = np.empty((0, 2))

# Array with results
result_array = np.zeros(((6+len(data_array)), 3))

# Save number of clusters, coordinates of cluster centers, row-, column-length of data
result_array[0][0] = cluster_number[0]
result_array[1][0] = center_1[0]
result_array[1][1] = center_1[1]
result_array[2][0] = center_2[0]
result_array[2][1] = center_2[1]
result_array[3][0] = center_3[0]
result_array[3][1] = center_3[1]
result_array[4][0] = loop_count
result_array[5][0] = data_array.shape[0]
result_array[5][1] = data_array.shape[1]

# Save the data points (x-, y-values)
for row in range(0, len(data_array)):
    result_array[(row + 6)][1] = data_array[row][0]
    result_array[(row + 6)][2] = data_array[row][1]

# Save the cluster number where the data points belongs to
for row in range(0, len(data_array)):
    for j in range(0, len(cluster_1)):
        if data_array[row][0] == cluster_1[j][0] and data_array[row][1] == cluster_1[j][1]:
            result_array[row+6][0] = 1
    for j in range(0, len(cluster_2)):
        if data_array[row][0] == cluster_2[j][0] and data_array[row][1] == cluster_2[j][1]:
            result_array[row+6][0] = 2
    for j in range(0, len(cluster_3)):
        if data_array[row][0] == cluster_3[j][0] and data_array[row][1] == cluster_3[j][1]:
            result_array[row+6][0] = 3

# Write the whole data into "output.csv"
index_count = 0
csvfile = open("output.csv", 'w', newline='')
for row in result_array:
    for column in row:
        if column == 0:
            csvfile.write("")
        elif column == 1 or column == 2 or column == 3:
            csvfile.write('%.0f;' % column)
        else:
            if index_count == 0 or index_count == 4 or index_count == 5:
                csvfile.write('%.0f;' % column)
            else:
                csvfile.write('%.9f;' % column)
    csvfile.write('\n')
    index_count += 1
csvfile.close()

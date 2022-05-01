#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Done by:
#   Andrey Ukhnalyov - me21m018
#   Yusifov Tamerlan - me21m006

import numpy as np
import pandas as pd
import random


def random_centroids(k, p):
    #Make random centroids and return it as list of tuples
    list_of_cent = []
    for i in range(k):
        x_list = [float(x_element) for x_element in p[0]]
        y_list = [float(y_element) for y_element in p[1]]
        x_rc = random.uniform(min(x_list), max(x_list))
        y_rc = random.uniform(min(y_list), max(y_list))
        list_of_cent.append((x_rc, y_rc))
    return list_of_cent


def k_means(centroids_on_start, data):
    data = data.astype('float')
    iterations = 0
    current_centroids = centroids_on_start
    while True:
        groups = []
        iterations += 1
        #Set the groups to the points
        for row in data.iterrows():
            dist1 = abs(row[1][0] - current_centroids[0][0]) + abs(row[1][1] - current_centroids[0][1])
            dist2 = abs(row[1][0] - current_centroids[1][0]) + abs(row[1][1] - current_centroids[1][1])
            dist3 = abs(row[1][0] - current_centroids[2][0]) + abs(row[1][1] - current_centroids[2][1])
            if dist1 < dist2 and dist1 < dist3:
                groups.append([0, row[1][0], row[1][1]])
            elif dist2 < dist1 and dist2 < dist3:
                groups.append([1, row[1][0], row[1][1]])
            elif dist3 < dist1 and dist3 < dist2:
                groups.append([2, row[1][0], row[1][1]])

        #Calculate new centroids as means of each group
        new_centroids = []
        for i in range(3):
            x_sum = 0
            y_sum = 0
            group_size = 0
            for point in groups:
                if point[0] == i:
                    group_size += 1
                    x_sum += point[1]
                    y_sum += point[2]
            x_new_cent = x_sum/group_size
            y_new_cent = y_sum/group_size
            new_centroids.append((x_new_cent, y_new_cent))

        #Stope iterations if centroids do not change any more or use new centroids for the next iteration
        if current_centroids == new_centroids:
            break
        else:
            current_centroids = new_centroids
    
    #Return a tuple of saved variables
    return iterations, current_centroids, groups


def main():
    #Preparing information from input file
    input_file = pd.read_csv('input.csv', sep=';', header=None)
    clusters = input_file.iloc[0, 0]
    clusters_int = int(clusters)
    rows = input_file.iloc[1, 0]
    columns = input_file.iloc[1, 1]
    points = input_file.iloc[2:, :].apply(lambda point: point.str.replace(',', '.'))
    
    #Create randon centroids
    centroids_initial = random_centroids(clusters_int, points)

    #Run k_means algorithm
    number_of_iterations, saved_centroids, groupped_points = k_means(centroids_initial, points)

    #Write output file
    with open('output.csv', mode='w') as file:
        file.write(f'{clusters};;\n')
        for cent in saved_centroids:
            file.write(f'{cent[0]};{cent[1]};\n')
        file.write(f'{number_of_iterations};;\n')
        file.write(f'{rows};{columns};\n')
        for point in groupped_points:
            file.write(f'{point[0]};{point[1]};{point[2]}\n')


if __name__ == '__main__':
    main()


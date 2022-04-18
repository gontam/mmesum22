# Import needed libraries:
import matplotlib.pyplot as plt
import os
import csv
import numpy as np
# import pandas as pd
# from sklearn.cluster import KMeans

# Read in Data (Cristina Soriano):
# ToDo: Read in Data with pandas (as Dataframe): see:
#  https://datatofish.com/import-csv-file-python-using-pandas/
# ToDo: Please delete the character (\n) from data set.
# ToDo: As you can see in, e.g. data[2] you have two values in one string please do two columns e.g.,
#  https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html
# ToDo: Make the data look like the table in excel that I uploaded in our folder.

# @Oliver Jovanović:
# get path of file and directory:
actFile = os.path.join(os.getcwd(), os.listdir(os.getcwd())[0])
csvFile = os.path.dirname(actFile)+'\input.csv'
direct = os.path.dirname(actFile)

# Read in csv file without Pandas (it does not work somehow...)
file = open(csvFile)
reader = csv.reader(file, delimiter=';')
rows = []
for row in reader:
    rows.append(row)
file.close()

# Number of Clusters
# Split clusters by "ï»¿"


def list_to_string(s):
    str1 = " "
    return str1.join(s)


cluster = rows[0]
a = list_to_string(cluster)
clHelp = list(a)
n_cluster = clHelp[3]

# Dimensions:
lst = rows[1]
splits = np.array_split(lst, 2)
data_length = int(splits[0])
dimensions = int(splits[1])

del a, clHelp, cluster, csvFile, actFile, direct, reader, row, rows[0], lst, splits
del rows[0]

# Carve out Data points:
data_points = rows

x = 0
while x <= data_length - 1:
    lst = data_points[x]
    splits = np.array_split(lst, dimensions)
    if x == 0:
        x_point = splits[dimensions - dimensions]
        y_point = splits[dimensions - 1]
    else:
        x_point = np.append(x_point, splits[dimensions - dimensions])
        y_point = np.append(y_point, splits[dimensions - 1])
    x += 1

del rows, splits, data_points, x, lst, file
# Implementation of Algorithm k-means (Oliver Jovanović):
# ToDo: Implement the algorithm using the given data.
# ToDo: Implement it from scratch:
#  https://www.analyticsvidhya.com/blog/2021/04/k-means-clustering-simplified-in-python/
# ToDo: Process the data as most efficient and usable for Marija Toshevska.

# Store Values in output_example.csv (new one) (Marija Toshevska).

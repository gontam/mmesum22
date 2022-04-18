# Import libraries:
import matplotlib.pyplot as plt
import os
import csv
# import pandas as pd
import numpy as np
from random import sample
# from sklearn.cluster import KMeans


# Read in Data (Cristina Soriano):
# Done by Oliver Jovanović:
#  Read in Data with pandas (as Dataframe): see:
#  https://datatofish.com/import-csv-file-python-using-pandas/
#  Please delete the character (\n) from data set.
#  As you can see in, e.g. data[2] you have two values in one string please do two columns e.g.,
#  https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html
#  Make the data look like the table in Excel that I uploaded in our folder.

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
n_cluster = int(clHelp[3])

# Dimensions:
lst = rows[1]
splits = np.array_split(lst, 2)
data_length = int(splits[0])
dimensions = int(splits[1])

# delete unnecessary variables:
del a, clHelp, cluster, csvFile, actFile, direct, reader, row, rows[0], lst, splits
del rows[0]

# Carve out Data points:
data_points = rows
x = 0
x_point = 0
y_point = 0
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

# Change "," in ".":
x = 0
lst = []
while x <= len(x_point) - 1:
    a = x_point[x]
    b = y_point[x]
    a = str(a).replace(',', '.')
    b = str(b).replace(',', '.')
    a = float(a)
    b = float(b)
    x_point[x] = a
    y_point[x] = b
    x += 1

# Cast every string into float:
x_point = list(map(float, x_point))
y_point = list(map(float, y_point))

# Put both points together:
data = np.vstack((np.array(x_point), np.array(y_point)))

# Delete unnecessary variables:
del rows, splits, data_points, x, a, b, lst, file

# Plot the points:
plt.scatter(data[0], data[1], alpha=.5)
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.xlabel('X-Values')
plt.ylabel('Y-Values')
plt.title('Have a look at the dat points')
plt.show()

# Implementation of Algorithm k-means (Oliver Jovanović):
# ToDo: Implement the algorithm using the given data.
# ToDo: Implement it from scratch:
#  https://www.analyticsvidhya.com/blog/2021/04/k-means-clustering-simplified-in-python/
# ToDo: Process the data as most efficient and usable for Marija Toshevska.

# Try it via from sklearn.cluster import KMeans (ModuleNotFoundError: No module named 'sklearn'):
# kmeans = KMeans(n_clusters=n_cluster, random_state=0).fit()

# Implementing it from scratch:

# Step 1: Select value of K:
# It is already selected by n_cluster.

# Step 2: Select random K points which will act as centroids.
np.random.seed(123)
# Select random seeds:
c1, c2, c3 = sample(range(0, 17), n_cluster)
# Centroid 1:
cent1 = np.empty([1, 2])
cent1[0, 0] = data[0, c1]
cent1[0, 1] = data[1, c1]

# Centroid 2:
cent2 = np.empty([1, 2])
cent2[0, 0] = data[0, c2]
cent2[0, 1] = data[1, c2]

# Centroid 3:
cent3 = np.empty([1, 2])
cent3[0, 0] = data[0, c3]
cent3[0, 1] = data[1, c3]

# Plot normal data and random centroids:
plt.scatter(data[0], data[1], c='black')
plt.scatter(cent1[0, 0], cent1[0, 1], c='red')
plt.scatter(cent2[0, 0], cent2[0, 1], c='red')
plt.scatter(cent3[0, 0], cent3[0, 1], c='red')
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.xlabel('X-Values')
plt.ylabel('Y-Values')
plt.title('Have a look at the dat points')
plt.show()


# Store Values in output_example.csv (new one) (Marija Toshevska).

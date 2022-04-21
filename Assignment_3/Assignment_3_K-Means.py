# Import libraries:
import matplotlib.pyplot as plt
import os
import csv
# import pandas as pd
import numpy as np
from random import sample
# from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist


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
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Have a look at the data points')
plt.show()

# Implementation of Algorithm k-means (Oliver Jovanović):
# ToDo: Implement the algorithm using the given data.
# ToDo: Implement it from scratch:
#  https://www.analyticsvidhya.com/blog/2021/04/k-means-clustering-simplified-in-python/
#  or https://anderfernandez.com/en/blog/kmeans-algorithm-python/
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
plt.scatter(cent2[0, 0], cent2[0, 1], c='blue')
plt.scatter(cent3[0, 0], cent3[0, 1], c='green')
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.xlabel('X-')
plt.ylabel('Y')
plt.title('Plot with Centroids')
plt.show()

del c1, c2, c3, sample
# Step 3: Assign all the points to the closest cluster centroid:
# Step 4: Recompute centroids of newly formed clusters:
# Step 5: Repeat Step 3 and 4:

#diff = 1
#j = 0

# Implement KMeans:
# Calculate distances:
x = 0
a = np.empty([1, 2])
tmp = np.empty([data_length, n_cluster])

while x <= data_length - 1:
    a[0, 0] = data[0, x]
    a[0, 1] = data[1, x]
    if(cent1 == a).any():
        dist1 = cdist(cent1, a, metric='cityblock')
        dist2 = cdist(cent2, a, metric='cityblock')
        dist3 = cdist(cent3, a, metric='cityblock')
        tmp[x, 0] = dist1
        tmp[x, 1] = dist2
        tmp[x, 2] = dist3
        x += 1
        pass
    elif(cent2 == a).any():
        dist1 = cdist(cent1, a, metric='cityblock')
        dist2 = cdist(cent2, a, metric='cityblock')
        dist3 = cdist(cent3, a, metric='cityblock')
        tmp[x, 0] = dist1
        tmp[x, 1] = dist2
        tmp[x, 2] = dist3
        x += 1
        pass
    elif(cent3 == a).any():
        dist1 = cdist(cent1, a, metric='cityblock')
        dist2 = cdist(cent2, a, metric='cityblock')
        dist3 = cdist(cent3, a, metric='cityblock')
        tmp[x, 0] = dist1
        tmp[x, 1] = dist2
        tmp[x, 2] = dist3
        x += 1
        pass
    else:
        dist1 = cdist(cent1, a, metric='cityblock')
        dist2 = cdist(cent2, a, metric='cityblock')
        dist3 = cdist(cent3, a, metric='cityblock')
        tmp[x, 0] = dist1
        tmp[x, 1] = dist2
        tmp[x, 2] = dist3
        x += 1

x = 0
cent1 = np.empty([data_length, dimensions])
cent2 = np.empty([data_length, dimensions])
cent3 = np.empty([data_length, dimensions])
# Rearrange the Centroids:
while x < data_length - 1:
    if(tmp[x, 0] == min(tmp[x])):
        cent1[x, dimensions - dimensions] = x_point[x]
        cent1[x, dimensions - 1] = y_point[x]
        cent2[x, dimensions - dimensions] = 0
        cent2[x, dimensions - 1] = 0
        cent3[x, dimensions - dimensions] = 0
        cent3[x, dimensions - 1] = 0
        x += 1
        pass
    elif(tmp[x, 1] == min(tmp[x])):
        cent1[x, dimensions - dimensions] = 0
        cent1[x, dimensions - 1] = 0
        cent2[x, dimensions - dimensions] = x_point[x]
        cent2[x, dimensions - 1] = y_point[x]
        cent3[x, dimensions - dimensions] = 0
        cent3[x, dimensions - 1] = 0
        x += 1
        pass
    elif(tmp[x, 2] == min(tmp[x])):
        cent1[x, dimensions - dimensions] = 0
        cent1[x, dimensions - 1] = 0
        cent2[x, dimensions - dimensions] = 0
        cent2[x, dimensions - 1] = 0
        cent3[x, dimensions - dimensions] = x_point[x]
        cent3[x, dimensions - 1] = y_point[x]
        x += 1
        pass
# Calculate mean:
# Centroid 1:
m1 = np.empty([1, 2])
m1[0, 0] = 0
m1[0, 1] = 0
n1 = 0
c1 = np.empty([1, 2])
c1[0, 0] = 0
c1[0, 1] = 0
i1 = 0
while i1 <= data_length - 1:
    if(cent1[i1] == 0).any():
        i1 += 1
        pass
    else:
        c1[0, 0] += float(cent1[i1, 0])
        c1[0, 1] += float(cent1[i1, 1])
        n1 += 1
        i1 += 1
        pass

# Centroid 2:
m2 = np.empty([1, 2])
m2[0, 0] = 0
m2[0, 1] = 0
n2 = 0
c2 = np.empty([1, 2])
c2[0, 0] = 0
c2[0, 1] = 0
i2 = 0
while i2 <= data_length - 1:
    if(cent2[i2] == 0).any():
        i2 += 1
        pass
    else:
        c2[0, 0] += float(cent2[i2, 0])
        c2[0, 1] += float(cent1[i2, 1])
        n2 += 1
        i2 += 1
        pass

# Centroid 3:
m3 = np.empty([1, 2])
m3[0, 0] = 0
m3[0, 1] = 0
n3 = 0
c3 = np.empty([1, 2])
c3[0, 0] = 0
c3[0, 1] = 0
i3 = 0
while i3 <= data_length - 1:
    if(cent3[i3] == 0).any():
        i3 += 1
        pass
    else:
        c3[0, 0] += float(cent3[i3, 0])
        c3[0, 1] += float(cent3[i3, 1])
        n3 += 1
        i3 += 1
        pass

# Calculate mean:
m1 = c1 / n1
m2 = c2 / n2
m3 = c3 / n3

# Store data points to centroids:
c1_x = 0
c1_y = 0
c2_x = 0
c2_y = 0
c3_x = 0
c3_y = 0
i1 = 0
i2 = 0
i3 = 0
i = 0
# Centroid 1:
while i1 <= data_length - 1:
    if(cent1[i1] == 0).any():
        i1 += 1
        pass
    else:
        c1_x = np.append(c1_x, float(cent1[i1, 0]))
        c1_y = np.append(c1_y, float(cent1[i1, 1]))
        i1 += 1
        i += 1
        pass


# Centroid 2:
i = 0
while i2 <= data_length - 1:
    if(cent2[i2] == 0).any():
        i2 += 1
        pass
    else:
        c2_x = np.append(c2_x, float(cent2[i2, 0]))
        c2_y = np.append(c2_y, float(cent2[i2, 1]))
        i2 += 1
        i += 1
        pass


# Centroid 3:
i = 0
while i3 <= data_length - 1:
    if(cent3[i3] == 0).any():
        i3 += 1
        pass
    else:
        c3_x = np.append(c3_x, float(cent3[i3, 0]))
        c3_y = np.append(c3_y, float(cent3[i3, 1]))
        i3 += 1
        i += 1
        pass

c1 = np.vstack((np.array(c1_x), np.array(c1_y)))
c2 = np.vstack((np.array(c2_x), np.array(c2_y)))
c3 = np.vstack((np.array(c3_x), np.array(c3_y)))

del c1_x, c1_y, c2_x, c2_y, c3_x, c3_y
# Plot new Centroids:
plt.scatter(data[0], data[1], c='black')
plt.scatter(m1[0, 0], m1[0, 1], c='red')
plt.scatter(c1[0], c1[1], c='red', alpha=.6)
plt.scatter(m2[0, 0], m2[0, 1], c='blue')
plt.scatter(c2[0], c2[1], c='blue', alpha=.6)
plt.scatter(m3[0, 0], m3[0, 1], c='green')
plt.scatter(c3[0], c3[1], c='green', alpha=.6)
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.xlabel('X-')
plt.ylabel('Y')
plt.title('Plot with new Centroids and data points')
plt.show()

# Store Values in output_example.csv (new one) (Marija Toshevska).

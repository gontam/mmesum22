# Project done and handed in by Cristiana Soriano, Marija Toshevksa, Oliver Jovanović
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
import math

# Read in Data:
reader = pd.read_csv('input.csv', header=None, delimiter=';')

# Change commas by dots and cast to numeric values type "float":
reader[0]=reader[0].str.replace(',','.')
reader[1]=reader[1].str.replace(',','.')
reader[0] = pd.to_numeric(reader[0], downcast="float")
reader[1] = pd.to_numeric(reader[1], downcast="float")

# Carve out Data points:
data = reader.iloc[2:len(reader)]

# Plot the points:
plt.scatter(data[0], data[1], alpha=.5)
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Have a look at the data points')
plt.show()

# Implementing KMeans from scratch:
# Step 1: Select value of K:
n_cluster = reader[0][0]
n_cluster = int(n_cluster)

# Step 2: Select random K points which will act as centroids.
i = 0
centroids = []

while i < n_cluster:
  np.random.seed(i)
  array = np.random.rand(2)
  centroids.append(array)
  i += +1

# Plot normal data and random centroids:
plt.scatter(data[0], data[1], c='black')
plt.scatter(centroids[0][0], centroids[0][1], c='red')
plt.scatter(centroids[1][0], centroids[1][1], c='blue')
plt.scatter(centroids[2][0], centroids[2][1], c='green')
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.xlabel('X-')
plt.ylabel('Y')
plt.title('Plot with Centroids')
plt.show()

# Calculate distances:
x = 0
a = np.empty([1, 2])
# Length of data:
data_length = data.shape[0]
dist = []
tmp = np.empty([data_length, n_cluster])
distances = []

while x <= data_length - 1:
  point = np.matrix(data.iloc[x],dtype='float64')
  dist1 = cdist(np.matrix(centroids[0]), point, metric='cityblock')
  dist2 = cdist(np.matrix(centroids[1]), point, metric='cityblock')
  dist3 = cdist(np.matrix(centroids[2]), point, metric='cityblock')
  dist.append(dist1)
  dist.append(dist2)
  dist.append(dist3)
  distances.append(dist)
  dist = []
  x += 1
  clusters = []
  for i in range(0, len(distances)):
      clusters.append(distances[i].index(min(distances[i])))

  # Plot clusters:
  for i in range(0, len(clusters)):
      if clusters[i] == 0:
          plt.scatter(data[0][i + 2], data[1][i + 2], c='red')
      elif clusters[i] == 1:
          plt.scatter(data[0][i + 2], data[1][i + 2], c='blue')
      elif clusters[i] == 2:
          plt.scatter(data[0][i + 2], data[1][i + 2], c='green')
  plt.scatter(centroids[0][0], centroids[0][1], c='red', marker='>')
  plt.scatter(centroids[1][0], centroids[1][1], c='blue', marker='>')
  plt.scatter(centroids[2][0], centroids[2][1], c='green', marker='>')
  plt.xlim([0, 1])
  plt.ylim([0, 1])
  plt.xlabel('X-')
  plt.ylabel('Y')
  plt.title('Cluster plots')
  plt.show()

clusters = []

for i in range(0,len(distances)):
  clusters.append(distances[i].index(min(distances[i])))

# Plot clusters:
for i in range(0, len(clusters)):
  if clusters[i] == 0:
    plt.scatter(data[0][i+2],data[1][i+2],c='red')
  elif clusters[i] == 1:
    plt.scatter(data[0][i+2],data[1][i+2],c='blue')
  elif clusters[i] == 2:
    plt.scatter(data[0][i+2],data[1][i+2],c='green')

plt.scatter(centroids[0][0],centroids[0][1], c='red', marker='>')
plt.scatter(centroids[1][0],centroids[1][1], c='blue', marker='>')
plt.scatter(centroids[2][0],centroids[2][1], c='green', marker='>')
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.xlabel('X-')
plt.ylabel('Y')
plt.title('Cluster plots')
plt.show()


maxerror = 10
iteration = 0
while maxerror > 0.5 and iteration < 10:
  clusters = []
  for i in range(0, len(distances)):
    clusters.append(distances[i].index(min(distances[i])))
  # Plot clusters:
  for i in range(0,len(clusters)):
    if clusters[i] == 0:
      plt.scatter(data[0][i+2],data[1][i+2], c='red')
    elif clusters[i] == 1:
      plt.scatter(data[0][i+2],data[1][i+2], c='blue')
    elif clusters[i] == 2:
      plt.scatter(data[0][i+2],data[1][i+2], c='green')
  plt.scatter(centroids[0][0],centroids[0][1], c='red', marker='>')
  plt.scatter(centroids[1][0],centroids[1][1], c='blue', marker='>')
  plt.scatter(centroids[2][0],centroids[2][1], c='green', marker='>')
  plt.xlim([0, 1])
  plt.ylim([0, 1])
  plt.xlabel('X-')
  plt.ylabel('Y')
  plt.title('Cluster plots')
  plt.show()

  # Recalculate centroids
  aux_01 = 0
  aux_11 = 0
  aux_02 = 0
  aux_12 = 0
  aux_03 = 0
  aux_13 = 0
  n_1 = 0
  n_2 = 0
  n_3 = 0
  for i in range(0, len(clusters)):
    if clusters[i] == 0:
      aux_01 = data[0][i+2]+aux_01
      aux_11 = data[1][i+2]+aux_11
      n_1 += 1
    elif clusters[i] == 1:
      aux_02 = data[0][i+2]+aux_02
      aux_12 = data[1][i+2]+aux_12
      n_2 += 1
    elif clusters[i] == 2:
      aux_03 = data[0][i+2]+aux_03
      aux_13 = data[1][i+2]+aux_13
      n_3 += 1
  oldcentroids = centroids
  centroids[0][0] = aux_01/n_1
  centroids[0][1] = aux_11/n_1
  centroids[1][0] = aux_02/n_2
  centroids[1][1] = aux_12/n_2
  centroids[2][0] = aux_03/n_3
  centroids[2][1] = aux_13/n_3
  errortotal = cdist(np.matrix(centroids), np.matrix(oldcentroids), metric='cityblock')

  for i in range(0,len(errortotal)):
    aux = max(errortotal[i])
    if aux < maxerror:
      maxerror = aux
      iteration = 0
    if aux == maxerror:
      iteration = iteration + 1
  print(maxerror)

# Inserting the cluster number before entries:
data.insert(0, "cluster", clusters)
df1 = pd.DataFrame({int(n_cluster)})
df2 = pd.DataFrame(centroids)
df3 = pd.DataFrame({math.floor(iteration)})
df4 = pd.DataFrame(reader.iloc[1:2])
df5 = pd.DataFrame(data)

# Concatenating dataframes:
vertical_concat = pd.concat([df1, df2, df3, df4, df5], axis=0)

# Writing in output.csv
vertical_concat.to_csv('output.csv')
# Import needed libraries:
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans

# Read in Data (Cristina Soriano):
# ToDo: Read in Data with pandas (as Dataframe): see:
#  https://datatofish.com/import-csv-file-python-using-pandas/
# ToDo: Please delete the character (\n) from data set.
# ToDo: As you can see in, e.g. data[2] you have two values in one string please do two columns e.g.,
#  https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html
# ToDo: Make the data look like the table in excel that I uploaded in our folder.
file = open('input.csv')
print(file)
data = []
for i in file:
    data.append(i)

cluster_num = data[0][3]

data[1] = data[1].split(';')
rows = data[1][0]
columns = data[1][1]

# Implementation of Algorithm k-means (Oliver Jovanović):
# ToDo: Implement the algorithm using the given data.
# ToDo: Process the data as most efficient and usable for Marija Toshevska.

# Store Values in output_example.csv (new one) (Marija Toshevska).

# Import needed libraries:
#import matplotlib.pyplot as plt
#import pandas as pd
import os
import csv
#from sklearn.cluster import KMeans

# Read in Data (Cristina Soriano):
# ToDo: Read in Data with pandas (as Dataframe): see:
#  https://datatofish.com/import-csv-file-python-using-pandas/
# ToDo: Please delete the character (\n) from data set.
# ToDo: As you can see in, e.g. data[2] you have two values in one string please do two columns e.g.,
#  https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html
# ToDo: Make the data look like the table in excel that I uploaded in our folder.

# @Oliver Jovanovic:
# get path of file and directory:
actFile = os.path.join(os.getcwd(), os.listdir(os.getcwd())[0])
csvFile = os.path.dirname(actFile)+'\input.csv'
direct = os.path.dirname(actFile)

# Read in csv file without Pandas (it does not work somehow...)
file = open(csvFile)
print(file)
reader = csv.reader(file)
rows = []
for row in reader:
    rows.append(row)

print(rows)


# Cristina Soriano:
#file = open('input.csv')
#print(file)
#data = []
#for i in file:
#    data.append(i)

#cluster_num = data[0][3]

#data[1] = data[1].split(';')
#rows = data[1][0]
#columns = data[1][1]

# Implementation of Algorithm k-means (Oliver Jovanović):
# ToDo: Implement the algorithm using the given data.
# ToDo: Implement it from scratch:
#  https://www.analyticsvidhya.com/blog/2021/04/k-means-clustering-simplified-in-python/
# ToDo: Process the data as most efficient and usable for Marija Toshevska.

# Store Values in output_example.csv (new one) (Marija Toshevska).

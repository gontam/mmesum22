Assignment 3 - K-means

Deadline: 24.04 23:59

For this assignment you are asked to implement the k-means clustering for a multi-dimensional data set.
The program should take an input file 'input.csv'.

The file is structured in the following way:

[number of clusters]
[number of data entries (rows)] [number of dimensions (columns)]
[matrix of the format: data entries * dimensions]

The code should create a file 'output.csv' with the following information in it:
[number of clusters] 
[list of the positions of the seed points]
[number of iterations]
[number of data entries] [number of dimensions]
[number of the corresponding cluster  + matrix of the size of entries * dimensions]

The seed points should be generated randomly within the max and min boundaries of each dimension. 
The number of iterations is the number of steps until the cluster centers do not change their position anymore. 

For the distance function use the Manhattan Distance.

Unless stated otherwise:
Every student must sumbit the 'A3_Lastname_Name.py'.
Please provide the names of both individuals who worked on the project in the first comment line of the file.
Additionally, please provide a link to your branch in the comments.

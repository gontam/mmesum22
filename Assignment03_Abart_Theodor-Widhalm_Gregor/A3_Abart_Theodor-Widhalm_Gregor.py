# Gregor Widhalm, me21m012
# Theodor Abart, me21m011

import csv
import numpy as np
import matplotlib.pyplot as plt
import random


def readInputFile(filename):

    # open csv file
    with open(filename, mode='r', newline='', encoding='utf-8-sig') as csvfile:
        # open reader, semicolon as delimiter
        csvreader = csv.reader(csvfile, delimiter=";")

        # first line --> number of clusters
        numberOfClusters = int(next(csvreader)[0])
        # second line --> number of rows + number of columns
        numberOfRows_numberOfColumns = next(csvreader)

        for x in range(numberOfRows_numberOfColumns.__len__()):
            numberOfRows_numberOfColumns[x] = int(numberOfRows_numberOfColumns[x])

        data = list(csvreader)

        # create float array
        floatArray = np.empty([data.__len__(), 2])

        lineCounter = 0
        # loop through all data lines
        for lineList in data:
            # loop through all data columns
            for y in range(len(lineList)):
                # append float to line of list --> each row should in the end
                # include as many columns of floats as the input file
                floatArray[lineCounter][y] = float(lineList[y].replace(',', '.'))
            lineCounter += 1
            # add a new line (list) as long as the line counter is smaller
            # than the number of rows of the input data
            # if lineCounter < data.__len__():
            # floatList.append(list())
    return numberOfClusters, np.array(numberOfRows_numberOfColumns), floatArray


def writeOutputFile(filename, numClusters, numRows_numColumns, centroids, iterationCounter, dataArray):
    with open(filename, mode='w', newline='', encoding='utf-8-sig') as file:
        writer = csv.writer(file, delimiter=";")
        # write number of clusters
        writer.writerow([str(numClusters)])

        # write coordinates of centroids
        for i in range(numClusters):
            writer.writerow([str(centroids[i][0]), str(centroids[i][1])])

        # write number of iterations
        writer.writerow(str(iterationCounter))

        # write number of entries and dimensions
        writer.writerow([str(numRows_numColumns[0]), str(numRows_numColumns[1])])

        for j in range(numRows_numColumns[0]):
            writer.writerow([str(int(dataArray[j, 2])), str(dataArray[j, 0]), str(dataArray[j, 1])])


def manhattenDist(a, b):
    # taken from https://www.statology.org/manhattan-distance-python/
    return sum(abs(val1 - val2) for val1, val2 in zip(a, b))


def plotResults(data, numClusters):
    # !!! PLOT SECTION !!!
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(dataArray[:, 0], data[:, 1], c='c')

    markerList = ['s', '*', 'X', 'p', 'P', 'H', 'D']
    colorList = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

    for j in range(numClusters):
        mask = (data[:, 2] == j)
        condArrayPlot = data[mask, :]

        ax1.scatter(condArrayPlot[:, 0], condArrayPlot[:, 1], c=colorList[j])
        ax1.scatter(centroids[j][0], centroids[j][1], c=colorList[j], marker=markerList[j])

    plt.show()


if __name__ == '__main__':

    # read the input file
    [numClusters, numRows_numColumns, dataArray] = readInputFile('input.csv')

    # define K random centroids
    centroids = np.empty([numClusters, 2])
    for i in range(numClusters):
        centroids[i] = np.array([random.uniform(min(dataArray[:, 0]), max(dataArray[:, 0])),
                                 random.uniform(min(dataArray[:, 1]), max(dataArray[:, 1]))])

    # add column to array and set it to default: -1
    dataArray = np.append(dataArray, np.ones([int(numRows_numColumns[0]), 1]), axis=1)
    dataArray[:, 2] = -1

    iterationCounter = 0
    while True:
        iterationCounter += 1

        # assignment to cluster according to minimal distance to centroid
        distToCentr = np.empty([numClusters])
        for j in range(dataArray.shape[0]):
            for k in range(centroids.shape[0]):
                distToCentr[k] = manhattenDist(dataArray[j, :], centroids[k, :])
            # assign to cluster where distance is minimal
            dataArray[j, 2] = np.argmin(distToCentr)

        centroids_prev = np.copy(centroids)
        # calculate new centroids
        for l in range(numClusters):
            # extract x and y values of specific cluster and save to condArray
            # in order to calc. mean of x and y for new centroid of cluster
            mask = (dataArray[:, 2] == l)
            condArray = dataArray[mask, :]
            for m in range(2):
                centroids[l, m] = np.mean(condArray[:, m])

        # in case of non-changing centroids --> DONE --> exit loop
        if np.array_equal(centroids_prev, centroids):
            break

    writeOutputFile('output.csv', numClusters, numRows_numColumns, centroids, iterationCounter, dataArray)

    plotResults(dataArray, numClusters)

    print('done.')

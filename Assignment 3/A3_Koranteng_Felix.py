#Collaborators: Raphael Hein & Felix Koranteng
#

import csv
import numpy as np
import random
from scipy.spatial import distance
import matplotlib.pyplot as plt


# ------------------------------- Class information to store data and run K-means  -------------------------------
class information:
    k = 0
    n = 0
    dim = 0
    datapoints = []
    cluster = [0] * 18
    centroids = []
    counter = 0

    # ------------------------------- Loading Input file CSV -------------------------------
    def parse_input(self):
        with open('input.csv', 'r', encoding='utf-8-sig', newline='') as file:
            data = csv.reader(file, delimiter=';')
            i = 0
            for line in data:
                if i == 0:
                    self.k = int(line[0])
                elif i == 1:
                    self.n = int(line[0])
                    self.dim = int(line[1])
                else:
                    x = float(line[0].replace(',', '.'))
                    y = float(line[1].replace(',', '.'))
                    self.datapoints.append(tuple((x, y)))
                i += 1
        return

    # ------------------------------- Writing results to output file (.csv) -------------------------------
    def parse_output(self):
        filename = 'output.csv'
        with open(filename, mode='w', newline='', encoding='utf-8-sig') as file:
            writer = csv.writer(file, delimiter=';')
            writer.writerow((str(int(self.k)), '', ''))

            for i in range(self.k):
                writer.writerow(
                    (str(round(float(self.centroids[i][0]), 9)), str(round(float(self.centroids[i][1]), 9)), ''))

            writer.writerow((str(int(self.counter)), '', ''))
            writer.writerow((str(int(self.n)), str(int(self.dim))))

            for i in range(self.n):
                writer.writerow(
                    (str(int(self.cluster[i])), str(float(self.datapoints[i][0])), str(float(self.datapoints[i][1]))))
        return

    # ------------------------------- Finding the closet centroid within the data (Manhattan distance = Cityblock) -------------------------------
    def findClosestCentroids(self):
        idx = 0
        for point in self.datapoints:
            dists = []
            for centr in self.centroids:
                dists.append(distance.cityblock(point, centr))
            self.cluster[idx] = (np.argmin(dists))
            idx += 1
        return

    # ------------------------------- Moving centroids based on the mean of each data point -------------------------------
    def calc_centroids(self):
        for clust in range(len(self.centroids)):
            sumx, sumy = 0, 0
            count = 0
            for idx in range(len(self.datapoints)):
                if self.cluster[idx] == clust:
                    count += 1
                    sumx += self.datapoints[idx][0]
                    sumy += self.datapoints[idx][1]
            if count == 0:
                count += 1
            new_cent = [sumx / count, sumy / count]
            self.centroids[clust] = np.array(new_cent)
        return

    # ------------------------------- K-means clustering -------------------------------
    def k_means(self):
        centr = []

        # ---------------------  initializing k random centroids  ---------------------
        for _ in range(self.k):
            x = random.uniform(min(self.datapoints[:][0]), max(self.datapoints[:][0]))
            y = random.uniform(min(self.datapoints[:][1]), max(self.datapoints[:][1]))
            centr.append(tuple((x, y)))
        self.centroids = np.array(centr)
        self.centroids.sort(axis=0)
        cnt = 0

        while (True):
            cnt += 1
            prev_cent = np.copy(self.centroids)
            self.findClosestCentroids()
            self.calc_centroids()
            comparison = prev_cent == self.centroids

            # ---------------------  termination condition no change in centroids   --------------------
            if comparison.all():
                print("Iterations: ", cnt)
                self.counter = cnt
                break
        return

    # ------------------------------- k-means most probable result (track result of n k-means runs) -------------------------------
    def k_means_prob(self, n):
        prob = [[0 for col in range(self.k)] for row in range(self.n)]
        for i in range(n):
            self.k_means()
            for i in range(self.n):
                prob[i][self.cluster[i]] += 1
        tracker = 0
        for point in prob:
            maxval = max(point)
            maxidx = point.index(maxval)
            self.cluster[tracker] = maxidx
            tracker += 1
            print("probs of ", tracker, " = ", point)

    # ------------------------------- Visual output (plotting) -------------------------------
    def plot_result(self):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        for j in range(self.n):
            if self.cluster[j] == 0:
                ax1.scatter(self.datapoints[j][0], self.datapoints[j][1], c='b')
            elif self.cluster[j] == 1:
                ax1.scatter(self.datapoints[j][0], self.datapoints[j][1], c='r')
            else:
                ax1.scatter(self.datapoints[j][0], self.datapoints[j][1], c='g')
        plt.savefig('scatter_plot')
        return

# ------------------------------- Initializing class   -------------------------------
trial = information()
trial.parse_input()

# ------------------------------- Run K-means  -------------------------------
trial.k_means()
trial.parse_output()

# ------------------------------- Most likely, n runs -------------------------------
trial.k_means_prob(100000)
trial.parse_output()

# ------------------------------- Visual output of the results -------------------------------
trial.plot_result()

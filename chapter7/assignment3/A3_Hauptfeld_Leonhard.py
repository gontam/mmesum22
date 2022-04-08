# Leonhard Hauptfeld, Christopher Schmidl, 
# https://github.com/gontam/mmesum22/tree/hauptfeld/chapter7/assignment3
import argparse
import csv
import numpy as np
import matplotlib.pyplot as plt

CLUSTER_NAMES = ["Hansi", "Schurli", "Peppi"]
CLUSTER_COLORS = ["red", "green", "blue"]

class Point:
    def __init__(self, position, cluster=None):
        self.position = position
        self.cluster = cluster
    def assign_cluster(self, cluster):
        # Set both the cluster to the point as well as the point to the cluster
        # for easier access
        self.cluster = cluster
        self.cluster.points.append(self)
    def get_dimension_count(self):
        return len(self.position)
class Cluster:
    def __init__(self, seed_point):
        self.id = id
        self.seed_point = seed_point
        self.points = []
        # Initialize the center point by the seed point
        self.point = Point(self.seed_point.position)
    def reset_points(self):
        self.points = []
    def recalculate(self):
        # Save the old point so we can calculate a difference to it later
        previous_point = Point(self.point.position)
        # Set the new center point position by calculating the mean of all the assigned points
        self.point.position = [np.mean([p.position[d] for p in self.points]) for d in range(0, self.point.get_dimension_count())]
        # Return the distance between the new and old center
        return manhattan_distance(previous_point, self.point, self.point.get_dimension_count())

def generate_clusters(points, cluster_count, dimension_count):
    # Calculate the maximum and minimum values per dimension, in a hugely inefficient way
    dimension_boundaries = [(np.min([p.position[d] for p in points]), np.max([p.position[d] for p in points])) for d in range(0, dimension_count)]
    # Create {cluster_count} number of cluster points, randomly between dimension min/max
    clusters = []
    for i in range(0, cluster_count):
        clusters.append(Cluster(Point([np.random.uniform(dimension_boundaries[d][0], dimension_boundaries[d][1]) for d in range(0, dimension_count)])))
    return clusters

def manhattan_distance(a, b, dimension_count=2):
    return np.sum([abs(a.position[d] - b.position[d]) for d in range(0, dimension_count)])

def k_means(points, clusters, dimension_count=2):
    position_diff = 0
    # First, reset all the point assignments
    for c in clusters:
        c.reset_points()
    # Then, assign all the points to their nearest clusters
    for p in points: 
        p.assign_cluster(clusters[np.argmin([manhattan_distance(p, c.point, dimension_count) for c in clusters])])
    # Finally, recalculate all cluster center points, sum up the positional differences
    for c in clusters: 
        position_diff += c.recalculate()
    return position_diff

def plot_clusters(clusters):
    c_idx = 0
    for c in clusters:
        # Plot all the cluster points
        plt.scatter([p.position[0] for p in c.points], [p.position[1] for p in c.points], color=CLUSTER_COLORS[c_idx])
        # Plot and annotate the cluster center point in different style
        plt.scatter([c.point.position[0]], [c.point.position[1]], color=CLUSTER_COLORS[c_idx], marker="D", edgecolors="black")
        plt.annotate("  " + CLUSTER_NAMES[c_idx], (c.point.position[0], c.point.position[1]))
        c_idx += 1
    plt.show()
def save_clusters(clusters, iteration_count, dimension_count, filename):
    with open(filename, mode='w', newline='', encoding='utf-8-sig') as csvfile:
        writer = csv.writer(csvfile, delimiter=';', quotechar='|')
        # Write the number of clusters
        writer.writerow([len(clusters)])
        # Write seed point positions
        for c in clusters: writer.writerow(c.seed_point.position)
        # Write the number of iterations
        writer.writerow([iteration_count])
        # Write the total number of points and dimensions
        writer.writerow([np.sum([len(c.points) for c in clusters]), dimension_count])
        # Write all the cluster's points
        c_idx = 0
        for c in clusters:
            for p in c.points:
                writer.writerow([c_idx] + p.position)
            c_idx += 1


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description='Run K-means on a CSV input file, output to CSV and/or show plot')
    parser.add_argument('--in-file', type=str, default="input.csv", help='CSV input file')
    parser.add_argument('--out-file', type=str, default="output.csv", help='CSV output file')
    parser.add_argument('--out-plot', action='store_true', help="Show plot")
    args = parser.parse_args()
    # Parse input file
    cluster_count = None
    dimension_count = None
    points = []
    with open(args.in_file, mode='r', newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile, delimiter=';', quotechar='|')
        row_index = 0
        for row in reader:
            if row_index == 0:
                cluster_count = int(row[0])
            elif row_index == 1:
                dimension_count = int(row[1])
            else:
                points.append(Point([float(row[d].replace(',', '.')) for d in range(0, dimension_count)]))
            row_index += 1
    # Show a summary
    print(f"{dimension_count} dimensions, {cluster_count} clusters, {len(points)} points")
    # Generate clusters with random seed points
    clusters = generate_clusters(points, cluster_count, dimension_count)
    # Run the loop until the cluster points don't change anymore
    diff = 1
    iteration_count = 0
    while diff != 0:
        diff = k_means(points, clusters, dimension_count)
        iteration_count += 1
    # Pipe to specified outputs
    if args.out_file: save_clusters(clusters, iteration_count, dimension_count, args.out_file)
    if args.out_plot: plot_clusters(clusters)
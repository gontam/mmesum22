import argparse
import csv
import numpy as np

CLUSTER_NAMES = ["Herbert", "Peter", "Schurli", "Peppi"]

class Point:
    def __init__(self, position, cluster=None):
        self.position = position
        self.cluster = cluster
class ClusterPoint(Point):
    def __init__(self, position, id):
        super().__init__(position)
        self.id = id

def generate_clusters():
    pass

def manhattan_distance(a, b, dimension_count=2):
    return np.sum([abs(a.position[d] - b.position[d]) for d in range(0, dimension_count)])

def k_means(points, clusters, dimension_count=2):
    for p in points:
        p.cluster = clusters[np.argmin([manhattan_distance(p, c, dimension_count) for c in clusters])]

if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description='Parse a FASTA formatted sequence and write the largest open reading frame to a file')
    parser.add_argument('--input', type=str, default="input.csv", help='CSV output file')
    parser.add_argument('--output', type=str, default="output.csv", help='CSV input file')
    args = parser.parse_args()
    # Parse input file
    cluster_count = None
    dimension_count = None
    points = []
    with open(args.input, mode='r', newline='', encoding='utf-8-sig') as csvfile:
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
    # Create {dimension_count} number of cluster points
    clusters = [ClusterPoint([0,0], i) for i in range(0, cluster_count)]
    # Run the loop
    k_means(points, clusters, dimension_count)
    print(", ".join([CLUSTER_NAMES[p.cluster.id] for p in points]))
    #for i in range(0, 100):
    #    pass
import argparse
import numpy as np
from sklearn import cluster, preprocessing

parser = argparse.ArgumentParser()
parser.add_argument('-i', help = 'filename');
parser.add_argument('-n', help = 'num_clusters', type = int);
parser.add_argument('-t', help = 'threshold', type = float, default = 0.0001);
parser.add_argument('-d', help = 'debug', action = 'store_true');

args = parser.parse_args()

data = np.loadtxt(args.i)
data = np.delete(data, 0, 1)

#data = preprocessing.normalize(data, norm = 'max', axis = 0);

clus = cluster.KMeans(n_clusters = args.n, verbose = args.d, tol = args.t, random_state = 42, precompute_distances = True, copy_x = False, algorithm = 'full')
#clus = cluster.Birch(n_clusters = args.n, threshold = 0.2, branching_factor = 25)
#clus = cluster.AgglomerativeClustering(n_clusters = args.n)

clus.fit(data)

np.savetxt(args.i + '.membership', clus.labels_, fmt = '%d')

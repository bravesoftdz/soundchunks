import argparse
import numpy as np
from sklearn import cluster, preprocessing

parser = argparse.ArgumentParser()
parser.add_argument('-i', help = 'filename');
parser.add_argument('-n', help = 'num_clusters', type = int);
parser.add_argument('-t', help = 'threshold', type = float, default = 0.0001);
parser.add_argument('-d', help = 'debug', action = 'store_true');
parser.add_argument('-a', help = 'alternate_method', action = 'store_true');

args = parser.parse_args()

data = np.loadtxt(args.i)
data = np.delete(data, 0, 1)

#data = preprocessing.normalize(data, norm = 'max', axis = 0);

if args.a:
    clus = cluster.AgglomerativeClustering(n_clusters = args.n)
else:
    clus = cluster.KMeans(n_clusters = args.n, verbose = args.d, tol = args.t, random_state = 42, precompute_distances = True, copy_x = False, algorithm = 'full')

#clus = cluster.SpectralClustering(n_clusters = args.n, assign_labels = 'discretize', gamma = 0.15, random_state = 42, eigen_tol = args.t)

clus.fit(data)

np.savetxt(args.i + '.membership', clus.labels_, fmt = '%d')

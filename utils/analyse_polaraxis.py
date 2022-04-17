import sys
import numpy as np
from sklearn.cluster import AgglomerativeClustering


data_1d = np.loadtxt("molecule_polaraxis.dat").T
data0 = np.array([np.array([x, 1]) for x in data_1d[0]])
data1 = np.array([np.array([x, 1]) for x in data_1d[1]])
data2 = np.array([np.array([x, 1]) for x in data_1d[2]])
data = np.array([data0, data1, data2])

caxis = int(np.loadtxt("caxis"))
dz = []
d1 = np.array([x[0] for x in data[caxis] if x[0] < 90])
d2 = np.array([x[0] for x in data[caxis] if x[0] > 90])
dz = [d1, d2]
# distance_threshold = 60
cluster = AgglomerativeClustering(n_clusters=3, affinity='euclidean', linkage='ward')
# cluster = AgglomerativeClustering(n_clusters=None, distance_threshold=distance_threshold, affinity='euclidean', linkage='complete')
dxy = []
for i in range(3):
    if i != caxis:
        cluster.fit(data[i])
        # print(cluster.labels_)
        d1 = np.array([data[i][k][0] for k in range(len(data[i])) if cluster.labels_[k] == 0])
        d2 = np.array([data[i][k][0] for k in range(len(data[i])) if cluster.labels_[k] == 1])
        d3 = np.array([data[i][k][0] for k in range(len(data[i])) if cluster.labels_[k] == 2])
        dxy.append([d1, d2, d3])

print("%10.3f %10.3f   %10.3f %10.3f         %10.3f %10.3f %10.3f   %10.3f %10.3f %10.3f"%(np.average(dz[0]), np.average(dz[1]), np.std(dz[0]), np.std(dz[1]), np.average(dxy[0][0]), np.average(dxy[0][1]), np.average(dxy[0][2]), len(dxy[0][0])/(len(dxy[0][0])+len(dxy[0][1])+len(dxy[0][2])), len(dxy[0][1])/(len(dxy[0][0])+len(dxy[0][1])+len(dxy[0][2])), len(dxy[0][2])/(len(dxy[0][0])+len(dxy[0][1])+len(dxy[0][2]))))

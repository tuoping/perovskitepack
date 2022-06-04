import sys
import numpy as np
from sklearn.cluster import AgglomerativeClustering
of1 = open("pbpbdistances.out", "w")
of2 = open("pbpbdistances-sorted.out", "w")
of1.write("T  average_135  average_45  std_135  std_45  size_domain_0 size_domain_1\n")
of2.write("T  average_short  average_long  std_short  std_long  size_domain_0 size_domain_1\n")
of3 = open("pbpbratio.out", "w")
of3.write("T  average_135/c  average_45/c  std_135/c  std_45/c  c\n")
for T in [i for i in range(20, 155, 5)]+[i for i in range(155, 440,10)]:
    dirname = None
    if T < 155:
        dirname = "annealT"+str(T)
    else:
        dirname = "T"+str(T)
    print(dirname)

    data1_1d = np.loadtxt(dirname+"/left_distances").T
    data2_1d = np.loadtxt(dirname+"/right_distances").T
    data3_1d = np.loadtxt(dirname+"/c_distances")
    # data = np.array([sorted([data1_1d[i], data2_1d[i]]) for i in range(len(data1_1d))])
    data = np.array([[data1_1d[i], data2_1d[i]] for i in range(len(data1_1d))]).T
    ratio = np.array([[data1_1d[i]/data3_1d[i]/np.sqrt(2), data2_1d[i]/data3_1d[i]/np.sqrt(2)] for i in range(len(data1_1d))]).T
    
    sorted_data = np.array([sorted([data1_1d[i], data2_1d[i]]) for i in range(len(data1_1d))]).T
    average_data = [np.average(data[0]), np.average(data[1])]
    std_data = [np.std(data[0]), np.std(data[1])]
    average_sorted_data = [np.average(sorted_data[0]), np.average(sorted_data[1])]
    std_sorted_data = [np.std(sorted_data[0]), np.std(sorted_data[1])]

    average_c = np.average(data3_1d)
    average_ratio = [np.average(ratio[0]), np.average(ratio[1])]
    std_ratio = [np.std(ratio[0]), np.std(ratio[1])]

    ofdomain = open(dirname+"/domain.dat", "w")
    idx_domain = np.array([int(data1_1d[i] > data2_1d[i]) for i in range(len(data1_1d))])
    size_domain = [np.sum(idx_domain == 0), np.sum(idx_domain == 1)]
    for idx in idx_domain:
        ofdomain.write("%d\n"%idx)
    ofdomain.close()

    # # distance_threshold = 60
    # cluster = AgglomerativeClustering(n_clusters=2, affinity='euclidean', linkage='ward')
    # # cluster = AgglomerativeClustering(n_clusters=None, distance_threshold=distance_threshold, affinity='euclidean', linkage='complete')
    # 
    # cluster.fit(data)
    # # print(cluster.labels_)
    # d1 = np.array([data[k][0] for k in range(len(data)) if cluster.labels_[k] == 0])
    # d2 = np.array([data[k][0] for k in range(len(data)) if cluster.labels_[k] == 1])

    of1.write("%d    %10.3f %10.3f   %10.3f %10.3f   %d %d\n"%(T, average_data[0]/np.sqrt(2), average_data[1]/np.sqrt(2), std_data[0], std_data[1], size_domain[0], size_domain[1]))
    of2.write("%d    %10.3f %10.3f   %10.3f %10.3f   %d %d\n"%(T, average_sorted_data[0]/np.sqrt(2), average_sorted_data[1]/np.sqrt(2), std_sorted_data[0], std_sorted_data[1], size_domain[0], size_domain[1]))
    of3.write("%d    %10.3f %10.3f   %10.3f %10.3f   %10.3f\n"%(T, average_ratio[0], average_ratio[1], std_ratio[0], std_ratio[1], average_c))

import sys
import numpy as np
from sklearn.cluster import AgglomerativeClustering
of1 = open("pbpbdistances.out", "w")
of2 = open("pbpbdistances-sorted.out", "w")
of1.write("T  average_135  average_45  std_135  std_45  size_domain_0 size_domain_1\n")
of2.write("T  average_short  average_long  std_short  std_long  size_domain_0 size_domain_1\n")
of3 = open("pbpbratio.out", "w")
of3.write("T  average_135/c  average_45/c  std_135/c  std_45/c  c\n")
for T in [50]+[i for i in range(52,64)]+[i for i in range(65,106)]+[i for i in range(110,136,5)]+[i for i in range(136,146)]+[i for i in range(147,151)]+[i for i in range(152,155)]+[170,190,195,200,205,210,215,220,225,230,250,255,260,265,270,275,280,285,290,295,300,305,310,320,340,350,400,450]:
    dirname = None
    if T >=50 & T <=155:
        dirname = "annealT"+str(T)
    if T in [i for i in range(170,330,20)]+[350]:
        dirname = "1-2ns-T"+str(T)
    if T in [195,200,205,215,220,225,295,300,305]:
        dirname = "0-2ns-T"+str(T)
    if T in [255,260,265,275,280,285,320,340]:
        dirname = "0-3ns-T"+str(T)
    if T in [400,450]:
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

    of1.write("%d    %10.3f %10.3f   %10.3f %10.3f   %d %d\n"%(T, average_data[0], average_data[1], std_data[0], std_data[1], size_domain[0], size_domain[1]))
    of2.write("%d    %10.3f %10.3f   %10.3f %10.3f   %d %d\n"%(T, average_sorted_data[0], average_sorted_data[1], std_sorted_data[0], std_sorted_data[1], size_domain[0], size_domain[1]))
    of3.write("%d    %10.3f %10.3f   %10.3f %10.3f   %10.3f\n"%(T, average_ratio[0], average_ratio[1], std_ratio[0], std_ratio[1], average_c))

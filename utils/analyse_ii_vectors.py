import sys
import numpy as np
import math
from sklearn.cluster import AgglomerativeClustering
of1 = open("ii_phi.out", "w")
# of2 = open("ii_theta.out", "w")
of1.write("T  average_phi  std_phi\n")
# of2.write("T  average_theta  std_theta\n")
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

    data = np.loadtxt(dirname+"/ii_vectors_caxis.dat").T
    phi = data[1]
    # theta = data[0]
    phiy = [np.abs(phi[i]) for i in range(len(phi)) if i % 3 == 1]
    average_phiy = np.average(phiy)
    std_phiy = np.std(phiy)
    # theta = [np.abs(theta[i]) for i in range(len(phi)) if i % 3 == 2]
    # average_theta = np.average(theta)
    # std_theta = np.std(theta)

    of1.write("%d    %10.3f %10.3f\n"%(T, average_phiy, std_phiy))
    # of2.write("%d    %10.3f %10.3f\n"%(T, average_theta, std_theta))

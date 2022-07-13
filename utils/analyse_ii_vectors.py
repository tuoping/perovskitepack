import sys
import numpy as np
import math
from sklearn.cluster import AgglomerativeClustering
of1 = open("ii_phi.out", "w")
# of2 = open("ii_theta.out", "w")
of1.write("T  average_phi  std_phi\n")
# of2.write("T  average_theta  std_theta\n")
for T in [100,110,120,130,140,155,175,195,215,235,255,275,295,315,335,355,375,395,415]:
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

import numpy as np
import math
import matplotlib.pyplot as plt
import sys

# filename = sys.argv[1]
f1 = "mesh_coords"

coord = np.loadtxt(f1, skiprows=0)[:512]
coord = coord.T

caxis = int(np.loadtxt("caxis"))
x = coord[(caxis+1)%3]
y = coord[(caxis+2)%3]
z = coord[caxis]
print(x, y)

# f2 = "ii_vectors_caxis.dat"
# data1 = np.loadtxt(f2)[:512*3].T[1]
# phiy = [data1[i] for i in range(512*3) if i%3 == 1]
# plt.scatter(x,z,c=phiy)
# plt.colorbar()
# plt.savefig("ii_phiy_xOz", dpi=1100, bbox_inches = "tight")

# data1 = np.loadtxt(f2)[:512*3].T[1]
# phix = [data1[i] for i in range(512*3) if i%3 == 0]
# plt.scatter(x,y,c=phix)
# plt.colorbar()
# plt.savefig("ii_phix", dpi=1100, bbox_inches = "tight")

# data1 = np.loadtxt(f2)[:512*3].T[1]
# phiz = [data1[i] for i in range(512*3) if i%3 == 2]
# plt.scatter(x,z,c=phiz,cmap="rainbow")
# plt.colorbar()
# plt.savefig("ii_phiz_xOz", dpi=1100, bbox_inches = "tight")

f2 = "mol_shortaxis_angle"
data1 = np.loadtxt(f2)[:512].T
# phix = [math.degrees(data1[0][i]) for i in range(len(data1[0]))]
# plt.scatter(x,y,c=phix)
# plt.colorbar()
# plt.savefig("molshortaxis_phix", dpi=1100, bbox_inches = "tight")

# phiy = [math.degrees(data1[1][i]) for i in range(len(data1[1]))]
# plt.scatter(x,y,c=phiy)
# plt.colorbar()
# plt.savefig("molshortaxis_phiy", dpi=1100, bbox_inches = "tight")

phiz = [math.degrees(data1[2][i]) for i in range(len(data1[2]))]
plt.scatter(x,y,c=phiz)
plt.colorbar()
plt.savefig("molshortaxis_phiz", dpi=1100, bbox_inches = "tight")
plt.show()

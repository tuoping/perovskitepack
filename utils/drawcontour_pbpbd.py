import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import sys

# filename = sys.argv[1]
f1 = "center_coords"

coord = np.loadtxt(f1, skiprows=0)
coord = coord.T

caxis = int(np.loadtxt("caxis"))
x = np.array(coord[(caxis+1)%3])
y = np.array(coord[(caxis+2)%3])
z = np.array(coord[caxis])


f2 = "left_distances"
data1 = np.loadtxt(f2)
f3 = "right_distances"
data2 = np.loadtxt(f3)

data = [data1[i]/data2[i] for i in range(len(data1))]
# data = [int(data1[i]/data2[i] > 1) for i in range(len(data1))]
plt.scatter(x,y,c=data, cmap="rainbow")
plt.colorbar()
plt.savefig("pbpbd_ratio", dpi=1100, bbox_inches = "tight")

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

# f2 = "mol_shortaxis_angle"
# data1 = np.loadtxt(f2)[:512].T
# phix = [math.degrees(data1[0][i]) for i in range(len(data1[0]))]
# plt.scatter(x,y,c=phix)
# plt.colorbar()
# plt.savefig("molshortaxis_phix", dpi=1100, bbox_inches = "tight")

# phiy = [math.degrees(data1[1][i]) for i in range(len(data1[1]))]
# plt.scatter(x,y,c=phiy)
# plt.colorbar()
# plt.savefig("molshortaxis_phiy", dpi=1100, bbox_inches = "tight")

# phiz = [math.degrees(data1[2][i]) for i in range(len(data1[2]))]
# plt.scatter(x,y,c=phiz)
# plt.colorbar()
# plt.savefig("molshortaxis_phiz", dpi=1100, bbox_inches = "tight")

'''
f2 = "mol_longaxis"
data1 = np.loadtxt(f2)[:512].T
for i in range(len(data1)):
    if data1[i][np.argmax(np.abs(data1[i]))] <0:
        data1[i] = -data1[i]
mol0 = [abs(data1[(caxis+1)%3][i]) for i in range(len(_z)) if _z[i] < 0]
mol1 = [abs(data1[(caxis+2)%3][i]) for i in range(len(_z)) if _z[i] < 0]
plt.quiver(x,y,mol0, mol1, headwidth=2, headlength=2, width = 0.005, scale=5, scale_units='inches')
plt.savefig("mollongaxis", dpi=1100, bbox_inches = "tight")

plt.figure()
# plt.rcParams["figure.figsize"] = (5,5)
plt.scatter(x_mesh, y_mesh)

segs1 = np.stack((x_mesh,y_mesh), axis=2)
segs2 = segs1.transpose(1,0,2)
plt.gca().add_collection(LineCollection(segs1))
plt.gca().add_collection(LineCollection(segs2))
mol0 = [abs(data1[(caxis+1)%3][i]) for i in range(len(_z)) if _y[i] < 0]
mol1 = [abs(data1[(caxis+3)%3][i]) for i in range(len(_z)) if _y[i] < 0]
plt.quiver(x,y,mol0, mol1, headwidth=2, headlength=2, width = 0.005, scale=5, scale_units='inches')
plt.savefig("mollongaxis_xOz", dpi=1100, bbox_inches = "tight")
# plt.quiver(x,z,data1[0], data1[2])
# plt.savefig("mollongaxis_xOz", dpi=1100, bbox_inches = "tight")
'''


plt.show()

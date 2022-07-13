import sys
import numpy as np
import math
import matplotlib.pyplot as plt
dirname = "./"
caxis = int(np.loadtxt(dirname+"/caxis"))
mesh_size = [30,30,2]

Totaltime = 1
Num_oct = mesh_size[0]*mesh_size[1]*mesh_size[2]
data = np.loadtxt(dirname+"/ii_vectors_caxis.dat")
data = np.reshape(data, [-1,Num_oct,3,2])
Totaltime = len(data)
print("Totaltime = ", Totaltime)
phiy= data[:,:,1,1]
# theta = data[0]
# theta = [np.abs(theta[i]) for i in range(len(phi)) if i % 3 == 2]
# average_theta = np.average(theta)
# std_theta = np.std(theta)

phiy_symbol = []
for m in range(Totaltime):
    idx = 0
    for i in range(mesh_size[0]):
        for j in range(mesh_size[1]):
            for k in range(mesh_size[2]):
                if caxis == 2:
                    phiy_symbol.append( phiy[m][idx]*math.pow(-1, i)*math.pow(-1, j))
                if caxis == 1:
                    phiy_symbol.append( phiy[m][idx]*math.pow(-1, i)*math.pow(-1, k))
                if caxis == 0:
                    phiy_symbol.append( phiy[m][idx]*math.pow(-1, j)*math.pow(-1, k))
                phiy_symbol[-1] = np.abs(phiy_symbol[-1])
                idx += 1
average_phiy = np.average(phiy_symbol)
std_phiy = np.std(phiy_symbol)

f1 = "mesh_coords"

coord = np.loadtxt(f1, skiprows=0)
coord = coord.T

caxis = int(np.loadtxt("caxis"))
x = np.array(coord[(caxis+1)%3])
y = np.array(coord[(caxis+2)%3])
z = np.array(coord[caxis])

fig = plt.figure()
plt.scatter(x,y,c=phiy_symbol, marker="s",s=60, cmap="rainbow")
plt.colorbar()
figfilename = "phiy_ii"
plt.savefig(figfilename, bbox_inches = "tight")
plt.show()

import sys
import numpy as np

data = np.loadtxt("molecule_longaxis.dat").T

caxis = int(np.loadtxt("caxis"))
dxy = []
for i in range(3):
    if i != caxis:
        dxy.append([])
        d1 = np.array([x for x in data[i] if x < 50])
        d2 = np.array([x for x in data[i] if x > 50])
        dxy[-1].append(d1)
        dxy[-1].append(d2)
print(np.average(data[caxis]), np.average(dxy[0][0]), np.average(dxy[0][1]), np.average(dxy[1][0]), np.average(dxy[1][1]), np.std(data[caxis]), np.std(dxy[0][0]), np.std(dxy[0][1]), np.std(dxy[1][0]), np.std(dxy[1][1]))

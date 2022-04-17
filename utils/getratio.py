import sys
import numpy as np
from copy import deepcopy
fname = "cellsize.out"
d = np.loadtxt(fname, skiprows=1)
of = open("ratio-cellsize.out", "w")
of.write("step   (a+b)/2/c/sqrt(2) a/c/sqrt(2) b/c/sqrt(2)\n")
for i in range(len(d)):
    c = sorted([d[i][1],d[i][2],d[i][3]])[0]
    a = sorted([d[i][1],d[i][2],d[i][3]])[1]
    b = sorted([d[i][1],d[i][2],d[i][3]])[2]
    r = (a+b)/2/c/np.sqrt(2)
    r1 = (a)/c/np.sqrt(2)
    r2 = (b)/c/np.sqrt(2)
    of.write("%5d    %11.5f  %11.5f  %11.5f\n" % (d[i][0],r,r1,r2))

of.close()

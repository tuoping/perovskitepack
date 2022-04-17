import dpdata
import numpy as np
import glob,os,re

def getorder(fname):
    f = re.match("traj/(.*?).lammpstrj", fname)
    return int(f[1])

def pbc_idx(idx):
    if idx > 2:
        return 0
    if idx < 0:
        return 2
    return idx

ifiles = sorted(glob.glob("traj/*.lammpstrj"), key=getorder)
of = open("cellsize.out", "w")
of.write("step    a  b  c\n")
of2 = open("ratio-cellsize.out", "w")
of2.write("step   a/c/sqrt(2)\n")
atot = 0.0
btot = 0.0
ctot = 0.0
caxis = int(np.loadtxt("caxis"))
for fname in ifiles:
    f = re.match("traj/(.*?).lammpstrj", fname)
    i = int(f[1])
    d = dpdata.System(fname, "lammps/dump")
    x_vec = d["cells"][0][pbc_idx(caxis+1)]
    y_vec = d["cells"][0][pbc_idx(caxis-1)]
    a_vec = y_vec-x_vec
    b_vec = x_vec+y_vec
    a = np.linalg.norm(a_vec)/8.0
    b = np.linalg.norm(b_vec)/8.0
    c = np.linalg.norm(d["cells"][0][caxis])/8.0
    atot += a
    btot += b
    ctot += c
    of.write("%5d    %11.5f  %11.5f  %11.5f\n" % (i/1000,a,b,c))
    # dim = sorted([a,b,c])
    of2.write("%5d    %11.5f\n" % (i/1000,(a+b)/2/c/np.sqrt(2)))

of.close()
of2.close()
atot /= float(len(ifiles))
btot /= float(len(ifiles))
ctot /= float(len(ifiles))
print(atot, btot, ctot)

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

#ifiles = sorted(glob.glob("traj/*.lammpstrj"), key=getorder)
ifiles = ["traj/"+str(i)+".lammpstrj" for i in range(0, 1210000, 1000)]
of = open("cellsize.out", "w")
of.write("step    a  b  c\n")
of2 = open("diag-cellsize.out", "w")
of2.write("step    a  b  c\n")
atot = 0.0
btot = 0.0
ctot = 0.0
# caxis = int(np.loadtxt("caxis"))
caxis = 2

for fname in ifiles:
    f = re.match("traj/(.*?).lammpstrj", fname)
    k = int(f[1])
    d = dpdata.System(fname, "lammps/dump")
    for i in range(1):
        len_a = np.linalg.norm(d["cells"][i][0])
        len_b = np.linalg.norm(d["cells"][i][1])
        len_c = np.linalg.norm(d["cells"][i][2])
        if k==0:
            a0=len_a
            b0=len_b
            c0=len_c
        x_vec = d["cells"][i][pbc_idx(caxis+1)]+d["cells"][i][pbc_idx(caxis-1)]
        y_vec = d["cells"][i][pbc_idx(caxis+1)]-d["cells"][i][pbc_idx(caxis-1)]
        z_vec = d["cells"][i][caxis]
        x = np.linalg.norm(x_vec)
        y = np.linalg.norm(y_vec)
        z = np.linalg.norm(z_vec)
        of.write("%10d    %11.5f  %11.5f  %11.5f\n" % (k,len_a/a0,len_b/b0,len_c/c0))
        of2.write("%5d    %11.5f  %11.5f  %11.5f\n" % (k,x,y,z))

of.close()
of2.close()

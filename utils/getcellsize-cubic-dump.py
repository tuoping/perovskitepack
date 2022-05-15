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

ifiles = ["md4.dump"]
of = open("cellsize.out", "w")
of.write("step    a  b  c\n")
of2 = open("ratio-cellsize.out", "w")
of2.write("step   a/c/sqrt(2)\n")
atot = 0.0
btot = 0.0
ctot = 0.0
caxis = int(np.loadtxt("caxis"))
k=0
for fname in ifiles:
    d = dpdata.System(fname, "lammps/dump")
    for i in range(d.get_nframes()):
        x_vec = d["cells"][i][pbc_idx(caxis+1)]
        y_vec = d["cells"][i][pbc_idx(caxis-1)]
        a_vec = y_vec-x_vec
        b_vec = x_vec+y_vec
        a = np.linalg.norm(a_vec)/8.0
        b = np.linalg.norm(b_vec)/8.0
        c = np.linalg.norm(d["cells"][i][caxis])/8.0
        atot += a
        btot += b
        ctot += c
        of.write("%5d    %11.5f  %11.5f  %11.5f\n" % (k,a,b,c))
        of2.write("%5d    %11.5f\n" % (k,(a+b)/2/c/np.sqrt(2)))
        k += 1

of.close()
of2.close()
atot /= float(d.get_nframes())
btot /= float(d.get_nframes())
ctot /= float(d.get_nframes())
print(atot, btot, ctot)

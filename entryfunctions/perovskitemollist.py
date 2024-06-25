from copy import deepcopy
import dpdata
from dpdata.lammps.lmp import box2lmpbox
# from pymatgen.core.operations import SymmOp
import numpy as np
from numpy.linalg import norm
import random, math
import sys, os
import datetime
from collections import Counter

from perovskitelattpack import *

if __name__ == "__main__":
    filename = sys.argv[1]
    # import cubic perovskite
    if filename.endswith(".lmp"):
        fmt = "lammps/lmp"
    else:
        if filename.endswith(".lammpstrj") or filename.endswith(".dump"):
            fmt = "lammps/dump"
        else:
            if filename.endswith(".vasp") or filename.endswith("POSCAR") or filename.endswith("CONTCAR"):
                fmt = "vasp/poscar"
            else:
                raise Exception("Unknown file format")
    cubic = FAPbI3(filename, fmt=fmt)
    
    # set axis according to "caxis"
    caxis = int(np.loadtxt("caxis"))
    ref_axis = np.eye(3)
    axis = np.eye(3)
    axis[2]=ref_axis[caxis]
    axis[0]=ref_axis[(caxis+1)%3]
    axis[1]=np.cross(axis[0], axis[2])
    cubic.set_axis(axis)
    
    # set cutoff of bonds
    cubic.setcutoff_I_Pb(5.0)
    cubic.setcutoff_CN_H(1.6)

    # start a mesh
    mesh_dim = [4,4,4]
    cubic.startmesh(mesh_dim, eps=0.0) 
    
    cubic.extract_mol()
    f = open("molecule_longaxis.dat", "w")
    for mol in cubic.molecules:
        f.write("%f %f %f\n"%(math.degrees(mol.angle_longaxis[0]), math.degrees(mol.angle_longaxis[1]), math.degrees(mol.angle_longaxis[2])))
    f.close()

    f = open("molecule_polaraxis.dat", "w")
    for mol in cubic.molecules:
        f.write("%f %f %f\n"%(math.degrees(mol.angle_polaraxis[0]), math.degrees(mol.angle_polaraxis[1]), math.degrees(mol.angle_polaraxis[2])))
    f.close()
    

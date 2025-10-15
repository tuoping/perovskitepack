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
import os,sys

sys.path.append("/home/tuoping/pkgs/perovskitepack")

from perovskitelattpack import *

if __name__ == "__main__":
    filename = sys.argv[1]
    # import cubic perovskite
    if filename.endswith(".lmp"):
        fmt = "lammps/lmp"
    elif filename.endswith(".lammpstrj") or filename.endswith(".dump"):
        fmt = "lammps/dump"
    elif filename.endswith(".vasp") or filename.endswith("POSCAR") or filename.endswith("CONTCAR"):
        fmt = "vasp/poscar"
    elif filename.endswith(".xyz") or filename.endswith(".extxyz"):
        fmt = "xyz"
        from ase.io import read
        traj = read(filename, format="extxyz", index=":")
    else:
        raise Exception("Unknown file format")
    if not ( filename.endswith(".xyz") or filename.endswith(".extxyz")) :
        traj = dpdata.System(filename, fmt=fmt)
        num_frames = traj.get_nframes()
    else:
        num_frames = len(traj)
    for i_frame in range(num_frames):
        print("Frame:: ", i_frame)
        frm = traj[i_frame]
        if filename.endswith(".xyz") or filename.endswith(".extxyz"):
            frm = dpdata.System(frm, fmt="ase/structure")
        cubic = FAPbI3(frm, fmt=fmt)
        
        # set axis according to "caxis"
        # caxis = int(np.loadtxt("caxis"))
        # ref_axis = np.eye(3)
        axis = np.eye(3)
        # axis[2]=ref_axis[caxis]
        # axis[0]=ref_axis[(caxis+1)%3]
        # axis[1]=np.cross(axis[0], axis[2])
        cubic.set_axis(axis)
        
        # set cutoff of bonds
        cubic.setcutoff_I_Pb(5.0)
        cubic.setcutoff_CN_H()
        
        try:
            indices_mol = np.load("indices_mol.npy")
            cubic.extract_mol_from_indices(indices_mol, moltype="FA")
        except:
            indices_mol = cubic.extract_mol(moltype="FA")
            np.save("indices_mol.npy", indices_mol)
        import os
        os.makedirs("unitaxis_mols", exist_ok=True)
        f = open(f"unitaxis_mols/molecule_unitlongaxis_frame{i_frame}.dat", "w")
        for mol in cubic.molecules:
            f.write("%f %f %f\n"%((mol.unitlongaxis[0]), (mol.unitlongaxis[1]), (mol.unitlongaxis[2])))
        f.close()
    
        f = open(f"unitaxis_mols/molecule_unitpolaraxis_frame{i_frame}.dat", "w")
        for mol in cubic.molecules:
            f.write("%f %f %f\n"%((mol.unitpolaraxis[0]), (mol.unitpolaraxis[1]), (mol.unitpolaraxis[2])))
        f.close()
        

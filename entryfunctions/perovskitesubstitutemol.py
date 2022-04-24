from copy import deepcopy
import dpdata
# from pymatgen.core.operations import SymmOp
import numpy as np
import sys
import time
from itertools import combinations

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
    # caxis = int(np.loadtxt("caxis"))
    # ref_axis = np.eye(3)
    # axis = np.eye(3)
    # axis[2]=ref_axis[caxis]
    # axis[0]=ref_axis[(caxis+1)%3]
    # axis[1]=np.cross(axis[0], axis[2])
    indices_mol = np.load("indices_mol.npy")

    traj = dpdata.System(filename, fmt=fmt)
    for i_frame in range(0,traj.get_nframes()):
        stime = time.time()
        frm = traj[i_frame]
        cubic = FAPbI3(frm, fmt=fmt)
        
        # # set axis according to "caxis"
        # cubic.set_axis(axis)
        
        # set cutoff of bonds
        cubic.setcutoff_I_Pb(5.0)
        cubic.setcutoff_CN_H(1.6)
    
        cubic.extract_mol_from_indices(indices_mol)
        # indices_mol = cubic.extract_mol()
        # np.save("indices_mol.npy", indices_mol)
        idx = 0
        for i_mol in combinations(range(len(cubic.molecules)), 7):
            print(i_mol[0], i_mol[1], i_mol[2], i_mol[3], i_mol[4], i_mol[5], i_mol[6])
            new_cubic = cubic.substitute_mol_by_idx(i_mol[0], "Cs")
            inext = i_mol[1] - 1 if i_mol[1] > i_mol[0] else i_mol[1]
            new_cubic = new_cubic.substitute_mol_by_idx(inext, "Cs")
            inext = i_mol[2] - 2
            new_cubic = new_cubic.substitute_mol_by_idx(inext, "Cs")
            inext = i_mol[3] - 3
            new_cubic = new_cubic.substitute_mol_by_idx(inext, "Cs")
            inext = i_mol[4] - 4
            new_cubic = new_cubic.substitute_mol_by_idx(inext, "Cs")
            inext = i_mol[5] - 5
            new_cubic = new_cubic.substitute_mol_by_idx(inext, "Cs")
            inext = i_mol[6] - 6
            new_cubic = new_cubic.substitute_mol_by_idx(inext, "Cs")
            new_cubic.cubic.to_vasp_poscar("POSCAR_substituted"+str(idx)+".vasp")
            idx += 1
        

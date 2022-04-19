from copy import deepcopy
import dpdata
# from pymatgen.core.operations import SymmOp
import numpy as np
import sys
import time

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

    traj = dpdata.System(filename, fmt=fmt)
    for i_frame in range(0,traj.get_nframes()):
        stime = time.time()
        frm = traj[i_frame]
        cubic = FAPbI3(frm, fmt=fmt)
        
        
        # set cutoff of bonds
        cubic.setcutoff_I_Pb(5.0)
        cubic.setcutoff_CN_H(1.6)
    
        # cubic.extract_mol_from_indices(indices_mol)
        indices_mol = cubic.extract_mol()
        np.save("indices_mol.npy", indices_mol)
        cubic.rotate_moleculelongaxis_by_idx(0, [0,0,1], 12*np.pi/12)
        cubic.cubic.to_vasp_poscar("POSCAR_rotated.vasp")
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
        cubic.cubic.to_lammps_lmp("POSCAR.lmp")
        
        # set cutoff of bonds
        cubic.setcutoff_I_Pb(5.0)
        cubic.setcutoff_CN_H(1.6)

        indices_mol = np.load("indices_mol.npy")
        cubic.extract_mol_from_indices(indices_mol)
        # indices_mol = cubic.extract_mol()
        # np.save("indices_mol.npy", indices_mol)
        '''
        # start a mesh
        mesh_dim = [30,30,2]
        cubic.startmesh(mesh_dim, eps=0.0) 
        ### map molecules to mesh
        Succeed = cubic.mesh.map_obj_mesh(cubic.molecules)
        if not Succeed:
            raise Exception("mapping failed")
        '''
        for idx, mol in enumerate(cubic.molecules):
            inter_center_coord = phys2Inter(mol.center_coord, cubic.cell)
            if inter_center_coord[0] > 0.5:
                # if mol.angle_longaxis[0] > mol.angle_longaxis[1]:
                #     cubic.rotate_moleculelongaxis_by_idx(idx, [0,0,1], math.radians(14))
                # else:
                #     cubic.rotate_moleculelongaxis_by_idx(idx, [0,0,1], math.radians(-14))
                pass
            else:
                if mol.angle_longaxis[0] > mol.angle_longaxis[1]:
                    cubic.rotate_moleculelongaxis_by_idx(idx, [0,0,1], math.radians(-28))
                else:
                    cubic.rotate_moleculelongaxis_by_idx(idx, [0,0,1], math.radians(28))
        cubic.cubic.to_lammps_lmp("POSCAR_rotated.lmp")

from copy import deepcopy
import dpdata
# from pymatgen.core.operations import SymmOp
import numpy as np
import sys,os
import time
from ase.geometry.geometry import get_distances

SCRIPT_DIR="/home/tuoping/pkgs/perovskitepack/"
sys.path.append(os.path.dirname(SCRIPT_DIR))
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
            # inter_center_coord = phys2Inter(mol.center_coord, cubic.cell)
            rotaxis = np.random.randn(3)  # Random vector in 3D
            rotaxis = rotaxis / np.linalg.norm(rotaxis)  # Normalize to make it a unit vector
            rotangle = np.random.uniform(0, 360)
            cubic.rotate_moleculelongaxis_by_idx(idx, rotaxis, math.radians(rotangle))
        # perturb I-Pb frame
        type_I = cubic.cubic["atom_names"].index("I")
        type_Pb = cubic.cubic["atom_names"].index("Pb")
        list_I = np.where(cubic.cubic["atom_types"]==type_I)
        list_Pb = np.where(cubic.cubic["atom_types"]==type_Pb)
        cubic.cubic["coords"][0][list_I[0]] += np.random.normal(scale=0.2, size=(list_I[0].shape[0], 3))
        cubic.cubic["coords"][0][list_Pb[0]] += np.random.normal(scale=0.2, size=(list_Pb[0].shape[0], 3))

        # cubic.cubic.to_lammps_lmp("POSCAR_rotated.lmp")
        cubic.cubic.to_vasp_poscar("POSCAR_rotated.vasp")

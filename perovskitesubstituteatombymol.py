from copy import deepcopy
import dpdata
# from pymatgen.core.operations import SymmOp
import numpy as np
import sys
import time
from itertools import combinations

from perovskitelattpack import *
from ase.geometry.geometry import get_distances
from pymatgen.core.operations import SymmOp

def rotate_long_axis(mol, axis, angle, angle_in_radians=True):
    op = SymmOp.from_origin_axis_angle(
                (0, 0, 0),
                axis=tuple(axis),
                angle=angle,
                angle_in_radians=angle_in_radians
            )
    m = op.rotation_matrix
    center_coord = np.mean(mol['coords'][0], axis=0)
    centered_coords = get_distances(mol['coords'][0], mol['coords'][0][0])[0][:,0,:]
    new_coords = np.zeros(mol['coords'][0].shape)

    for idx, c in enumerate(centered_coords):
        new_coords[idx] = np.dot(m, c.T).T
    for i in range(len(new_coords)):
        new_coords[i] += center_coord

    for i in range(len(new_coords)):
        mol['coords'][0][i] = new_coords[i]

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
    # indices_mol = np.load("indices_mol.npy")
    sample_mol = dpdata.System("mol.vasp", fmt="vasp/poscar")
    traj = dpdata.System(filename, fmt=fmt)
    for i_frame in range(0,traj.get_nframes()):
        stime = time.time()
        frm = traj[i_frame]
        cubic = FAPbI3(frm, fmt=fmt, type_map=["I", "Pb", 'C'])
        
        # # set axis according to "caxis"
        # cubic.set_axis(axis)
        
        # set cutoff of bonds
        cubic.setcutoff_I_Pb(5.0)
        cubic.setcutoff_CN_H(1.6)
    
        # cubic.extract_mol_from_indices(indices_mol)
        indices_mol = np.where(cubic.cubic['atom_types'] == cubic.cubic['atom_names'].index("C"))[0]
        print(indices_mol)
        for i_substituted, i_mol in enumerate(indices_mol[:]):
            print(i_mol)
            while(True):
                rotaxis = np.random.randn(3)  # Random vector in 3D
                rotaxis = rotaxis / np.linalg.norm(rotaxis)  # Normalize to make it a unit vector
                rotangle = np.random.uniform(0, 360)
                rotate_long_axis(sample_mol, rotaxis, math.radians(rotangle))
                new_cubic = cubic.substitute_idx_by_molecule(i_mol-i_substituted, sample_mol)
                distmat = get_distances(new_cubic.cubic['coords'][0], new_cubic.cubic['coords'][0])[1]
                distances = np.array([distmat[j][i] for j in range(distmat.shape[0]) for i in range(j+1, distmat.shape[1])])
                print(">>> min dist::", distances.min())
                if np.all(distances > 1):
                    break
            cubic = new_cubic
        new_cubic.cubic.to_vasp_poscar("POSCAR_substituted.vasp")
        

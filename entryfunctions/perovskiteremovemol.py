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

    #indices_mol = np.load("indices_mol.npy")
    traj = dpdata.System(filename, fmt=fmt)
    caxis = int(np.loadtxt("caxis"))
    ref_axis = np.eye(3)
    axis = np.eye(3)
    axis[2]=ref_axis[caxis]
    axis[0]=ref_axis[(caxis+1)%3]
    axis[1]=np.cross(axis[0], axis[2])
    molecules = []
    for i_frame in range(0,min(traj.get_nframes(), 1000),10):
        stime = time.time()
        frm = traj[i_frame]
        cubic = FAPbI3(frm, fmt=fmt)
        
        # # set axis according to "caxis"
        cubic.set_axis(axis)
        
        # set cutoff of bonds
        cubic.setcutoff_I_Pb(5.0)
        cubic.setcutoff_CN_H(1.6)
    
        if i_frame == 0:
            indices_mol = cubic.extract_mol(moltype = "MA")
            np.save("indices_mol.npy", indices_mol)
        else:
            cubic.extract_mol_from_indices(indices_mol, moltype="MA")

        mesh_dim = [4,4,4]
        cubic.startmesh(mesh_dim, eps=0.0) 
        Succeed = cubic.mesh.map_obj_mesh(cubic.molecules)
        if not Succeed:
            cubic.startmesh(mesh_dim, eps=-1.0) 
            Succeed = cubic.mesh.map_obj_mesh(cubic.molecules)
            if not Succeed:
                cubic.startmesh(mesh_dim, eps=-1.5) 
                Succeed = cubic.mesh.map_obj_mesh(cubic.molecules)
                if not Succeed:
                    cubic.startmesh(mesh_dim, eps=1.5) 
                    Succeed = cubic.mesh.map_obj_mesh(cubic.molecules)
                    if not Succeed:
                        cubic.startmesh(mesh_dim, eps=1.0) 
                        Succeed = cubic.mesh.map_obj_mesh(cubic.molecules)
                        if not Succeed:
                            cubic.startmesh(mesh_dim, eps=-2.0) 
                            Succeed = cubic.mesh.map_obj_mesh(cubic.molecules)
                            if not Succeed:
                                cubic.startmesh(mesh_dim, eps=-2.5) 
                                Succeed = cubic.mesh.map_obj_mesh(cubic.molecules)
                                if not Succeed:
                                    cubic.startmesh(mesh_dim, eps=2.5) 
                                    Succeed = cubic.mesh.map_obj_mesh(cubic.molecules)
                                    if not Succeed:
                                        cubic.startmesh(mesh_dim, eps=2.0) 
                                        Succeed = cubic.mesh.map_obj_mesh(cubic.molecules)
                                        if not Succeed:
                                            raise Exception("mapping failed")
        mesh = cubic.mesh
        center0 = mesh.get_mesh_point(2,2,2).obj
        center1 = mesh.get_mesh_point(3,2,2).obj
        center2 = mesh.get_mesh_point(2,3,2).obj
        center3 = mesh.get_mesh_point(3,3,2).obj
        center4 = mesh.get_mesh_point(2,2,3).obj
        center5 = mesh.get_mesh_point(3,2,3).obj
        center6 = mesh.get_mesh_point(2,3,3).obj
        center7 = mesh.get_mesh_point(3,3,3).obj
        molecules.append([center0, center1, center2, center3, center4, center5, center6, center7])
        # indices_mol = cubic.extract_mol()
        # np.save("indices_mol.npy", indices_mol)
    cubic = FAPbI3(traj[-1],fmt=fmt)
    cubic.extract_mol_from_indices(indices_mol)
    print(cubic.cubic)
    for i_mol in range(len(cubic.molecules)):
        new_cubic = cubic.remove_mol(0)
        cubic = new_cubic
    print(cubic.cubic)
    for mollist in molecules:
        for mol in mollist:
            cubic.add_atom("C", mol.coords_mol[0])
            cubic.add_atom("N", mol.coords_mol[1])
            cubic.add_atom("N", mol.coords_mol[2])
            cubic.add_atom("H", mol.coords_mol[3])
            cubic.add_atom("H", mol.coords_mol[4])
            cubic.add_atom("H", mol.coords_mol[5])
            cubic.add_atom("H", mol.coords_mol[6])
            cubic.add_atom("H", mol.coords_mol[7])
    cubic.cubic.to_vasp_poscar("1mol.vasp")

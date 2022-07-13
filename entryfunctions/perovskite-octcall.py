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
    # read "caxis"
    caxis = int(np.loadtxt("caxis"))
    ref_axis = np.eye(3)
    axis = np.eye(3)
    axis[2]=ref_axis[caxis]
    axis[0]=ref_axis[(caxis+1)%3]
    axis[1]=np.cross(axis[0], axis[2])
    # Output files
    ofr = open("right_distances", "w")
    ofl = open("left_distances", "w")
    ofc = open("center_coords", "w")
    ofm = open("mesh_coords", "w")
    f = open("ii_vectors_caxis.dat", "w")
    ofcd = open("c_distances", "w")
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
    for i in range(1):
        print("frame "+str(i))
        frm = traj[i]
        cubic = FAPbI3(frm, fmt=fmt)
        # set axis according to "caxis"
        cubic.set_axis(axis)
        
        # set cutoff of bonds
        cubic.setcutoff_I_Pb(5.0)
    
        # start a mesh
        mesh_dim = [30,30,2]
        cubic.startmesh(mesh_dim, eps=0.0)
        try:
            indices_oct = np.load("indices_oct.npy")
            cubic.extract_octahedron_from_indices(indices_oct)
        except:
            indices_oct = cubic.extract_octahedron()
            np.save("indices_oct.npy", indices_oct)

    
        ### map octahedra to mesh
        Succeed = cubic.mesh.map_obj_mesh(cubic.octahedra)
        if not Succeed:
            print("eps = 2.0")
            cubic.startmesh(mesh_dim, eps=2.0) 
            Succeed = cubic.mesh.map_obj_mesh(cubic.octahedra)
            if not Succeed:
                print("eps = -2.0")
                cubic.startmesh(mesh_dim, eps=-2.0)
                Succeed = cubic.mesh.map_obj_mesh(cubic.octahedra)
                if not Succeed:
                    print("skippinng frame "+str(i))
                    continue

        num = cubic.mesh.mesh_size
        mesh = cubic.mesh

        right_distances = []
        left_distances = []
        c_distances = []
        for i in range(num[0]):
            for j in range(num[1]):
                for k in range(num[2]):
                    center = mesh.get_mesh_point(i,j,k)
                    ofc.write("%f %f %f\n"%(center.obj.center_coord[0], center.obj.center_coord[1], center.obj.center_coord[2]))
                    ofm.write("%d %d %d\n"%(i,j,k))
                    if caxis == 0:
                        leftlow = mesh.get_mesh_point(i,j-1,k-1)
                        leftup = mesh.get_mesh_point(i,j-1,k+1)
                        c_d = mesh.get_mesh_point(i-1,j,k)
                    if caxis == 1:
                        leftlow = mesh.get_mesh_point(i-1,j,k-1)
                        leftup = mesh.get_mesh_point(i+1,j,k-1)
                        c_d = mesh.get_mesh_point(i,j-1,k)
                    if caxis == 2:
                        leftlow = mesh.get_mesh_point(i-1,j-1,k)
                        leftup = mesh.get_mesh_point(i-1,j+1,k)
                        c_d = mesh.get_mesh_point(i,j,k-1)
                    
                    right_distances.append(np.linalg.norm(distance(center.obj.center_coord, leftlow.obj.center_coord, cubic.cell)))
                    left_distances.append(np.linalg.norm(distance(center.obj.center_coord, leftup.obj.center_coord, cubic.cell)))
                    c_distances.append(np.linalg.norm(distance(center.obj.center_coord, c_d.obj.center_coord, cubic.cell)))
                    ofr.write("%f\n"%(right_distances[-1]))
                    ofl.write("%f\n"%(left_distances[-1]))
                    ofcd.write("%f\n"%(c_distances[-1]))
                  
        for oct in cubic.octahedra:
            f.write("%f %f\n"%(math.degrees(oct.theta_II_vectors[0]), math.degrees(oct.phi_II_vectors[0])))
            f.write("%f %f\n"%(math.degrees(oct.theta_II_vectors[1]), math.degrees(oct.phi_II_vectors[1])))
            f.write("%f %f\n"%(math.degrees(oct.theta_II_vectors[2]), math.degrees(oct.phi_II_vectors[2])))
            f.write("\n")
    f.close()
    ofr.close()
    ofl.close()
    ofc.close()
    ofm.close()
    ofcd.close()

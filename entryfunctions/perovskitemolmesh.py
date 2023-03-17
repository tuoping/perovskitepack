from copy import deepcopy
import dpdata
# from pymatgen.core.operations import SymmOp
import numpy as np
import sys
import time
lib_path="/home/tuoping/MAPbI3/phasediagram/cubic-natom24576sc16162-confparallel/perovskitepack"
sys.path.append(lib_path)
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
    caxis = int(np.loadtxt("caxis"))
    ref_axis = np.eye(3)
    axis = np.eye(3)
    axis[2]=ref_axis[caxis]
    axis[0]=ref_axis[(caxis+1)%3]
    axis[1]=np.cross(axis[0], axis[2])

    traj = dpdata.System(filename, fmt=fmt)
    for i_frame in range(0,traj.get_nframes()):
        print("frame"+str(i_frame))
        stime = time.time()
        frm = traj[i_frame]
        cubic = FAPbI3(frm, fmt=fmt)
        
        # set axis according to "caxis"
        cubic.set_axis(axis)
        
        # set cutoff of bonds
        cubic.setcutoff_I_Pb(5.0)
        cubic.setcutoff_CN_H(1.6)
    
        try:
            indices_mol = np.load("indices_mol.npy")
            cubic.extract_mol_from_indices(indices_mol, moltype="MA")
        except:
            indices_mol = cubic.extract_mol(moltype="MA")
            np.save("indices_mol.npy", indices_mol)
        
        # start a mesh
        ### map molecules to mesh
        mesh_dim = [2,32,32]
        try:
            molmesh = np.loadtxt("objmesh.dat")   
            cubic.startmesh(mesh_dim, eps=0.0) 
            cubic.mesh.read_obj_mesh(molmesh, cubic.molecules)
        except:
            cubic.startmesh(mesh_dim, eps=0.0) 
            Succeed = cubic.mesh.map_obj_mesh(cubic.molecules)
            if not Succeed:
                cubic.startmesh(mesh_dim, eps=-2.0) 
                Succeed = cubic.mesh.map_obj_mesh(cubic.molecules)
                if not Succeed:
                    cubic.startmesh(mesh_dim, eps=2.0) 
                    Succeed = cubic.mesh.map_obj_mesh(cubic.molecules)
                    if not Succeed:
                        raise Exception("mapping failed")

        ### Output coordinates & molecule vectors
        num = cubic.mesh.mesh_size
        mesh = cubic.mesh
        # ofl = open("mol_longaxis_frame"+str(i_frame), "w")
        ofs = open("mol_polaraxis_frame"+str(i_frame), "w")
        # ofla = open("mol_longaxis_angle_frame"+str(i_frame), "w")
        ofsa = open("mol_polaraxis_angle_frame"+str(i_frame), "w")
        ofc = open("mol_center_coords_frame"+str(i_frame), "w")
       
        for i in range(num[0]):
            for j in range(num[1]):
                for k in range(num[2]):
                    center = mesh.get_mesh_point(i,j,k)
                    # center = mesh.mesh[i][j][k]
                    ofc.write("%f %f %f\n"%(center.obj.center_coord[0], center.obj.center_coord[1], center.obj.center_coord[2]))
                    
                    # ofl.write("%f %f %f\n"%((center.obj.unitlongaxis[0], center.obj.unitlongaxis[1], center.obj.unitlongaxis[2])))
                    ofs.write("%f %f %f\n"%(center.obj.unitpolaraxis[0], center.obj.unitpolaraxis[1], center.obj.unitpolaraxis[2]))
    
                    # ofla.write("%f %f %f\n"%((math.degrees(center.obj.angle_longaxis[0]), math.degrees(center.obj.angle_longaxis[1]), math.degrees(center.obj.angle_longaxis[2]))))
                    ofsa.write("%f %f %f\n"%(math.degrees(center.obj.angle_polaraxis[0]), math.degrees(center.obj.angle_polaraxis[1]), math.degrees(center.obj.angle_polaraxis[2])))
        
        # ofl.close()
        ofs.close()
        # ofla.close()
        ofsa.close()
        ofc.close()
        etime = time.time()
        print(etime-stime)
        stime = etime

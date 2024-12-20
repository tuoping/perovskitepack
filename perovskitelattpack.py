from copy import deepcopy
from token import TYPE_COMMENT
import dpdata
from dpdata.lammps.lmp import box2lmpbox
from dpdata.lammps.dump import box2dumpbox
from pymatgen.core.operations import SymmOp
import numpy as np
from numpy.linalg import norm
import random, math
import sys, os
import datetime
from collections import Counter
import time
from ase.geometry.geometry import get_distances

def bubble_sortby_ii(alist,ii):
    length = len(alist)
    order = list(range(length))
    for i in range(length-1):
        for j in range(1, length-i):
            if alist[j-1][ii] > alist[j][ii]:
                t = deepcopy(alist[j-1])
                alist[j-1] = deepcopy(alist[j])
                alist[j] = deepcopy(t)
                intt = deepcopy(order[j-1])
                order[j-1] = deepcopy(order[j])
                order[j] = deepcopy(intt)
    return order


def apply_pbc(dx,cell):
    rdx=[np.mod(dx[0],1.0),np.mod(dx[1],1.0),np.mod(dx[2],1.0)]
    if(rdx[0]<-0.5):
        rdx[0]=rdx[0]+1.0
    if(rdx[1]<-0.5):
        rdx[1]=rdx[1]+1.0
    if(rdx[2]<-0.5):
        rdx[2]=rdx[2]+1.0
    if(rdx[0]>0.5):
        rdx[0]=rdx[0]-1.0
    if(rdx[1]>0.5):
        rdx[1]=rdx[1]-1.0
    if(rdx[2]>0.5):
        rdx[2]=rdx[2]-1.0
    return np.array(rdx)

def phys2Inter(dx, cell):
    invcell = np.linalg.inv(cell)
    return np.matmul(invcell, dx)

def Inter2phys(rdx, cell):
    return np.matmul(cell, rdx)
    

def distance(x1, x2, cell):
    physd = x2-x1
    interd = phys2Inter(physd, cell)
    rinterd = apply_pbc(interd, cell)
    d = Inter2phys(rinterd, cell)
    return d

def center(x1, x2, cell):
    d = distance(x1, x2, cell)
    center = x1+d/2
    return center

def write_vasp_poscar(cell, elem_type, elem_num, coords, filename="POSCAR"):
    fconf = open(filename, "w")
    fconf.write("{:<5s}".format("VASP") + "\n")
    fconf.write("{:<5d}".format(1) + "\n")
    fconf.write("%f  %f  %f"%(cell[0][0], cell[0][1], cell[0][2]) + "\n")
    fconf.write("%f  %f  %f"%(cell[1][0], cell[1][1], cell[1][2]) + "\n")
    fconf.write("%f  %f  %f"%(cell[2][0], cell[2][1], cell[2][2]) + "\n")
    type_str = ""
    for e in elem_type:
        type_str += "{:<5s}".format(e)
    type_str += "\n"
    num_str = ""
    for n in elem_num:
        num_str += "{:<5d}".format(n)
    num_str += "\n"
    fconf.write(type_str)
    fconf.write(num_str)
    fconf.write("C" + "\n")
    for c in coords:
        fconf.write("%f  %f  %f"%(c[0], c[1], c[2]) + "\n")
    fconf.close()

class Octahedron(object):
    """An octahedron is defined by a Pb atom and 6 I atoms."""

    """
    Attributes:
        cell: the perovskite supercell.
        coords_oct: the coordinates of Pb and I atoms.
        indices_oct: the indices of Pb and I atoms.
        mesh_point: the mesh position of Pb atom.

    Typical usage example:
        oct = Octahedron(indices_oct, coords_oct, cell)
        oct.setII_vectors(axis)
        oct.set_mesh_point(mesh)

    """
    
    def __init__(self, indices, coords, cell):
        self.cell = cell
        self.coords_oct = coords
        self.indices_oct = indices
        self.set_center_postion(coords[0])
        self.setPbI_vectors()
        
    def set_center_postion(self, coord):
        assert len(coord) == 3, "Pb coord must be a list of 3 numbers"
        _c = phys2Inter(np.array(coord), self.cell)
        pc = apply_pbc(_c,self.cell)+0.5
        self.center_coord = Inter2phys(pc,self.cell)


    def setPbI_vectors(self):
        self.PbI_vectors = []
        # step1: calculate all I-Pb vectors
        for idx in range(1, len(self.indices_oct)):
            v = distance(self.coords_oct[idx], self.coords_oct[0], self.cell)
            self.PbI_vectors.append(v)


    def setII_vectors(self, axis):
        """
        Set II vectors, and create class attribute "self.II_vectors" with II vectors in the order [parallel to x, parallel to y, parallel to z]

        Args:
            axis: cell. For example: [[1,0,0],[0,1,0],[0,0,1]]
        """
        normal_axis = np.array([axis[0]/np.linalg.norm(axis[0]), axis[1]/np.linalg.norm(axis[1]), axis[2]/np.linalg.norm(axis[2])])
        
        # find pb_I_vectors in z direction
        z_vectors = np.zeros([2,3])
        for idx in range(len(self.PbI_vectors)):
            v = self.PbI_vectors[idx]
            normal_v = v/np.linalg.norm(v)
            if np.degrees(np.abs(np.arccos(np.dot(normal_v, normal_axis[2])))) < 30:
                z_vectors[0] = v
            if np.degrees(np.abs(np.arccos(np.dot(normal_v, normal_axis[2])))) > 150:
                z_vectors[1] = v
        
        x_vectors = np.zeros([2,3])
        for idx in range(len(self.PbI_vectors)):
            v = self.PbI_vectors[idx]
            normal_v = v/np.linalg.norm(v)
            if np.degrees(np.abs(np.arccos(np.dot(normal_v, normal_axis[0])))) < 30:
                x_vectors[0] = v
            if np.degrees(np.abs(np.arccos(np.dot(normal_v, normal_axis[0])))) > 150:
                x_vectors[1] = v
        
        y_vectors = np.zeros([2,3])
        for idx in range(len(self.PbI_vectors)):
            v = self.PbI_vectors[idx]
            normal_v = v/np.linalg.norm(v)
            if np.degrees(np.abs(np.arccos(np.dot(normal_v, normal_axis[1])))) < 30:
                y_vectors[0] = v
            if np.degrees(np.abs(np.arccos(np.dot(normal_v, normal_axis[1])))) > 150:
                y_vectors[1] = v
        
        self.II_vectors = []
        self.II_vectors.append(x_vectors[0]-x_vectors[1])
        self.II_vectors.append(y_vectors[0]-y_vectors[1])
        self.II_vectors.append(z_vectors[0]-z_vectors[1])

    def angle_PbI_vectors(self,axis):
        """
        Calculate the angle between Pb-I vectors and the axis, and create class attribute "self.theta_PbI_vectors" & "self.phi_PbI_vectors"

        Args:
            axis: cell. For example: [[1,0,0],[0,1,0],[0,0,1]]
        """
        normal_axis = np.array([axis[0]/np.linalg.norm(axis[0]), axis[1]/np.linalg.norm(axis[1]), axis[2]/np.linalg.norm(axis[2])])
        self.theta_PbI_vectors = []
        self.phi_PbI_vectors = []
        for v in self.PbI_vectors:
            normal_v = v/np.linalg.norm(v)
            self.theta_PbI_vectors.append(np.arccos(np.dot(normal_v, normal_axis[2])))
            proj_v_a2 = np.dot(normal_v, normal_axis[2]) * normal_axis[2] 
            self.phi_PbI_vectors.append(np.arccos(np.dot((normal_v-proj_v_a2), normal_axis[0]))-0.5*np.pi)

    def angle_II_vectors(self, axis):
        """
        Calculate the angle between I-I vectors and the axis, and create class attribute "self.theta_II_vectors" & "self.phi_II_vectors"

        Args:
            axis: cell. For example: [[1,0,0],[0,1,0],[0,0,1]]
        """
        normal_axis = np.array([axis[0]/np.linalg.norm(axis[0]), axis[1]/np.linalg.norm(axis[1]), axis[2]/np.linalg.norm(axis[2])])
        self.theta_II_vectors = []
        self.phi_II_vectors = []
        for v in self.II_vectors:
            normal_v = v/np.linalg.norm(v)
            self.theta_II_vectors.append(np.arccos(np.dot(normal_v, normal_axis[2])))
            proj_v_a2 = np.dot(normal_v, normal_axis[2]) * normal_axis[2] 
            self.phi_II_vectors.append(np.arccos(np.dot((normal_v-proj_v_a2), normal_axis[0]))-0.5*np.pi)




class mesh_point(object):

    def __init__(self, index, coord, obj = None, i_obj = None):
        self.index = index
        self.coord = np.array(coord)
        # self.boundaries = np.array(boundaries)
        self.obj = obj
        self.i_obj = i_obj

    def set_obj(self, obj:object, c, i):
        if self.obj is not None:
            print(self.coord, self.obj.center_coord, i, c)
            # raise Exception("mesh_point already occupied")
            return False
        self.obj = obj
        self.i_obj = i
        return True

class Mesh(object):

    def __init__(self, cell, num, eps=0.0):
        self.mesh_size = num
        self.eps = eps
        self.mesh = []
        index = 0
        for i in range(num[0]):
            self.mesh.append([])
            for j in range(num[1]):
                self.mesh[i].append([])
                for k in range(num[2]):
                    self.mesh[i][j].append([])
                    mesh = mesh_point(index, [i,j,k])
                    self.mesh[i][j][k] = mesh
                    index += 1
        self.cell = cell/num
        self.supercell = cell
        
    def _pbc_mesh_coord(self, x, y, z):
        x = (x+self.mesh_size[0])%self.mesh_size[0]
        y = (y+self.mesh_size[1])%self.mesh_size[1]
        z = (z+self.mesh_size[2])%self.mesh_size[2]
        return x,y,z

    def get_mesh_point(self, _x, _y, _z):
        x,y,z = self._pbc_mesh_coord(_x,_y,_z)
        return self.mesh[x][y][z]
        
    def map_coord_mesh(self, coord):
        mesh_point_list = []
        c_list = np.vstack((coord[np.newaxis,:] +[[0, 0, 0]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[-self.eps, 0, 0]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[0, self.eps, 0]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[ 0,-self.eps, 0]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[0, 0, self.eps]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[0, 0,-self.eps]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[0,  self.eps, self.eps]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[0, -self.eps,-self.eps]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[0,  self.eps,-self.eps]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[0, -self.eps, self.eps]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[self.eps, -self.eps, 0]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[-self.eps, self.eps, 0]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[self.eps, self.eps, 0]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[-self.eps, -self.eps, 0]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[self.eps, self.eps, self.eps]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[-self.eps, -self.eps, -self.eps]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[self.eps, -self.eps, -self.eps]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[-self.eps, -self.eps, self.eps]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[-self.eps, self.eps, -self.eps]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[-self.eps, self.eps, self.eps]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[self.eps, -self.eps, self.eps]]))
        c_list = np.vstack((c_list,coord[np.newaxis,:] +[[self.eps, self.eps, -self.eps]]))
        for c in c_list:
            icoord = phys2Inter(c, self.cell)
            x = int(np.floor(icoord[0]))
            y = int(np.floor(icoord[1]))
            z = int(np.floor(icoord[2]))
            _x,_y,_z = self._pbc_mesh_coord(x,y,z)
            mesh_point_list.append(self.mesh[_x][_y][_z])
        return mesh_point_list, c_list

    def map_obj_mesh(self, obj):
        coord_list = []
        cart_coord_list = []
        for i in range(len(obj)):
            direct_coord = phys2Inter(obj[i].center_coord, self.cell)
            x = np.floor(direct_coord[0])
            y = np.floor(direct_coord[1])
            z = np.floor(direct_coord[2])
            _x,_y,_z = self._pbc_mesh_coord(x,y,z)
            # print(obj[i].center_coord, direct_coord, [x,y,z], [_x,_y,_z])
            icoord = np.array([_x,_y,_z])
            coord_list.append(icoord)
            cart_coord_list.append(obj[i].center_coord)
        coord_list = np.array(coord_list)
        order = np.lexsort([coord_list[:,2], coord_list[:,1], coord_list[:,0]])
        order_mesh = np.reshape(order, self.mesh_size)
        print(self.mesh_size)
        for i in range(self.mesh_size[0]):
            for j in range(self.mesh_size[1]):
                for k in range(self.mesh_size[2]):
                    print(i,j,k,order_mesh[i][j][k], cart_coord_list[order_mesh[i][j][k]], coord_list[order_mesh[i][j][k]])
                    mesh_point = self.mesh[i][j][k]
                    Succeed = mesh_point.set_obj(obj[order_mesh[i][j][k]], cart_coord_list[order_mesh[i][j][k]], order_mesh[i][j][k])
                    if not Succeed:
                        return False
        self.output_mesh()
        return True

    def read_obj_mesh(self, objmesh, obj):
        for m in range(objmesh.shape[0]):
            i = int(objmesh[m][0])
            j = int(objmesh[m][1])
            k = int(objmesh[m][2])
            i_obj = int(objmesh[m][3])
            Succeed = self.mesh[i][j][k].set_obj(obj[i_obj], obj[i_obj].center_coord, i_obj)
            if not Succeed:
                raise Exception("Fail to set_obj:: ", i_obj, i,j,k)

    def output_mesh(self):
        foutput = open("objmesh.dat", "w")
        for i in range(self.mesh_size[0]):
            for j in range(self.mesh_size[1]):
                for k in range(self.mesh_size[2]):
                    foutput.write("%d  %d  %d   %d\n"%(i,j,k,self.mesh[i][j][k].i_obj))
        foutput.close()
        


class Molecule(object):

    def __init__(self, indices, coords, cell):
        self.cell = cell
        self.coords_mol = np.array(coords)
        self.indices_mol = indices
        self.set_longaxis()
        self.set_polaraxis()
        self.set_center_postion(coords[0])

    def set_longaxis(self, vec=np.array([0.0,0.0,1.0]), length=2.32):
        self.unitlongaxis = np.array(vec)
        self.length_longaxis = length
        self.longaxis = length * vec
    
    def set_polaraxis(self, vec=np.array([0.0,1.0,0.0]), length=0.61):
        self.unitpolaraxis = np.array(vec)
        self.length_polaraxis = length
        self.polaraxis = length * vec

    def set_center_postion(self, coord):
        assert len(coord) == 3, "Pb coord must be a list of 3 numbers"
        self.center_coord = coord

    def set_angle_longaxis(self, axis=[0,0,1]):
        self.angle_longaxis = []
        self.angle_longaxis.append( math.acos(np.abs(np.dot(self.unitlongaxis, axis[0])/(np.linalg.norm(axis[0])))))
        self.angle_longaxis.append( math.acos(np.abs(np.dot(self.unitlongaxis, axis[1])/(np.linalg.norm(axis[1])))))
        self.angle_longaxis.append( math.acos(np.abs(np.dot(self.unitlongaxis, axis[2])/(np.linalg.norm(axis[2])))))

    def set_angle_polaraxis(self, axis=[0,0,1]):
        self.angle_polaraxis = []
        self.angle_polaraxis.append( math.acos(np.dot(self.unitpolaraxis, axis[0])/(np.linalg.norm(axis[0]))))
        self.angle_polaraxis.append( math.acos(np.dot(self.unitpolaraxis, axis[1])/(np.linalg.norm(axis[1]))))
        self.angle_polaraxis.append( math.acos(np.dot(self.unitpolaraxis, axis[2])/(np.linalg.norm(axis[2]))))

    def set_mesh_point(self, mesh):
        _C = phys2Inter(self.center_coord,self.cell)
        pC = apply_pbc(_C,self.cell)+0.5
        C = Inter2phys(pC,self.cell)
        self.mesh_point = mesh.map_coord_mesh(C)

    def rotate_long_axis(self, axis, angle, angle_in_radians=True):
        op = SymmOp.from_origin_axis_angle(
                    (0, 0, 0),
                    axis=tuple(axis),
                    angle=angle,
                    angle_in_radians=angle_in_radians
                )
        m = op.rotation_matrix
        centered_coords = np.zeros(self.coords_mol.shape)
        for i in range(len(self.coords_mol)):
            centered_coords[i] = get_distances(self.center_coord, self.coords_mol[i], self.cell, pbc=True)[0][0][0]
        new_coords = np.zeros(self.coords_mol.shape)
        for idx, c in enumerate(centered_coords):
            new_coords[idx] = np.dot(m, c.T).T
        for i in range(len(new_coords)):
            new_coords[i] += self.center_coord
        return new_coords

    def set_normal_vec(self):
        vec1 = distance(self.coords_mol[0], self.coords_mol[1], self.cell)
        vec2 = distance(self.coords_mol[0], self.coords_mol[2], self.cell)
        if sum(vec1<vec2) > 1:
            self.normal = np.cross(vec1, vec2)
        else:
            self.normal = np.cross(vec2, vec1)


class FAPbI3(object):

    def __init__(self, filename, fmt="lammps/dump", type_map=["I","Pb","C","N","H"], _name=None):
        self._name = _name
        self.setcutoff_CN_H()
        self.setcutoff_I_Pb()
        self.mols = []
        self.octahedra = []
        self.types = {}
        self.coords = {}
        self.indices_elem = {}
        
        if isinstance(filename, str):
            self._from_file(filename, fmt, type_map)
        else:
            self._from_obj(filename, type_map, fmt)
    
        self.cell = np.array(self.cubic["cells"][0])
        #invcell = np.linalg.inv(cell)
        #dcoord = np.matmul(coord, invcell)
        for elem_type in type_map:
            self.indices_elem[elem_type]= np.where(self.cubic["atom_types"]==self.types[elem_type])
            self.coords[elem_type] = self.cubic["coords"][0][self.indices_elem[elem_type][0]]
        
        self.cubic_6types = deepcopy(self.cubic)
        self.cubic_6types["atom_names"].append("Ha")
        self.cubic_6types["atom_numbs"].append(0)
        self.types_6types = deepcopy(self.types)
        self.types_6types["Ha"] = 5
        self.set_axis()


    def _from_file(self, filename, fmt, type_map):
        if fmt == "lammps/lmp":
            self.cubic = dpdata.System(filename, "lammps/lmp")
            for idx,elem_type in enumerate(type_map):
                self.cubic["atom_names"][idx] = elem_type
        else:
            if fmt == "lammps/dump":
                self.cubic = dpdata.System(filename, "lammps/dump")
                for idx,elem_type in enumerate(type_map):
                    self.cubic["atom_names"][idx] = elem_type    
            else:
                if fmt == "vasp/poscar":
                    self.cubic = dpdata.System(filename, "vasp/poscar")
                else:
                    raise Exception("Unknown format")
        for elem_type in type_map:
            self.types[elem_type] =  self.cubic["atom_names"].index(elem_type)
        assert self.cubic.get_nframes() == 1, print(self.cubic.get_nframes())

    def _from_obj(self, filename, type_map, fmt="vasp/poscar"):
        self.cubic = filename
        if fmt == "lammps/lmp":
            for idx,elem_type in enumerate(type_map):
                self.cubic["atom_names"][idx] = elem_type
        elif fmt == "lammps/dump":
            for idx,elem_type in enumerate(type_map):
                    self.cubic["atom_names"][idx] = elem_type
        elif fmt == "vasp/poscar":
            pass
        else:
            raise Exception("Unknown format")
        for elem_type in type_map:
            self.types[elem_type] =  self.cubic["atom_names"].index(elem_type)
        assert self.cubic.get_nframes() == 1, print(self.cubic.get_nframes())


    def set_axis(self, axis = np.eye(3)):
        self.axis = axis
        
    def startmesh(self, num, eps=0.0):
        # size_x = np.linalg.norm(self.cell[0])
        # size_y = np.linalg.norm(self.cell[1])
        # size_z = np.linalg.norm(self.cell[2])
        # self.mesh = Mesh([size_x, size_y, size_z], num, eps)
        self.mesh = Mesh(self.cell, num, eps)

    def assigh_molecule_to_mesh(self):
        for molecule in self.molecules:
            molecule.set_mesh_point(self.mesh)

    def assigh_oct_to_mesh(self):
        for oct in self.octahedra:
            oct.set_mesh_point(self.mesh)

    def rotate_moleculelongaxis_by_idx(self, idx_mol, axis, angle):
        new_coords = self.molecules[idx_mol].rotate_long_axis(axis, angle)
        for i in range(len(new_coords)):
            idx_atom = self.molecules[idx_mol].indices_mol[i]
            self.molecules[idx_mol].coords_mol[i] = new_coords[i]
            self.cubic["coords"][0][idx_atom] = new_coords[i]

    def extract_mol_from_indices(self, indices_molecules, moltype="FA"):
        cell = self.cell

        self.molecules = []
        for indices in indices_molecules:
            if moltype == "FA":
                coords = self.cubic["coords"][0][indices]
                vecNN = distance(coords[1], coords[2], cell)
                dNN = np.linalg.norm(vecNN)
                nvecNN = vecNN/dNN

                cNN = center(coords[1], coords[2], cell)
                vecCN = distance(coords[0], cNN, cell)
                dCN = np.linalg.norm(vecCN)
                nvecCN = vecCN/dCN

                molecule = Molecule(indices, coords, cell)
                molecule.set_longaxis(nvecNN, dNN)
                molecule.set_polaraxis(nvecCN, dCN)
                molecule.set_angle_longaxis(self.axis)
                molecule.set_angle_polaraxis(self.axis)
                self.molecules.append(molecule)
            elif moltype == "MA":
                coords = self.cubic["coords"][0][indices]
                vecCN = distance(coords[0], coords[1], cell)
                dCN = np.linalg.norm(vecCN)
                nvecCN = vecCN/dCN

                cCN = center(coords[0], coords[1], cell)

                molecule = Molecule(indices, coords, cell)
                #molecule.set_longaxis(nvecNN, dNN)
                molecule.set_polaraxis(nvecCN, dCN)
                #molecule.set_angle_longaxis(self.axis)
                molecule.set_angle_polaraxis(self.axis)
                molecule.set_center_postion(cCN)
                self.molecules.append(molecule)
            else:
                raise Exception("Unknown molecule type")


    def extract_mol(self, indices_molecules = None, moltype="FA"):
        if indices_molecules is not None:
            self.extract_mol_from_indices(indices_molecules)
            return indices_molecules
        cell = self.cell
        coords_C = self.coords["C"]
        coords_N = self.coords["N"]
        coords_H = self.coords["H"]
        list_C = self.indices_elem["C"]
        list_N = self.indices_elem["N"]
        list_H = self.indices_elem["H"]
        cubic_6types = self.cubic_6types

        idx = 0
        self.molecules = []
        indices_molecules = []
        if moltype == "FA":
            for idx_C,C in enumerate(coords_C):
                # stime = time.time()
                indices = []
                indices.append(list_C[0][idx_C])
                coords = []
                coords.append(C)
                for idx_N,N in enumerate(coords_N):
                    if len(indices) == 3:
                        break
                    d = get_distances(N, C, cell, pbc=True)[1][0][0]
                    if d < self.cutoff_CN_H["CN"]:
                        coords.append(N)
                        indices.append(list_N[0][idx_N])
                # Find Ha: H atoms bonded to C
                for idx_H,H in enumerate(coords_H):
                    if len(indices) == 4:
                        break
                    d = get_distances(H, C, cell, pbc=True)[1][0][0]
                    if d < self.cutoff_CN_H["CH"]:
                        coords.append(H)
                        indices.append(list_H[0][idx_H])
                        self.cubic_6types["atom_types"][indices[-1]] = self.types_6types["Ha"]
                        self.cubic_6types["atom_numbs"][self.types_6types["H"]] -= 1
                        self.cubic_6types["atom_numbs"][self.types_6types["Ha"]] += 1
                # Find H atoms bonded to N
                for idx_N in indices[1:3]:
                    N = self.cubic["coords"][0][idx_N]
                    if len(indices) == 8:
                        break
                    for idx_H,H in enumerate(coords_H):
                        d = get_distances(H, N, cell, pbc=True)[1][0][0]
                        if d < self.cutoff_CN_H["NH"]:
                            coords.append(H)
                            indices.append(list_H[0][idx_H])

                idx += 1
                # etime = time.time()
                # print("find molecule: ", etime-stime)
                # stime = etime
                assert len(indices) == 8, print(np.array(indices)+1) 

                vecNN = distance(coords[1], coords[2], cell)
                dNN = np.linalg.norm(vecNN)
                nvecNN = vecNN/dNN

                cNN = center(coords[1], coords[2], cell)
                vecCN = distance(coords[0], cNN, cell)
                dCN = np.linalg.norm(vecCN)
                nvecCN = vecCN/dCN

                molecule = Molecule(indices, coords, cell)
                molecule.set_longaxis(nvecNN, dNN)
                molecule.set_polaraxis(nvecCN, dCN)
                molecule.set_angle_longaxis(self.axis)
                molecule.set_angle_polaraxis(self.axis)
                self.molecules.append(molecule)
                # etime = time.time()
                # print("calculate vector: ", etime-stime)
                # stime = etime
                indices_molecules.append(indices)       
            assert(len(self.cubic_6types["coords"][0]) == sum(self.cubic_6types["atom_numbs"]))
        elif moltype == "MA":
            for idx_C,C in enumerate(coords_C):
                # stime = time.time()
                indices = []
                indices.append(list_C[0][idx_C])
                coords = []
                coords.append(C)
                for idx_N,N in enumerate(coords_N):
                    if len(indices) == 2:
                        break
                    d = get_distances(N, C, cell, pbc=True)[1][0][0]
                    if d < self.cutoff_CN_H["CN"]:
                        coords.append(N)
                        indices.append(list_N[0][idx_N])
                assert len(indices) == 2, print("C+N:", np.array(indices)+1) 
                # Find Ha: H atoms bonded to C
                for idx_H,H in enumerate(coords_H):
                    if len(indices) == 5:
                        break
                    d = get_distances(H, C, cell, pbc=True)[1][0][0]
                    if d < self.cutoff_CN_H["CH"]:
                        coords.append(H)
                        indices.append(list_H[0][idx_H])
                        self.cubic_6types["atom_types"][indices[-1]] = self.types_6types["Ha"]
                        self.cubic_6types["atom_numbs"][self.types_6types["H"]] -= 1
                        self.cubic_6types["atom_numbs"][self.types_6types["Ha"]] += 1
                assert len(indices) == 5, print("C+N+Ha:", np.array(indices)+1) 
                # Find H atoms bonded to N
                for idx_N in indices[1:2]:
                    N = self.cubic["coords"][0][idx_N]
                    if len(indices) == 8:
                        break
                    for idx_H,H in enumerate(coords_H):
                        if len(indices) == 8:
                            break
                        d = get_distances(H, N, cell, pbc=True)[1][0][0]
                        if d < self.cutoff_CN_H["NH"]:
                            coords.append(H)
                            indices.append(list_H[0][idx_H])

                idx += 1
                # etime = time.time()
                # print("find molecule: ", etime-stime)
                # stime = etime
                assert len(indices) == 8, print(np.array(indices)+1) 

                #vecNN = distance(coords[1], coords[2], cell)
                #dNN = np.linalg.norm(vecNN)
                #nvecNN = vecNN/dNN
                vecCN = distance(coords[0], coords[1], cell)
                dCN = np.linalg.norm(vecCN)
                nvecCN = vecCN/dCN

                cCN = center(coords[0], coords[1], cell)

                molecule = Molecule(indices, coords, cell)
                #molecule.set_longaxis(nvecNN, dNN)
                molecule.set_polaraxis(nvecCN, dCN)
                #molecule.set_angle_longaxis(self.axis)
                molecule.set_angle_polaraxis(self.axis)
                molecule.set_center_postion(cCN)
                self.molecules.append(molecule)
                # etime = time.time()
                # print("calculate vector: ", etime-stime)
                # stime = etime
                indices_molecules.append(indices)       
            assert(len(self.cubic_6types["coords"][0]) == sum(self.cubic_6types["atom_numbs"]))
        else:
            raise Exception("Unknown moltype")
        return indices_molecules

    def substitute_mol_by_idx(self,idx, new_atom_name="Cs"):
        coord = self.molecules[idx].center_coord
        new_FAPbI3 = self.remove_mol(idx)
        if new_atom_name in new_FAPbI3.cubic["atom_names"]:
            new_FAPbI3.cubic["atom_numbs"][list(new_FAPbI3.cubic["atom_names"]).index(new_atom_name)] += 1
            new_FAPbI3.cubic.data["atom_types"] = np.append(new_FAPbI3.cubic["atom_types"], list(new_FAPbI3.cubic["atom_names"]).index(new_atom_name))
            new_FAPbI3.cubic.data["coords"] = np.vstack((new_FAPbI3.cubic["coords"][0], np.array([coord])))[np.newaxis, :]
        else:
            new_FAPbI3.cubic.add_atom_names([new_atom_name])
            new_FAPbI3.cubic["atom_numbs"][-1] = 1
            new_FAPbI3.cubic.data["atom_types"] = np.append(new_FAPbI3.cubic["atom_types"], list(new_FAPbI3.cubic["atom_names"]).index(new_atom_name))
            new_FAPbI3.cubic.data["coords"] = np.vstack((new_FAPbI3.cubic["coords"][0], np.array([coord])))[np.newaxis, :]
        return new_FAPbI3

    def add_atom_names(self, atom_names):
        self.cubic.add_atom_names(atom_names)
        for a in atom_names:
            self.types[a] = self.cubic["atom_names"].index(a)


    def add_atom(self, atom_name, atom_coord):
        '''
        self.data = {}
        self.data['atom_numbs'] = []
        self.data['atom_names'] = []
        self.data['atom_types'] = []
        self.data['orig'] = np.array([0, 0, 0])
        self.data['cells'] = []
        self.data['coords'] = []
        '''
        if atom_name not in self.cubic["atom_names"]:
            self.add_atom_names([atom_name])
        atom_type = self.types[atom_name]
        self.cubic.data["atom_numbs"][atom_type] += 1
        self.cubic.data["atom_types"] = np.append(self.cubic.data["atom_types"], atom_type)
        self.cubic.data["coords"] = np.vstack((self.cubic.data["coords"][0], np.array([atom_coord])))[np.newaxis, :]

    def remove_mol(self, idx):
        removed_atom_idx = self.molecules[idx].indices_mol
        molecules = deepcopy(self.molecules)
        for i_mol, mol in enumerate(molecules):
            if i_mol == idx:
                continue
            for ridx in sorted(removed_atom_idx, reverse=True):
                mol.indices_mol = [idx-1 if idx > ridx else idx for idx in mol.indices_mol]
        picked_atom_idx = np.delete(np.arange(sum(self.cubic["atom_numbs"])), removed_atom_idx)
        new_sys = self.cubic.pick_atom_idx(picked_atom_idx)

        new_molecules = molecules
        new_molecules.pop(idx)

        new_FAPbI3 = self.copy()
        new_FAPbI3.cubic = new_sys
        new_FAPbI3.molecules = new_molecules
        return new_FAPbI3

    def remove_atom(self, idx):
        if idx > sum(self.cubic["atom_numbs"]):
            print(idx, sum(self.cubic["atom_numbs"]))
            raise Exception("Idx out of range")
        removed_atom_idx = idx
        picked_atom_idx = np.delete(np.arange(sum(self.cubic["atom_numbs"])), removed_atom_idx)
        new_sys = self.cubic.pick_atom_idx(picked_atom_idx)
        new_FAPbI3 = self.copy()
        new_FAPbI3.cubic = new_sys
        new_FAPbI3.types = self.types
        return new_FAPbI3

    def copy(self):
        """Returns a copy of the system.  """
        return self.__class__(deepcopy(self.cubic))


    def setcutoff_I_Pb(self, cutoff = 3.5):
        self.cutoff_I_Pb = cutoff

    def setcutoff_CN_H(self, cutoff ={"CH": 1.9, "NH": 1.6, "CN": 2.0}):
        self.cutoff_CN_H = cutoff

    def extract_octahedron_from_indices(self, indices_octahedra):
        cell = self.cell
        coords_I = self.coords["I"]
        coords_Pb = self.coords["Pb"]
        list_I = self.indices_elem["I"]
        list_Pb = self.indices_elem["Pb"]
        cubic = self.cubic

        idx = 0
        # coords_oct = []
        # nvec_oct = []
        # indices_oct = []
        self.octahedra = []
        for indices in indices_octahedra:
            coords = self.cubic["coords"][0][indices]
            oct = Octahedron(indices, coords, cell)
            oct.angle_PbI_vectors(self.axis)
            oct.setII_vectors(self.axis)
            oct.angle_II_vectors(self.axis)
            self.octahedra.append(oct)


    def extract_octahedron(self):
        # coords_I, list_I, coords_Pb, list_Pb, cubic
        cell = self.cell
        coords_I = self.coords["I"]
        coords_Pb = self.coords["Pb"]
        list_I = self.indices_elem["I"]
        list_Pb = self.indices_elem["Pb"]
        cubic = self.cubic

        idx = 0
        # coords_oct = []
        # nvec_oct = []
        # indices_oct = []
        self.octahedra = []
        indices_octahedra = []
        distances_PbI = get_distances(coords_Pb, coords_I, cell, pbc=True)

        for idx_Pb,Pb in enumerate(coords_Pb):
            indices = []
            indices.append(list_Pb[0][idx_Pb])
            coords = []
            coords.append(Pb)
            # Find I bonded to Pb
            for idx_I,d in enumerate(distances_PbI[1][idx_Pb]):
                I = coords_I[idx_I]
                # if len(indices) == 6:
                #     break
                if d < self.cutoff_I_Pb:
                    coords.append(I)
                    indices.append(list_I[0][idx_I])
            # assert len(indices) == 7, print(np.array(indices)+1) 
            idx += 1

            oct = Octahedron(indices, coords, cell)
            oct.angle_PbI_vectors(self.axis)
            oct.setII_vectors(self.axis)
            oct.angle_II_vectors(self.axis)
            self.octahedra.append(oct)
            indices_octahedra.append(indices)
        return indices_octahedra
    
    def dump_lammps_lmp(self, outfilename = "formated.lmp", formlist = None):
        fdump = open(outfilename, "w")
        fdump.write("Formated\n")
        fdump.write("%d atoms\n"%(self.cubic.get_natoms()))
        fdump.write("%d atom types\n"%(len(self.cubic["atom_names"])))
        lohi,tilt=box2lmpbox(self.cubic["orig"],self.cubic["cells"][0])
        fdump.write("%f %f xlo xhi\n"%(lohi[0][0], lohi[0][1]))
        fdump.write("%f %f ylo yhi\n"%(lohi[1][0], lohi[1][1]))
        fdump.write("%f %f zlo zhi\n"%(lohi[2][0], lohi[2][1]))
        fdump.write("%f %f %f xy xz yz\n"%(tilt[0], tilt[1], tilt[2]))
        fdump.write("\n")
        fdump.write("Atoms\n")
        fdump.write("\n")
        if formlist is None:
            for i in range(sum(self.cubic["atom_numbs"])):
                fdump.write("%d %d  %f %f %f\n"%(i+1,self.cubic["atom_types"][i]+1, self.cubic["coords"][0][i][0], self.cubic["coords"][0][i][1], self.cubic["coords"][0][i][2]))
        else:
            for i in range(sum(self.cubic["atom_numbs"])):
                fdump.write("%d %d  %f %f %f %f\n"%(i+1,self.cubic["atom_types"][i]+1, self.cubic["coords"][0][i][0], self.cubic["coords"][0][i][1], self.cubic["coords"][0][i][2], formlist[i]))
        fdump.close()

    
    def dump_lammps_dump(self, outfilename = "formated.dump", formlist = None, append="w", step = 0):
        fdump = open(outfilename, append)
        fdump.write("ITEM: TIMESTEP\n")
        fdump.write("%d\n"%(step))
        fdump.write("ITEM: NUMBER OF ATOMS\n")
        fdump.write("%d\n"%(self.cubic.get_natoms()))
        fdump.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
        lohi,tilt=box2dumpbox(self.cubic["orig"],self.cubic["cells"][0])
        fdump.write("%f %f %f  xlo xhi xy\n"%(lohi[0][0], lohi[0][1], tilt[0]))
        fdump.write("%f %f %f  ylo yhi xz\n"%(lohi[1][0], lohi[1][1], tilt[1]))
        fdump.write("%f %f %f  zlo zhi yz\n"%(lohi[2][0], lohi[2][1], tilt[2]))
        
        if formlist is None:
            fdump.write("ITEM: ATOMS id type x y z\n")
            for i in range(sum(self.cubic["atom_numbs"])):
                fdump.write("%d %d  %f %f %f\n"%(i+1,self.cubic["atom_types"][i]+1, self.cubic["coords"][0][i][0], self.cubic["coords"][0][i][1], self.cubic["coords"][0][i][2]))
        else:
            fdump.write("ITEM: ATOMS id type x y z format\n")
            for i in range(sum(self.cubic["atom_numbs"])):
                fdump.write("%d %d  %f %f %f %f\n"%(i+1,self.cubic["atom_types"][i]+1, self.cubic["coords"][0][i][0], self.cubic["coords"][0][i][1], self.cubic["coords"][0][i][2], formlist[i]))
        fdump.close()




def molecular_order_parameter_by_mesh(mesh, mollist, cell):
    mesh = mesh.T
    vecdot_x = 0.0
    vecdot_y = 0.0
    vecdot_z = 0.0
    corr_x = 0.0
    corr_y = 0.0
    corr_z = 0.0
    for idx_mol in range(len(mollist)):
        x_mol = mesh[0][idx_mol]
        y_mol = mesh[1][idx_mol]
        z_mol = mesh[2][idx_mol]
        
        if x_mol - 1 < 0:
            x_mol_left = max(mesh[0])
        else:
            x_mol_left = x_mol - 1
        if y_mol - 1 < 0:
            y_mol_left = max(mesh[1])
        else:
            y_mol_left = y_mol - 1
        if z_mol - 1 < 0:
            z_mol_left = max(mesh[2])
        else:
            z_mol_left = z_mol - 1

        if x_mol + 1 > max(mesh[0]):
            x_mol_right = 0
        else:
            x_mol_right = x_mol + 1
        if y_mol + 1 > max(mesh[1]):
            y_mol_right = 0
        else:
            y_mol_right = y_mol + 1
        if z_mol + 1 > max(mesh[2]):
            z_mol_right = 0
        else:
            z_mol_right = z_mol + 1

        xlimit = np.where(mesh[0]==x_mol)[0]
        ylimit = np.where(mesh[1]==y_mol)[0]
        zlimit = np.where(mesh[2]==z_mol)[0]
        
        xplane = np.concatenate((np.array(np.where(mesh[0]==x_mol_left)[0]), np.array(np.where(mesh[0]==x_mol_right)[0])))
        xnlist = np.intersect1d(np.intersect1d(xplane, ylimit), zlimit)
        
        yplane = np.concatenate((np.array(np.where(mesh[1]==y_mol_left)[0]), np.array(np.where(mesh[1]==y_mol_right)[0])))
        ynlist = np.intersect1d(np.intersect1d(yplane, xlimit), zlimit)
        
        zplane = np.concatenate((np.array(np.where(mesh[2]==z_mol_left)[0]), np.array(np.where(mesh[2]==z_mol_right)[0])))
        znlist = np.intersect1d(np.intersect1d(zplane, xlimit), ylimit)

        assert(len(xnlist)==2)
        assert(len(ynlist)==2)
        assert(len(znlist)==2)
        
        assert(np.abs(np.dot(mollist[idx_mol], mollist[xnlist[0]]))-1.0 < 1e-6)
        assert(np.abs(np.dot(mollist[idx_mol], mollist[ynlist[0]]))-1.0 < 1e-6)
        assert(np.abs(np.dot(mollist[idx_mol], mollist[znlist[0]]))-1.0 < 1e-6)
        assert(np.abs(np.dot(mollist[idx_mol], mollist[xnlist[1]]))-1.0 < 1e-6)
        assert(np.abs(np.dot(mollist[idx_mol], mollist[ynlist[1]]))-1.0 < 1e-6)
        assert(np.abs(np.dot(mollist[idx_mol], mollist[znlist[1]]))-1.0 < 1e-6)
        corr_x += math.acos(min(np.abs(np.dot(mollist[idx_mol], mollist[xnlist[0]])),1.0))
        corr_y += math.acos(min(np.abs(np.dot(mollist[idx_mol], mollist[ynlist[0]])),1.0))
        corr_z += math.acos(min(np.abs(np.dot(mollist[idx_mol], mollist[znlist[0]])),1.0))
        corr_x += math.acos(min(np.abs(np.dot(mollist[idx_mol], mollist[xnlist[1]])),1.0))
        corr_y += math.acos(min(np.abs(np.dot(mollist[idx_mol], mollist[ynlist[1]])),1.0))
        corr_z += math.acos(min(np.abs(np.dot(mollist[idx_mol], mollist[znlist[1]])),1.0))
        vecdot_x += (min(np.abs(np.dot(mollist[idx_mol], mollist[xnlist[0]])),1.0))
        vecdot_y += (min(np.abs(np.dot(mollist[idx_mol], mollist[ynlist[0]])),1.0))
        vecdot_z += (min(np.abs(np.dot(mollist[idx_mol], mollist[znlist[0]])),1.0))
        vecdot_x += (min(np.abs(np.dot(mollist[idx_mol], mollist[xnlist[1]])),1.0))
        vecdot_y += (min(np.abs(np.dot(mollist[idx_mol], mollist[ynlist[1]])),1.0))
        vecdot_z += (min(np.abs(np.dot(mollist[idx_mol], mollist[znlist[1]])),1.0))
          
    print(vecdot_x, vecdot_y, vecdot_z)
    return corr_x/2/math.pi/len(mollist), corr_y/2/math.pi/len(mollist), corr_z/2/math.pi/len(mollist)

def frame_order_parameter_by_mesh(mesh, octlist, cell):
    mesh = mesh.T

    corr_x = 1.0
    corr_y = 1.0
    corr_z = 1.0
    for idx_mol in range(len(octlist)):
        x_mol = mesh[0][idx_mol]
        y_mol = mesh[1][idx_mol]
        z_mol = mesh[2][idx_mol]
        
        if x_mol - 1 < 0:
            x_mol_left = max(mesh[0])
        else:
            x_mol_left = x_mol - 1
        if y_mol - 1 < 0:
            y_mol_left = max(mesh[1])
        else:
            y_mol_left = y_mol - 1
        if z_mol - 1 < 0:
            z_mol_left = max(mesh[2])
        else:
            z_mol_left = z_mol - 1

        if x_mol + 1 > max(mesh[0]):
            x_mol_right = 0
        else:
            x_mol_right = x_mol + 1
        if y_mol + 1 > max(mesh[1]):
            y_mol_right = 0
        else:
            y_mol_right = y_mol + 1
        if z_mol + 1 > max(mesh[2]):
            z_mol_right = 0
        else:
            z_mol_right = z_mol + 1

        xlimit = np.where(mesh[0]==x_mol)[0]
        ylimit = np.where(mesh[1]==y_mol)[0]
        zlimit = np.where(mesh[2]==z_mol)[0]
        
        xplane = np.concatenate((np.array(np.where(mesh[0]==x_mol_left)[0]), np.array(np.where(mesh[0]==x_mol_right)[0])))
        xnlist = np.intersect1d(np.intersect1d(xplane, ylimit), zlimit)
        
        yplane = np.concatenate((np.array(np.where(mesh[1]==y_mol_left)[0]), np.array(np.where(mesh[1]==y_mol_right)[0])))
        ynlist = np.intersect1d(np.intersect1d(yplane, xlimit), zlimit)
        
        zplane = np.concatenate((np.array(np.where(mesh[2]==z_mol_left)[0]), np.array(np.where(mesh[2]==z_mol_right)[0])))
        znlist = np.intersect1d(np.intersect1d(zplane, xlimit), ylimit)

        assert(len(xnlist)==2)
        assert(len(ynlist)==2)
        assert(len(znlist)==2)
        
        assert(np.dot(octlist[idx_mol][0], octlist[xnlist[0]][0])-1.0 < 1e-6)
        assert(np.dot(octlist[idx_mol][1], octlist[xnlist[0]][1])-1.0 < 1e-6)
        assert(np.dot(octlist[idx_mol][2], octlist[xnlist[0]][2])-1.0 < 1e-6)
        assert(np.dot(octlist[idx_mol][0], octlist[xnlist[1]][0])-1.0 < 1e-6)
        assert(np.dot(octlist[idx_mol][1], octlist[xnlist[1]][1])-1.0 < 1e-6)
        assert(np.dot(octlist[idx_mol][2], octlist[xnlist[1]][2])-1.0 < 1e-6)

        corr_x *= np.dot(octlist[idx_mol][1], octlist[xnlist[0]][1])
        corr_x *= np.dot(octlist[idx_mol][2], octlist[xnlist[0]][2])
        corr_x *= np.dot(octlist[idx_mol][1], octlist[xnlist[1]][1])
        corr_x *= np.dot(octlist[idx_mol][2], octlist[xnlist[1]][2])
        corr_y *= np.dot(octlist[idx_mol][0], octlist[ynlist[0]][0])
        corr_y *= np.dot(octlist[idx_mol][2], octlist[ynlist[0]][2])
        corr_y *= np.dot(octlist[idx_mol][0], octlist[ynlist[1]][0])
        corr_y *= np.dot(octlist[idx_mol][2], octlist[ynlist[1]][2])
        corr_z *= np.dot(octlist[idx_mol][0], octlist[znlist[0]][0])
        corr_z *= np.dot(octlist[idx_mol][1], octlist[znlist[0]][1])
        corr_z *= np.dot(octlist[idx_mol][0], octlist[znlist[1]][0])
        corr_z *= np.dot(octlist[idx_mol][1], octlist[znlist[1]][1])
          
    return math.pow(corr_x, 2/len(octlist)), math.pow(corr_y, 2/len(octlist)), math.pow(corr_z, 2/len(octlist))

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
    cubic.setcutoff_CN_H({"CH": 1.9, "NH": 1.6, "CN": 2.0})

    # start a mesh
    mesh_dim = [8,8,8]
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
    

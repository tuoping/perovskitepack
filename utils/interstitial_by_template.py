import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
import ase, ase.io
from pymatgen.core.operations import SymmOp
from ase.geometry.geometry import get_distances

def rotate_long_axis(mol, axis, angle, angle_in_radians=True):
    op = SymmOp.from_origin_axis_angle(
                (0, 0, 0),
                axis=tuple(axis),
                angle=angle,
                angle_in_radians=angle_in_radians
            )
    m = op.rotation_matrix
    new_coords = np.zeros(mol.shape)

    for idx, c in enumerate(mol):
        new_coords[idx] = np.dot(m, c.T).T

    for i in range(len(new_coords)):
        mol[i] = new_coords[i]

from ase import Atoms
def read_cp2k_restart(fname):
    lines = open(fname, 'r').readlines()
    flag_reading_coords = False
    for i in range(len(lines)):
        if "&CELL" == lines[i].strip():
            cell = []
            strings_line = lines[i+1].split()
            cell.append([float(x) for x in strings_line[1:]])
            strings_line = lines[i+2].split()
            cell.append([float(x) for x in strings_line[1:]])
            strings_line = lines[i+3].split()
            cell.append([float(x) for x in strings_line[1:]])  
        if "&COORD" in lines[i]:
            coords = []
            symbols = []
            flag_reading_coords = True
            continue
        if "&END COORD" in lines[i]:
            flag_reading_coords = False
        if "NUMBER_OF_ATOMS" in lines[i]:
            num_atoms = int(lines[i].split()[-1])
            break
        if  flag_reading_coords:
            strings_line = lines[i].split()
            symbols.append(strings_line[0])
            coords.append([float(x) for x in strings_line[1:]])

    assert num_atoms == len(coords), "  ".join([str(num_atoms), str(len(coords))])
    atoms = Atoms("".join(symbols),
                  positions = coords,
                  cell = cell,
                  pbc = [1,1,1])
    return atoms

frame_Pb = {
    'cubic': [],
    '2d': [],
    'hex': []
}

frame_Pb['cubic'] = np.array([
    [0, 0, 0],
])
frame_Pb['2d'] = np.array([
    [0, 0, 0],
])
frame_Pb['hex'] = np.array([
    [0, 0, 0],
    [0.46975, 0, 0],
    [0.22000, 0.5, 0.],
    [0.72067, 0.5, 0.],
    [0., 0., 0.5],
    [0.46545, 0., 0.5],
    [0.22423, 0.5, 0.5],
    [0.71121, 0.5, 0.5],

])


frame_interstitial = {
    'cubic': [],
    '2d': [],
    'hex': []
}

frame_interstitial['cubic'] = np.array([
    [0.5, 0.5, 0],
    [0, 0.5, 0.5],
    [0.5, 0, 0.5]
])
frame_interstitial['2d'] = np.array([
    [0.5, 0, 0],
])
frame_interstitial['hex'] = np.mod(np.array([
    [0.97717,    0.31731,    0.20065],
    [0.45542,    0.31407,    0.32450],
    [0.24907,    0.86111,    0.24642],
    [0.72762,    0.86185,    0.30436],
    [0.23894,    0.15948,    0.70733],
    [0.71689,    0.14218,    0.83866],
    [0.48616,    0.66276,    0.80475],
    [0.95445,    0.62964,    0.87989],      
]) + np.array([0, 0, 0.5]), 1.)

from copy import deepcopy
def shrink_frac_coord(frame, key, sc):
    frac_c = np.zeros([len(frame[key])*np.prod(sc),3])
    frac_c[:len(frame[key]), 0] = frame[key][:, 0]/sc[0]
    frac_c[:len(frame[key]), 1] = frame[key][:, 1]/sc[1]
    frac_c[:len(frame[key]), 2] = frame[key][:, 2]/sc[2]

    inv_sc = 1./sc
    idx = len(frame[key])
    for i in range(sc[0]):
        for j in range(sc[1]):
            for k in range(sc[2]):
                if i == 0 and j == 0 and k == 0:
                    continue
                frac_c[idx:idx+len(frame[key]), 0] = frame[key][:, 0]/sc[0] + inv_sc[0]*i
                frac_c[idx:idx+len(frame[key]), 1] = frame[key][:, 1]/sc[1] + inv_sc[1]*j
                frac_c[idx:idx+len(frame[key]), 2] = frame[key][:, 2]/sc[2] + inv_sc[2]*k
                idx += len(frame[key])
    return frac_c


crystal = read_cp2k_restart('_cubic-cell3/cubic-cell3-1.restart')
key = 'cubic'
sc = np.array([2, 2, 2])
target_num_inter = 6

pos = crystal.positions
cell = crystal.cell
atomic_numbers = np.array(crystal.get_atomic_numbers())
idx_Pb = np.where(atomic_numbers == 82)[0]

# build a grid (rows x cols)
frac_coord_lead_latt = shrink_frac_coord(frame_Pb, key, sc)
coord_lead_latt = frac_coord_lead_latt@crystal.cell.T

# lead_latt = Atoms(
#     "".join(['Pb']*len(frac_coord_lead_latt)),
#     positions = frac_coord_lead_latt@crystal.cell.T,
#     cell = crystal.cell,
#     pbc = [1,1,1]
#     )
# ase.io.write('lead_latt.xyz', lead_latt, format='extxyz')
origin = np.array( [crystal.positions[idx_Pb][:,0].min(), crystal.positions[idx_Pb][:,1].min(), crystal.positions[idx_Pb][:,2].min(), ])
print('origin = ', origin)
frac_coord_inter_latt = shrink_frac_coord(frame_interstitial, key, sc)
coord_inter_latt = frac_coord_inter_latt@crystal.cell.T + origin
lead_inter_latt = Atoms(
    "".join(['Pb']*len(frac_coord_lead_latt)+['O']*len(frac_coord_inter_latt)),
    positions = np.concatenate([frac_coord_lead_latt@crystal.cell.T, frac_coord_inter_latt@crystal.cell.T], axis=0),
    cell = crystal.cell,
    pbc = [1,1,1]
    )
ase.io.write('lead_inter_latt.xyz', lead_inter_latt, format='extxyz')
'''
# cost matrix: points x grid (use 'euclidean', 'sqeuclidean', or 'cityblock')
C = cdist(coord_lead_latt, pos[idx_Pb], metric="euclidean")  # or 'sqeuclidean' / 'cityblock'
# Hungarian algorithm (min-cost assignment)
row_ind, col_ind = linear_sum_assignment(C)  # use maximize=True to maximize scores
total_cost = C[row_ind, col_ind].sum()
print("Assignments (point_i -> grid_j):")
for i, j in zip(row_ind, col_ind):
    print(f"point {j} ({pos[idx_Pb][i]}) -> grid {i} ({coord_lead_latt[j]})  cost={C[i, j]:.3f}")
print("Total cost:", total_cost)
'''
import math
water = ase.io.read('water.vasp', format='vasp')
centered_water_pos = water.positions - np.mean(water.positions, axis=0)
coord_inter_water = []
symbols_inter = []
idx_inter_pos = np.arange(len(coord_inter_latt))
np.random.shuffle(idx_inter_pos) 
for i in idx_inter_pos:
    for k in range(10):
        rotaxis = np.random.randn(3)  # Random vector in 3D
        rotaxis = rotaxis / np.linalg.norm(rotaxis)  # Normalize to make it a unit vector
        rotangle = np.random.uniform(0, 360)
        rotate_long_axis(centered_water_pos, rotaxis, math.radians(rotangle))
        pos_ = coord_inter_latt[i] + centered_water_pos
        distmat = get_distances(crystal.positions, pos_, cell=crystal.cell, pbc=True)[1]
        distances = np.array([distmat[l2][l1] for l2 in range(distmat.shape[0]) for l1 in range(distmat.shape[1])])
        # print(pos_)
        # print(coord_inter_latt[i])
        if np.any(distances < 0.5):
            # print("Distance between water and crystal < 1.5", k, i)
            continue
        if len(coord_inter_water) < 1:
            break
        distmat = get_distances(coord_inter_water, pos_, cell=crystal.cell, pbc=True)[1]
        distances = np.array([distmat[l2][l1] for l2 in range(distmat.shape[0]) for l1 in range(l2+1, distmat.shape[1])])
        if np.any(distances < 0.5):
            # print("Distance between water and water < 1.5", k, i)
            continue
        # print("Found pos for water", k, i)
        break
    if k == 9:
        continue
    for j in range(len(pos_)):
        coord_inter_water.append(pos_[j])
        symbols_inter.append(water.get_chemical_symbols()[j])
    print(i,len(coord_inter_water), k)
    if len(coord_inter_water) >= target_num_inter*3:
        break
    print(i,len(coord_inter_water), k)
coord_inter_water = np.array(coord_inter_water)

dist_waters = get_distances(coord_inter_water[0], coord_inter_water[3], cell=crystal.cell, pbc=True)[1]
print("Distance between two waters = ", dist_waters)
formula = "".join(crystal.get_chemical_symbols() + symbols_inter)
out_atoms = Atoms(
    formula, 
    positions = np.concatenate([crystal.positions, coord_inter_water], axis=0),
    cell = crystal.cell,
    pbc = [1,1,1]
)

ase.io.write("output.xyz", out_atoms, format='extxyz')

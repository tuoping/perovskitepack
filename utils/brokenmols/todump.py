import dpdata 
from dpdata.lammps.lmp import box2lmpbox
from dpdata.lammps.dump import box2dumpbox

def dump_lammps_dump(outfilename, conf, ids, ve, types, step):
    fdump = open(outfilename, "w")
    fdump.write("ITEM: TIMESTEP\n")
    fdump.write("%d\n"%(step))
    fdump.write("ITEM: NUMBER OF ATOMS\n")
    fdump.write("%d\n"%(conf.get_natoms()))
    fdump.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
    lohi,tilt=box2dumpbox(conf["orig"],conf["cells"][0])
    fdump.write("%f %f %f  xlo xhi xy\n"%(lohi[0][0], lohi[0][1], tilt[0]))
    fdump.write("%f %f %f  ylo yhi xz\n"%(lohi[1][0], lohi[1][1], tilt[1]))
    fdump.write("%f %f %f  zlo zhi yz\n"%(lohi[2][0], lohi[2][1], tilt[2]))

    fdump.write("ITEM: ATOMS id type x y z vx vy vz\n")
    for i in range(conf.get_natoms()):
        idx_ve = ids.index(i+1)
        assert types[idx_ve] == conf["atom_types"][i]+1, print("fail at ", i, types[idx_ve], conf["atom_types"][i]+1)
        fdump.write("%d %d  %f %f %f  %f %f %f\n"%(i+1,conf["atom_types"][i]+1, conf["coords"][0][i][0], conf["coords"][0][i][1], conf["coords"][0][i][2], ve[idx_ve][0], ve[idx_ve][1], ve[idx_ve][2]))
    fdump.close()

def read_ve_dump(infilename):
    lines = open(infilename).readlines()
    step = int(lines[1].split()[0])
    ids = []
    ve = []
    types = []
    for line in lines[9:]:
        l = line.split()
        ids.append(int(l[0]))
        types.append(int(l[1]))
        ve.append([float(l[5]), float(l[6]), float(l[7])])
    return ids, ve, types, step


conf_new = dpdata.System("end_new.lmp", "lammps/lmp")
ids, ve, types, step = read_ve_dump("end.dump")
assert len(ids) == conf_new.get_natoms()
dump_lammps_dump("end_new_velocity.dump", conf_new, ids, ve, types, step)

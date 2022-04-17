import glob,os

logfile = open("log.lammps","r")
lines = logfile.readlines()
for idx in range(len(lines)-1, 0, -1):
    try:
        data.append(lines[idx])
    except:
        print("Skipping line "+str(idx))
    if "Step Temp PotEng" in lines[idx]:
        break
    if "Loop time" in lines[idx]:
        data = []
data.pop(-1)
of = open("thermo.dat", "w")
of.write("".join(data))

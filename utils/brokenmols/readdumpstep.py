import sys
lines = open(sys.argv[1], "r").readlines()
timestep = int(sys.argv[2])
outputfile = open(f"step{timestep}.dump", "w")
writing = False
for idx_line,line in enumerate(lines):
    if "TIMESTEP" in line:
        if writing: break
        current_timestep = int(lines[idx_line+1])
        if current_timestep == timestep: writing = True
    if writing: outputfile.write(line)
            
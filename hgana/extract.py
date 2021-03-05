import os
import sys


# Convert trr to gro at specified time
def convert_trr(time, orient):
    os.system("gmx_mpi trjconv -f run.trr -s run.tpr -o bound_o"+str(orient)+"_"+str(time).zfill(7)+"ps.gro -dump "+str(time)+"  >> extract.log 2>&1 <<EOF\n0\nEOF\n")


# Main methods
if __name__ == "__main__":
    # Define parameters
    dt = 2000
    com = 0.05
    orient = [0.35, 0.65]
    counter = [3, 3]
    convert_count = [0, 0]
    
    # Run through COLVAR
    with open("COLVAR", "r") as file_in:
        for line in file_in:
            if not "#" in line and not "@" in line:
                line_data = line.split()
                time = int(line_data[0].split(".")[0])
                if time % dt == 0:
                    if float(line_data[1]) < com:
                        if float(line_data[2]) < orient[0] and convert_count[0] < counter[0]:
                            convert_trr(time, 1)
                            convert_count[0] += 1
                        elif float(line_data[2]) > orient[1] and convert_count[1] < counter[1]:
                            convert_trr(time, 2)
                            convert_count[1] += 1

                if time % 100000 == 0:
                    sys.stdout.write("Finished frame "+"%7i"%time+"...\r")
                    sys.stdout.flush()

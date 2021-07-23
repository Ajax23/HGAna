################################################################################
# Extract Module                                                               #
#                                                                              #
"""Extract structures from COLAVR file for free energy simulations."""
################################################################################


import os
import sys
import shutil


def convert_trr(time, orient, out_link):
    """Convert trr to gro at specified time using gromacs.

    Parameters
    ----------
    time : integer
        Time in ps
    orient : string
        Orientation identifier
    """
    # Process output
    if not out_link[-1] == "/":
        out_link += "/"

    # Run gromacs
    if shutil.which("gmx_mpi"):
        os.system("gmx_mpi trjconv -f run.trr -s run.tpr -o "+out_link+"bound_o"+str(orient)+"_"+str(time).zfill(7)+"ps.gro -dump "+str(time)+"  >> "+out_link+"extract.log 2>&1 <<EOF\n0\nEOF\n")


def extract(file_link, out_link, dt=2000, com=0.05, orient=[0.35, 0.65], num_out=[3, 3]):
    """Extract complex structure from COLVAR file.

    Parameters
    ----------
    file_link : string
        File link to COLVAR file
    out_link : string
        Link to output folder
    dt : integer, optional
        Time step between possible outputs - trr output frequency in ps
    com : float, optional
        Distance from host com to guest com to consider
    orient : list, optional
        Orientation cutoffs to consider
    num_out : list, optional
        Number of output structures to generate
    """
    # Initialize
    convert_count = [0, 0]

    # Run through COLVAR
    with open(file_link, "r") as file_in:
        for line in file_in:
            if not "#" in line and not "@" in line:
                line_data = line.split()
                time = int(line_data[0].split(".")[0])
                if time % dt == 0:
                    if float(line_data[1]) < com:
                        if float(line_data[2]) < orient[0] and convert_count[0] < num_out[0]:
                            convert_trr(time, 1, out_link)
                            convert_count[0] += 1
                        elif float(line_data[2]) > orient[1] and convert_count[1] < num_out[1]:
                            convert_trr(time, 2, out_link)
                            convert_count[1] += 1

                if time % 100000 == 0:
                    sys.stdout.write("Finished frame "+"%7i"%time+"...\r")
                    sys.stdout.flush()

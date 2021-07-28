################################################################################
# Extract Module                                                               #
#                                                                              #
"""Extract structures from COLAVR file for free energy simulations."""
################################################################################


import os
import sys
import shutil
import numpy as np
import pandas as pd


def convert_gro(time, orient, out_link):
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
                            convert_gro(time, 1, out_link)
                            convert_count[0] += 1
                        elif float(line_data[2]) > orient[1] and convert_count[1] < num_out[1]:
                            convert_gro(time, 2, out_link)
                            convert_count[1] += 1

                if time % 100000 == 0:
                    sys.stdout.write("Finished frame "+"%7i"%time+"...\r")
                    sys.stdout.flush()


def restraints(file_link, out_link, conditions={1: [0.0, 0.7]}, atom_list={"a": 22, "b": 85, "c": 127, "A": 148, "B": 150, "C": 151}):
    """Extract restraints for decoupling simulations of the specified state
    within a specific error.

    Parameters
    ----------
    file_link : string
        File link to COLVAR file
    out_link : string
        Link to output folder
    conditions : dictionary, optional
        Dictionary conditions to check - key: [min, max]
    optional : dictionary, optional
        Atom ids for the sampled atoms in the COLVAR file as described by Boresch and Karplus

    Returns
    -------
    restraints : dictionary
        dictionary containing the restraints with names and values
    """
    # Read data
    with open(file_link, "r") as file_in:
        data = []
        for line in file_in:
            if not "#" in line:
                data.append([float(x) for x in line.split()])

    # Filter data
    filter = []
    for time_step in data:
        is_add = 0
        for col, cond in conditions.items():
            if time_step[col] >= cond[0] and time_step[col] <= cond[1]:
                is_add += 1
        if is_add==len(conditions):
            filter.append(time_step)

    filter_df = pd.DataFrame(filter)
    filter_df.columns = ["time", "d1", "d2", "d3", "a1", "r_aA", "theta_A", "theta_B", "phi_A", "phi_B", "phi_C"]

    # Extract distances and angles
    restraints = {"r_aA": {}, "theta_A": {}, "theta_B": {}, "phi_A": {}, "phi_B": {}, "phi_C": {}}

    for restraint in restraints:
        rest_data = filter_df[restraint] if restraint=="r_aA" else [x*180/np.pi for x in filter_df[restraint]]
        restraints[restraint] = round(np.mean(rest_data), 2)

    # Write to output file
    a = {atom_name: "%3i"%atom_id for atom_name, atom_id in atom_list.items()}
    d = {name: "%6.2f"%val for name, val in restraints.items()}

    with open(out_link, "w") as file_out:
        file_out.write("[ intermolecular_interactions ]\n")
        file_out.write("[ bonds ]\n")
        file_out.write("; ai   aj   type  bA    kA   bB    kB\n")
        file_out.write("   "+a["a"]+"  "+a["A"]+"  6     "+d["r_aA"]+"  0.0  "+d["r_aA"]+"  500\n")
        file_out.write("\n")

        file_out.write("[ angles ]\n")
        file_out.write("; ai   aj   ak   type  thA    fcA  thB    fcB\n")
        file_out.write("   "+a["b"]+"  "+a["a"]+"  "+a["A"]+"  1     "+d["theta_A"]+"  0.0  "+d["theta_A"]+"  50.0\n")
        file_out.write("   "+a["a"]+"  "+a["A"]+"  "+a["B"]+"  1     "+d["theta_B"]+"  0.0  "+d["theta_B"]+"  50.0\n")
        file_out.write("\n")

        file_out.write("[ dihedrals ]\n")
        file_out.write("; ai    aj   ak   al   type  thA     fcA  thB     fcB\n")
        file_out.write("   "+a["c"]+"  "+a["b"]+"  "+a["a"]+"  "+a["A"]+"  2     "+d["phi_A"]+"  0.0  "+d["phi_A"]+"  50.0\n")
        file_out.write("   "+a["b"]+"  "+a["a"]+"  "+a["A"]+"  "+a["B"]+"  2     "+d["phi_B"]+"  0.0  "+d["phi_B"]+"  50.0\n")
        file_out.write("   "+a["a"]+"  "+a["A"]+"  "+a["B"]+"  "+a["C"]+"  2     "+d["phi_C"]+"  0.0  "+d["phi_C"]+"  50.0\n\n\n")

    return restraints

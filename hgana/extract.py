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


def convert_gro(time, out_link):
    """Convert trr to gro at specified time using gromacs.

    Parameters
    ----------
    time : integer
        Time in ps
    out_link : string
        Link to output file
    """
    if shutil.which("gmx_mpi"):
        os.system("gmx_mpi trjconv -f run.trr -s run.tpr -o "+out_link+" -dump "+str(time)+" >> /dev/null 2>&1 <<EOF\n0\nEOF\n")


def extract(file_link, out_link=".", dt=2000, conditions={1: [0.0, 0.5], 2: [0.2, 0.4]}, num=0):
    """Extract complex structure from COLVAR file.

    Parameters
    ----------
    file_link : string
        File link to COLVAR file
    out_link : string, optional
        Link to output folder
    dt : integer, optional
        Time step between possible outputs - trr output frequency in ps
    conditions : dictionary, optional
        Dictionary conditions to check - key: [min, max]
    num : list, optional
        Number of output structures to generate

    Returns
    -------
    struct : dictionary
        Dictionary of all output structures with values for the conditions - time: [cond 1, ...]
    """
    # Run through COLVAR
    structs = {}
    num_convert = 0
    with open(file_link, "r") as file_in:
        for line in file_in:
            if not "#" in line and not "@" in line:
                time = int(line.split()[0].split(".")[0])
                line_data = [float(x) for x in line.split()]
                # Check time
                if time % dt == 0:
                    # Apply conditions
                    is_convert = 0
                    for col, cond in conditions.items():
                        if line_data[col] >= cond[0] and line_data[col] <= cond[1]:
                            is_convert += 1
                    if is_convert==len(conditions):
                        # Set output name
                        out_name = out_link
                        out_name += "/" if not out_link[-1]=="/" else ""
                        out_name += "complex_"
                        for col in conditions.keys():
                            out_name += "%.2f"%line_data[col]+"_"
                        out_name += str(time).zfill(7)+"ps.gro"

                        # Convert
                        convert_gro(time, out_name)
                        structs[time] = []
                        for col in conditions.keys():
                            structs[time].append(line_data[col])

                        # Check number of conversions
                        num_convert += 1
                        if num and num_convert==num:
                            break

                # Progress
                if time % 100000 == 0:
                    sys.stdout.write("Finished frame "+"%7i"%time+"...\r")
                    sys.stdout.flush()
        print()

    return structs


def restraints(file_link, out_link="restraints.top", conditions={1: [0.0, 0.7]}, atom_list={"a": 22, "b": 85, "c": 127, "A": 148, "B": 150, "C": 151}):
    """Extract restraints for decoupling simulations of the specified state
    within a specific error.

    Parameters
    ----------
    file_link : string
        File link to COLVAR file
    out_link : string, optional
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
    # Run through COLVAR data
    with open(file_link, "r") as file_in:
        data = []
        for line in file_in:
            if not "#" in line and not "@" in line:
                # Covert data
                line_data = [float(x) for x in line.split()]

                # Filter data
                is_add = 0
                for col, cond in conditions.items():
                    if line_data[col] >= cond[0] and line_data[col] <= cond[1]:
                        is_add += 1
                if is_add==len(conditions):
                    data.append(line_data)

    df = pd.DataFrame(data)
    df.columns = ["time", "d1", "d2", "d3", "a1", "r_aA", "theta_A", "theta_B", "phi_A", "phi_B", "phi_C"]

    # Extract distances and angles
    restraints = {"r_aA": {}, "theta_A": {}, "theta_B": {}, "phi_A": {}, "phi_B": {}, "phi_C": {}}

    for restraint in restraints:
        rest_data = df[restraint] if restraint=="r_aA" else [x*180/np.pi for x in df[restraint]]
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

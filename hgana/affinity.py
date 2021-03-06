################################################################################
# Analyze Class                                                                #
#                                                                              #
"""All necessary function for analysing the resulted data."""
################################################################################


import os
import sys
import math
import random

import numpy as np

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import hgana.utils as utils


def sample(link_colvar, link_out, com, orient, is_force=False):
    """This function samples the bound and unbound instance of a molecule and
    cyclodextrin. This function is to be run on the cluster due to a high time
    consumption. The output, a data object, is then used to calculate the
    binding addinity using functions :func:`number` and :func:`time`.

    The data structure is a dictionary containing the counted instances over
    the time periods.

    .. math::
        \\begin{bmatrix}
        \\text{u}&\\begin{pmatrix}t_{u,1}&\\dots&t_{u,N_u}\\end{pmatrix}\\\\
        \\text{b}&\\begin{pmatrix}t_{b,1}&\\dots&t_{b,N_b}\\end{pmatrix}\\\\
        \\text{b_1}&\\begin{pmatrix}t_{b,1}^1&\\dots&t_{b_1,N_b^1}\\end{pmatrix}\\\\
        \\text{b_2}&\\begin{pmatrix}t_{b,1}^2&\\dots&t_{b_2,N_b^2}\\end{pmatrix}
        \\end{bmatrix}

    with bound time length :math:`t_{b,i}` and unbound time length
    :math:`t_{u,i}` with number of bound instances :math:`N_b` and unbound ones
    :math:`N_u`. The bound instances are furthermore separated in orientations
    :math:`1` and :math:`2`.

    Parameters
    ----------
    link_colvar : string
        Link to COLVAR file generated by PLUMED
    link_out : string
        Link to output object file
    com : list
        List of index and cut-off value for the centre of mass entry in the
        COLVAR file
    orient : list
        List of index and cut-off value for the orientational entry in the
        COLVAR file
    is_force : bool, optional
        True to force re-extraction of data
    """

    # Check if already sampled
    if not os.path.exists(link_out) or is_force:

        # Extract data
        with open(link_colvar, "r") as file_in:
            # Initialize data structure
            # count_num = {"b": 0, "b_1": 0, "b_2": 0, "u": 0}
            count = {"b": [], "b_1": [], "b_2": [], "u": []}
            counter = {"b_u": 1, "o_1_2": 1}
            data_old = None
            state_old = None

            # Run through frames
            for frame_id, frame in enumerate(file_in):
                if not "#" in frame and not "@" in frame:
                    data = [float(x) if i>0 else int(x.split(".")[0]) for i, x in enumerate(frame.split())]

                    # Remove duplicates
                    if not data == data_old:
                        is_in = data[com[0]] < com[1]
                        is_or_1 = data[orient[0]] < orient[1] if is_in else None

                        # Ignore first frame
                        if state_old is not None:
                            # State is still the same
                            if is_in == state_old[0]:
                                # Add to bound/unbound counter
                                counter["b_u"] += 1
                                # Check if orientation is still the same in case molecule is still bound
                                if is_in:
                                    # Same orientation
                                    if is_or_1 == state_old[1]:
                                        # Add to orientation counter
                                        counter["o_1_2"] += 1
                                    # Orientation changed
                                    else:
                                        # Add instance to list of last orientation
                                        if state_old[1]:
                                            count["b_1"].append(counter["o_1_2"])
                                        else:
                                            count["b_2"].append(counter["o_1_2"])
                                        # Reset counter
                                        counter["o_1_2"] = 1
                            # Molecule state changed
                            else:
                                # Molecule was bound
                                if state_old[0]:
                                    # Add to bound list
                                    count["b"].append(counter["b_u"])
                                    # Add instance to list of last orientation
                                    if state_old[1]:
                                        count["b_1"].append(counter["o_1_2"])
                                    else:
                                        count["b_2"].append(counter["o_1_2"])
                                    # Reset counter
                                    counter["o_1_2"] = 1
                                # Molecule was unbound
                                else:
                                    # Add to unbound list
                                    count["u"].append(counter["b_u"])
                                # Reset counter
                                counter["b_u"] = 1

                        # Save last frame
                        data_old = data
                        state_old = [is_in, is_or_1]

                # Progress
                if frame_id%100000 == 0:
                    sys.stdout.write("Finished frame "+"%7i" % (frame_id)+"...\r")
                    sys.stdout.flush()

            print()

            # Add last counter to list
            if state_old[0]:
                # Add to bound list
                count["b"].append(counter["b_u"])
                # Add instance to list of last orientation
                if state_old[1]:
                    count["b_1"].append(counter["o_1_2"])
                else:
                    count["b_2"].append(counter["o_1_2"])
            # Molecule was unbound
            else:
                # Add to unbound list
                count["u"].append(counter["b_u"])


        # Define output dictionary
        inp = {"com": com, "orient": orient}
        output = {"count": count, "inp": inp}

        # Save data
        utils.save(output, link_out)

    # File already exists
    else:
        print("Object file already exists. If you wish to overwrite the file set the input *is_force* to True.")


def number(data_link, temp, volume, num_mol=1):
    """This function calculates the binding affinity :math:`\\Delta G_N`. This
    is done through a brute force summation of all instances the cyclodextrin
    is occupied :math:`N_b` and how often it is not :math:`N_u`.
    These numbers are then used to calculate the binding free energy

    .. math::

        \\Delta G_N=-RT\\ln\\frac{N_b}{N_u}-RT\\ln\\frac{V}{NV_0}

    with number of solute molecules :math:`N`, standard gas constant
    :math:`R`, temperature :math:`T`, simulation box volume :math:`V` and
    standard state volume :math:`V_0=1.661 \\text{nm}^3`.

    Parameters
    ----------
    data_link : string
        Sampled affinity object
    temp : float
        Simulated temperature in :math:`\\text{K}`
    volume : float
        Simulation box volume in :math:`\\text{m}^3`
    num_mol : integer, optional
        Number of molecules

    Returns
    -------
    table : DataFrame
        Pandas DataFrame containing the binding affinity calculated for the
        total number of instances and for both orientations. Values are given
        in :math:`\\frac{\\text{kJ}}{\\text{mol}}` and :math:`\\frac{\\text{kcal}}{\\text{mol}}`
    """
    # Load data object
    sample = utils.load(data_link)

    # Get data
    data = sample["count"]

    # Calculate variables
    RT = [-8.314e-3*temp, -1.986e-3*temp]  # kJ/mol, kcal/mol
    box_g = [x*np.log(volume/1.661e-27/num_mol) for x in RT]  # kJ/mol, kcal/mol

    # Calculate binding affinities
    delta_g = {}
    delta_g["dG"]    = [RT[i]*np.log(sum(data["b"])/  sum(data["u"]))+box_g[i] if sum(data["u"])>0 else box_g[i] for i in range(len(RT))]
    delta_g["dG_O1"] = [RT[i]*np.log(sum(data["b_1"])/sum(data["u"]))+box_g[i] if sum(data["u"])>0 else box_g[i] for i in range(len(RT))]
    delta_g["dG_O2"] = [RT[i]*np.log(sum(data["b_2"])/sum(data["u"]))+box_g[i] if sum(data["u"])>0 else box_g[i] for i in range(len(RT))]

    # Output
    return pd.DataFrame.from_dict(delta_g, orient="index", columns=["kJ/mol", "kcal/mol"])


def time(data_link, cutoff, temp, volume, num_mol=1, is_std=False):
    """This function calculates the binding affinity :math:`\\Delta G_T`.
    This is done by determining the association rate constant :math:`k_\\text{On}`
    and dissociation rate constant :math:`k_\\text{Off}`

    .. math::

        \\begin{array}{cc}
            k_\\text{On}=\\dfrac{1}{\\langle t_u\\rangle\\cdot C_g},&
            k_\\text{Off}=\\dfrac{1}{\\langle t_b\\rangle}
        \\end{array}

    with solute concentration :math:`C_g` of the free state in the complex

    .. math::

        C_g=\\frac{N}{V},

    number of solute molecules :math:`N`, box volume :math:`V`,
    average bound time :math:`\\langle t_b\\rangle` and average unbound time
    :math:`\\langle t_u\\rangle`. The two averages are calculated by
    determining the time of the host molecule being in a bound state
    before changing to an unbound one and vice versa. The resulting bound
    :math:`t_b` and unbound :math:`t_u` time instances, with instance
    quantities :math:`M_b` and :math:`M_u`, are then filtered by
    a minimal dwelling time :math:`c_{min}` resulting in the time
    average by normalizing the sum of all instances by their count

    .. math::

        \\begin{array}{cc}
            \\langle t_b\\rangle=\\dfrac1{M_b}\\sum_{i=1}^{M_b}t_{b,i},&
            \\langle t_u\\rangle=\\dfrac1{M_u}\\sum_{i=1}^{M_u}t_{u,i}.
        \\end{array}

    Finally the binding affinity results from

    .. math::

        \\Delta G_T&=-RT\\ln\\frac{k_\\text{On}C_0}{k_\\text{Off}}
        =-RT\\ln\\frac{\\langle t_b\\rangle}{\\langle t_u\\rangle}-RT\\ln\\frac{C_0}{C_g}\\\\
        &=-RT\\ln\\frac{\\langle t_b\\rangle}{\\langle t_u\\rangle}-RT\\ln\\frac{V}{NV_0}

    with standard gas constant :math:`R`, temperature :math:`T`, standard
    state concentration :math:`C_0=1\\frac{\\text{mol}}{\\text{l}}=V_0^{-1}` and standard
    state volume :math:`V_0=1.661nm`.

    The standard deviation can be determined by creating permutations
    containing a percentage of elements from the time arrays :math:`t_{b,i}`
    and :math:`t_{u,j}`

    .. math::
        \\begin{array}{cc}
            \\boldsymbol{P}_u=\\sum_k^{N_p}p_{u,k},&
            p_{u,k}=[t_{u,0},\\dots,t_{u,N_u-x}]
        \\end{array}

    with permutation matrix :math:`\\boldsymbol{P}`, permutation :math:`p`,
    number of permutations :math:`N_p` and element percentage :math:`1-x`.
    This permutation matrix is then used to calculate multiple association
    and dissociation rates creating a pool of values, of which the
    standard deviation can be determined using

    .. math::

        \\text{std}(k_\\text{On})=\\sqrt{\\frac{\\sum_{k=1}^{N_p}\\left(k_{\\text{On},k}-\\bar k_\\text{On}\\right)^2}{N_p-1}}

    with mean value :math:`\\bar k_\\text{On}` of the association rates
    calculated from the permutation matrix. Similarly, the standard deviation of
    the binding affinity can be determined by

    .. math::

        \\text{std}(\\Delta G_T)=\\sqrt{\\frac{\\sum_{k=1}^{N_p}\\left(\\Delta G_T^k-\\Delta\\bar G_T\\right)^2}{N_p-1}}

    with mean value :math:`\\Delta\\bar G_T`.

    Parameters
    ----------
    data_link : string
        Sampled affinity object
    cutoff : integer
        Time calculation cut-off in :math:`\\text{ps}`
    temp : float
        Simulated temperature in :math:`\\text{K}`
    volume : float
        Simulation box volume in :math:`\\text{m}^3`
    num_mol : integer, optional
        Number of molecules
    is_std : bool, optional
        True to calculate standard deviation

    Returns
    -------
    table : DataFrame
        Pandas DataFrame of binding affinity in :math:`\\frac{\\text{kJ}}{\\text{mol}}`
        and :math:`\\frac{\\text{kcal}}{\\text{mol}}`,
        :math:`k_\\text{On}` in :math:`\\frac{1}{\\text{s}}` and
        :math:`k_\\text{Off}` in
        :math:`\\frac{\\text{dm}^3}{\\text{mol}\\cdot\\text{s}}`, and optionally
        the standard deviations of these values
    """
    # Load data object
    sample = utils.load(data_link)

    # Get data
    data = sample["count"]

    # Process input
    N = num_mol
    V = volume

    # Calculate variables
    NA = 6.022e23   # #/mol
    V0 = 1.661e-27  # m^3
    C0 = 1/V0       # #/m^3
    Cg = N/V        # #/m^3
    RT = [-8.314e-3*temp, -1.986e-3*temp]  # kJ/mol, kcal/mol
    box_g = [x*np.log(C0/Cg) for x in RT]

    # Apply cuttoff
    data = {x: [y for y in data[x] if y >= cutoff] for x in data}

    # Calculate mean time value
    mean_b = sum(data["b"])/len(data["b"])*10e-12 if len(data["b"]) > 0 else 0
    mean_u = sum(data["u"])/len(data["u"])*10e-12 if len(data["u"]) > 0 else 0

    # Calculate the standard deviation
    if is_std:
        percent = 0.8
        set_size = 500
        ids = [[y for y in range(len(data[x]))] for x in ["b", "u"]]
        rand = [[], []]

        for i in range(set_size):
            random.shuffle(ids[0])
            random.shuffle(ids[1])

            rand[0].append([data["b"][ids[0][k]] for k in range(len(ids[0])) if k < math.floor(len(ids[0])*percent)])
            rand[1].append([data["u"][ids[1][k]] for k in range(len(ids[1])) if k < math.floor(len(ids[1])*percent)])

        # Calculate mean time
        mean_rand_b = [sum(x)/len(x)*10e-12 if len(x) > 0 else 0 for x in rand[0]]
        mean_rand_u = [sum(x)/len(x)*10e-12 if len(x) > 0 else 0 for x in rand[1]]

        # Append rate
        k_rand_b = [1/x/Cg*NA*1e3 if x > 0 else 0 for x in mean_rand_u]  # dm^3/mol/s
        k_rand_u = [1/x if x > 0 else 0 for x in mean_rand_b]
        d_g_rand_kj = [RT[0]*np.log(mean_rand_b[i]/mean_rand_u[i])+box_g[0] if mean_rand_u[i] > 0 else 0 for i in range(set_size)]
        d_g_rand_kc = [RT[1]*np.log(mean_rand_b[i]/mean_rand_u[i])+box_g[1] if mean_rand_u[i] > 0 else 0 for i in range(set_size)]

    # Initialize output variable
    delta_g = [cutoff]
    idx = ["Cutoff (ps)"]

    # Calculate delta g
    delta_g.append(RT[0]*np.log(mean_b/mean_u)+box_g[0] if mean_u > 0 else 0)  # kJ/mol
    idx.append("dG (kJ/mol)")
    if is_std:
        delta_g.append(np.std(d_g_rand_kj))
        idx.append("dG_std (kJ/mol)")

    delta_g.append(RT[1]*np.log(mean_b/mean_u)+box_g[1] if mean_u > 0 else 0)  # kJ/mol
    idx.append("dG (kcal/mol)")
    if is_std:
        delta_g.append(np.std(d_g_rand_kc))
        idx.append("dG_std (kcal/mol)")

    # Calculate association rate k_on
    delta_g.append(1/mean_u/Cg*NA*1e3 if mean_u > 0 else 0)  # dm^3/mol/s
    idx.append("k_on (dm^3/mol/s)")
    if is_std:
        delta_g.append(np.std(k_rand_b))
        idx.append("k_on_std (dm^3/mol/s)")

    # Calculate dissociation rate k_off
    delta_g.append(1/mean_b if mean_b > 0 else 0)  # 1/s
    idx.append("k_off (1/s)")
    if is_std:
        delta_g.append(np.std(k_rand_u))
        idx.append("k_off_std (1/s)")

    # Return table
    return pd.DataFrame([delta_g], columns=idx)


def time_series(link_colvar, com, dt=100):
    """This function plots the bound and unbound instances over time for a given
    spatial bound cutoff.

    Parameters
    ----------
    link_colvar : string
        Link to COLVAR file generated by PLUMED
    com : list
        List of index and cut-off value for the centre of mass entry in the
        COLVAR file
    dt : integer, optional
        Step size in ps between frames to plot
    """
    # Extract data
    with open(link_colvar, "r") as file_in:
        state = []
        data_old = None

        # Run through frames
        for frame_id, frame in enumerate(file_in):
            if not "#" in frame and not "@" in frame:
                data = [float(x) if i>0 else int(x.split(".")[0]) for i, x in enumerate(frame.split())]

                # Remove duplicates
                if not data == data_old:
                    if frame_id%dt == 0:
                        state.append(int(data[com[0]] < com[1]))

                    # Save last frame
                    data_old = data

            # Progress
            if frame_id%100000 == 0:
                sys.stdout.write("Finished frame "+"%7i" % (frame_id)+"...\r")
                sys.stdout.flush()

        print()

        sns.lineplot(x=[x*dt/1000000 for x in range(len(state))], y=state)

        # Set labels
        plt.xlabel(r"Time ($\mu$s)")
        plt.ylabel("State (Bound - 1, Unbound - 0)")


def hist(link_colvar, column, legend=None, conditions=None):
    """This function plots a histogram of the provided COLVAR file. For given
    columns.

    Conditions are dictionaries with the key being the column to which the
    condition is applied. The value of each entry is a list containing on the
    first entry the column id and the second entry, the cutoff value.

    Parameters
    ----------
    link_colvar : string
        Link to COLVAR file generated by PLUMED
    column : integer, list
        column id or list of ids to plot the histogram for
    legend : list, string, None, optional
        Legend entry or a list of those
    conditions : dictionary, None, optional
        Conditions for plotting

    Examples
    --------
    Histogram with condition

    .. code-block:: python

        hist("COLVAR", [1, 2], conditions={2: [1, 0.5])
    """
    # Process user input
    columns = column if isinstance(column, list) else [column]
    if legend is not None:
        legend = legend if isinstance(legend, list) else [legend]

    # Extract data
    data = {x: [] for x in columns}
    with open(link_colvar, "r") as file_in:
        # Run through file
        for line in file_in:
            # Remove comments
            if not "#" in line:
                line_data = [float(x) for x in line.split()]
                # Run through requested columns
                for column in columns:
                    is_append = True
                    # Check if condition is fulfilled if given
                    if conditions is not None and column in conditions:
                        condition = conditions[column]
                        is_append = line_data[condition[0]] < condition[1]
                    if is_append:
                        data[column].append(line_data[column])

    # Plot
    plot_data = {legend[col_id]: data[column] for col_id, column in enumerate(columns)}
    sns.histplot(plot_data, stat="density", bins="rice", kde=True)

    # Set labels
    plt.xlabel("Distance (nm)")
    plt.ylabel("Normalized quantity")

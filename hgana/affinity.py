################################################################################
# Addinity Module                                                              #
#                                                                              #
"""All necessary function for analysing the host-guest simulation."""
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


def sample(file_link, out_link, conditions={1: [0, 0.7]}):
    """This function samples the bound and unbound instance of a molecule and
    cyclodextrin. The output, a data object, is then used to calculate the
    binding affinity using functions :func:`number` and :func:`time`.

    The data structure is a dictionary containing the counted instances over
    the time periods.

    Parameters
    ----------
    file_link : string
        Link to parameter file generated by PLUMED
    out_link : string
        Link to output object file
    conditions : dictionary, optional
        Dictionary conditions to check - key: [min, max]
    """
    # Initialize
    counter = 0
    series = {"b": [], "u": []}

    data_old = []
    state_old = None

    # Extract data
    with open(file_link, "r") as file_in:
        # Run through frames
        for frame_id, frame in enumerate(file_in):
            # Remove comments
            if not "#" in frame and not "@" in frame:
                data = [float(x) if i>0 else int(x.split(".")[0]) for i, x in enumerate(frame.split())]
                # Remove duplicates
                if not data == data_old:
                    # Process conditions
                    is_in = 0
                    for col, cond in conditions.items():
                        if data[col] >= cond[0] and data[col] <= cond[1]:
                            is_in += 1
                    is_in = is_in==len(conditions)

                    # Check state - remove first frame
                    if state_old is not None:
                        # Same state - (still bound/still unbound)
                        if is_in == state_old:
                            counter += 1
                        # State changed
                        else:
                            # Add counter to time series
                            series["b" if state_old else "u"].append(counter)
                            # Reset counter
                            counter = 1

                # Save last frame
                data_old = data
                state_old = is_in

            # Progress
            if frame_id%100000 == 0:
                sys.stdout.write("Finished frame "+"%7i"%(frame_id)+"...\r")
                sys.stdout.flush()
        print()

        # Add last counter to list
        series["b" if state_old else "u"].append(counter)

    # Save data
    utils.save({"inp": conditions, "series": series}, out_link)


def number(data_link, T, V):
    """This function calculates the binding affinity :math:`\\Delta G_N`. This
    is done through a brute force summation of all instances the cyclodextrin
    is occupied :math:`N_b` and how often it is not :math:`N_u`.
    These numbers are then used to calculate the binding free energy

    .. math::

        \\Delta G_N=-RT\\ln\\frac{N_b}{N_u}-RT\\ln\\frac{V}{V_0}

    with standard gas constant :math:`R`, temperature :math:`T`, simulation box
    volume :math:`V` and standard state volume :math:`V_0=1.661 \\text{nm}^3`.

    Parameters
    ----------
    data_link : string
        Sampled affinity object
    T : float
        Simulated temperature in :math:`\\text{K}`
    V : float
        Simulation box volume in :math:`\\text{m}^3`

    Returns
    -------
    dG : dictionary
        Dictionary containing the binding affinity calculated for the
        total number of instances. Values are given in
        :math:`\\frac{\\text{kJ}}{\\text{mol}}` **kJ/mol** and
        :math:`\\frac{\\text{kcal}}{\\text{mol}}` **kcal/mol**
    """
    # Load data object
    sample = utils.load(data_link)

    # Get data
    data = sample["series"]

    # Calculate variables
    RTs = {"kJ/mol": -8.314e-3*T, "kcal/mol": -1.986e-3*T}  # kJ/mol, kcal/mol
    V0 = 1.661e-27  # m^3

    N_b = sum(data["b"])
    N_u = sum(data["u"])

    log_V = np.log(V/V0)
    log_N = np.log(N_b/N_u) if N_u>0 else 0

    # Calculate binding affinity
    dG = {key: RT*log_N+RT*log_V for key, RT in RTs.items()}

    # Output
    return dG


def time(data_link, T, V, dt, len_frame=1, is_std=False):
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
    T : float
        Simulated temperature in :math:`\\text{K}`
    V : float
        Simulation box volume in :math:`\\text{m}^3`
    dt : integer
        Time calculation cut-off in ps
    len_frame : float, optional
        Length of a frame in ps
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
    data = sample["series"]

    # Set constants
    NA = 6.022e23   # #/mol
    V0 = 1.661e-27  # m^3
    RTs = {"kJ/mol": -8.314e-3*T, "kcal/mol": -1.986e-3*T}  # kJ/mol, kcal/mol
    log_V = np.log(V/V0)

    # Apply cuttoff
    d_frames = round(dt/len_frame)
    data = {x: [y for y in data[x] if y >= d_frames] for x in data}

    # Calculate mean time value
    t_b = np.mean(data["b"])*10e-12
    t_u = np.mean(data["u"])*10e-12

    # Calculate binding affinity
    log_t = np.log(t_b/t_u) if t_u>0 else 0
    k_on = V/t_u*NA*1e3 if t_u>0 else 0
    k_off = 1/t_b
    dG = {key: RT*log_t+RT*log_V for key, RT in RTs.items()}

    # Calculate the standard deviation
    if is_std:
        # Initialize
        ratio = 0.8
        num_std = 300
        num_entry = {key: math.floor(ratio*len(val)) for key, val in data.items()}

        # Generate permutations
        t_std = {key: [np.mean(random.sample(val, num_entry[key])) for i in range(num_std)] for key, val in data.items()}

        # Calculate standard deviation
        log_t_std = [np.log(t_std["b"][i]/t_std["u"][i]) for i in range(num_std) if t_std["u"][i]>0]
        k_on_std = np.std([V/t_std["u"][i]*NA*1e3 for i in range(num_std) if t_std["u"][i]>0])
        k_off_std = np.std([1/t_std["b"][i] for i in range(num_std)])
        dG_std = {key: np.std([RT*log_t_std[i]+RT*log_V for i in range(num_std)]) for key, RT in RTs.items()}

    # Initialize output variable
    output = [dt]
    idx = ["Cutoff (ps)"]
    vals = [dG["kJ/mol"], dG["kcal/mol"], k_on, k_off]
    vals_std = [dG_std["kJ/mol"], dG_std["kcal/mol"], k_on_std, k_off_std] if is_std else []
    units = ["kJ/mol", "kcal/mol", "dm^3/mol/s", "1/s"]
    names = ["dG", "dG", "k_on", "k_off"]

    # Process table
    for i in range(len(vals)):
        output.append(vals[i])
        idx.append(names[i]+" ("+units[i]+")")
        if is_std:
            output.append(vals_std[i])
            idx.append(names[i]+"_std ("+units[i]+")")

    # Return table
    return pd.DataFrame([output], columns=idx)


def double_decoupling(temp, dG_hydr, dG_host, restraints, forces={"K_r": 500, "K_theta_A" : 50, "K_theta_B" : 50, "K_phi_A" : 50, "K_phi_B" : 50, "K_phi_C" : 50}):
    """Calculate the double decoupling free enthalpy difference according to

    .. math::

        \\Delta G^\\text{DD}=\\Delta G_{\\mathrm{u\\rightarrow b}}^\\mathrm{M}=-\\Delta G_\\mathrm{hyd}^\\mathrm{M}-\\Delta G_{\\mathrm{b\\rightarrow tor}}^{\\mathrm{M\\rightarrow M'}}-\\Delta G_{\\mathrm{tor}}^{\\mathrm{M'}}.

    Hereby :math:`\\Delta G_\\mathrm{hyd}^\\mathrm{M}` is the hydration enthalpy
    energy from the unbound state, :math:`\\Delta G_{\\mathrm{b\\rightarrow tor}}^{\\mathrm{M\\rightarrow M'}}`
    the free enthalpy difference from turning on orientational and translational
    restrants (tor) and deactivating the intramolecular interactions of the
    guest molecule with its surroundings. The free enthalpy differnece due to
    restraints :math:`\\Delta G_{\\mathrm{tor}}^{\\mathrm{M'}}` is then
    subtracted. This contribution can be calculated analytaically as proposed
    by Boresch and Karplus

    .. math::

            \\Delta G_{\\mathrm{tor}}^{\\mathrm{M'}}=-RT\\ln\\left[\\frac{8\\pi V^0(K_rK_{\\theta_A}K_{\\theta_B}K_{\\phi_A}K_{\\phi_B}K_{\\phi_C})^\\frac{1}{2}}{r_{aA,0}^2\\sin\\theta_{A,0}\\sin\\theta_{B,0}(2\\pi RT)^3}\\right].

    Parameters
    ----------
    temp : float
        System temperature in K
    dG_hydr : float
        Hydration free enthalpy :math:`\\Delta G_\\mathrm{hyd}^\\mathrm{M}` in kJ/mol
    dG_host : float
        Free enthalpy difference from turning on restraints and deactivating
        intramolecular interactions :math:`\\Delta G_{\\mathrm{b\\rightarrow tor}}^{\\mathrm{M\\rightarrow M'}}`
        in kJ/mol
    restraints : dictionary
        Restraints applied during simulation
    forces : dictionary, optional
        Restraint forces applied during simulation

    Returns
    -------
    dG : dictionary
        Dictionary containing the binding free enthalpy **dG_tot** and the
        analytical free enthalpy difference due to restraints **dG_rest**
    """
    # Initialize
    V = 1.661  # nm^3
    R = 8.314e-3   # kJ/K/mol
    T = temp  # K
    pi = np.pi
    to_rad = pi/180

    # Set restraints
    r_aA_0    = restraints["r_aA"]         # nm
    theta_A_0 = restraints["theta_A"]*to_rad  # rad
    theta_B_0 = restraints["theta_B"]*to_rad  # rad

    # Set forces
    K_r = forces["K_r"]  # kJ/mol/nm^2
    K_theta_A = forces["K_theta_A"]  # kJ/mol/rad^2
    K_theta_B = forces["K_theta_B"]  # kJ/mol/rad^2
    K_phi_A = forces["K_phi_A"]  # kJ/mol/rad^2
    K_phi_B = forces["K_phi_B"]  # kJ/mol/rad^2
    K_phi_C = forces["K_phi_C"]  # kJ/mol/rad^2

    # Calculate restraint
    dG_rest = -R*T*np.log((8*pi**2*V*(K_r*K_theta_A*K_theta_B*K_phi_A*K_phi_B*K_phi_C)**0.5)/(r_aA_0**2*np.sin(theta_A_0)*np.sin(theta_B_0)*(2*pi*R*T)**3))

    # Calculate binding free enthalpy
    dG_tot = -dG_hydr-dG_host-dG_rest

    return {"dG_tot": dG_tot, "dG_rest": dG_rest}


def plot_time(link_colvar, com, dt=100, kwargs={}):
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
        Step size between frames
    kwargs: dict, optional
        Dictionary with plotting parameters (only for given intent)
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

        sns.lineplot(x=[x*dt/1000000 for x in range(len(state))], y=state, **kwargs)

        # Set labels
        plt.xlabel(r"Time ($\mu$s)")
        plt.ylabel("State (Bound - 1, Unbound - 0)")

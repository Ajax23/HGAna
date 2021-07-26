################################################################################
# Adsorption Class                                                             #
#                                                                              #
"""Parallelization of the MC algorithm to calcuilate adsorption isotherms."""
################################################################################


import math
import numpy as np
import seaborn as sns
import itertools as it
import multiprocessing as mp
import matplotlib.pyplot as plt

import hgana.utils as utils

from hgana.box import Box
from hgana.mc import MC


class Adsorption(Box):
    """This class run calculates adsorption isotherms using the MC algorithm.

    Parameters
    ----------
    size : list, optional
        Three dimensional cell count in each dimension - leave empty if only
        plotting is intended.
    """
    def __init__(self, size=[0, 0, 0]):
        super(Adsorption, self).__init__(size)

    def add_mol(self, num, is_move=True, name="", struct=""):
        """Add a molecule to the system. The interaction matrix is set up and
        filled with empty values.

        Parameters
        ----------
        num : integer, list
            Number of molecules to be added
        is_move : bool, optional
            True to move molecule during simulation
        name : string, optional
            Molecule name
        struct : string, optional
            File link to structure file
        """
        # Process input
        num = num if isinstance(num, list) else [num]

        if max(num) > self._free:
            print("Number of molecules is larger than remaining free cells")
            return

        # Remove number of molecules from number of cells
        self._free -= max(num)

        # Add interaction matrix entry to existing molecules
        for val in self._im.values():
            val[len(self._mols)] = 0

        # Add new molecule
        self._mols[len(self._mols)] = {"num": num, "is_move": is_move, "name": name, "struct": struct}

        # Add new interaction matrix entry
        self._im[len(self._im)] = {idx: 0 for idx in self._mols.keys()}

    def _run_helper(self, systems, temp, steps_equi, steps_prod, binding, pb_f, n_print, out, traj):
        """Run Monte Carlo algorithm.

        Parameters
        ----------
        systems : list
            List of system touples containing number of molecules of each type
        temp : float
            Simulation temperature in Kelvin
        steps_equi : integer
            Number of MC steps in the equilibration phase
        steps_prod : integer
            Number of MC steps in the production phase
        binding : list, optional
            Systems to calculate the binding probability for
        pb_f : list, optional
            Length and frequency of probability calculation list
        n_print : integer, optional
            Number of steps to print output for
        out : list, optional
            Probability output file link and frequency in steps
        traj : list, optional
            Trajectory output file link and frequency in steps

        Returns
        -------
        res_dict : dictionary
            Dictionary of system and corresponding results
        """
        # Initialize
        res_dict = {}

        for system in systems:
            # Set up system
            box = Box(self._size)
            for mol_id, mol in self._mols.items():
                box.add_mol(system[mol_id], mol["is_move"], mol["name"], mol["struct"])
            box.set_im(self._im)

            # Run MC
            mc = MC(box, temp)
            res_dict[system] = mc.run(steps_equi, steps_prod, binding, pb_f, n_print, out, traj)

        return res_dict

    def run(self, temp, steps_equi, steps_prod, out_link, binding=[{"host": 0, "guest": 1}], pb_f=[1000, 100], n_print=1000, out=["", 0], traj=["", 0], np=0, is_parallel=True):
        """Run Monte Carlo algorithm.

        Parameters
        ----------
        temp : float
            Simulation temperature in Kelvin
        steps_equi : integer
            Number of MC steps in the equilibration phase
        steps_prod : integer
            Number of MC steps in the production phase
        out_link : string
            File link for output
        binding : list, optional
            Systems to calculate the binding probability for
        pb_f : list, optional
            Length and frequency of probability calculation list
        n_print : integer, optional
            Number of steps to print output for
        out : list, optional
            Probability output file link and frequency in steps
        traj : list, optional
            Trajectory output file link and frequency in steps
        np : integer, optional
            Number of cores to use
        is_parallel : bool, optional
            True to run parallelization

        Returns
        -------
        res_dict : dictionary
            Dictionary of system and corresponding results
        """
        # Create all possible configurations
        systems = list(it.product(*[x["num"] for x in self._mols.values()]))

        # Run helper function
        if is_parallel:
            # Get number of cores
            np = np if np and np<=mp.cpu_count() else mp.cpu_count()

            # Calculate number of systems per processor
            sys_num = math.floor(len(systems)/np)

            # Divide systems on processors
            if sys_num==0:
                sys_np = [[system] for system in systems]
            else:
                sys_np = [systems[sys_num*i:] if i == np-1 else systems[sys_num*i:sys_num*(i+1)] for i in range(np)]

            # Run parallel search
            pool = mp.Pool(processes=len(sys_np) if len(sys_np)<np else np)
            results = [pool.apply_async(self._run_helper, args=(x, temp, steps_equi, steps_prod, binding, pb_f, n_print, out, traj)) for x in sys_np]
            res_dict = {}
            for res in results:
                res_dict.update(res.get())

            # Destroy object
            del results
        else:
            res_dict = self._run_helper(systems, temp, steps_equi, steps_prod, binding, pb_f, n_print, out, traj)

        # Save results
        utils.save([self, res_dict], out_link)

        # Return results
        return res_dict

    def plot_pb(self, results_link, mol_x, mol_y, p_b_id):
        """Visualize adsorption results.

        Parameters
        ----------
        results : string
            File link to results object
        mol_x : integer
            Molecule id to show on x-axis
        mol_y : integer
            Molecule id to show in legend
        p_b_id : touple
            Reults index of the p_b list
        """
        # Load results
        res = utils.load(results_link)

        # Get data
        mols = res[0].get_mols()
        results = {key: val["p_b"] for key, val in res[1].items()}

        y = mols[mol_y]["num"]
        y.sort()
        for num_y in y:
            x = mols[mol_x]["num"]
            x.sort()
            p_b = []
            for num_x in x:
                for system, result in results.items():
                    if system[mol_x] == num_x and system[mol_y] == num_y:
                        p_b.append(np.mean(result[p_b_id]))

            sns.lineplot(x=x, y=p_b)

        plt.legend(y)

    def plot_ads(self, results_link, system, val_x, val_y, is_plot=True):
        """Visualize adsorption numbers for results.

        Parameters
        ----------
        results : string
            File link to results object
        system : dictionary
            host **host** and guest **guest** molecule for system search
        val_x : dictionary
            Molecule id **mol_id**, p_b value **p_b**, and bound/unbound **bu** to show on x-axis
        val_y : dictionary
            Molecule id **mol_id**, p_b value **p_b**, and bound/unbound **bu** to show on y-axis
        is_plot : bool, optional
            True to plot lines

        Returns
        -------
        plot_dict : dictionary
            Dictionary with host number as keys and x, y lists as values
        """
        # Load results
        res = utils.load(results_link)

        # Get data
        mols = res[0].get_mols()
        results = {key: val["p_b"] for key, val in res[1].items()}

        # Get molecule numbers
        mol_x = val_x["mol_id"]
        mol_y = val_x["mol_id"]
        mol_h = system["host"]
        mol_g = system["guest"]
        num_h = mols[mol_h]["num"]
        num_g = mols[mol_g]["num"]

        # Run through systems
        plot_dict = {}
        for i in num_h:
            x = []
            y = []
            for j in num_g:
                for sys, result in results.items():
                    if sys[mol_h] == i and sys[mol_g] == j:
                        num_x = i if mol_x==mol_h else j
                        num_y = i if mol_y==mol_h else j
                        x.append(np.mean(result[val_x["p_b"]])*num_x if val_x["bu"]=="b" else ((1-np.mean(result[val_x["p_b"]]))*num_x))
                        y.append(np.mean(result[val_y["p_b"]])*num_y if val_y["bu"]=="b" else ((1-np.mean(result[val_y["p_b"]]))*num_y))

            # Sort
            x.sort()
            y.sort()
            plot_dict[i] = [x, y]

        # Plot
        if is_plot:
            legend = []
            for leg, xy in plot_dict.items():
                sns.lineplot(x=xy[0], y=xy[1])
                legend.append(leg)
            plt.legend(legend)

        return plot_dict

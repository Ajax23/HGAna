################################################################################
# MC Class                                                                     #
#                                                                              #
"""Monte Carlo algorithm for calculating occupation probability."""
################################################################################


import sys
import random
import numpy as np


class MC:
    """This class run a Monte Carlo algorithm in order to calculate the
    probability of complex creating between two molecules.

    Parameters
    ----------
    box : Box
        Theoretical box system
    temp : float
        Simulation temperature in Kelvin
    """
    def __init__(self, box, temp):
        self._box = box
        self._mols = box.get_mols()
        self._im = box.get_im()

        self._temp = temp  # K
        self._beta = 1/8.314e-3/temp  # kJ/mol

        self._move_list = [key for key, mol in self._mols.items() if mol["is_move"]]
        self._lattice = {x: [] for x in range(box.get_cells())}
        self._occupied = {x: [] for x in self._mols.keys()}

        # Fill lattice
        cell_id = 0
        for mol_id, props in self._mols.items():
            for i in range(props["num"]):
                self._lattice[cell_id].append(mol_id)
                self._occupied[mol_id].append(cell_id)
                cell_id += 1

    def _move(self, mol_id, old, new):
        """Move molecule from an old cell to a new cell.

        Parameters
        ----------
        mol_id : integer
            Molecule index
        old : integer
            Old cell index
        new : integer
            New cell index
        """
        # Update lattice
        self._lattice[old].remove(mol_id)
        self._lattice[new].append(mol_id)

        # Update occupancy
        self._occupied[mol_id].remove(old)
        self._occupied[mol_id].append(new)

    def _metropolis(self, dE):
        """Performs the acceptance criterion of the Metropolis–Hastings
        algorithm.

        Parameters
        ----------
        dE : float
            The difference in energy between the new and old state

        Returns
        -------
        accepted : bool
            True if move is accepted
        """
        # Choose random number
        rand = random.uniform(0, 1)

        # Check metropolis
        if dE <= 0:
            return True
        else:
            return rand < min(1, np.exp(-self._beta*dE))

    def _run_phase(self, steps, binding, pb_f, n_print, out, traj):
        """Run Monte Carlo algorithm for a number of steps.

        Parameters
        ----------
        steps : integer
            Number of MC steps
        binding : list
            Systems to calculate the binding probability for
        pb_f : list
            Length and frequency of probability calculation list
        n_print : integer
            Number of steps to print output for
        out : list
            Probability output file link and frequency in steps
        traj : list
            Trajectory output file link and frequency in steps
        """
        # Initialize
        im = self._im
        mols = self._mols
        lattice = self._lattice
        occupied = self._occupied

        # Output format
        out_form = "%"+str(len(str(steps)))+"i"

        # Acceptance numbers
        n_acc = 0
        n_rej = 0

        # Binding probability
        p_b = {(x["host"], x["guest"]): [] for x in binding}

        # Run through MC steps
        for step_id in range(steps):
            # Choose random molecule
            mol_id = random.choice(self._move_list)

            # Choose random old and new position
            pos_old = random.choice(occupied[mol_id])
            pos_new = random.choice(list(lattice.keys()))

            # Get occupancy
            cell_old = lattice[pos_old]
            cell_new = lattice[pos_new]

            # Run through states
            ## New cell is empty
            if not cell_new:
                ### Old unbound
                if len(cell_old) == 1:
                    is_accept = True
                ### Old bound
                else:
                    #### Check unbinding
                    is_accept = self._metropolis(-im[cell_old[0]][cell_old[1]])
            ## New Cell filled
            else:
                ### Only one molecule in new cell
                if len(cell_new) == 1:
                    #### Old unbound
                    if len(cell_old) == 1:
                        ##### Metropolis check binding, and interaction matrix if binding is possible (im > 0)
                        if im[mol_id][cell_new[0]]:
                            is_accept = self._metropolis(im[mol_id][cell_new[0]])
                        else:
                            is_accept = False
                    #### Old bound
                    else:
                        ##### Metropolis check unbinding and then binding
                        if self._metropolis(-im[cell_old[0]][cell_old[1]]):
                            is_accept = self._metropolis(im[mol_id][cell_new[0]])
                        else:
                            is_accept = False
                ### Complex in new cell
                else:
                    is_accept = False

            # Proccess acceptance
            if is_accept:
                n_acc += 1
                self._move(mol_id, pos_old, pos_new)
            else:
                n_rej += 1

            # Calculate average binding probability
            if pb_f[0] and ((step_id+1)%pb_f[1]==0 or step_id==0 or step_id==steps-1):
                for host, guest in list(p_b.keys()):
                    if len(p_b[(host, guest)])==pb_f[0]:
                        p_b[(host, guest)].pop(0)
                    p_b[(host, guest)].append(0)

                    if is_accept or step_id<=2:
                        for cell in list(set(occupied[host])):
                            if host in lattice[cell] and guest in lattice[cell]:
                                p_b[(host, guest)][-1] += 1
                        p_b[(host, guest)][-1] /= mols[host]["num"]
                    else:
                        p_b[(host, guest)][-1] += p_b[(host, guest)][-2]

            # Progress
            if n_print and ((step_id+1)%n_print==0 or step_id==0 or step_id==steps-1):
                out_string = out_form%(step_id+1)+"/"+out_form%steps
                out_string += "  - acc/rej="+"%.3f"%(n_acc/n_rej if n_rej>0 else 0) # +" "+str(n_acc)+" "+str(n_rej)
                for host, guest in list(p_b.keys()):
                    out_string += ", p_b("+str(host)+","+str(guest)+")="+"%.4f"%np.mean(p_b[(host, guest)])+"+-"+"%.4f"%np.std(p_b[(host, guest)])
                sys.stdout.write(out_string+"\r")
                sys.stdout.flush()

        if n_print:
            print()

        return {"p_b": p_b, "n_acc": n_acc, "n_rej": n_rej}

    def run(self, steps_equi, steps_prod, binding=[{"host": 0, "guest": 1}], pb_f=[1000, 100], n_print=1000, out=["", 0], traj=["", 0]):
        """Run Monte Carlo algorithm.

        Parameters
        ----------
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
        """
        # Run equilibration phase
        if n_print:
            print("Running equilibration phase ...")
        self._run_phase(steps_equi, binding, [0, 0], 0, ["", 0], ["", 0])

        # Run Production phase
        if n_print:
            print("Running production phase ...")
        return self._run_phase(steps_prod, binding, pb_f, n_print, out, traj)

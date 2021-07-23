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

    def _run_phase(self, steps, binding, dt):
        """Run Monte Carlo algorithm for a number of steps.

        Parameters
        ----------
        steps : integer
            Number of MC steps
        binding : list
            Systems to calculate the binding probability for
        dt : integer
            Output frequency in steps
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

            # Calculation and output
            if dt:
                # Calculate average binding probability
                for host, guest in list(p_b.keys()):
                    p_b[(host, guest)].append(0)
                    for cell in list(set(occupied[host])):
                        if host in lattice[cell] and guest in lattice[cell]:
                            p_b[(host, guest)][-1] += 1
                    p_b[(host, guest)][-1] /= mols[host]["num"]

            # Progress
                if (step_id+1)%dt==0 or step_id==0 or step_id==steps-1:
                    sys.stdout.write(out_form%(step_id+1)+"/"+out_form%steps+" - "+"n_acc = "+"%5i"%n_acc+", n_rej = "+"%5i"%n_rej+", mean(0, 1) = "+"%.5f"%np.mean(p_b[(0, 1)])+", std(0, 1) = "+"%.5f"%np.std(p_b[(0, 1)])+"\r")
                    sys.stdout.flush()

        if dt:
            print()

        return {"n_acc": n_acc, "n_rej": n_rej}

    def run(self, steps_equi, steps_prod, binding=[{"host": 0, "guest": 1}], dt=1000):
        """Run Monte Carlo algorithm.

        Parameters
        ----------
        steps_equi : integer
            Number of MC steps in the equilibration phase
        steps_prod : integer
            Number of MC steps in the production phase
        binding : list
            Systems to calculate the binding probability for
        dt : integer, optional
            Output frequency in steps
        """
        # Run equilibration phase
        self._run_phase(steps_equi, binding, 0)

        # Run Production phase
        self._run_phase(steps_prod, binding, dt)

################################################################################
# Box Class                                                                    #
#                                                                              #
"""Theoretical simulation box for the Monte Carlo algorithm."""
################################################################################


class Box:
    """This class creates a hypothetical simulation box, to which molecules are
    added. Interaction of molecules are defined using an interaction matrix.

    Parameters
    ----------
    size : list
        Three dimensional cell count in each dimension
    """
    def __init__(self, size):
        self._size = size
        self._cells = size[0]*size[1]*size[2]
        self._mols = {}
        self._im = {}
        self._free = size[0]*size[1]*size[2]

    def add_mol(self, num, is_move=True, name="", struct=""):
        """Add a molecule to the system. The interaction matrix is set up and
        filled with empty values.

        Parameters
        ----------
        num : integer
            Number of molecules to be added
        is_move : bool, optional
            True to move molecule during simulation
        name : string, optional
            Molecule name
        struct : string, optional
            File link to structure file
        """
        # Process input
        if num > self._free:
            print("Number of molecules is larger than remaining free cells")

        # Remove number of molecules from number of cells
        self._free -= num

        # Add interaction matrix entry to existing molecules
        for val in self._im.values():
            val[len(self._mols)] = 0

        # Add new molecule
        self._mols[len(self._mols)] = {"num": num, "is_move": is_move, "name": name, "struct": struct}

        # Add new interaction matrix entry
        self._im[len(self._im)] = {idx: 0 for idx in self._mols.keys()}


    ##################
    # Setter Methods #
    ##################
    def set_interaction(self, mol_a, mol_b, val):
        """Add interaction parameter between two molecules.

        Parameters
        ----------
        mol_a : integer
            Index of first molecule
        mol_b : integer
            Index of second molecule
        val : float
            Interaction parameter
        """
        self._im[mol_a][mol_b] = val
        self._im[mol_b][mol_a] = val


    ##################
    # Getter Methods #
    ##################
    def get_size(self):
        """Return the lattice size.

        Returns
        -------
        size : list
            Number of cells in each dimension
        """
        return self._size

    def get_cells(self):
        """Return the number of cells.

        Returns
        -------
        cells : integer
            Number of cells
        """
        return self._cells

    def get_mols(self):
        """Return the molecule dictionary.

        Returns
        -------
        mols : dictionary
            Molecule dictionary containing all information
        """
        return self._mols

    def get_num(self, idx):
        """Return the molecule number.

        Parameters
        ----------
        idx : int
            Molecule index

        Returns
        -------
        num : integer
            Molecule number
        """
        return self._mols[idx]["num"]

    def get_name(self, idx):
        """Return the molecule name.

        Parameters
        ----------
        idx : int
            Molecule index

        Returns
        -------
        name : string
            Molecule name
        """
        return self._mols[idx]["name"]

    def get_im(self):
        """Return the interaction matrix.

        Returns
        -------
        im : dictionary
            Interaction matrix
        """
        return self._im

    def get_interaction(self, mol_a, mol_b):
        """Return the interaction parameter between to molecules.

        Parameters
        ----------
        mol_a : integer
            Index of first molecule
        mol_b : integer
            Index of second molecule

        Returns
        -------
        val : float
            Interaction parameter
        """
        return self._im[mol_a][mol_b]

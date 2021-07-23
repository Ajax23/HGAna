import os
import sys

import shutil
import unittest

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import hgana as hga


class UserModelCase(unittest.TestCase):
    #################
    # Remove Output #
    #################
    @classmethod
    def setUpClass(self):
        if os.path.isdir("tests"):
            os.chdir("tests")

        folder = 'output'
        hga.utils.mkdirp(folder)
        hga.utils.mkdirp(folder+"/temp")
        open(folder+"/temp.txt", 'a').close()

        for filename in os.listdir(folder):
            file_path = os.path.join(folder, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

    #########
    # Utils #
    #########
    def test_utils(self):
        file_link = "output/test/test.txt"

        hga.utils.mkdirp("output/test")

        self.assertEqual(hga.utils.column([[1, 1, 1], [2, 2, 2]]), [[1, 2], [1, 2], [1, 2]])

        hga.utils.save([1, 1, 1], file_link)
        self.assertEqual(hga.utils.load(file_link), [1, 1, 1])

        self.assertEqual(round(hga.utils.mumol_m2_to_mols(3, 100), 4), 180.66)
        self.assertEqual(round(hga.utils.mols_to_mumol_m2(180, 100), 4), 2.989)
        self.assertEqual(round(hga.utils.mmol_g_to_mumol_m2(0.072, 512), 2), 0.14)
        self.assertEqual(round(hga.utils.mmol_l_to_mols(30, 1000), 4), 18.066)
        self.assertEqual(round(hga.utils.mols_to_mmol_l(18, 1000), 4), 29.8904)

        print()
        hga.utils.toc(hga.utils.tic(), message="Test", is_print=True)
        self.assertEqual(round(hga.utils.toc(hga.utils.tic(), is_print=True)), 0)


    ################
    # Cyclodextrin #
    ################
    def test_cd(self):
        self.skipTest("Temporary")

        # Set style
        sns.set_style("white",{"xtick.bottom": True,'ytick.left': True})
        sns.set_context("paper")
        sns.set_palette(sns.color_palette("deep"))

        # Count bound and unbound instances
        print()
        hga.affinity.sample("data/COLVAR", "output/count.obj", [1, 0.5], [2, 0.46], is_force=True)
        hga.affinity.sample("data/COLVAR", "output/count.obj", [1, 0.5], [2, 0.46], is_force=False)

        # Calculate binding affinity through brute-force summation
        print()
        table = hga.affinity.number("output/count.obj", 298.15, 31.3707e-27)
        print(table)

        # Calculate binding affinity through association and dissociation rates
        print()
        tables = [hga.affinity.time("output/count.obj", cutoff, 298.15, 31.3707e-27) for cutoff in [100*x for x in range(11)]]
        table = pd.concat(tables)
        print(table)

        # Test standard deviation
        print()
        tables = [hga.affinity.time("output/count.obj", cutoff, 298.15, 31.3707e-27, is_std=True) for cutoff in [100*x for x in range(11)]]
        table = pd.concat(tables)
        print(table)

        # Plot histogram
        plt.figure(figsize=(6, 4))
        hga.affinity.hist("data/COLVAR", [1, 2], ["Centers of Mass", "Oxygenes"], conditions={2: [1, 0.5]})
        plt.savefig("output/affinity.pdf", format="pdf", dpi=1000)
        # plt.show()

        # Plot time series
        plt.figure(figsize=(6, 4))
        hga.affinity.time_series("data/COLVAR", [1, 2])
        plt.savefig("output/time_series.pdf", format="pdf", dpi=1000)
        # plt.show()

    ######
    # MC #
    ######
    def test_box(self):
        # Initialize
        box = hga.Box([10, 10, 10])
        print()

        # Add molecules
        box.add_mol(10)
        box.add_mol(10)
        box.add_mol(10)
        print(box.get_mols())

        # Set interaction matrix
        box.set_interaction(0, 1, 10)
        print(box.get_interaction(1, 0))
        print(box.get_im())

    def test_mc(self):
        # Set up box
        box = hga.Box([10, 10, 10])
        box.add_mol(10, is_move=False)
        box.add_mol(10)
        box.set_interaction(0, 1, -15)
        box.set_interaction(0, 0, 0)

        # Initialize
        mc = hga.MC(box, 298)
        print()
        print(mc._occupied)

        # Move
        mc._move(0, 0, 999)
        print(mc._occupied)

        # Run
        mc.run(10000000, 100000000, dt=1000)
        for occ in mc._occupied.values():
            occ.sort()
            print(occ)



if __name__ == '__main__':
    unittest.main(verbosity=2)

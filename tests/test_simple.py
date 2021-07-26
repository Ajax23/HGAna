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

        # Set style
        sns.set_style("white",{"xtick.bottom": True,'ytick.left': True})
        sns.set_context("paper")
        sns.set_palette(sns.color_palette("deep"))

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


    ###########
    # Extract #
    ###########
    def test_extract(self):
        # self.skipTest("Temporary")

        hga.extract.extract("data/COLVAR", "output", com=0.9)


    ############
    # Affinity #
    ############
    def test_affinity(self):
        # self.skipTest("Temporary")

        # Count bound and unbound instances
        print()
        hga.affinity.sample("data/COLVAR", "output/count.obj", [1, 0.5], [2, 0.46], is_force=True)
        hga.affinity.sample("data/COLVAR", "output/count.obj", [1, 0.9], [2, 0.20], is_force=False)

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

        # Add molecules
        box.add_mol(10)
        box.add_mol(10)
        box.add_mol(10)
        self.assertEqual(box.get_mols(), {0: {'num': 10, 'is_move': True, 'name': '', 'struct': ''}, 1: {'num': 10, 'is_move': True, 'name': '', 'struct': ''}, 2: {'num': 10, 'is_move': True, 'name': '', 'struct': ''}})

        # Set interaction matrix
        box.set_interaction(0, 1, 10)
        self.assertEqual(box.get_interaction(1, 0), 10)
        self.assertEqual(box.get_im(), {0: {0: 0, 1: 10, 2: 0}, 1: {0: 10, 1: 0, 2: 0}, 2: {0: 0, 1: 0, 2: 0}})
        box.set_im({0: {0: 0, 1: 15, 2: 0}, 1: {0: 15, 1: 0, 2: 0}, 2: {0: 0, 1: 0, 2: 0}})
        self.assertEqual(box.get_im(), {0: {0: 0, 1: 15, 2: 0}, 1: {0: 15, 1: 0, 2: 0}, 2: {0: 0, 1: 0, 2: 0}})

        # Getter functions
        self.assertEqual(box.get_size(), [10, 10, 10])
        self.assertEqual(box.get_num(0), 10)
        self.assertEqual(box.get_name(0), "")
        self.assertEqual(box.get_struct(0), "")

        # Test error
        print()
        self.assertIsNone(box.add_mol(100000))

    def test_mc(self):
        # self.skipTest("Temporary")
        print()

        # Set up box
        box = hga.Box([10, 10, 10])
        box.add_mol(10, is_move=False)
        box.add_mol(10)
        box.set_interaction(0, 1, -15)
        box.set_interaction(0, 0, 0)

        # Initialize
        mc = hga.MC(box, 298)
        self.assertEqual(mc._occupied, {0: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 1: [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]})

        # Move
        mc._move(0, 0, 999)
        self.assertEqual(mc._occupied, {0: [1, 2, 3, 4, 5, 6, 7, 8, 9, 999], 1: [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]})

        # Run
        mc.run(100000, 100000, binding=[{"host": 0, "guest": 1}], pb_f=[1000, 50], n_print=1000)
        for occ in mc._occupied.values():
            occ.sort()
            print(occ)

        # Add Inhibitor
        box.add_mol(10)
        box.set_interaction(0, 2, -10)
        mc = hga.MC(box, 298)
        mc.run(100000, 100000, binding=[{"host": 0, "guest": 1}, {"host": 0, "guest": 2}], pb_f=[1000, 50], n_print=1000)
        for occ in mc._occupied.values():
            occ.sort()
            print(occ)

    def test_ads(self):
        # self.skipTest("Temporary")
        print()

        # Test simulation
        ads = hga.Adsorption([10, 10, 10])
        ads.add_mol([1, 10], is_move=False)
        # ads.add_mol([x for x in range(1, 20+1, 5)])
        ads.add_mol([x for x in range(1, 100+1, 20)])
        ads.set_interaction(0, 1, -15)
        # ads.set_interaction(0, 2, -10)

        # Run single
        ads.run(298, 100000, 10000, "output/ads.obj", binding=[{"host": 0, "guest": 1}, {"host": 1, "guest": 0}], pb_f=[1000, 50], n_print=1000, is_parallel=True)
        ads.run(298, 100000, 10000, "output/ads.obj", binding=[{"host": 0, "guest": 1}, {"host": 1, "guest": 0}], pb_f=[1000, 50], n_print=1000, is_parallel=False)
        # print(results)

        # Plot pb
        plt.figure(figsize=(6, 4))
        ads.plot_pb("output/ads.obj", 1, 0, (0, 1))
        plt.savefig("output/ads_pb.pdf", format="pdf", dpi=1000)
        # plt.show()

        # Plot number
        plt.figure(figsize=(6, 4))
        ads.plot_ads("output/ads.obj", {"host": 0, "guest": 1}, {"mol_id": 1, "p_b": (1, 0), "bu": "u"}, {"mol_id": 1, "p_b": (1, 0), "bu": "b"})
        plt.savefig("output/ads_num.pdf", format="pdf", dpi=1000)
        # plt.show()

        # Test special cases
        ads = hga.Adsorption([10, 10, 10])
        ads.add_mol(1, is_move=False)
        ads.add_mol(1)
        ads.set_interaction(0, 1, -15)
        ads.run(298, 100000, 10000, "output/ads.obj", binding=[{"host": 0, "guest": 1}], pb_f=[1000, 50], n_print=1000, is_parallel=True)

        # Test error
        ads.add_mol(1000000000)


if __name__ == '__main__':
    unittest.main(verbosity=2)

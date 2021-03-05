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


    ################
    # Cyclodextrin #
    ################
    def test_cd(self):
        # Set style
        sns.set_style("white",{"xtick.bottom": True,'ytick.left': True})
        sns.set_context("paper")
        sns.set_palette(sns.color_palette("deep"))

        # Plot histogram
        hga.affinity.hist("data/COLVAR", [1, 2], ["Centers of Mass", "Oxygenes"], conditions={2: [1, 0.5]})
        plt.savefig("output/affinity.pdf", format="pdf", dpi=1000)
        # plt.show()

        # Count bound and unbound instances
        print()
        hga.affinity.sample("data/COLVAR", "output/count.obj", [1, 0.5], [2, 0.46], is_force=True)

        # Calculate binding affinity through brute-force summation
        print()
        table = hga.affinity.number("output/count.obj", 298.15, 31.3707e-27)
        print(table)

        # Calculate binding affinity through association and dissociation rates
        print()
        tables = [hga.affinity.time("output/count.obj", cutoff, 298.15, 31.3707e-27) for cutoff in [100*x for x in range(11)]]
        table = pd.concat(tables)
        print(table)


if __name__ == '__main__':
    unittest.main(verbosity=2)

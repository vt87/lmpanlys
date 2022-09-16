import sys
import os
from paircoeffs_martini import paircoeffs_martini

m = paircoeffs_martini(fname = "martini.inp",
                       database = "../utils/database_martini/martini_v3.0.0.itp",
                       cutoff = 1.2,
                       inner_cutoff = 0.9,
                       fftype = "lj/gromacs",
                       ofile = "paircoeffs.txt")

import pyshtools as pysh
import numpy as np

cs_coef = pysh.SHCoeffs.from_file("CS_20_0309_LS.txt", dtype=np.float64)
fig, ax = cs_coef.plot_spectrum2d(show=False)

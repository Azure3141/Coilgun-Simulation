import numpy as np
import classes

# Constants
mu0 = 4 * np.pi * pow(10, -7)
ambient_temp = 300

# Density, heat capacity, resistivity, alpha
copper = classes.Material(8960, 385, 1.68 * pow(10, -8), 0.00382 * pow(10, -8))

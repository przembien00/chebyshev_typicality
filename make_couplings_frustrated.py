import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

# Define the lattice row sizes
Lx = 4
Ly = 5
N = Lx*Ly


# Create the coupling matrix with periodic boundary conditions
couplings = np.zeros((N, N))

def site_index(x, y):
    return x + y * Lx


for x in range(Lx):
    for y in range(Ly):
        i = site_index(x, y)
        # Neighbour to the right (periodic)
        j = site_index((x + 1) % Lx, y)
        couplings[i, j] = 1
        couplings[j, i] = 1
        # Neighbour above (periodic)
        j = site_index(x, (y + 1) % Ly)
        couplings[i, j] = 1
        couplings[j, i] = 1
        # Next-nearest neighbours (diagonals, periodic)
        j = site_index((x + 1) % Lx, (y + 1) % Ly)
        couplings[i, j] = 0.5
        couplings[j, i] = 0.5
        j = site_index((x - 1) % Lx, (y + 1) % Ly)
        couplings[i, j] = 0.5
        couplings[j, i] = 0.5

print(couplings)

# Save as HDF5
with h5.File(f"Couplings/Square_Frustrated_PBC_N={N}.hdf5", "w") as f:
    all = f.create_group("all")
    all.create_dataset("J_ij", data=couplings, dtype='float64')
    all.attrs["num_Spins"] = N


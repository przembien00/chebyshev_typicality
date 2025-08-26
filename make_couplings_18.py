import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

# Define the diamond lattice row sizes
row_sizes = [2, 4, 6, 4, 2]
coords = []
y = 0
for n in row_sizes:
    x0 = (max(row_sizes) - n) // 2
    for x in range(n):
        coords.append((y, x0 + x))
    y += 1

def build_couplings_matrix(coords):
    """Build adjacency matrix for nearest neighbors on a square lattice."""
    N = len(coords)
    mat = np.zeros((N, N), dtype=int)
    for i, (yi, xi) in enumerate(coords):
        for j, (yj, xj) in enumerate(coords):
            if i != j and abs(xi - xj) + abs(yi - yj) == 1:
                mat[i, j] = 1
    return mat

couplings = build_couplings_matrix(coords)

# Add PBC couplings between opposite sides of the diamond
# Manually specify the pairs as per your example
pbc_pairs = [(0, 11), (0, 15), (2, 15), (2, 17), (6, 11), (6, 17),
             (1, 6), (1, 12), (5, 12), (5, 16), (11, 16) ]
for i, j in pbc_pairs:  
    couplings[i, j] = 1
    couplings[j, i] = 1

# Save as HDF5
with h5.File("Couplings/Square_NN_PBC_N=18.hdf5", "w") as f:
    all = f.create_group("all")
    all.create_dataset("J_ij", data=couplings, dtype='float64')
    all.attrs["num_Spins"] = 18

# Optional: visualize the lattice
plt.figure(figsize=(5, 5))
for idx, (y, x) in enumerate(coords):
    plt.plot(x, -y, 'o', color='k')
    plt.text(x, -y, str(idx), ha='center', va='center', color='r')
    # Draw lines to neighbors
    for jdx, (y2, x2) in enumerate(coords):
        if couplings[idx, jdx]:
            plt.plot([x, x2], [-y, -y2], 'b-', alpha=0.3)
plt.axis('off')
plt.title('Diamond-shaped lattice (18 sites)')
plt.show()
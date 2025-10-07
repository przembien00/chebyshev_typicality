import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

J = 1
N = 20

couplings = np.random.uniform(0.9, 1.1, (N, N))  # Random couplings
couplings = np.triu(couplings, 1) + np.triu(couplings, 1).T  # Make it symmetric
couplings -= np.diag(np.diag(couplings))  # Set diagonal to zero
print(couplings)
JQ = np.sum(couplings**2, axis=0)[0]**0.5
couplings *= J / JQ  # Normalize
print(JQ)
# Save as HDF5
with h5.File(f"Couplings/Square_Random_N={N}.hdf5", "w") as f:
    all = f.create_group("all")
    all.create_dataset("J_ij", data=couplings, dtype='float64')
    all.attrs["num_Spins"] = N

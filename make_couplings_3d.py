import numpy as np
import h5py as h5

Lx = 3
Ly = 3
Lz = 2
N = Lx * Ly * Lz


def site_index(x, y, z):
    return x + y * Lx + z * Lx * Ly


couplings = np.zeros((N, N))

for x in range(Lx):
    for y in range(Ly):
        for z in range(Lz):
            i = site_index(x, y, z)
            for dx, dy, dz in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
                nx, ny, nz = (x + dx) % Lx, (y + dy) % Ly, (z + dz) % Lz
                j = site_index(nx, ny, nz)
                couplings[i, j] = 1
                couplings[j, i] = 1

with h5.File(f"Couplings/Cuboid_NN_PBC_Lx={Lx}_Ly={Ly}_Lz={Lz}_N={N}.hdf5", "w") as f:
    grp = f.create_group("all")
    grp.create_dataset("J_ij", data=couplings, dtype="float64")
    grp.attrs["num_Spins"] = N
    grp.attrs["Lx"] = Lx
    grp.attrs["Ly"] = Ly
    grp.attrs["Lz"] = Lz

print(f"Saved Lx={Lx} Ly={Ly} Lz={Lz} N={N} cuboid coupling file.")
print(couplings)

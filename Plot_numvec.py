import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from scipy.optimize import curve_fit

def F(x, a):
    return a/x**0.5

def ImportData(physical_data, project_name = ""):
    # process the inserted data:
    root_folder = "Data"

    # determine the folder:
    foldername = root_folder + "/"
    if project_name != "":
        foldername += project_name + "/"

    # determine file and return data:
    filename = foldername + physical_data + ".hdf5"
    all = h5.File( filename, 'r' )

    # discretization
    params =  all['parameters']
    disc = np.linspace(0., 5, params.attrs['num_TimePoints'])

    return all, disc

def ImportData_ED(physical_data, project_name = ""):
    # process the inserted data:
    root_folder = "Data/ED_new"

    # determine the folder:
    foldername = root_folder + "/"
    if project_name != "":
        foldername += project_name + "/"

    # determine file and return data:
    filename = foldername + physical_data + ".hdf5"
    all = h5.File( filename, 'r' )

    # discretization
    params =  all['parameters']
    disc = np.linspace(0., 1, params.attrs['num_TimePoints'])

    return all, disc

N_cores = 8
N_spins = 20
N_array = np.array( [1, 3, 9, 12, 15, 20, 30 ] )
plt.figure(constrained_layout=True)
all_conv, times = ImportData(f"NumVec_Re/ISO__Square_NN_PBC_N={N_spins}__beta=0__rescale=0.5__numVecPerCore={N_array[-1]}")
N_array = N_array[0:-1]
G_conv = np.array( [ gab for gab in all_conv['results']['Re_correlation']][0] )
# G_conv = np.concatenate((G_conv,np.flip(G_conv)))

sqsums = np.array([])

for n in N_array:

    all, times = ImportData(f"NumVec_Re/ISO__Square_NN_PBC_N={N_spins}__beta=0__rescale=0.5__numVecPerCore={n}")

    # all_ed, times_ed =  ImportData_ED("ISO_Square_NN_PBC_N=16__ISO_Disordered_Blockwise__rescale=0.5_1_5")

    G = np.array( [ gab for gab in all['results']['Re_correlation']][0] )
    # print(G_1[-1])
    # G_ed = np.array( [ gab for gab in all_ed['results'][f'1.00']['fluctuation']][0] )
    # print(np.abs(G_1[-1]))
    # G_mirror = np.concatenate((G_1,np.flip(G_1)))
    plt.plot(times, G, label = rf'NumVec={N_cores*n}')
    sqsums = np.append(sqsums, np.sum(np.abs(G-G_conv)**2)**0.5/len(G))


# plt.plot(times_ed, G_ed, 'b--', label = rf'ED')

# plt.ylim(0, 0.5)
plt.xlabel(r'$\tau$/$\beta$')
plt.ylabel(r'difference')
plt.legend()

plt.savefig(f"Plots/Test_NumVec_beta=0_{N_spins}.pdf")
plt.clf()

par, cov = curve_fit(F, N_cores*N_array, sqsums)
print(par[0], cov)

plt.xlabel('N')
plt.ylabel('error conv')
# plt.yscale('log')
plt.loglog(N_cores*N_array, sqsums, 'o')
plt.loglog(N_cores*N_array, F(N_cores*N_array, par[0]))
plt.savefig(f"Plots/NumVec_errors_beta=0_{N_spins}.pdf")
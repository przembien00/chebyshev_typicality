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
    disc = np.linspace(0., 1, params.attrs['num_TimePoints'])

    return all, disc

def ImportData_ED(physical_data, project_name = ""):
    # process the inserted data:
    root_folder = "Data/ED"

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

N_array = np.array( [1, 3, 5, 10, 25, 30] )
plt.figure(constrained_layout=True)

all_conv, times = ImportData(f"NumVec_DOUBLE/ISO__Square_NN_PBC_N=16__beta=0.5__rescale=0.5__numVecPerCore=30")
G_conv = np.array( [ 4 * gab for gab in all_conv['results']['g_zz']] )
# print(G_conv[-1])

sqsums = np.array([])

for n in N_array:

    all_1, times = ImportData(f"NumVec_DOUBLE/ISO__Square_NN_PBC_N=16__beta=0.5__rescale=0.5__numVecPerCore={n}")

    all_ed, times_ed =  ImportData_ED("ISO_Square_NN_PBC_N=16__ISO_Disordered_Blockwise__rescale=0.5")

    G_1 = np.array( [ gab for gab in all_1['results']['g_zz']] )
    # print(G_1[-1])
    G_ed = np.array( [ gab for gab in all_ed['results'][f'0.50']['fluctuation']][0] )
    # print(np.abs(G_ed[-1]-G_1[-1]))
    # G_mirror = np.concatenate((G[:int(len(G)/2)],np.flip(G[:int(len(G)/2)])))
    plt.plot(times, np.abs(G_1), label = rf'NumVec={4*n}')
    sqsums = np.append(sqsums, np.sum(np.abs(G_1-G_conv)**2)**0.5/201)

plt.plot(times_ed, G_ed, '--', label = rf'ED')

plt.xlabel(r'$\tau$/$\beta$')
plt.ylabel(r'difference')
plt.legend()

plt.savefig("Plots/Test_NumVec_beta=0.5.pdf")
plt.clf()

par, cov = curve_fit(F, 4*N_array, sqsums)
print(par[0], cov)

plt.xlabel('N')
plt.ylabel('error conv')
# plt.yscale('log')
plt.plot(4*N_array, sqsums, 'o')
plt.plot(4*N_array, F(4*N_array, par[0]))
plt.savefig("Plots/NumVec_errors_beta=0.5.pdf")
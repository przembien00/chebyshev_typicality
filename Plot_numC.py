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
    disc = np.linspace(0., 1., 2*params.attrs['num_TimePoints'])

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

N = 12
# N_array = np.array( [ 3, 6, 9, 12, 15, 18 ] )
N_array = np.array( [50, 100, 200, 250, 300 ] )
plt.figure(constrained_layout=True)
all_conv, times =  ImportData(f"numCouplings_Random/ISO__Random__N={N}__beta=1__h_z=2__numConfigs={N_array[-1]}")
G_conv = np.array( [ gab for gab in all_conv['results']['Re_correlation']][3] )
G_conv = np.concatenate((G_conv,np.flip(G_conv)))
plt.plot(times, G_conv, label = rf'numCouplings={N_array[-1]}')
N_array = N_array[0:-1]


sqsums = np.array([])

for n in N_array:

    all, times =  ImportData(f"NumCouplings_Random/ISO__Random__N={N}__beta=1__h_z=2__numConfigs={n}")

    # all_ed, times_ed =  ImportData_ED("ISO_Square_NN_PBC_N=16__ISO_Disordered_Blockwise__rescale=0.5_1_5")

    G = np.array( [ gab for gab in all['results']['Re_correlation']][3] )
    # print(G_1[-1])
    # G_ed = np.array( [ gab for gab in all_ed['results'][f'1.00']['fluctuation']][0] )
    # print(np.abs(G_1[-1]))
    G = np.concatenate((G,np.flip(G)))
    plt.plot(times, G, label = rf'numCouplings={n}')
    sqsums = np.append(sqsums, np.sum(np.abs(G-G_conv)**2)**0.5/len(G))


# plt.plot(times_ed, G_ed, 'b--', label = rf'ED')

plt.xlim(0, 1)
plt.ylim(0.2425, 0.2504)
plt.xlabel(r'$\tau$/$\beta$')
plt.ylabel(r'$g_{zz}$($\tau$)')
plt.legend()

plt.savefig(f"Plots/Test_numCouplings_{N}.pdf")
plt.clf()

par, cov = curve_fit(F, N_array, sqsums)
print(par[0], cov)

plt.xlabel('N')
plt.ylabel('error conv')
# plt.yscale('log')
plt.loglog(N_array, sqsums, 'o')
plt.loglog(N_array, F(N_array, par[0]))
plt.savefig(f"Plots/numCouplings_{N}_errors.pdf")
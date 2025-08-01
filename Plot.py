import numpy as np
import matplotlib.pyplot as plt
import h5py as h5


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

beta_array = [1, 2, 3, 4, 5]

sqsums = np.array([])

for beta in beta_array:

    all_1, times = ImportData(f"ISO__Square_NN_PBC_N=16__beta={beta:.2g}__rescale=0.5")

    all_ed, times_ed =  ImportData_ED("ISO_Square_NN_PBC_N=16__ISO_Disordered_Blockwise__rescale=0.5_1_5")

    G_1 = np.array( [ gab for gab in all_1['results']['g_zz']] )
    G_ed = np.array( [ gab for gab in all_ed['results'][f'{beta:.2f}']['fluctuation']][0] )
    # print(np.abs(G_ed[-1]-G_1[-1]))
    G_mirror = np.concatenate((G_1,np.flip(G_1)))
    plt.plot(times_ed, G_mirror, label = rf'Chebyshev, $\beta$={beta:.2g}')
    plt.plot(times_ed, G_ed, '--', label = rf'ED, $\beta$={beta:.2g}')
    sqsums = np.append(sqsums, np.sum(np.abs(G_mirror-G_ed)**2)**0.5/len(G_mirror))

plt.xlabel(r'$\tau$/$\beta$')
plt.ylabel(r'$g_{xx}$($\tau$)')
# plt.ylim(0.24, 0.26)
plt.legend()

plt.savefig("Plots/Test_lowT.pdf")
plt.clf()
plt.plot(beta_array, sqsums, 'o')
plt.yscale('log')
plt.savefig("Plots/errors_beta_lowT.pdf")
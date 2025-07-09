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

beta_array = [0.1, 0.2, 0.3, 0.4, 0.5]

for beta in beta_array:

    all, times = ImportData(f"ISO__Square_NN_PBC_N=16__beta={beta:.2g}__rescale=0.5")
    all_ed, times_ed =  ImportData_ED("ISO_Square_NN_PBC_N=16__ISO_Disordered_Blockwise__rescale=0.5")

    G = [ gab for gab in all['results']['g_zz']]
    G_ed = [ gab for gab in all_ed['results'][f'{beta:.2f}']['fluctuation']]
    print(G_ed[0][-1])
    print(G[-1])
    G_mirror = np.concatenate((G[:int(len(G)/2)],np.flip(G[:int(len(G)/2)])))
    plt.plot(times, G_mirror, label = 'Chebyshev')
    plt.plot(times_ed, G_ed[0], '--', label = 'ED')

plt.legend()

plt.savefig("Plots/Test.pdf")
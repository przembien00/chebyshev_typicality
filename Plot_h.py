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
    disc = np.linspace(0., 1, 2*params.attrs['num_TimePoints'])

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

beta_array = [0.2, 0.5, 1, 1.5, 2, 2.5]

for beta in beta_array:

    all_20, times = ImportData(f"ISO__Square_NN_PBC_N=20__beta={beta:.2g}__rescale=0.5")
    all_24, times = ImportData(f"ISO__Square_NN_PBC_N=24__beta={beta:.2g}__rescale=0.5")

    G_20 = np.array( [ gab for gab in all_20['results']['Re_correlation']][0] )
    G_24 = np.array( [ gab for gab in all_24['results']['Re_correlation']][0] )
    G_err_20 = np.array( [ gab for gab in all_20['results']['Re_stds']][0] )
    G_err_24 = np.array( [ gab for gab in all_24['results']['Re_stds']][0] )

    G_20 = np.concatenate((G_20,np.flip(G_20)))
    G_err_20 = np.concatenate((G_err_20,np.flip(G_err_20)))
    plt.errorbar(times, G_20, yerr=G_err_20, label=rf'CET N=20, $\beta={beta:.2g}$')

    G_24 = np.concatenate((G_24,np.flip(G_24)))
    G_err_24 = np.concatenate((G_err_24,np.flip(G_err_24)))
    plt.errorbar(times, G_24, yerr=G_err_24, linestyle=':', label=rf'CET N=24, $\beta={beta:.2g}$')

plt.xlabel(r'$\tau$/$\beta$')
plt.ylabel(r'$g_{xx}$($\tau$)')
# plt.ylim(0.24, 0.26)
plt.legend(fontsize=8)

plt.savefig("Plots/Plot_20v24.pdf")

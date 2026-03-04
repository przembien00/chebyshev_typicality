import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
plt.style.use('ggplot')


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

beta_array = [0.5, 1., 1.5, 2.]

for beta in beta_array:

    all, times = ImportData(f"ISO__Square_NN_PBC_N=20__beta={beta:.2g}__h_z=0.5__rescale=-0.5")

    G = np.array( [ gab for gab in all['results']['Re_correlation']][0] )
    G_err = np.array( [ gab for gab in all['results']['Re_stds']][0] )
    G = np.concatenate((G,np.flip(G)))
    G_err = np.concatenate((G_err,np.flip(G_err)))
    plt.errorbar(times, G, yerr=G_err, label=rf'CET N=20, $\beta={beta:.2g}$')

plt.xlabel(r'$\tau$/$\beta$')
plt.ylabel(r'$g_{xx}$($\tau$)')
# plt.ylim(0.24, 0.26)
plt.legend(fontsize=8)

plt.savefig("Plots/Test_h.pdf")

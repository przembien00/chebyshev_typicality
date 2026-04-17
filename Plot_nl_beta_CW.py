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

N_array = [18, 19, 20]

beta_array = [0.2, 0.4, 0.6, 0.8, 1]
beta_array_0 = [0., 0.2, 0.4, 0.6, 0.8, 1]

for N in N_array:
    G0 = [0.]
    G0_err = [0.]
    for beta in beta_array:

        all, times = ImportData(f"ISO__CW_N={N}__site=1__beta={beta}", project_name="Pair_Correlations")

        G = np.array( [ gab for gab in all['results']['Re_correlation']][0] )
        G_err = np.array( [ gab for gab in all['results']['Re_stds']][0] )
        G0.append(G[0])
        G0_err.append(G_err[0])
    plt.errorbar(beta_array_0, G0, yerr=G0_err, label = f'N={N}')

plt.xlabel(r'$\beta J_Q$')
plt.ylabel(r'$g_{12}(0)$')
plt.legend()

plt.savefig(f"Plots/Pair_Correlations/pcorr_beta_CW_AFM.pdf")

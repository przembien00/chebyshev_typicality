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

N_array = np.array([9, 10, 11, 12, 13])
nC_array = [15000, 5000, 3000, 1000, 250]

sqsums = np.array([])

C = []
C_err = []

for N, Nc in zip(N_array, nC_array):

    all, times = ImportData(f"ISO__Random__N={N}__site=1__beta=1__numConfigs={Nc}", project_name="Pair_Correlations")

    G = np.array( [ gab for gab in all['results']['Re_correlation']][0] )
    G_err = np.array( [ gab for gab in all['results']['Re_stddev']][0] )
    G_err = np.concatenate([G_err, np.flip(G_err)])
    G = np.concatenate([G, np.flip(G)])
    # plt.errorbar(times, G, yerr=G_err, label = rf'N={N}')
    # plt.plot(times, G '--', label = rf'N={N}')
    C.append(-G[0])
    C_err.append(G_err[0])

plt.ylabel(r'$g^{xx}_{12}$(0)')
plt.xlabel(r'$N$')
plt.errorbar(N_array, C, yerr=C_err)
# plt.ylim(0.00001, 0.005)
# plt.plot(N_array, 0.005/N_array)
# plt.xlabel(r'$\tau /\beta$')
# plt.ylabel(r'$g^{xx}$($\tau$)')
# plt.ylim(-0.02, 0.26)
# plt.xlim(0,1)
# plt.legend()

plt.savefig("Plots/Plot_nlcorr.pdf")
plt.clf()
# plt.plot(beta_array, sqsums, 'o')
# plt.yscale('log')
# plt.savefig("Plots/errors_beta.pdf")
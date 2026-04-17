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

site_array = [1, 2, 6]

N = 9
# system = f"Random_N={N}"
system = "FM"

if system == "AFM":
    r = 0.5
else:
    r = -0.5

beta_array = [0.2, 0.4, 0.6, 1]
beta_array_0 = [0., 0.2, 0.4, 0.6, 1.]
N_array = [18, 19, 20]


for N in N_array:
    VV0=[0.25*(N-2)/(N-1)]
    for beta in beta_array:
        # all_NN, times = ImportData(f"ISO__Square_NN_PBC_N=20__site=1__beta={beta}__rescale={r}", project_name="Pair_Correlations")
        # all_NNN, times = ImportData(f"ISO__Square_NN_PBC_N=20__site=6__beta={beta}__rescale={r}", project_name="Pair_Correlations")
        # all_NNNN, times = ImportData(f"ISO__Square_NN_PBC_N=20__site=2__beta={beta}__rescale={r}", project_name="Pair_Correlations")
        all, times = ImportData(f"ISO__CW_N={N}__site=1__beta={beta}", project_name="Pair_Correlations")

        G = np.array( [ gab for gab in all['results']['Re_correlation']][0] )
        G = np.concatenate([G, np.flip(G)])
        VV_err = G
        # VV0.append(VV_err[0]*(N-2)+0.25)
        VV0.append(0.25*(N-2)/(N-1)  + VV_err[0]*(N-2 +1/(N-1)) )
    plt.plot(beta_array_0, VV0, label=rf'N={N}')
# plt.plot(beta_array_0, [0.25]*len(beta_array_0), label = 'auto', linestyle='--')


plt.xlabel(r'$\beta J_Q$')
plt.ylabel(r'$\overline{V^x_0(0)V^x_1(0)}$')
# plt.xlim(0, 1)
# plt.ylim(0.2, 0.28)
plt.legend()

plt.savefig(f"Plots/Pair_Correlations/VV_ij_beta_CW_AFM.pdf")
plt.clf()
# plt.plot(beta_array, sqsums, 'o')
# plt.yscale('log')
# plt.savefig("Plots/errors_beta.pdf")
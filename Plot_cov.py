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
system = "AFM"

if system == "AFM":
    r = 0.5
else:
    r = -0.5

beta_array = [0.2, 0.4, 0.6, 0.8, 1]

V0 = []
C_NN = []
C_NNN = []

for beta in beta_array:


    all_NN, times = ImportData(f"ISO__Square_NN_PBC_N=20__site=1__beta={beta}__rescale={r}", project_name="Pair_Correlations")
    all_NNN, times = ImportData(f"ISO__Square_NN_PBC_N=20__site=6__beta={beta}__rescale={r}", project_name="Pair_Correlations")
    all_NNNN, times = ImportData(f"ISO__Square_NN_PBC_N=20__site=2__beta={beta}__rescale={r}", project_name="Pair_Correlations")
    # all, times = ImportData(f"ISO__Random__N=9__site=1__beta=1__numConfigs=8000", project_name="Pair_Correlations")

    G_NN = np.array( [ gab for gab in all_NN['results']['Re_correlation']][0] )
    G_3N = np.array( [ gab for gab in all_NNN['results']['Re_correlation']][0] )
    G_4N = np.array( [ gab for gab in all_NNNN['results']['Re_correlation']][0] )
    VV_err = 2 * G_3N + G_4N

    # G = np.array( [ gab for gab in all['results']['Re_correlation']][0] )
    # G = np.concatenate([G, np.flip(G)])
    # VV_err = (N-1)*(N-2)/N*G

    VV = np.concatenate([VV_err, np.flip(VV_err)])
    # G = np.concatenate([G, np.flip(G)])


    # plt.plot(times, VV, label=rf'$\beta J_Q$={beta}')
    # plt.plot(times, G, '--', label=rf'$\beta J_Q$={beta}, auto')

    V0.append(VV[0]+0.25)
    C_NN.append(G_NN[0])
    C_NNN.append(G_3N[0])

# plt.plot(beta_array, V0, label = 'full')
# plt.plot(beta_array, [0.25]*len(beta_array), label = 'auto', linestyle='--')
plt.plot(beta_array, C_NN, label = 'NN')
plt.plot(beta_array, C_NNN, label = 'NNN')

plt.xlabel(r'$\beta J_Q$')
plt.ylabel(r'$<S_iS_j>$')
# plt.ylabel(r'$\overline{V^x(\tau)V^x(0)}$')
# plt.xlim(0, 1)
# plt.ylim(0.2, 0.28)
plt.legend()

plt.savefig(f"Plots/Plot_NN_{system}.pdf")
plt.clf()
# plt.plot(beta_array, sqsums, 'o')
# plt.yscale('log')
# plt.savefig("Plots/errors_beta.pdf")
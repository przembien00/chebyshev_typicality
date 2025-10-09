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

def ImportData_spinDMFT( spin_model, physical_data = "", project = "", selfcons = True, extension = "", postfix=None ):

    # process the inserted data:
    root_folder = "Data/spinDMFT/"
    if physical_data != "":
        physical_data = "__" + physical_data
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = root_folder
    if not selfcons:
        foldername += "noselfcons/"
    if project != "":
        foldername += project + "/"

    # determine file and return data:
    filename = foldername + "spinmodel=" + spin_model + physical_data + extension + ".hdf5"
    all = h5.File( filename, 'r' )

    # discretization
    params =  all['parameters']
    disc = np.linspace(0., 1, params.attrs['num_TimePoints'])
    
    return all, disc 

def InfLim( N_list, data_list ):
    inf_data = []
    for t in range(len(data_list[0])):
        data_t = np.array( [data[t] for data in data_list] )
        coeffs = np.polyfit( 1/np.array(N_list), data_t, 1 )
        inf_lim = coeffs[-1]
        inf_data.append( inf_lim )
    return np.array(inf_data)

beta_array = [0.2, 0.6, 0.8, 1]
N_array = [16, 18, 20, 24]
plt.style.use('ggplot')

markers = ['v', '^', 's', 'x', 'D', '+']

i=0
for beta in beta_array:
    corr_list= []
    for N in N_array:
        all, times = ImportData(f"FM/ISO__Square_NN_PBC_N={N}__beta={beta:.2g}__rescale=-0.5")
        G = np.array( [ gab for gab in all['results']['Re_correlation']][0] )
        corr_list.append(G)
    G_inf = InfLim(N_array, corr_list)

    all_sdmft, times_sdmft = ImportData_spinDMFT("ISO", physical_data=f"beta={beta:.2g}", project="", extension="")
    G_sdmft = np.array( [ gab for gab in all_sdmft['results']['Re_correlation']][0] )
    G = np.concatenate((G_inf,np.flip(G_inf)))
    
    plt.plot(times_sdmft, G_sdmft, label = rf'spinDMFT, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=20)
    plt.plot(times, G, '--', label = rf'Chebyshev N$\rightarrow\infty$, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=20)
    i+=1

plt.xlabel(r'$\tau$/$\beta$')
plt.ylabel(r'$g_{zz}$($\tau$)')
plt.xlim(0, 1)
plt.legend(fontsize=7)
plt.savefig("Plots/spinDMFT_vs_FM_inf.pdf", dpi=1000)
plt.clf()
# plt.plot(beta_array, sqsums, 'o')
# plt.yscale('log')
# plt.savefig("Plots/errors_beta_18.pdf")
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
    root_folder = "Data/"
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

def ImportData_CspinDMFT( spin_model, config, other_physical_data = None, project = None, extension = None, warn_ETV = True, postfix = None ):
    # determine the folder:
    root_folder = "Data/"
    foldername = root_folder
    if project != None:
        foldername += project + "/"

    # interpret input:
    if isinstance(spin_model,str):
        spin_model = "spinmodel=" + spin_model
    else: # spin_model is a list of two elements
        spin_model = "spinmodel=" + spin_model[0] + "__mfspinmodel=" + spin_model[1]
    if other_physical_data != None:
        other_physical_data = "__" + other_physical_data
    else:
        other_physical_data = ""
    if extension != None:
        extension = "_" + extension
    else:
        extension = ""

    # create file str:
    filename = foldername + spin_model + "__config=" + config + other_physical_data + extension

    # check whether the parameters file exists and if required add suffix ETV (eigenvalue threshold violated):
    filename += ".hdf5"
    all = h5.File( filename, 'r' )


    # discretization
    params =  all['parameters']
    disc = np.linspace(0., 1, params.attrs['num_TimePoints'])

    return all, disc 

sqsums = np.array([])

# plt.style.use('ggplot')

markers = ['v', '^', 's', 'x', 'D', 'p']

N=12
beta_array = [0.5]

fig, ax = plt.subplots(figsize=(8, 3.5))

i=0
for beta in beta_array:
    all_FM, times_FM = ImportData(f"ISO__Square_NN_PBC_N=20__beta=0.5__rescale=-0.5")
    all_AFM, times_AFM = ImportData(f"ISO__Square_NN_PBC_N=20__beta=0.5__rescale=0.5")
    all_sdmft_FM, times_sdmft_FM = ImportData_CspinDMFT("ISO", f"Square_2D_N=2_NN_J=-0.5", other_physical_data=f"beta={beta:.2g}", project="CspinDMFT")
    all_sdmft_AFM, times_sdmft_AFM = ImportData_CspinDMFT("ISO", f"Square_2D_N=2_NN_J=0.5", other_physical_data=f"beta={beta:.2g}", project="CspinDMFT")
    G_FM = np.array( [ gab for gab in all_FM['results']['Re_correlation'][0] ] )
    G_AFM = np.array( [ gab for gab in all_AFM['results']['Re_correlation'][0] ] )
    G_sdmft = [gab for gab in all_sdmft_FM['results']['Re_correlation']['1-1'][0] ]

    G_sdmft_FM = np.array( [ gab for gab in all_sdmft_FM['results']['Re_correlation']['1-1'][0] ] )
    G_sdmft_FM_err = np.array( [ gab for gab in all_sdmft_FM['runtimedata']['Re_correlation_sample_stds']['2-2'][0] ] )
    G_sdmft_AFM = np.array( [ gab for gab in all_sdmft_AFM['results']['Re_correlation']['1-1'][0] ] )
    G_sdmft_AFM_err = np.array( [ gab for gab in all_sdmft_AFM['runtimedata']['Re_correlation_sample_stds']['2-2'][0] ] )
    G_FM = np.concatenate((G_FM,np.flip(G_FM)))
    G_AFM = np.concatenate((G_AFM,np.flip(G_AFM)))

    ax.errorbar(times_sdmft_FM, G_sdmft_FM, yerr=G_sdmft_FM_err, label = rf'spinDMFT, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, color='blue')
    ax.plot(times_FM, G_FM, ls=':', label = rf'Chebyshev, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, color='darkblue')
    ax.errorbar(times_sdmft_AFM, G_sdmft_AFM, yerr=G_sdmft_AFM_err, label = rf'spinDMFT, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, color='red')
    ax.plot(times_AFM, G_AFM, ls=':', label = rf'Chebyshev, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, color='darkred')
    i+=1

marker_p = []
for m in markers:
    sup_fig, sup_ax = plt.subplots()
    picture, = sup_ax.plot([1], marker=m, color="black")
    marker_p.append(picture)
label_list = []
for beta in beta_array:
    label_list.append(rf'$\beta J_Q$={beta:.2g}')

ax.legend(marker_p, label_list, loc=4)
ax.set_xlabel(r'$\tau$/$\beta$')
ax.set_ylabel(r'$g_{xx}$($\tau$)')
ax.set_xlim(0, 1)
# ax.set_ylim(0.165, 0.252)

fig.savefig(f"Plots/CspinDMFT.pdf", dpi=1000)
# fig.clf()
# plt.plot(beta_array, sqsums, 'o')
# plt.yscale('log')
# plt.savefig("Plots/errors_beta_18.pdf")
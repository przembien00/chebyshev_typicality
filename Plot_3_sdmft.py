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

beta_array = [0.2, 0.6, 0.8, 1]

sqsums = np.array([])

# plt.style.use('ggplot')

markers = ['v', '^', 's', 'x']

b_list =  ['lightskyblue', 'cornflowerblue', 'dodgerblue', 'blue', 'mediumblue']
r_list =  ['lightcoral', 'indianred', 'red', 'firebrick', 'darkred']
g_list =  ['lightgreen', 'limegreen', 'green', 'forestgreen', 'darkgreen']
fig, ax = plt.subplots()
i=0
for beta in beta_array:

    all_FM, times = ImportData(f"FM/ISO__Square_NN_PBC_N=24__beta={beta:.2g}__rescale=-0.5")
    all_AFM, times = ImportData(f"AFM/ISO__Square_NN_PBC_N=24__beta={beta:.2g}__rescale=0.5")
    all_sdmft, times_sdmft = ImportData_spinDMFT("ISO", physical_data=f"beta={beta:.2g}", project="spinDMFT", extension="")
    G_FM = np.array( [ gab for gab in all_FM['results']['Re_correlation']][0] )
    G_AFM = np.array( [ gab for gab in all_AFM['results']['Re_correlation']][0] )
    G_sdmft = np.array( [ gab for gab in all_sdmft['results']['Re_correlation']][0] )
    G_FM = np.concatenate((G_FM,np.flip(G_FM)))
    G_AFM = np.concatenate((G_AFM,np.flip(G_AFM)))
    
    ax.plot(times_sdmft, G_sdmft, label = rf'spinDMFT, $\beta J_Q$={beta:.2g}', marker=markers[i], color='green', markevery=25)
    ax.plot(times, G_FM, label = rf'Chebyshev, $\beta J_Q$={beta:.2g}', marker=markers[i], color='dodgerblue', markevery=25)
    ax.plot(times, G_AFM, label = rf'Chebyshev AFM, $\beta J_Q$={beta:.2g}', marker=markers[i], color='crimson', markevery=25)
    i+=1

ax.set_xlabel(r'$\tau$/$\beta$')
ax.set_ylabel(r'$g_{xx}$($\tau$)')
ax.set_xlim(0, 1)
marker_p = []
for m in markers:
    sup_fig, sup_ax = plt.subplots()
    picture, = sup_ax.plot([1], marker=m, color="black")
    marker_p.append(picture)
label_list = []
for beta in beta_array:
    label_list.append(rf'$\beta J_Q$={beta:.2g}')
ax.legend(marker_p, label_list, loc=4)

fig.savefig("Plots/SDMFT_AFM_FM.pdf", dpi=1000)

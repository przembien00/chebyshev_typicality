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

sqsums = np.array([])

# plt.style.use('ggplot')

markers = ['v', '^', 's', 'x', 'D', 'p']

N=12
beta_array = [0.2, 0.5, 1, 1.5, 2, 2.5]

fig, ax = plt.subplots(figsize=(8, 3.5))

i=0
for beta in beta_array:
    if beta > 1.8:
        all, times = ImportData(f"Random_Couplings/ISO__Random__N={N}__beta={beta:.2g}__numConfigs=800")
    else:
        all, times = ImportData(f"Random_Couplings/ISO__Random__N={N}__beta={beta:.2g}__numConfigs=1000")
    all_sdmft, times_sdmft = ImportData_spinDMFT("ISO", physical_data=f"beta={beta:.2g}", project="spinDMFT", extension="")
    G = np.array( [ gab for gab in all['results']['Re_correlation']][0] )
    G_err = np.array( [ gab for gab in all['results']['Re_stddev']][0] )
    G_sdmft = np.array( [ gab for gab in all_sdmft['results']['Re_correlation']][0] )
    G = np.concatenate((G,np.flip(G)))
    G_err = np.concatenate((G_err,np.flip(G_err)))

    ax.plot(times_sdmft, G_sdmft, label = rf'spinDMFT, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, color='limegreen')
    ax.errorbar(times, G, yerr=G_err, ls=':', label = rf'Chebyshev, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, color='darkviolet', errorevery=5)
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
ax.set_ylim(0.165, 0.252)

fig.savefig(f"Plots/spinDMFT_Random_N={N}.pdf", dpi=1000)
# fig.clf()
# plt.plot(beta_array, sqsums, 'o')
# plt.yscale('log')
# plt.savefig("Plots/errors_beta_18.pdf")
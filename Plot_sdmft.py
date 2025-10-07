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

beta_array = [0.2, 0.6, 0.8, 1]

sqsums = np.array([])

plt.style.use('ggplot')



for beta in beta_array:

    # all_16, times = ImportData(f"ISO__Square_NN_PBC_N=16__beta={beta:.2g}__rescale=0.5")
    all, times = ImportData(f"FM/ISO__Square_NN_PBC_N=24__beta={beta:.2g}__rescale=-0.5")
    all_sdmft, times_sdmft = ImportData_spinDMFT("ISO", physical_data=f"beta={beta:.2g}", project="", extension="")

    # all_ed, times_ed =  ImportData_ED("ISO_Square_NN_PBC_N=16__ISO_Disordered_Blockwise__rescale=0.5_1_5")

    # G_16 = np.array( [ gab for gab in all_16['results']['g_zz']] )
    G = np.array( [ gab for gab in all['results']['Re_correlation']][0] )
    # G_ed = np.array( [ gab for gab in all_ed['results'][f'{beta:.2f}']['fluctuation']][0] )
    G_sdmft = np.array( [ gab for gab in all_sdmft['results']['Re_correlation']][0] )
    # print(np.abs(G_ed[-1]-G_1[-1]))
    # G_mirror_16 = np.concatenate((G_16,np.flip(G_16)))
    G = np.concatenate((G,np.flip(G)))
    # plt.plot(times, G_mirror_16, label = rf'Chebyshev 16, $\beta$={beta:.2g}')
    # plt.plot(times_ed, G_ed, '--', label = rf'ED, $\beta$={beta:.2g}')
    
    plt.plot(times_sdmft, G_sdmft, label = rf'spinDMFT, $\beta J_Q$={beta:.2g}')
    plt.plot(times, G, '--', label = rf'Chebyshev N=24, $\beta J_Q$={beta:.2g}')

    # sqsums = np.append(sqsums, np.sum(np.abs(G_mirror-G_ed)**2)**0.5/len(G_mirror))

plt.xlabel(r'$\tau$/$\beta$')
plt.ylabel(r'Im $g_{xy}$($\tau$)')
plt.xlim(0, 1)
plt.legend(fontsize=7)
plt.savefig("Plots/spinDMFT_vs_AFM.pdf", dpi=1000)
plt.clf()
# plt.plot(beta_array, sqsums, 'o')
# plt.yscale('log')
# plt.savefig("Plots/errors_beta_18.pdf")
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5


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

beta_array = [1, 2, 3, 4]

sqsums = np.array([])

for beta in beta_array:

    all_1, times = ImportData(f"Glass/ISO__Square_Glass_N=16__beta={beta:.2g}")
    all_sdmft, times_sdmft = ImportData_spinDMFT("ISO", physical_data=f"beta={beta:.2g}", project="", extension="")


    # all_ed, times_ed =  ImportData_ED("ISO_Square_NN_PBC_N=16__ISO_Disordered_Blockwise__rescale=0.5_02_08")

    G_1 = np.array( [ gab for gab in all_1['results']['g_zz']] )
    G_sdmft = np.array( [ gab for gab in all_sdmft['results']['correlation']] )
    # G_ed = np.array( [ gab for gab in all_ed['results'][f'{beta:.2f}']['fluctuation']][0] )

    G_mirror = np.concatenate((G_1,np.flip(G_1)))
    plt.plot(times, G_mirror, label = rf'Chebyshev, $\beta$={beta:.2g}')
    # plt.plot(times_ed, G_ed, '--', label = rf'ED, $\beta$={beta:.2g}')
    plt.plot(times_sdmft, G_sdmft[0], '--', label = rf'spinDMFT, $\beta$={beta:.2g}')
    # sqsums = np.append(sqsums, np.sum(np.abs(G_mirror-G_ed)**2)**0.5/len(G_mirror))

plt.xlabel(r'$\tau$/$\beta$')
plt.ylabel(r'$g_{xx}$($\tau$)')
# plt.ylim(0.24, 0.26)
plt.legend()

plt.savefig("Plots/Test_glass.pdf")
plt.clf()

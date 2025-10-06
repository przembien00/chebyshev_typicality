import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
plt.style.use('ggplot')


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

beta_array = [0.2]

for beta in beta_array:

    all_1, times = ImportData(f"Mag_Field/ISO__Square_NN_PBC_N=24__beta={beta:.2g}__h_z=2__rescale=-0.5")

    G_1 = np.array( [ gab for gab in all_1['results']['Im_correlation']][1] )
    # print(np.abs(G_ed[-1]-G_1[-1]))
    G_1 = np.concatenate((G_1,-np.flip(G_1)))
    plt.plot(times, G_1, label=rf'CET N=24, $\beta={beta:.2g}$')

plt.xlabel(r'$\tau$/$\beta$')
plt.ylabel(r'$g_{xx}$($\tau$)')
# plt.ylim(0.24, 0.26)
plt.legend(fontsize=8)

plt.savefig("Plots/Test_h.pdf")

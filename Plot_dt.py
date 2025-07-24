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
    disc = np.linspace(0., 1, params.attrs['num_TimePoints'])

    return all, disc

def ImportData_ED(physical_data, project_name = ""):
    # process the inserted data:
    root_folder = "Data/ED"

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

dt_array = [0.001, 0.0005]

all_conv, times = ImportData(f"dt_therm/ISO__Square_NN_PBC_N=16__beta=2__rescale=0.5__dt=0.0005")
G_conv = np.array( [ gab for gab in all_conv['results']['g_zz']] )

for dt in dt_array:

    all_1, times = ImportData(f"dt_therm/ISO__Square_NN_PBC_N=16__beta=2__rescale=0.5__dt={dt}")

    all_ed, times_ed =  ImportData_ED("ISO_Square_NN_PBC_N=16__ISO_Disordered_Blockwise__rescale=0.5X")

    G_1 = np.array( [ gab for gab in all_1['results']['g_zz']] )
    print(G_1)
    G_ed = np.array( [ gab for gab in all_ed['results'][f'2.00']['fluctuation']][0] )
    # print(np.abs(G_ed[-1]-G_1[-1]))
    # G_mirror = np.concatenate((G[:int(len(G)/2)],np.flip(G[:int(len(G)/2)])))
    plt.plot(times, G_1, label = rf'dt={dt}')
    # plt.plot(times_ed, G_ed, '--', label = rf'ED')

# plt.ylim(0,0.5)
plt.xlabel(r'$\tau$/$\beta$')
plt.ylabel(r'$g_{xx}$($\tau$)')
plt.legend()

plt.savefig("Plots/Test_dt.pdf")
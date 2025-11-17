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

beta_array = [0.6, 0.8, 1]
J_L = 0.0
h_z = 2.
foldername = f"Plots/Random_h_z={h_z:.1g}/"

markers = ['v', '^', 's', 'x', 'D', '+']
fig_xx, ax_xx = plt.subplots()
fig_xy, ax_xy = plt.subplots()
fig_zz, ax_zz = plt.subplots()

i=0
for beta in beta_array:

    # all, times = ImportData(f"Mag_Field/ISO__Square_NN_PBC_N=20__beta={beta:.2g}__h_z=2__rescale=0.5")
    all_sdmft, times_sdmft = ImportData_spinDMFT("ISO", physical_data=f"beta={beta:.2g}__h=z_h_abs={h_z:.1g}", project="spinDMFT", extension="")
    all, times =  ImportData(f"Random_Couplings/ISO__Random__beta={beta:.2g}__h_z={h_z:.1g}__numConfigs=30")
    Re_G = np.array( [ gab for gab in all['results']['Re_correlation']] )
    Im_G = np.array( [ gab for gab in all['results']['Im_correlation']] )
    Re_G_sdmft = np.array( [ gab for gab in all_sdmft['results']['Re_correlation']] )
    Im_G_sdmft = np.array( [ gab for gab in all_sdmft['results']['Im_correlation']] )
    G_xx = np.concatenate((Re_G[0],np.flip(Re_G[0])))
    G_zz = np.concatenate((Re_G[3],np.flip(Re_G[3])))
    G_xy = np.concatenate((Im_G[1],-np.flip(Im_G[1])))
    
    ax_xx.plot(times_sdmft, Re_G_sdmft[0], label = rf'spinDMFT, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25)
    ax_xx.plot(times, G_xx, ':', label = rf'Chebyshev, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25)
    ax_xy.plot(times_sdmft, Im_G_sdmft[1], label = rf'spinDMFT, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25)
    ax_xy.plot(times, G_xy, ':', label = rf'Chebyshev, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25)
    ax_zz.plot(times_sdmft, Re_G_sdmft[3], label = rf'spinDMFT, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25)
    ax_zz.plot(times, G_zz, ':', label = rf'Chebyshev, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25)
    i+=1

for ax in [ax_xx, ax_xy, ax_zz]:
    ax.set_xlabel(r'$\tau$/$\beta$')
    ax.set_xlim(0, 1)
    ax.legend(fontsize=7)

ax_xx.set_ylabel(r'$g_{xx}$($\tau$)')
ax_xy.set_ylabel(r'Im $g_{xy}$($\tau$)')
ax_zz.set_ylabel(r'$g_{zz}$($\tau$)')

fig_xx.savefig(foldername+"g_xx.pdf", dpi=1000)
fig_xy.savefig(foldername+"g_xy.pdf", dpi=1000)
fig_zz.savefig(foldername+"g_zz.pdf", dpi=1000)
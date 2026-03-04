import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import scipy.optimize



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

def f(N, a, b):
    return a / N + b

def f(N, a, b):
    return a / N + b

def extrapolate( N_list, data_list, error_list ):
    inf_data = []
    inf_error = []
    for t in range(len(data_list[0])):
        data_t = np.array( [data[t] for data in data_list] )
        err_t = np.array( [error[t] for error in error_list] )
        par, cov = scipy.optimize.curve_fit(f, N_list, data_t, sigma=err_t, absolute_sigma=True )
        inf_lim = par[-1]
        err = np.sqrt( cov[-1,-1] )
        inf_data.append( inf_lim )
        inf_error.append( err )
    return np.array(inf_data), np.array(inf_error)

beta_array = [0.2, 0.5, 1, 2, 2.5]
N_array = [11, 12, 13]
h_z = 0.5
foldername = f"Plots/Random_h_z={h_z:.1g}_N=inf/"

markers = ['v', '^', 's', 'x', 'D', '+']
fig_xx, ax_xx = plt.subplots()
fig_xy, ax_xy = plt.subplots()
fig_zz, ax_zz = plt.subplots()



i=0
for beta in beta_array:
    all_11, times = ImportData(f"Random_Couplings/ISO__Random__N=11__beta={beta:.2g}__h_z={h_z:.1g}__numConfigs=1000")
    all_12, times = ImportData(f"Random_Couplings/ISO__Random__N=12__beta={beta:.2g}__h_z={h_z:.1g}__numConfigs=250")
    all_13, times = ImportData(f"Random_Couplings/ISO__Random__N=13__beta={beta:.2g}__h_z={h_z:.1g}__numConfigs=50")
    
    G_xx_list = []
    G_xx_err = []
    G_xy_list = []
    G_xy_err = []
    G_zz_list = []
    G_zz_err = []

    for all in [all_11, all_12, all_13]:
        Re_G = np.array( [ gab for gab in all['results']['Re_correlation']] )
        Im_G = np.array( [ gab for gab in all['results']['Im_correlation']] )
        Re_G_err = np.array( [ gab for gab in all['results']['Re_stddev']] )
        Im_G_err = np.array( [ gab for gab in all['results']['Im_stddev']] )
        G_xx_list.append( Re_G[0] )
        G_xy_list.append( Im_G[1] )
        G_zz_list.append( Re_G[3] )
        G_xx_err.append( Re_G_err[0] )
        G_xy_err.append( Im_G_err[1] )
        G_zz_err.append( Re_G_err[3] )

    G_xx_inf, G_xx_inf_err = extrapolate( N_array, G_xx_list, G_xx_err )
    G_xy_inf, G_xy_inf_err = extrapolate( N_array, G_xy_list, G_xy_err )
    G_zz_inf, G_zz_inf_err = extrapolate( N_array, G_zz_list, G_zz_err )
    G_xx_inf[0] = 0.25
    G_zz_inf[0] = 0.25
    G_xx = np.concatenate((G_xx_inf,np.flip(G_xx_inf)))
    G_zz = np.concatenate((G_zz_inf,np.flip(G_zz_inf)))
    G_xy = np.concatenate((G_xy_inf,np.flip(-G_xy_inf)))
    G_xx_err = np.concatenate((G_xx_inf_err,np.flip(G_xx_inf_err)))
    G_zz_err = np.concatenate((G_zz_inf_err,np.flip(G_zz_inf_err)))
    G_xy_err = np.concatenate((G_xy_inf_err,np.flip(G_xy_inf_err)))

    all_sdmft, times_sdmft = ImportData_spinDMFT("ISO", physical_data=f"beta={beta:.2g}__h=z_h_abs={h_z:.1g}", project="spinDMFT", extension="")

    Re_G_sdmft = np.array( [ gab for gab in all_sdmft['results']['Re_correlation']] )
    Im_G_sdmft = np.array( [ gab for gab in all_sdmft['results']['Im_correlation']] )

    ax_xx.plot(times_sdmft, Re_G_sdmft[0], label = rf'spinDMFT, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, color='forestgreen')
    ax_xx.errorbar(times, G_xx, yerr=G_xx_err, ls=':', label = rf'Chebyshev, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, errorevery=5, color='darkviolet')
    ax_xy.plot(times_sdmft, Im_G_sdmft[1], label = rf'spinDMFT, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, color='forestgreen')
    ax_xy.errorbar(times, G_xy, yerr=G_xy_err, ls=':', label = rf'Chebyshev, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, errorevery=5, color='darkviolet')
    ax_zz.plot(times_sdmft, Re_G_sdmft[3], label = rf'spinDMFT, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, color='forestgreen')
    ax_zz.errorbar(times, G_zz, yerr=G_zz_err, ls=':', label = rf'Chebyshev, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, errorevery=5, color='darkviolet')
    i+=1

marker_p = []
for m in markers:
    sup_fig, sup_ax = plt.subplots()
    picture, = sup_ax.plot([1], marker=m, color="black")
    marker_p.append(picture)
label_list = []
for beta in beta_array:
    label_list.append(rf'$\beta J_Q$={beta:.2g}')

for ax in [ax_xx, ax_xy, ax_zz]:
    ax.set_xlabel(r'$\tau$/$\beta$')
    ax.set_xlim(0, 1)
    ax.legend(marker_p, label_list, loc=4)
    ax.legend(fontsize=7)

ax_xx.set_ylabel(r'$g_{xx}$($\tau$)')
ax_xy.set_ylabel(r'Im $g_{xy}$($\tau$)')
ax_zz.set_ylabel(r'$g_{zz}$($\tau$)')

fig_xx.savefig(foldername+"g_xx.pdf", dpi=1000)
fig_xy.savefig(foldername+"g_xy.pdf", dpi=1000)
fig_zz.savefig(foldername+"g_zz.pdf", dpi=1000)
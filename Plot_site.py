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
    if params.attrs["evol_type"]=="imaginary":
        disc = np.linspace(0., 1., 2*params.attrs['num_TimePoints'])
    else:
        disc = np.linspace(0., params.attrs['Tmax'], params.attrs['num_TimePoints'])

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

site_array = [0, 1, 2, 6]

sqsums = np.array([])

for site in site_array:

    all_1, times = ImportData(f"ISO__Square_NN_PBC_N=20__site={site}__beta=0__rescale=0.5", project_name="Real_Time")

    G = np.array( [ gab for gab in all_1['results']['Re_correlation']][0] )
    plt.plot(times, G, label = rf'Chebyshev, site={site}')
    # plt.plot(times_ed, G_ed, '--', label = rf'ED, $\beta$={beta:.2g}')
    # sqsums = np.append(sqsums, np.sum(np.abs(G_mirror-G_ed)**2)**0.5/len(G_mirror))

plt.xlabel(r'$tJ_Q$')
plt.ylabel(r'$g_{xx}$($t$)')
# plt.ylim(0.24, 0.26)
plt.legend()

plt.savefig("Plots/Test.pdf")
plt.clf()
# plt.plot(beta_array, sqsums, 'o')
# plt.yscale('log')
# plt.savefig("Plots/errors_beta.pdf")
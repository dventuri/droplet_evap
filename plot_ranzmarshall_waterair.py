import numpy as np
import matplotlib.pyplot as plt

# Diameter plot
exp_data = np.loadtxt("ranzmarshall_dp.csv", delimiter=",")
dvt_data = np.loadtxt("water_air-dp.dat", delimiter=',')
mfs_data = np.loadtxt("MFSim-dp.dat", delimiter=',')
abg_data = np.loadtxt("ranzmarshall_abgail_AS.dat", delimiter=",")
jam_data = np.loadtxt("dados_raio_mixAS_d2_com_epsilon.dat")

fig, ax = plt.subplots()
ax.set_xlabel('Tempo [s]')
ax.set_ylabel(r'$(D_d/D_{d,0})^2$')
ax.axis([-0.0001, 900, 0, 1.2])
ax.xaxis.set_major_locator(plt.MultipleLocator(100))
ax.xaxis.set_minor_locator(plt.MultipleLocator(50))
ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
ax.plot(dvt_data[:,0], (dvt_data[:,1]/dvt_data[0,1])**2,
    ls='-', color='k', #marker='o',
    label='Diego'
)
ax.plot(mfs_data[:,0], (mfs_data[:,1]/mfs_data[0,1])**2,
    ls='-', color='b', #marker='o',
    label='MFSim'
)
ax.plot(jam_data[:,0], jam_data[:,1]/jam_data[0,1],
    ls='-', color='r', #marker='s',
    label='Jamille'
)
ax.plot(abg_data[:,0], abg_data[:,1],
    ls='--', color='gray',
    label='Abgail'
)
ax.scatter(exp_data[:,0], exp_data[:,1],
    marker='*',
    label='Ranz-Marshall'
)
ax.grid(color='lightgrey',ls='-.')
ax.legend(facecolor="white", framealpha=1, frameon=1)
fig.tight_layout(pad=0.15)
fig.savefig('water_air_droplet_diameter.png',dpi=1200,format='png')


# Temperature plot
dvt_data = np.loadtxt("water_air-Temp.dat", delimiter=',')
mfs_data = np.loadtxt("MFSim-Temp.dat", delimiter=',')
abg_data = np.loadtxt("temperature_abgail.dat", delimiter=',')
jam_data = np.loadtxt("dados_temp_med_mixAS_d2_com_epsilon.dat")

fig, ax = plt.subplots()
ax.set_xlabel('Tempo [s]')
ax.set_ylabel('Temperatura [K]')
ax.axis([-0.0001, 900, 281, 285])
ax.xaxis.set_major_locator(plt.MultipleLocator(100))
ax.xaxis.set_minor_locator(plt.MultipleLocator(50))
ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
ax.plot(dvt_data[:,0], dvt_data[:,1],
    ls='-', color='k', #marker='o',
    label='Diego'
)
ax.plot(mfs_data[:,0], mfs_data[:,1],
    ls='-', color='b', #marker='o',
    label='MFSim'
)
ax.plot(jam_data[:,0], jam_data[:,1],
    ls='-', color='r', #marker='s',
    label='Jamille'
)
ax.plot(abg_data[:,0], abg_data[:,1],
    ls='--', color='gray',
    label='Abgail'
)
ax.grid(color='lightgrey',ls='-.')
ax.legend(facecolor="white", framealpha=1, frameon=1)
fig.tight_layout(pad=0.15)
fig.savefig('water_air_droplet_temperature.png',dpi=1200,format='png')

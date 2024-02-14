import numpy as np
import matplotlib.pyplot as plt

# Diameter plot
dvt_data_10 = np.loadtxt("water_hcs-dp-10micra.dat", delimiter=',')
dvt_data_100 = np.loadtxt("water_hcs-dp-100micra.dat", delimiter=',')
dvt_data_1000 = np.loadtxt("water_hcs-dp-1000micra.dat", delimiter=',')

fig, ax = plt.subplots()
ax.set_xlabel('Distância [m]')
ax.set_ylabel(r'$D_d/D_{d,0}$')
ax.axis([-0.5, 30, 0.84, 1.05])
ax.xaxis.set_major_locator(plt.MultipleLocator(5))
ax.xaxis.set_minor_locator(plt.MultipleLocator(2.5))
ax.yaxis.set_major_locator(plt.MultipleLocator(0.05))
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.025))
ax.plot(dvt_data_10[:,0]*20, (dvt_data_10[:,1]/dvt_data_10[0,1])**2,
    ls='--', color='m',
    label=r'$D_{d,0} = 1\times10^{-5}$ m'
)
ax.plot(dvt_data_100[:,0]*20, (dvt_data_100[:,1]/dvt_data_100[0,1])**2,
    ls='--', color='b', marker='x', markevery=100000,
    label=r'$D_{d,0} = 1\times10^{-4}$ m'
)
ax.plot(dvt_data_1000[:,0]*20, (dvt_data_1000[:,1]/dvt_data_1000[0,1])**2,
    ls=':', color='r', marker='s', mfc='w', markevery=100000,
    label=r'$D_{d,0} = 1\times10^{-3}$ m'
)
# ax.vlines(0.332, 0.84, 1.035,
#           colors='gray', linestyles='-.')
# props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# ax.text(0.05, 0.15, '~ 33 cm', transform=ax.transAxes, fontsize=10,
#         verticalalignment='top', bbox=props)
ax.grid(color='lightgrey',ls='-.')
ax.legend(facecolor="white", framealpha=1, frameon=1)
fig.tight_layout(pad=0.15)
fig.savefig('water_hcs_droplet_diameter.png',dpi=1200,format='png')


# Temperature plot
dvt_data_10 = np.loadtxt("water_hcs-Temp-10micra.dat", delimiter=',')
dvt_data_100 = np.loadtxt("water_hcs-Temp-100micra.dat", delimiter=',')
dvt_data_1000 = np.loadtxt("water_hcs-Temp-1000micra.dat", delimiter=',')

fig, ax = plt.subplots()
ax.set_xlabel('Distância [m]')
ax.set_ylabel('Temperatura [K]')
ax.axis([-0.5, 30, 333.15, 367.5])
ax.xaxis.set_major_locator(plt.MultipleLocator(5))
ax.xaxis.set_minor_locator(plt.MultipleLocator(2.5))
ax.yaxis.set_major_locator(plt.MultipleLocator(5))
ax.yaxis.set_minor_locator(plt.MultipleLocator(2.5))
ax.plot(dvt_data_10[:,0]*20, dvt_data_10[:,1],
    ls='--', color='m',
    label=r'$D_{d,0} = 1\times10^{-5}$ m'
)
ax.plot(dvt_data_100[:,0]*20, dvt_data_100[:,1],
    ls='--', color='b', marker='x', markevery=100000,
    label=r'$D_{d,0} = 1\times10^{-4}$ m'
)
ax.plot(dvt_data_1000[:,0]*20, dvt_data_1000[:,1],
    ls=':', color='r', marker='s', mfc='w', markevery=100000,
    label=r'$D_{d,0} = 1\times10^{-3}$ m'
)
# ax.vlines(0.332, 333.15, 364.9,
#           colors='gray', linestyles='-.')
# props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# ax.text(0.05, 0.15, '~ 33 cm', transform=ax.transAxes, fontsize=10,
#         verticalalignment='top', bbox=props)
ax.grid(color='lightgrey',ls='-.')
ax.legend(facecolor="white", framealpha=1, frameon=1)
fig.tight_layout(pad=0.15)
fig.savefig('water_hcs_droplet_temperature.png',dpi=1200,format='png')

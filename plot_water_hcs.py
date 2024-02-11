import numpy as np
import matplotlib.pyplot as plt

# Diameter plot
dvt_data = np.loadtxt("water_hcs-dp.dat", delimiter=',')

fig, ax = plt.subplots()
ax.set_xlabel('Distância [m]')
ax.set_ylabel(r'$D_d/D_{d,0}$')
ax.axis([-1, 20, 0.9, 1.05])
ax.xaxis.set_major_locator(plt.MultipleLocator(5))
ax.xaxis.set_minor_locator(plt.MultipleLocator(2.5))
ax.yaxis.set_major_locator(plt.MultipleLocator(0.02))
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.01))
ax.plot(dvt_data[:,0]*20, (dvt_data[:,1]/dvt_data[0,1])**2,
    ls='-', color='k', #marker='o',
)
ax.grid(color='lightgrey',ls='-.')
ax.legend(facecolor="white", framealpha=1, frameon=1)
fig.tight_layout(pad=0.15)
fig.savefig('water_hcs_droplet_diameter.png',dpi=1200,format='png')


# Temperature plot
dvt_data = np.loadtxt("water_hcs-Temp.dat", delimiter=',')

fig, ax = plt.subplots()
ax.set_xlabel('Distância [m]')
ax.set_ylabel('Temperatura [K]')
ax.axis([-1, 20, 333.15, 367.5])
ax.xaxis.set_major_locator(plt.MultipleLocator(5))
ax.xaxis.set_minor_locator(plt.MultipleLocator(2.5))
ax.yaxis.set_major_locator(plt.MultipleLocator(5))
ax.yaxis.set_minor_locator(plt.MultipleLocator(2.5))
ax.plot(dvt_data[:,0]*20, dvt_data[:,1],
    ls='-', color='k', #marker='o',
)
ax.grid(color='lightgrey',ls='-.')
ax.legend(facecolor="white", framealpha=1, frameon=1)
fig.tight_layout(pad=0.15)
fig.savefig('water_hcs_droplet_temperature.png',dpi=1200,format='png')

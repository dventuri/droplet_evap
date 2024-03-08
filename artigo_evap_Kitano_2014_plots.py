import numpy as np
import matplotlib.pyplot as plt

# Diameter plot
d0 = 700E-6
kit_data_00 = np.loadtxt("Kitano-00-dp.dat", delimiter=',')
kit_data_01 = np.loadtxt("Kitano-01-dp.dat", delimiter=',')
kit_data_02 = np.loadtxt("Kitano-02-dp.dat", delimiter=',')
kit_data_03 = np.loadtxt("Kitano-03-dp.dat", delimiter=',')

fig, ax = plt.subplots()
ax.set_xlabel(r'$t/d_0^2$ [s/mm²]')
ax.set_ylabel(r'$(d/d_0)^2$ [-]')
ax.axis([0, 16, 0.0, 1.2])
ax.xaxis.set_major_locator(plt.MultipleLocator(2))
ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
ax.plot(kit_data_00[:,0]/(d0*1000)**2, (kit_data_00[:,1]**2)/d0**2,
    ls='-', color='k',
    label='0.1MPa,471K'
)
ax.plot(kit_data_01[:,0]/(d0*1000)**2, (kit_data_01[:,1]**2)/d0**2,
    ls='--', color='orange',
    label='0.5MPa,468K'
)
ax.plot(kit_data_02[:,0]/(d0*1000)**2, (kit_data_02[:,1]**2)/d0**2,
    ls=':', color='darkblue',
    label='1.0MPa,466K'
)
ax.plot(kit_data_03[:,0]/(d0*1000)**2, (kit_data_03[:,1]**2)/d0**2,
    ls='-.', color='m',
    label='2.0MPa,452K'
)
ax.grid(color='lightgrey',ls='-.')
ax.legend(facecolor="white", framealpha=1, frameon=1)
fig.tight_layout(pad=0.15)
fig.savefig('Kitano-low-dp.png',dpi=1200,format='png')


# Temperature plot
kit_data_00 = np.loadtxt("Kitano-00-Temp.dat", delimiter=',')
kit_data_01 = np.loadtxt("Kitano-01-Temp.dat", delimiter=',')
kit_data_02 = np.loadtxt("Kitano-02-Temp.dat", delimiter=',')
kit_data_03 = np.loadtxt("Kitano-03-Temp.dat", delimiter=',')

fig, ax = plt.subplots()
ax.set_xlabel(r'$t/d_0^2$ [s/mm²]')
ax.set_ylabel(r'$T$ [K]')
ax.axis([0, 16, 300, 440])
ax.xaxis.set_major_locator(plt.MultipleLocator(2))
ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
ax.yaxis.set_major_locator(plt.MultipleLocator(20))
ax.yaxis.set_minor_locator(plt.MultipleLocator(10))
ax.plot(kit_data_00[:,0]/(d0*1000)**2, kit_data_00[:,1],
    ls='-', color='k',
    label='0.1MPa,471K'
)
ax.plot(kit_data_01[:,0]/(d0*1000)**2, kit_data_01[:,1],
    ls='--', color='orange',
    label='0.5MPa,468K'
)
ax.plot(kit_data_02[:,0]/(d0*1000)**2, kit_data_02[:,1],
    ls=':', color='darkblue',
    label='1.0MPa,466K'
)
ax.plot(kit_data_03[:,0]/(d0*1000)**2, kit_data_03[:,1],
    ls='-.', color='m',
    label='2.0MPa,452K'
)
ax.grid(color='lightgrey',ls='-.')
ax.legend(facecolor="white", framealpha=1, frameon=1)
fig.tight_layout(pad=0.15)
fig.savefig('Kitano-low-Temp.png',dpi=1200,format='png')

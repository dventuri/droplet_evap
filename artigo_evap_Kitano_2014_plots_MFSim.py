import numpy as np
import matplotlib.pyplot as plt

# Diameter plot
d0 = 700E-6
kit_data_00 = np.loadtxt("Kitano-00-dp.dat", delimiter=',')
kit_mfsim_itc = np.loadtxt("Kitano-00-dp-MFSim-itc.txt", delimiter=' ')
kit_mfsim_ftc = np.loadtxt("Kitano-00-dp-MFSim-ftc.txt", delimiter=' ')
kit_mfsim_etc = np.loadtxt("Kitano-00-dp-MFSim-etc.txt", delimiter=' ')
kit_mfsim_iago = np.loadtxt("Kitano-00-dp-MFSim-iago.txt", delimiter=' ')

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
    label='Python,ITC'
)
ax.scatter(kit_mfsim_itc[:,0]/(d0*1000)**2, (kit_mfsim_itc[:,1]**2)/d0**2,
    c='white', marker='o', edgecolors='black',
    label='MFSim,ITC'
)
ax.scatter(kit_mfsim_ftc[:,0]/(d0*1000)**2, (kit_mfsim_ftc[:,1]**2)/d0**2,
    c='white', marker='s', edgecolors='black',
    label='MFSim,FTC'
)
ax.scatter(kit_mfsim_iago[:,0]/(d0*1000)**2, (kit_mfsim_iago[:,1]**2)/d0**2,
    c='white', marker='*', edgecolors='black',
    label='MFSim,FTC,100'
)
ax.scatter(kit_mfsim_etc[:,0]/(d0*1000)**2, (kit_mfsim_etc[:,1]**2)/d0**2,
    c='white', marker='^', edgecolors='black',
    label='MFSim,ETC'
)
ax.set_title('0.1MPa,471K')
ax.grid(color='lightgrey',ls='-.')
ax.legend(facecolor="white", framealpha=1, frameon=1)
fig.tight_layout(pad=0.15)
# fig.savefig('Kitano-low-dp.png',dpi=1200,format='png')


# Temperature plot
kit_data_00 = np.loadtxt("Kitano-00-Temp.dat", delimiter=',')
kit_mfsim_itc = np.loadtxt("Kitano-00-Temp-MFSim-itc.txt", delimiter=' ')
kit_mfsim_ftc = np.loadtxt("Kitano-00-Temp-MFSim-ftc.txt", delimiter=' ')
kit_mfsim_etc = np.loadtxt("Kitano-00-Temp-MFSim-etc.txt", delimiter=' ')
kit_mfsim_iago = np.loadtxt("Kitano-00-Temp-MFSim-iago.txt", delimiter=' ')

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
    label='Python,ITC'
)
ax.scatter(kit_mfsim_itc[:,0]/(d0*1000)**2, kit_mfsim_itc[:,1],
    c='white', marker='o', edgecolors='black',
    label='MFSim,ITC'
)
ax.scatter(kit_mfsim_ftc[:,0]/(d0*1000)**2, kit_mfsim_ftc[:,1],
    c='white', marker='s', edgecolors='black',
    label='MFSim,FTC'
)
ax.scatter(kit_mfsim_iago[:,0]/(d0*1000)**2, kit_mfsim_iago[:,1],
    c='white', marker='*', edgecolors='black',
    label='MFSim,FTC,100'
)
ax.scatter(kit_mfsim_etc[:,0]/(d0*1000)**2, kit_mfsim_etc[:,1],
    c='white', marker='^', edgecolors='black',
    label='MFSim,ETC'
)
ax.grid(color='lightgrey',ls='-.')
ax.legend(facecolor="white", framealpha=1, frameon=1)
fig.tight_layout(pad=0.15)
# fig.savefig('Kitano-low-Temp.png',dpi=1200,format='png')

import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI


def Lv_Perry(T,Tc,M,A,B,C,D):
    Tr = T/Tc
    return A*(1 - Tr)**(B - C*Tr + D*Tr**2)/M

def pv_Perry():
    pass


# Funcoes MFSim

def Lv(T,Tc,M,A,B):
    Tr = T/Tc
    return A*(1.0 - Tr)**B/M*1.0E6

def rho(T, Tc, A, B, C):
    D = -(1 - T/Tc)**C
    return 1000*A*B**D

def cp(T, M, A, B, C, D):
    return 1000*(A + B*T + C*T**2 + D*T**3)/M

def k(T, A, B, C):
    return A + B*T + C*T**2

def mu(T, A, B, C, D):
    return 1E-3*10.0**(A + B/T + C*T + D*T**2);


#Funcoes Coolprop

def Lv_CP(T, species):
    Hv = PropsSI('H','T',T,'Q',1,species)
    Hl = PropsSI('H','T',T,'Q',0,species)
    return Hv - Hl

def rho_CP(T, P, species):
    return PropsSI('D','T',T,'P',P,species)

def cp_CP(T, P, species):
    return PropsSI('C','T',T,'P',P,species)

def k_CP(T, P, species):
    return PropsSI('L','T',T,'P',P,species)

def mu_CP(T, P, species):
    return PropsSI('V','T',T,'P',P,species)



if __name__ == "__main__":

    T_droplet = np.linspace(298.15, 398.15, 100)
    P = 0.1E6

    # Água
    fig, ax = plt.subplots()
    ax.set_xlabel(r'$T$ [$\degree$C]')
    ax.set_ylabel(r'$L_v$ [J/kg]')
    ax.axis([25, 125, 2.1E6, 2.5E6])
    ax.xaxis.set_major_locator(plt.MultipleLocator(25))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(12.5))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.1E6))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.05E6))
    ax.plot(T_droplet-273.15, Lv(T_droplet,647.13,18.015,54,0.34),
        ls='-', color='k',
        label='MFSim'
    )
    ax.plot(T_droplet-273.15, Lv_CP(T_droplet, "Water"),
        ls='--', color='gray',
        label='Coolprop'
    )
    ax.set_title('Água')
    ax.grid(color='lightgrey',ls='-.')
    ax.legend(facecolor="white", framealpha=1, frameon=1)
    fig.tight_layout(pad=0.15)

    fig, ax = plt.subplots()
    ax.set_xlabel(r'$T$ [$\degree$C]')
    ax.set_ylabel(r'$\rho_l$ [kg/m³]')
    ax.axis([25, 125, 925, 1025])
    ax.xaxis.set_major_locator(plt.MultipleLocator(25))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(12.5))
    ax.yaxis.set_major_locator(plt.MultipleLocator(25))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(12.5))
    ax.plot(T_droplet-273.15, rho(T_droplet,647.13,0.325,0.27,0.23),
        ls='-', color='k',
        label='MFSim'
    )
    ax.plot(T_droplet-273.15, rho_CP(T_droplet, P, "Water"),
        ls='--', color='gray',
        label='Coolprop'
    )
    ax.set_title('Água')
    ax.grid(color='lightgrey',ls='-.')
    ax.legend(facecolor="white", framealpha=1, frameon=1)
    fig.tight_layout(pad=0.15)

    fig, ax = plt.subplots()
    ax.set_xlabel(r'$T$ [$\degree$C]')
    ax.set_ylabel(r'$C_{p,l}$ [J/kg.K]')
    ax.axis([25, 125, 4150, 4250])
    ax.xaxis.set_major_locator(plt.MultipleLocator(25))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(12.5))
    ax.yaxis.set_major_locator(plt.MultipleLocator(25))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(12.5))
    ax.plot(T_droplet-273.15, cp(T_droplet,18.015,-22.417,0.87697,-0.0025704,2.4838E-06),
        ls='-', color='k',
        label='MFSim'
    )
    ax.plot(T_droplet-273.15, cp_CP(T_droplet, P, "Water"),
        ls='--', color='gray',
        label='Coolprop'
    )
    ax.set_title('Água')
    ax.grid(color='lightgrey',ls='-.')
    ax.legend(facecolor="white", framealpha=1, frameon=1)
    fig.tight_layout(pad=0.15)

    fig, ax = plt.subplots()
    ax.set_xlabel(r'$T$ [$\degree$C]')
    ax.set_ylabel(r'$k_l$ [W/m.K]')
    ax.axis([25, 125, 0.6, 0.7])
    ax.xaxis.set_major_locator(plt.MultipleLocator(25))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(12.5))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.025))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.0125))
    ax.plot(T_droplet-273.15, k(T_droplet,-0.35667,0.005057,-6.1071E-06),
        ls='-', color='k',
        label='MFSim'
    )
    ax.plot(T_droplet-273.15, k_CP(T_droplet, P, "Water"),
        ls='--', color='gray',
        label='Coolprop'
    )
    ax.set_title('Água')
    ax.grid(color='lightgrey',ls='-.')
    ax.legend(facecolor="white", framealpha=1, frameon=1)
    fig.tight_layout(pad=0.15)

    fig, ax = plt.subplots()
    ax.set_xlabel(r'$T$ [$\degree$C]')
    ax.set_ylabel(r'$\mu_l$ [Pa.s]')
    ax.axis([25, 125, 2E-4, 9E-4])
    ax.xaxis.set_major_locator(plt.MultipleLocator(25))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(12.5))
    ax.yaxis.set_major_locator(plt.MultipleLocator(1E-4))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5E-4))
    ax.plot(T_droplet-273.15, mu(T_droplet,-11.6225,1.9490E03,2.1641E-02,-1.5990E-05),
        ls='-', color='k',
        label='MFSim'
    )
    ax.plot(T_droplet-273.15, mu_CP(T_droplet, P, "Water"),
        ls='--', color='gray',
        label='Coolprop'
    )
    ax.set_title('Água')
    ax.grid(color='lightgrey',ls='-.')
    ax.legend(facecolor="white", framealpha=1, frameon=1)
    fig.tight_layout(pad=0.15)


    # n-Heptano
    fig, ax = plt.subplots()
    ax.set_xlabel(r'$T$ [$\degree$C]')
    ax.set_ylabel(r'$L_v$ [J/kg]')
    ax.axis([25, 125, 2.9E5, 6.8E5])
    ax.xaxis.set_major_locator(plt.MultipleLocator(25))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(12.5))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.5E5))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25E5))
    ax.plot(T_droplet-273.15, Lv(T_droplet,540.13,100.21,50.28826,-0.217),
        ls='-', color='k',
        label='MFSim'
    )
    ax.plot(T_droplet-273.15, Lv_CP(T_droplet, "n-Heptane"),
        ls='--', color='gray',
        label='Coolprop'
    )
    ax.set_title('n-Heptano')
    ax.grid(color='lightgrey',ls='-.')
    ax.legend(facecolor="white", framealpha=1, frameon=1)
    fig.tight_layout(pad=0.15)

    fig, ax = plt.subplots()
    ax.set_xlabel(r'$T$ [$\degree$C]')
    ax.set_ylabel(r'$\rho_l$ [kg/m³]')
    ax.axis([25, 125, 575, 700])
    ax.xaxis.set_major_locator(plt.MultipleLocator(25))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(12.5))
    ax.yaxis.set_major_locator(plt.MultipleLocator(25))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(12.5))
    ax.plot(T_droplet-273.15, rho(T_droplet,540.13,0.2326,0.26042,0.281),
        ls='-', color='k',
        label='MFSim'
    )
    ax.plot(T_droplet-273.15, rho_CP(T_droplet, P, "n-Heptane"),
        ls='--', color='gray',
        label='Coolprop'
    )
    ax.set_title('n-Heptano')
    ax.grid(color='lightgrey',ls='-.')
    ax.legend(facecolor="white", framealpha=1, frameon=1)
    fig.tight_layout(pad=0.15)

    fig, ax = plt.subplots()
    ax.set_xlabel(r'$T$ [$\degree$C]')
    ax.set_ylabel(r'$C_{p,l}$ [J/kg.K]')
    ax.axis([25, 125, 250, 2600])
    ax.xaxis.set_major_locator(plt.MultipleLocator(25))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(12.5))
    ax.yaxis.set_major_locator(plt.MultipleLocator(500))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(250))
    ax.plot(T_droplet-273.15, cp(T_droplet,100.21,-98.7132,1.02021,-0.0031805,4.2735E-06),
        ls='-', color='k',
        label='MFSim'
    )
    ax.plot(T_droplet-273.15, cp_CP(T_droplet, P, "n-Heptane"),
        ls='--', color='gray',
        label='Coolprop'
    )
    ax.set_title('n-Heptano')
    ax.grid(color='lightgrey',ls='-.')
    ax.legend(facecolor="white", framealpha=1, frameon=1)
    fig.tight_layout(pad=0.15)

    fig, ax = plt.subplots()
    ax.set_xlabel(r'$T$ [$\degree$C]')
    ax.set_ylabel(r'$k_l$ [W/m.K]')
    ax.axis([25, 125, -0.4, 0.2])
    ax.xaxis.set_major_locator(plt.MultipleLocator(25))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(12.5))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.05))
    ax.plot(T_droplet-273.15, k(T_droplet,-0.2167,-2.9108E-4,-7.6582E-08),
        ls='-', color='k',
        label='MFSim'
    )
    ax.plot(T_droplet-273.15, k_CP(T_droplet, P, "n-Heptane"),
        ls='--', color='gray',
        label='Coolprop'
    )
    ax.set_title('n-Heptano')
    ax.grid(color='lightgrey',ls='-.')
    ax.legend(facecolor="white", framealpha=1, frameon=1)
    fig.tight_layout(pad=0.15)

    fig, ax = plt.subplots()
    ax.set_xlabel(r'$T$ [$\degree$C]')
    ax.set_ylabel(r'$\mu_l$ [Pa.s]')
    ax.axis([25, 125, 1.5E-4, 4E-4])
    ax.xaxis.set_major_locator(plt.MultipleLocator(25))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(12.5))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.5E-4))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25E-4))
    ax.plot(T_droplet-273.15, mu(T_droplet,-5.2620,7.4043E02,1.2279E-02,-1.4466E-05),
        ls='-', color='k',
        label='MFSim'
    )
    ax.plot(T_droplet-273.15, mu_CP(T_droplet, P, "n-Heptane"),
        ls='--', color='gray',
        label='Coolprop'
    )
    ax.set_title('n-Heptano')
    ax.grid(color='lightgrey',ls='-.')
    ax.legend(facecolor="white", framealpha=1, frameon=1)
    fig.tight_layout(pad=0.15)

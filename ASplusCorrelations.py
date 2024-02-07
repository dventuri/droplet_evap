import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from hamopy import ham_library as ham
from aux import *
from scipy.integrate import odeint

def dpm_model(y, t):
    dp, T_droplet_inlet = y

    # Input data
    T_gas_inlet = 418.15            # in [K], equals 145 °C
    P_system_init = 124642.54753    # in [Pa], equals 1.271 kgf/cm²
    # T_droplet_inlet = 333.15        # in [K], equals 60 °C
    u_dpm = 20
    u_gas = 20


    # Number of components
    N_gas = 6
    N_liq = 1


    # Bulk gas-phase
    gas = ct.Solution('hc_water_mod.yaml')
    gcomp1 = "C2H6:1.37E-02, NC6H14:3.45E-02, NC7H16:0.135360243, C9H19-1:0.361590943,"
    gcomp2 = " C11H21:0.243069983, NC13H28:3.14E-02, H2O:0.180422265"
    gas.TPY = T_gas_inlet, P_system_init, gcomp1+gcomp2
    gas()   # prints a summary of the cantera gas object

    # Bulk gas-phase physical properties
    gas_rho = gas.density_mass
    gas_mu = gas.viscosity
    gas_k = gas.thermal_conductivity
    gas_cp = gas.cp_mass
    gas_D = gas.mix_diff_coeffs[-1]
    gas_Y_vapor = gas.Y[-1]


    # Molar composition at droplet surface
    x_liq_surf = np.empty(N_gas + N_liq)
    liquid_vapor_pressure = ham.p_sat(T_droplet_inlet)  # ! No caso atual somente água liq
    x_liq_surf[-1] = liquid_vapor_pressure/P_system_init
    complementary_gas = 1-x_liq_surf[-1]
    sum_molar_gas_sup = sum(gas.X[:-1])
    for n_comp in range(N_gas):
        x_liq_surf[n_comp] = gas.X[n_comp]/sum_molar_gas_sup*complementary_gas
    x_liq_surf[0] += (1 - sum(x_liq_surf))  # ensure sum = 1 using first component as complementary

    # Mass composition at droplet surface
    Y_liq_surf = np.empty(N_gas + N_liq)
    sum_xi_Wi = sum(x_liq_surf*gas.molecular_weights)
    for n_comp in range(N_gas+N_liq):
        Y_liq_surf[n_comp] = x_liq_surf[n_comp]*gas.molecular_weights[n_comp]/sum_xi_Wi
    Y_liq_surf[0] += (1 - sum(Y_liq_surf))  # ensure sum = 1 using first component as complementary


    # Film around droplet
    T_film = calc_prop_at_droplet_vapor_layer(T_droplet_inlet, T_gas_inlet)
    Y_film = np.zeros(N_gas + N_liq)
    Y_film[-1] = calc_prop_at_droplet_vapor_layer(Y_liq_surf[-1], gas.X[-1])
    complementary_gas = 1-Y_film[-1]
    for n_comp in range(N_gas):
        Y_film[n_comp] = gas.Y[n_comp]/sum(gas.Y[:-1])*complementary_gas
    gas.TPY = T_film, P_system_init, Y_film

    # Film vapor-phase physical properties
    film_D = gas.mix_diff_coeffs[-1]
    film_rho = gas.density_mass


    # Film, only vapor
    gas.TPY = T_film, P_system_init, np.concatenate((np.zeros(N_gas),Y_film[-1]), axis=None)

    # Film vapor-phase physical properties
    film_vapor_mu = gas.viscosity
    film_vapor_k = gas.thermal_conductivity
    film_vapor_cp = gas.cp_mass


    # Film, remaining gas
    gas.TPY = T_film, P_system_init, np.concatenate((Y_film[:-1], 0), axis=None)

    # Film gas-phase physical properties
    film_gas_rho = gas.density_mass
    film_gas_mu = gas.viscosity
    film_gas_k = gas.thermal_conductivity
    film_gas_cp = gas.cp_mass
    film_gas_D = gas.mix_diff_coeffs[-1]


    # Correlations
    Re = gas_rho*abs(u_dpm - u_gas)*dp/gas_mu
    f_Re = calc_fRe(Re)
    Le = film_gas_k/(film_gas_rho*film_D*film_gas_cp)
    Pr = film_gas_cp*film_gas_mu/film_gas_k
    Sc = film_gas_mu/(film_gas_rho*film_D)
    Nu_0 = 1 + ((1 + Re*Pr)**(1/3))*f_Re
    Sh_0 = 1 + ((1 + Re*Sc)**(1/3))*f_Re


    # Mass transfer quantities
    Bm = (Y_liq_surf[-1] - gas_Y_vapor)/(1 - Y_liq_surf[-1])
    Fm = ((1 + Bm)**0.7)*(np.log(1+Bm)/Bm)
    Sh_star = 2 + (Sh_0 - 2)/Fm
    mdot = np.pi*film_gas_rho*film_gas_D*dp*Sh_star*np.log(1+Bm)


    # Thermal quantities
    Ql = 0
    Hl = PropsSI('H','P',P_system_init,'Q',1,'Water') - PropsSI('H','P',P_system_init,'Q',0,'Water')
    for k in range(50):
        Bt = film_vapor_cp*(T_gas_inlet - T_droplet_inlet)/(Hl + Ql/mdot)
        Ft = ((1 + Bt)**0.7)*(np.log(1+Bt)/Bt)
        Nu_star = 2 + (Nu_0 - 2)/Ft
        phi = (film_vapor_cp/film_gas_cp)*(Sh_star/Nu_star)*(1/Le)
        Bt_n = (1+Bm)**phi - 1
        Ql = mdot*((film_vapor_cp*(T_gas_inlet - T_droplet_inlet))/Bt_n - Hl)
        if(abs(Bt_n - Bt) < 1E-12): break

    drdt = - 2*mdot/(np.pi*980*dp**2)
    mp = np.pi/6*dp**3*980
    dTdt = Ql/(mp*4180)

    return [drdt, dTdt]


if __name__=='__main__':
    y0 = [1E-4, 333.15]
    t = np.linspace(0,1,100001)

    sol = odeint(dpm_model, y0, t)

    plt.figure()
    plt.plot(t*20, sol[:, 0]*10**6)
    plt.xlabel('Distance [m]')
    plt.ylabel('Dp [micron]')

    plt.figure()
    plt.plot(t*20, sol[:, 1]-273.15)
    plt.xlabel('Distance [m] ')
    plt.ylabel('Temperature (degree Celsius)')

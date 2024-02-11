import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from hamopy import ham_library as ham
from aux import *
from scipy.integrate import odeint

def dpm_model(y, t, Ambient_Gas, Liquid, Rel_vel):

    """
    Defines a model to solve the Abramzon-Sirignano evaporation
    of a single droplet.

    Receives a vector 'y' with [droplet_diameter, droplet_temperature]
    and a vector of the integration variable (in this instance, time 't')

    In the current implementation, the gas can be defined by any Cantera
    solution with any number of components.
    However, the droplet must be monocomponent (uses CoolProp for physical
    properties).

    'Ambient_Gas' is a Cantera Solution Object, where the evaporating species
    that comes from the droplet must always be last in the vector order.

    'Liquid' is the name of the liquid species in the Coolprop database.
    """

    # Unpack vars
    dp, T_droplet = y
    T_gas = Ambient_Gas.T
    P_amb = Ambient_Gas.P
    Gas_composition = Ambient_Gas.X


    # Number of components
    N_spt = len(Ambient_Gas.X)
    N_gas = N_spt - 1
    N_liq = 1                   # ! Hard-coded for now


    # Bulk gas-phase physical properties
    gas_rho = Ambient_Gas.density_mass
    gas_mu = Ambient_Gas.viscosity
    gas_k = Ambient_Gas.thermal_conductivity
    gas_cp = Ambient_Gas.cp_mass
    gas_D = Ambient_Gas.mix_diff_coeffs[-1]
    gas_Y_vapor = Ambient_Gas.Y[-1]


    # Molar composition at droplet surface
    Liquid_name = Liquid.strip().lower()
    x_liq_surf = np.empty(N_spt)
    if(Liquid_name == 'water'):
        liquid_vapor_pressure = ham.p_sat(T_droplet)
    else:
        liquid_vapor_pressure = PropsSI('P','T',T_droplet,'Q',1,Liquid_name)
    x_liq_surf[-1] = liquid_vapor_pressure/P_amb
    complementary_gas = 1-x_liq_surf[-1]
    sum_molar_gas_sup = np.sum(Ambient_Gas.X[:-1])
    for n_comp in range(N_gas):
        x_liq_surf[n_comp] = Ambient_Gas.X[n_comp]/sum_molar_gas_sup*complementary_gas
    x_liq_surf[0] += (1 - sum(x_liq_surf))  # ensure sum = 1 using first component as complementary


    # Mass composition at droplet surface
    Y_liq_surf = np.empty(N_gas + N_liq)
    sum_xi_Wi = sum(x_liq_surf*Ambient_Gas.molecular_weights)
    for n_comp in range(N_gas+N_liq):
        Y_liq_surf[n_comp] = x_liq_surf[n_comp]*Ambient_Gas.molecular_weights[n_comp]/sum_xi_Wi
    Y_liq_surf[0] += (1 - sum(Y_liq_surf))  # ensure sum = 1 using first component as complementary


    # Film around droplet
    T_film = calc_prop_at_droplet_vapor_layer(T_droplet, T_gas)
    Y_film = np.zeros(N_spt)
    Y_film[-1] = calc_prop_at_droplet_vapor_layer(Y_liq_surf[-1], Ambient_Gas.X[-1])
    complementary_gas = 1-Y_film[-1]
    for n_comp in range(N_gas):
        Y_film[n_comp] = Ambient_Gas.Y[n_comp]/sum(Ambient_Gas.Y[:-1])*complementary_gas
    Ambient_Gas.TPY = T_film, P_amb, Y_film

    # Film vapor-phase physical properties
    film_D = Ambient_Gas.mix_diff_coeffs[-1]
    film_rho = Ambient_Gas.density_mass


    # Film, only vapor
    Ambient_Gas.TPY = T_film, P_amb, np.concatenate((np.zeros(N_gas),Y_film[-1]), axis=None)

    # Film vapor-phase physical properties
    film_vapor_mu = Ambient_Gas.viscosity
    film_vapor_k = Ambient_Gas.thermal_conductivity
    film_vapor_cp = Ambient_Gas.cp_mass


    # Film, remaining gas
    Ambient_Gas.TPY = T_film, P_amb, np.concatenate((Y_film[:-1], 0), axis=None)

    # Film gas-phase physical properties
    film_gas_rho = Ambient_Gas.density_mass
    film_gas_mu = Ambient_Gas.viscosity
    film_gas_k = Ambient_Gas.thermal_conductivity
    film_gas_cp = Ambient_Gas.cp_mass
    film_gas_D = Ambient_Gas.mix_diff_coeffs[-1]


    # Restore original Ambient Gas for next loop
    Ambient_Gas.TPX = T_gas, P_amb, Gas_composition


    # Correlations
    Re = gas_rho*abs(Rel_vel)*dp/gas_mu
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
    Hl = PropsSI('H','P',P_amb,'Q',1,Liquid_name) - PropsSI('H','P',P_amb,'Q',0,Liquid_name)
    for k in range(50):
        Bt = film_vapor_cp*(T_gas - T_droplet)/(Hl + Ql/mdot)
        Ft = ((1 + Bt)**0.7)*(np.log(1+Bt)/Bt)
        Nu_star = 2 + (Nu_0 - 2)/Ft
        phi = (film_vapor_cp/film_gas_cp)*(Sh_star/Nu_star)*(1/Le)
        Bt_n = (1+Bm)**phi - 1
        Ql = mdot*((film_vapor_cp*(T_gas - T_droplet))/Bt_n - Hl)
        if(abs(Bt_n - Bt) < 1E-12): break


    rho_liq = PropsSI('D','T',T_droplet,'Q',0,Liquid_name)
    cp_liq = PropsSI('C','T',T_droplet,'Q',0,Liquid_name)
    ddpdt = - 2*mdot/(np.pi*rho_liq*dp**2)
    mp = np.pi/6*dp**3*rho_liq
    dTdt = Ql/(mp*cp_liq)


    return [ddpdt, dTdt]


if __name__=='__main__':

    # Ranz-Marshall: Water droplet in stagnat air
    gas = ct.Solution('air.yaml')
    gcomp = "O2:0.21, N2:0.78, AR:0.01, H2O:0"
    gas.TPX = 298, 101325, gcomp

    y0 = [1.048E-3, 282]
    t = np.linspace(0,800,800001)

    sol = odeint(dpm_model, y0, t, args=(gas, 'Water', 0))

    np.savetxt("water_air-dp.dat",
               np.column_stack((t, sol[:,0])),
               fmt=['%.3f','%10.5e'],
               delimiter=',')

    np.savetxt("water_air-Temp.dat",
               np.column_stack((t, sol[:,1])),
               fmt=['%.3f','%10.5e'],
               delimiter=',')


    # Water droplet in hot hydrocarbon gas
    gas = ct.Solution('hc_water_mod.yaml')
    gcomp1 = "C2H6:1.37E-02, NC6H14:3.45E-02, NC7H16:0.135360243, C9H19-1:0.361590943,"
    gcomp2 = " C11H21:0.243069983, NC13H28:3.14E-02, H2O:0.180422265"
    gas.TPY = 418.15, 124642.54753, gcomp1+gcomp2

    y0 = [1E-4, 333.15]
    t = np.linspace(0,1,100001)

    sol = odeint(dpm_model, y0, t, args=(gas, 'Water', 0))

    np.savetxt("water_hcs-dp.dat",
               np.column_stack((t, sol[:,0])),
               fmt=['%.3f','%10.5e'],
               delimiter=',')

    np.savetxt("water_hcs-Temp.dat",
               np.column_stack((t, sol[:,1])),
               fmt=['%.3f','%10.5e'],
               delimiter=',')

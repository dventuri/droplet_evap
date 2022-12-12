from math import pi, sqrt, log
from CoolProp.CoolProp import PropsSI
from hamopy import ham_library as ham

# water on air
Wwater = 18.01528   # molecular weight
Wair = 28.9647

# ambient conditions, far from droplet
Patm = 101325       # Pa
P = Patm            # Pa
T_gas = 417.5         # K
Y_vapor_inf = 0.2
x_vapor_inf = (Y_vapor_inf/Wwater)/(
    (Y_vapor_inf/Wwater) + ((1-Y_vapor_inf)/Wair)
)

# droplet conditions, surface
T_droplet = 333   # K
d_droplet = 1e-3
P_Fs = ham.p_sat(T_droplet)
x_Fs = P_Fs/P
Y_Fs = x_Fs*Wwater/(x_Fs*Wwater + (1-x_Fs)*Wair)

# convection
Urel = abs(20-15)

# near-droplet conditions, 1/3 rule
T_mix = T_droplet + 1/3*(T_gas - T_droplet)
Y_vapor_mix = Y_Fs + 1/3*(Y_vapor_inf - Y_Fs)
Y_air_mix = 1 - Y_vapor_mix
x_vapor_mix = x_Fs + 1/3*(x_vapor_inf - x_Fs)
x_air_mix = 1 - x_vapor_mix

#mu_mix Wilke 1950 rule
mu_vapor_Tmix = PropsSI('V', 'T', T_mix, 'Q', 1, 'water')
mu_air_Tmix   = PropsSI('V', 'T', T_mix, 'P', P, 'air')

mu_mix = (
    x_vapor_mix*mu_vapor_Tmix/(
        (
            x_vapor_mix*(1 + sqrt(mu_vapor_Tmix/mu_vapor_Tmix)*(Wwater/Wwater)**(1/4))**2 / (sqrt(8*(1 + Wwater/Wwater)))
        ) + 
        (
            x_air_mix*(1 + sqrt(mu_vapor_Tmix/mu_air_Tmix)*(Wair/Wwater)**(1/4))**2 / (sqrt(8*(1 + Wwater/Wair)))
        )
    ) + 
    x_air_mix*mu_air_Tmix/(
        (
            x_vapor_mix*(1 + sqrt(mu_air_Tmix/mu_vapor_Tmix)*(Wwater/Wair)**(1/4))**2 / (sqrt(8*(1 + Wair/Wwater)))
        ) +
        (
            x_air_mix*(1 + sqrt(mu_air_Tmix/mu_air_Tmix)*(Wair/Wair)**(1/4))**2 / (sqrt(8*(1 + Wair/Wair)))
        )
    )
)

#k_mix Wilke 1950 rule
k_vapor_Tmix = PropsSI('L', 'T', T_mix, 'Q', 1, 'water')
k_air_Tmix   = PropsSI('L', 'T', T_mix, 'P', P, 'air')

k_mix = (
    x_vapor_mix*k_vapor_Tmix/(
        (
            x_vapor_mix*(1 + sqrt(k_vapor_Tmix/k_vapor_Tmix)*(Wwater/Wwater)**(1/4))**2 / (sqrt(8*(1 + Wwater/Wwater)))
        ) + 
        (
            x_air_mix*(1 + sqrt(k_vapor_Tmix/k_air_Tmix)*(Wair/Wwater)**(1/4))**2 / (sqrt(8*(1 + Wwater/Wair)))
        )
    ) + 
    x_air_mix*k_air_Tmix/(
        (
            x_vapor_mix*(1 + sqrt(k_air_Tmix/k_vapor_Tmix)*(Wwater/Wair)**(1/4))**2 / (sqrt(8*(1 + Wair/Wwater)))
        ) +
        (
            x_air_mix*(1 + sqrt(k_air_Tmix/k_air_Tmix)*(Wair/Wair)**(1/4))**2 / (sqrt(8*(1 + Wair/Wair)))
        )
    )
)

#rho_mix
rho_vapor_Tmix = PropsSI('D', 'T', T_mix, 'Q', 1, 'water')
rho_air_Tmix   = PropsSI('D', 'T', T_mix, 'P', P, 'air')
rho_mix = Y_vapor_mix*rho_vapor_Tmix + Y_air_mix*rho_air_Tmix

#cp_mix
cp_vapor_Tmix = PropsSI('C', 'T', T_mix, 'Q', 1, 'water')
cp_air_Tmix   = PropsSI('C', 'T', T_mix, 'P', P, 'air')
cp_mix = Y_vapor_mix*cp_vapor_Tmix + Y_air_mix*cp_air_Tmix

#Diffusivity @Tmix - melhorar, só considera água em ar ??
Diff_mix = ham.D_va(T_mix)

#Lewis mix
Le_mix = k_mix/(rho_mix*Diff_mix*cp_mix)

#Prandtl mix
Pr_mix = cp_mix*mu_mix/k_mix

#Schmidt mix
Sc_mix = mu_mix/(rho_mix*Diff_mix)

#Re
rho_air_Tinf = PropsSI('D', 'T', T_gas, 'P', P, 'air')
Re_mix = rho_air_Tinf*Urel*d_droplet/mu_mix

#Nussel
Nu = 2 + 0.6*sqrt(Re_mix)*(Pr_mix**(1/3))

# Sherwood
Sh = 2 + 0.6*sqrt(Re_mix)*(Sc_mix**(1/3))

#mass transfer
Bm = (Y_Fs - Y_vapor_inf)/(1 - Y_Fs)
Fm = ((1 + Bm)**0.7)*log(1 + Bm)/Bm
Sh_star = 2 + (Sh - 2)/Fm
mdot = pi*rho_mix*Diff_mix*d_droplet*Sh_star*log(1 + Bm)

# Saturated vapor enthalpy of Water at 1 atm in J/kg
H_V = PropsSI('H','P',101325,'Q',1,'Water')

# Saturated liquid enthalpy of Water at 1 atm in J/kg
H_L = PropsSI('H','P',101325,'Q',0,'Water')

# Latent heat of vaporization of Water at 1 atm in J/kg
dH_water = H_V - H_L


Bt = cp_vapor_Tmix*(T_gas - T_droplet)/(dH_water + 0/mdot)
Ft = ((1 + Bt)**0.7)*log(1 + Bt)/Bt
Nu_star = 2 + (Nu -2)/Ft
phi = (cp_vapor_Tmix/cp_mix)*(Sh_star/Nu_star)*1/Le_mix
Bt = (1 + Bm)**phi - 1

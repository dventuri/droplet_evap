from math import pi, sqrt, log
from CoolProp.CoolProp import PropsSI
from hamopy import ham_library as ham

# water on air

Patm = 101325       # Pa
P = Patm            # Pa
T_droplet = 332.5   # K
T_gas = 418         # K
T_mix = T_droplet + 1/3*(T_gas - T_droplet)

#########

P_Fs = ham.p_sat(T_droplet)
x_Fs = P_Fs/P
Y_Fs = x_Fs*18.01528/(x_Fs*18.01528 + (1-x_Fs)*28.9647)

#mu_mix Wilke 1950 rule
mu_vapor_Tmix = PropsSI('V', 'T', T_mix, 'Q', 1, 'water')
mu_air_Tmix   = PropsSI('V', 'T', T_mix, 'P', P, 'air')

mu_mix = (
    x_Fs*mu_vapor_Tmix/(
        (
            x_Fs*(1 + sqrt(mu_vapor_Tmix/mu_vapor_Tmix)*(18.01528/18.01528)**(1/4))**2 / (sqrt(8*(1 + 18.01528/18.01528)))
        ) + 
        (
            (1-x_Fs)*(1 + sqrt(mu_vapor_Tmix/mu_air_Tmix)*(28.9647/18.01528)**(1/4))**2 / (sqrt(8*(1 + 18.01528/28.9647)))
        )
    ) + 
    (1 - x_Fs)*mu_air_Tmix/(
        (
            x_Fs*(1 + sqrt(mu_air_Tmix/mu_vapor_Tmix)*(18.01528/28.9647)**(1/4))**2 / (sqrt(8*(1 + 28.9647/18.01528)))
        ) +
        (
            (1-x_Fs)*(1 + sqrt(mu_air_Tmix/mu_air_Tmix)*(28.9647/28.9647)**(1/4))**2 / (sqrt(8*(1 + 28.9647/28.9647)))
        )
    )
)

#k_mix Wilke 1950 rule
k_vapor_Tmix = PropsSI('L', 'T', T_mix, 'Q', 1, 'water')
k_air_Tmix   = PropsSI('L', 'T', T_mix, 'P', P, 'air')

k_mix = (
    x_Fs*k_vapor_Tmix/(
        (
            x_Fs*(1 + sqrt(k_vapor_Tmix/k_vapor_Tmix)*(18.01528/18.01528)**(1/4))**2 / (sqrt(8*(1 + 18.01528/18.01528)))
        ) + 
        (
            (1-x_Fs)*(1 + sqrt(k_vapor_Tmix/k_air_Tmix)*(28.9647/18.01528)**(1/4))**2 / (sqrt(8*(1 + 18.01528/28.9647)))
        )
    ) + 
    (1 - x_Fs)*k_air_Tmix/(
        (
            x_Fs*(1 + sqrt(k_air_Tmix/k_vapor_Tmix)*(18.01528/28.9647)**(1/4))**2 / (sqrt(8*(1 + 28.9647/18.01528)))
        ) +
        (
            (1-x_Fs)*(1 + sqrt(k_air_Tmix/k_air_Tmix)*(28.9647/28.9647)**(1/4))**2 / (sqrt(8*(1 + 28.9647/28.9647)))
        )
    )
)

#rho_mix
rho_vapor_Tmix = PropsSI('D', 'T', T_mix, 'Q', 1, 'water')
rho_air_Tmix   = PropsSI('D', 'T', T_mix, 'P', P, 'air')
rho_mix = Y_Fs*rho_vapor_Tmix + (1-Y_Fs)*rho_air_Tmix

#Cp_vapor_mix
cp_vapor_Tmix = PropsSI('C', 'T', T_mix, 'Q', 1, 'water')

#Cp_mix
cp_air_Tmix   = PropsSI('C', 'T', T_mix, 'P', P, 'air')
cp_mix = Y_Fs*cp_vapor_Tmix + (1-Y_Fs)*cp_air_Tmix

#Diffusivity @Tmix - melhorar, só considera água em ar ??
Diff_mix = ham.D_va(T_mix)

#Lewis mix
Le_mix = k_mix/(rho_mix*Diff_mix*cp_mix)

#Prandtl mix
Pr_mix = cp_mix*mu_mix/k_mix

#Schmidt mix
Sc_mix = mu_mix/(rho_mix*Diff_mix)

#Re
Urel = 5
d_droplet = 1e-4
rho_air_Tinf = PropsSI('D', 'T', T_gas, 'P', P, 'air')
Re_mix = rho_air_Tinf*Urel*d_droplet/mu_mix


#Nussel
Nu = 2 + 0.6*sqrt(Re_mix)*(Pr_mix**(1/3))

# Sherwood
Sh = 2 + 0.6*sqrt(Re_mix)*(Sc_mix**(1/3))


#mass transfer
Bm = (Y_Fs - 0)/(1 - Y_Fs)
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

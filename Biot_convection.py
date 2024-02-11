# From Incropera et al., 2011, 6. ed.
# Biot < 0.1

def water_thermal_conductivity(Temperature):
    '''
    Returns liquid water thermal conductivity for given temperature in Kelvin
    from Yaws, 2015, Handbook of physical properties
    '''
    if (Temperature<273.15 or Temperature>635):
        print("WARNING: Temperature out of range for water thermal conductivity estimation")

    return -3.5667E-1 + 5.057E-3*Temperature - 6.1071E-6*Temperature**2


def water_density(Temperature):
    '''
    Returns liquid water density for given temperature in Kelvin
    from Yaws, 2015, Handbook of physical properties
    '''
    if (Temperature<290 or Temperature>647.13):
        print("WARNING: Temperature out of range for water density estimation")

    return 0.325 * 0.27**(-(1.0-Temperature/647.13)**0.23)*1000


def water_viscosity(Temperature):
    '''
    Returns liquid water viscosity for given temperature in Kelvin
    from Yaws, 2015, Handbook of physical properties
    '''
    if (Temperature<273.15 or Temperature>646.15):
        print("WARNING: Temperature out of range for water viscosity estimation")

    water_visc = 10**(
        -11.6225 + 
        1.9490E+03/Temperature +
        2.1641E-02*Temperature +
        -1.5990E-05*Temperature**2
    )/1000

    return water_visc


def air_thermal_conductivity(Temperature):
    '''
    Returns air thermal conductivity for given temperature in Kelvin
    from Yaws, 2015, Handbook of physical properties
    '''
    if (Temperature<100 or Temperature>1500):
        print("WARNING: Temperature out of range for air thermal conductivity estimation")

    return -3.8603E-04 + 1.0311E-04*Temperature - 5.4199E-08*Temperature**2 + 1.7429E-11*Temperature**3


def air_viscosity(Temperature):
    '''
    Returns air viscosity for given temperature in Kelvin
    from Yaws, 2015, Handbook of physical properties
    '''
    if (Temperature<100 or Temperature>1500):
        print("WARNING: Temperature out of range for air viscosity estimation")
    
    visc = 4.5608 + 7.0077E-1*Temperature - 3.7287E-4*Temperature**2 + 9.0437E-8

    return visc*1E-7


def Reynolds(density, velocity, characteristic_lenght, viscosity):
    '''
    Returns Reynolds number
    '''
    return density*velocity*characteristic_lenght/viscosity


def nusselt_sphere(Reynolds, Prandtl, viscosity_ratio):
    '''
    Returns Nusselt number for a sphere according to Whitaker's correlation (1972)
    Reference: Incropera et al., 2011, 6. ed.
    '''
    if (Reynolds<3.5 or Reynolds>7.6E4):
        print("WARNING: Reynolds out of range for Nusselt calculation")
    if (Prandtl<0.71 or Prandtl>380):
        print("WARNING: Prandtl out of range for Nusselt calculation")
    if (viscosity_ratio<1.0 or viscosity_ratio>3.2):
        print("WARNING: Viscosity ratio out of range for Nusselt calculation")
    
    Nusselt = 2.0 + (0.4*Reynolds**0.5 + 0.06*Reynolds**(2/3))*(Prandtl**0.4)*(viscosity_ratio**0.25)

    return Nusselt


def convective_coefficient(Nusselt, thermal_conductivity, characteristic_lenght):
    
    return Nusselt*thermal_conductivity/characteristic_lenght


def biot_sphere(convective_coefficient, radius, thermal_conductivity):

    return convective_coefficient*(radius/3)/thermal_conductivity


D = 5e-4

Pr_air = 0.708
mu = air_viscosity(418)
mu_s = air_viscosity(332.5)
rho_air = 1.2
k_air = air_thermal_conductivity(418)

Re = Reynolds(rho_air,10.0,D,mu)
Nu = nusselt_sphere(Re, Pr_air, mu/mu_s)
h = convective_coefficient(Nu, k_air, D)

k_water = water_thermal_conductivity(332.5)
Bi = biot_sphere(h,D/2,k_water)
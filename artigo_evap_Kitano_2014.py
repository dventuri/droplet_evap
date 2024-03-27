import numpy as np
import cantera as ct
from scipy.integrate import odeint
from AbramzonSirignano import as_evap

# Kitano/Nomura: n-heptane droplet in stagnat air
gas = ct.Solution('N2_n-Heptane.yaml')
gas.transport_model = 'UnityLewis'
# gcomp = "O2:0.21, N2:0.78, AR:0.01, H2O:0"
gcomp = "N2:1"
Ts = [471, 468, 466, 452]
Ps = [0.1E6, 0.5E6, 1.0E6, 2.0E6]

y0 = [700E-6, 300]
t = np.linspace(0,8,8000001)

for count, (T, P) in enumerate(zip(Ts, Ps)):
    gas.TPX = T, P, gcomp

    sol = odeint(as_evap, y0, t, args=(gas, 'n-Heptane', 0))

    np.savetxt("Kitano-0"+str(count)+"-dp.dat",
                np.column_stack((t[::1000], sol[::1000,0])),
                fmt=['%.3f','%e'],
                delimiter=',')

    np.savetxt("Kitano-0"+str(count)+"-Temp.dat",
                np.column_stack((t[::1000], sol[::1000,1])),
                fmt=['%.3f','%e'],
                delimiter=',')

import numpy as np

T_film = 300

sigma = np.array([3.621, 5.644])
sigma_va = np.sum(sigma)/2

epsilon = np.array([97.53, 537.597])
epsilon_va = np.sqrt(epsilon[0]*epsilon[1])

Mv = 28.013
Ma = 100.21
Mva = 2/(1/Mv + 1/Ma)

T_mod = T_film/epsilon_va
Fc_D = 1.06036/(T_mod**0.15610) + 0.19300/np.exp(0.47635*T_mod) + 1.03587/np.exp(1.52996*T_mod) + 1.76474/np.exp(3.89411*T_mod)
film_D = (3.03 - 0.98/np.sqrt(Mva))*(T_film**1.5)/(1.01325*sigma_va**2*Fc_D*np.sqrt(Mva))*1E-3
print("D_prausnitz = ", film_D/100**2)

film_D = (3.03 - 0.98/np.sqrt(Mva))*(T_film**1.5)/(101325*sigma_va**2*Fc_D*np.sqrt(Mva))*1E-2
print("D_prausnitz_mod = ", film_D)

film_D = (3.03 - 0.98/np.sqrt(Mva))*(T_film**1.5)/(10.1325*sigma_va**2*Fc_D*np.sqrt(Mva))*1E-7
print("D_mfsim = ", film_D)

film_D = (2.17 - 1/np.sqrt(2*Mva))*(np.sqrt(2)*T_film**1.5)/(1*sigma_va**2*Fc_D*np.sqrt(Mva))*1E-3
print("D_cremasco = ", film_D/100**2)

import numpy as np
import matplotlib.pyplot as plt
import scipy

#Follow https://arxiv.org/pdf/0805.2708

gamma = 1 #to be given by Steck for Rb and Cs
C = 1 #indpendednt parameters to be modelled
delta = 1
t= 1 
phi = 0


def L_n(n,Delta,omega_m):
    a = (Delta-n*omega_m)
    return gamma**2/(gamma**2 + a**2)

def D_n(n,Delta,omega_m):
    a = (Delta-n*omega_m)
    return gamma*a/(gamma**2 + a**2)

def Beat(omega_m, Delta): #Consider only carrier and first order sidebands
    Bess_0 = scipy.special.jv(0,delta)
    Bess_1 = scipy.special.jv(1,delta)
    term_0 = (C/np.sqrt(gamma**2 + omega_m**2))*Bess_0*Bess_1
    term_1 = (L_n(-1,Delta,omega_m) - L_n(-1/2,Delta,omega_m) + 
              L_n(1/2,Delta,omega_m) - L_n(1,Delta,omega_m))*np.cos(omega_m*t+phi)
    term_2 = (D_n(-1,Delta,omega_m) - D_n(-1/2,Delta,omega_m) - 
              D_n(1/2,Delta,omega_m) + D_n(1,Delta,omega_m))*np.sin(omega_m*t+phi)
    return term_1, term_2, term_0*(term_1 + term_2)

print(Beat(1,1))

Delta_array = np.linspace(-10,10,1000)
plot = plt.figure(figsize=(10, 6),dpi =300)
plt.plot(Delta_array, Beat(1,Delta_array)[0], label='Beat Signal', color='blue')
plt.show()
# Single particle width calculations using the wave function from FRESCO CRC  
#
# Author: Eli Temanson
# Date: April, 2023
# 
# requirements: matplotlib, pandas, scipy, mpmath, numpy

import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate, integrate
import mpmath as mp
import numpy as np

HBARC=197.3269804  #MeV*fm
# mu=901.909/931.494 # amu 
mu=1007825.031898*24990428.31/(1007825.031898 + 24990428.31) * 1E-6 # amu
Zp=1
Zt=13

E=0.4154 # MeV 0.4154(8)
#E=0.4146 # MeV low limit
#E=0.4162 # MeV upper limit
k=0.218735*np.sqrt(mu*E)
print('wave number, k = ', "{:6f}".format(k))
eta=0.1574854*Zp*Zt*np.sqrt(mu/E)
print('sommerfeld parameter, coulomb eta = ', "{:6f}".format(eta))

r0=1.25
r0=r0-r0*0.04 # 5% change
ChannelRadius=r0*(1.0 + 25.0**(1.0/3.0))
print('channel radius = ', "{:6f}".format(ChannelRadius), ' (fm)')

print('wavelength to sie ratio = ', "{:6f}".format(k*ChannelRadius))

# regular coulomb wave fuction
def regCWF(l,eta,z):
    func = np.frompyfunc(mp.coulombf, 3, 1)
    return func(l,eta,z)
# irregular coulomb wave function
def irregCWF(l,eta,z): 
    func = np.frompyfunc(mp.coulombg, 3, 1)
    return func(l,eta,z)

wave_function = pd.read_csv('fort.58',sep='\s+',header=None,skiprows=1,skipfooter=1,engine='python')

# Normalize the u^2 wavefunction
norm = interpolate.UnivariateSpline(wave_function[0],wave_function[1]**2).integral(0,30)
# norm = interpolate.UnivariateSpline(wave_function[0],wave_function[1]**2).integral(0,ChannelRadius)
wave_function[1] = abs(wave_function[1])

# print('norm = ',norm)

wave_function[1] = wave_function[1]/np.sqrt(norm)
wf = interpolate.CubicSpline(wave_function[0], wave_function[1])

reduced_width_sqr = 0.5*ChannelRadius*wf(ChannelRadius)**2
print('Dimenisonless single-particle reduced width = ', "{:6f}".format(reduced_width_sqr))

penetrability = k*ChannelRadius/(regCWF(0,eta,k*ChannelRadius)**2 + irregCWF(0,eta,k*ChannelRadius)**2)
print('penetrability = ' + mp.nstr(penetrability, 8))

print('Single-Particle Width = ' \
    + mp.nstr(reduced_width_sqr*penetrability*2*HBARC**2/(931.494*mu*ChannelRadius**2)*1E6)+ ' (eV)')

fig, ax = plt.subplots()
# plt.style.use('fivethirtyeight')

ax.plot(wave_function[0], wf(wave_function[0])**2)
ax.axvline(x=ChannelRadius,color='r',linestyle='--',label='Channel Radius = {:.2f}'.format(ChannelRadius))
# ax.axvline(x=4.31,color='k',linestyle='--')
# ax.axvline(x=6.37,color='k',linestyle='--')
ax.set_xlim([0,20])
ax.set_ylabel('$u(r)^{2}$',fontsize=14)
ax.set_xlabel('Radius (fm)',fontsize=14)
ax.legend()

# eq1 = (r"$\Gamma_p = C^{2}S_{l=0} \frac{\hbar^2}{\mu}$")
# ax.text(10, 0.2, eq1, color="C2", fontsize=18)

plt.show()



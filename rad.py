# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 19:16:16 2022

@author: trubin
"""

import astropy.constants as const
import astropy.units as u
import xraylib as xrl

import matplotlib.pyplot as plt
import numpy as np
import scipy.special
from plasmapy.formulary.radiation import thermal_bremsstrahlung

def integrate(integrand,r):
    res = np.zeros(integrand.shape)
    for i in range(res.shape[0]):
        res[i] = np.trapz(integrand[0:i],r[0:i])
    return res




# Shield element
element = "U"
Element_num = xrl.SymbolToAtomicNumber(element)
density = xrl.ElementDensity(Element_num)  #g/cm^3


#looking at a monoenergetic photon beam

Energy = 200 #kev
Intensity = 100 #optical power per unit area W/m^2


sigma_tot = xrl.CS_Total(Element_num,Energy) # cm^2/g

# slab
I = lambda x,x0,mu: np.exp(-mu*(x-x0)) # x in cm
# cylinder
Ir = lambda x,x0,mu: np.exp(-mu*(x-x0))/x # x in cm

# =============================================================================
# x = np.linspace(0,10,1000)+1
# plt.title("intensity [W/cm^2]")
# plt.loglog(x,Intensity*I(x,x[0],sigma_tot*density),label = "Slab, "+element)
# plt.loglog(x,Intensity*Ir(x,x[0],sigma_tot*density),label = "Cylinder, "+element)
# 
# element = "Fe"
# Element_num = xrl.SymbolToAtomicNumber(element)
# density = xrl.ElementDensity(Element_num)  #g/cm^3
# sigma_tot = xrl.CS_Total(Element_num,Energy) # cm^2/g
# plt.loglog(x,Intensity*I(x,x[0],sigma_tot*density),label = "Slab, "+element)
# plt.loglog(x,Intensity*Ir(x,x[0],sigma_tot*density),label = "Cylinder, "+element)
# 
# 
# element = "Pb"
# Element_num = xrl.SymbolToAtomicNumber(element)
# density = xrl.ElementDensity(Element_num)  #g/cm^3
# sigma_tot = xrl.CS_Total(Element_num,Energy) # cm^2/g
# plt.loglog(x,Intensity*I(x,x[0],sigma_tot*density),label = "Slab, "+element)
# plt.loglog(x,Intensity*Ir(x,x[0],sigma_tot*density),label = "Cylinder, "+element)
# 
# element = "U"
# Element_num = xrl.SymbolToAtomicNumber(element)
# density = xrl.ElementDensity(Element_num)  #g/cm^3
# sigma_tot = xrl.CS_Total(Element_num,Energy) # cm^2/g
# plt.loglog(x,Intensity*I(x,x[0],sigma_tot*density),label = "Slab, "+element)
# plt.loglog(x,Intensity*Ir(x,x[0],sigma_tot*density),label = "Cylinder, "+element)
# 
# plt.legend()
# plt.tight_layout()
# plt.show()
# =============================================================================

# Now, looking at the continuous spectrum:
# The power spectral density of the plasma:

# electron temperature
Te = 6e2 * u.keV
ion_species_B = "B-11 5+"
ion_species_p = "H-1 1+"

# B11 density
nB = 3/7*1e21 * u.cm ** -3
# proton density
nH = 4/7*1e21 * u.cm ** -3

ne = 5*nB+nH

device_radius = 1 * u.m


frequencies = np.arange(15, 19.9, 0.001)
frequencies = (10 ** frequencies) / u.s  #Hz

energies = (frequencies * const.h.si).to(u.keV) #keV

omega = 2*np.pi*frequencies # rad/s


spectrum_B = thermal_bremsstrahlung(frequencies, ne, Te, nB, ion_species=ion_species_B).to(u.W * u.s / u.m ** 3)
spectrum_p = thermal_bremsstrahlung(frequencies, ne, Te, nH, ion_species=ion_species_p).to(u.W * u.s / u.m ** 3)
spectrum = spectrum_p+spectrum_B

# =============================================================================
# 
# lbl = "$T_e$ = {:.1e} eV,  ".format(Te.value) + "$n_e$ = {:.1e} 1/cm^3".format(ne.value)
# plt.loglog(energies, spectrum_B, label="B")
# plt.loglog(energies, spectrum_p, label="p")
# plt.loglog(energies, spectrum, label="B+p")
# 
# plt.title(f"Thermal Bremsstrahlung Spectrum"+lbl)
# plt.xlabel("Energy (keV)")
# plt.ylabel("Power Spectral Density (W s/m^3)")
# plt.legend()
# plt.show()
# =============================================================================

 
vol = np.pi * device_radius ** 2

power_spectrum = 0.5*spectrum*device_radius ** 2 # W / (m rad/s)


sigma_tot = np.zeros(energies.shape)
for i in range(energies.shape[0]):
    try: sigma_tot[i] = xrl.CS_Total(Element_num,energies[i].to_value())  # cm^2/g
    except: 
        sigma_tot[i] = np.nan
sigma_tot[np.isnan(sigma_tot)] = sigma_tot[np.isfinite(sigma_tot)][0]

mu = sigma_tot*density/ u.cm
        
# =============================================================================
# plt.semilogy(energies, spectrum, label = "spectrum x=a")
# 
# plt.semilogy(energies, spectrum*Ir((1+0.0001)*device_radius.to_value(u.cm),device_radius.to_value(u.cm),mu.value), label = "spectrum x=1.00001*a")
# plt.legend()
# plt.show()
# =============================================================================

dw = np.gradient(omega)  # Frequency step size

#T = lambda r, r0, mu: np.trapz(power_spectrum * (np.exp(mu*r0)*(scipy.special.expi(-mu*r)-scipy.special.expi(-mu*r0))-np.log(r/r0)),omega)


kT_prime = lambda r, r0, mu: np.sum(power_spectrum * (np.exp(-mu*(r-r0))-1)*dw)/r

radius = np.linspace(device_radius.to(u.cm), 1.0001*device_radius.to(u.cm),10000)
kT_prime_vals = np.zeros(radius.shape)*(u.W/u.m**2)


for i in range(len(radius)):
    kT_prime_vals[i]= kT_prime(radius[i],radius[0],mu).to(u.W/u.m**2)

plt.plot(radius,kT_prime_vals)
plt.show()


plt.plot(radius.to_value(u.m), integrate(kT_prime_vals.to_value(u.W/u.m**2),radius.to_value(u.m)))
plt.show()

total_energy = (np.trapz(spectrum, omega) * vol).to(u.W/u.m)
print("Total Power: {:.2e} W/m".format(total_energy.value))

















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


# Shield element
element = "C"
Element_num = xrl.SymbolToAtomicNumber(element)
density = xrl.ElementDensity(Element_num)  #g/cm^3


#looking at a monoenergetic photon beam

Energy = 200 #kev
Intensity = 100 #optical power per unit area W/m^2


sigma_tot = xrl.CS_Total(Element_num,Energy) # cm^2/g

# slab
I = lambda x,x0,mu: Intensity*np.exp(-mu*(x-x0)) # x in cm
# cylinder
Ir = lambda x,x0,mu: Intensity*np.exp(-mu*(x-x0))*x0/x # x in cm

x = np.linspace(0,10,1000)+1
plt.title("intensity [W/cm^2]")
plt.loglog(x,I(x,x[0],sigma_tot*density),label = "Slab, "+element)
plt.loglog(x,Ir(x,x[0],sigma_tot*density),label = "Cylinder, "+element)

element = "Fe"
Element_num = xrl.SymbolToAtomicNumber(element)
density = xrl.ElementDensity(Element_num)  #g/cm^3
sigma_tot = xrl.CS_Total(Element_num,Energy) # cm^2/g
plt.loglog(x,I(x,x[0],sigma_tot*density),label = "Slab, "+element)
plt.loglog(x,Ir(x,x[0],sigma_tot*density),label = "Cylinder, "+element)


element = "Pb"
Element_num = xrl.SymbolToAtomicNumber(element)
density = xrl.ElementDensity(Element_num)  #g/cm^3
sigma_tot = xrl.CS_Total(Element_num,Energy) # cm^2/g
plt.loglog(x,I(x,x[0],sigma_tot*density),label = "Slab, "+element)
plt.loglog(x,Ir(x,x[0],sigma_tot*density),label = "Cylinder, "+element)

element = "U"
Element_num = xrl.SymbolToAtomicNumber(element)
density = xrl.ElementDensity(Element_num)  #g/cm^3
sigma_tot = xrl.CS_Total(Element_num,Energy) # cm^2/g
plt.loglog(x,I(x,x[0],sigma_tot*density),label = "Slab, "+element)
plt.loglog(x,Ir(x,x[0],sigma_tot*density),label = "Cylinder, "+element)

plt.legend()
plt.tight_layout()
plt.show()


























from plasmapy.formulary.radiation import thermal_bremsstrahlung

frequencies = np.arange(15, 19.9, 0.01)
frequencies = (10 ** frequencies) / u.s  #Hz

energies = (frequencies * const.h.si).to(u.eV) #eV

ne = 1e22 * u.cm ** -3
Te = 6e2*1000 * u.eV
ion_species = "B-11 5+"


spectrum = thermal_bremsstrahlung(frequencies, ne, Te, ion_species=ion_species)

spectrum = spectrum.to(u.W * u.s / u.m ** 3)

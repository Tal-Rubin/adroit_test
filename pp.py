import astropy.constants as const
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from plasmapy.formulary.radiation import thermal_bremsstrahlung

frequencies = np.arange(15, 19.9, 0.001)
frequencies = (10 ** frequencies) / u.s  #Hz

energies = (frequencies * const.h.si).to(u.keV) #keV

omega = 2*np.pi*frequencies

Te = 6e2*1000 * u.eV
ion_species_B = "B-11 5+"
ion_species_p = "H-1 1+"

nB = 3/7*1e21 * u.cm ** -3
nH = 4/7*1e21 * u.cm ** -3

ne = 5*nB+nH


spectrum_B = thermal_bremsstrahlung(frequencies, ne, Te, nB, ion_species=ion_species_B).to(u.W * u.s / u.m ** 3)
spectrum_p = thermal_bremsstrahlung(frequencies, ne, Te, nH, ion_species=ion_species_p).to(u.W * u.s / u.m ** 3)


lbl = "$T_e$ = {:.1e} eV,  ".format(Te.value) + "$n_e$ = {:.1e} 1/cm^3".format(ne.value)
plt.loglog(energies, spectrum_B, label="B")
plt.loglog(energies, spectrum_p, label="p")
plt.loglog(energies, spectrum_p+spectrum_B, label="B+p")

plt.title(f"Thermal Bremsstrahlung Spectrum"+lbl)
plt.xlabel("Energy (keV)")
plt.ylabel("Power Spectral Density (W s/m^3)")
plt.legend()
plt.show()

spectrum = spectrum_p+spectrum_B
 
vol = 0.5 * u.cm ** 3
dw = 2 * np.pi * np.gradient(frequencies)  # Frequency step size
total_energy = (np.sum(spectrum * dw) * vol).to(u.W)
print("Total Power: {:.2e} W".format(total_energy.value))
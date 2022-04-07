import astropy.constants as const
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from plasmapy.formulary.radiation import thermal_bremsstrahlung

frequencies = np.arange(15, 19.9, 0.01)
frequencies = (10 ** frequencies) / u.s

energies = (frequencies * const.h.si).to(u.eV)

ne = 1e22 * u.cm ** -3
Te = 6e2*1000 * u.eV
#ion_species = "B-11 5+"
ion_species = "C-12 4+"
spectrum = thermal_bremsstrahlung(frequencies, ne, Te, ion_species=ion_species)

print(spectrum.unit)

lbl = "$T_e$ = {:.1e} eV,\n".format(Te.value) + "$n_e$ = {:.1e} 1/cm^3".format(ne.value)
plt.plot(energies, spectrum, label=lbl)
plt.title(f"Thermal Bremsstrahlung Spectrum")
plt.xlabel("Energy (eV)")
plt.ylabel("Power Spectral Density (W s/m^3)")
plt.legend()
plt.show()

spectrum = spectrum.to(u.W * u.s / u.m ** 3)
spectrum.unit

t = 5 * u.ns
vol = 0.5 * u.cm ** 3
dw = 2 * np.pi * np.gradient(frequencies)  # Frequency step size
total_energy = (np.sum(spectrum * dw) * t * vol).to(u.J)
print("Total Energy: {:.2e} J".format(total_energy.value))


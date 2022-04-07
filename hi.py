import numpy as np
import matplotlib.pyplot as plt
import xraylib as xrl


element = "U"
Element_num = xrl.SymbolToAtomicNumber(element)


num_bins = 1000
min_E = 0
max_E = 600



bin_lims = np.linspace(min_E, max_E ,num_bins +1) #kev
bin_E = (bin_lims[1:]+bin_lims[:-1]) * 0.5 #kev


def get_sigma_tot(element, bin_lims, order=2):
    bin_sigma_tot = np.zeros(num_bins)
    if order==1:
        for i in range(num_bins):
    	    bin_sigma_tot[i] = xrl.CS_Total(element,0.5*(bin_lims[i]+bin_lims[i+1])) # cm^2/g
    elif order==2:
        for i in range(num_bins):
            loc = np.array([0.5*(bin_lims[i]+bin_lims[i+1]) + 0.5*(bin_lims[i+1]-bin_lims[i])*3**-0.5,
						0.5*(bin_lims[i]+bin_lims[i+1]) - 0.5*(bin_lims[i+1]-bin_lims[i])*3**-0.5])
            bin_sigma_tot[i] = 0.5*xrl.CS_Total(element,loc[0]) + 0.5*xrl.CS_Total(element,loc[1])# cm^2/g


    return bin_sigma_tot

bin_sigma_tot = get_sigma_tot(Element_num, bin_lims, 1)
density = xrl.ElementDensity(Element_num)  #g/cm^3
atomic_weight = xrl.AtomicWeight(Element_num)  #g/mol

plt.figure()
plt.plot(bin_E,bin_sigma_tot*density)
plt.xlabel("E [kev]")
plt.ylabel("$\\sigma \\rho [cm^-1]$")
plt.tight_layout()

plt.figure()
plt.semilogy(bin_E,bin_sigma_tot*density*bin_E)
plt.xlabel("E [kev]")
plt.ylabel("$\\sigma \\rho E [kev cm^-1]$")
plt.tight_layout()


plt.show()

print("element "+element)
print("Total energy deposited {:.2e}".format(sum(bin_sigma_tot*density*bin_E)*(bin_lims[1]-bin_lims[0])))

print("density {} g/cm^3".format(density))



















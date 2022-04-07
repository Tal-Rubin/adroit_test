import numpy as np
import matplotlib.pyplot as plt
import xraylib as xrl

num_bins = 200
min_E = 0
max_E = 600



bin_lims = np.linspace(min_E, max_E ,num_bins +1) #kev
bin_E = (bin_lims[1:]+bin_lims[:-1]) * 0.5


def get_sigma_tot(element, bin_lims, order):
    bin_sigma_tot = np.zeros(num_bins)
    for i in range(num_bins):
    	bin_sigma_tot[i] = xrl.CS_Total(element,0.5*(bin_lims[i]+bin_lims[i+1]))
    return bin_sigma_tot

bin_sigma_tot = get_sigma_tot(6, bin_lims, 1)



plt.figure()
plt.plot(bin_E,bin_sigma_tot)


plt.figure()
plt.plot(bin_E,bin_sigma_tot*bin_E)

plt.show()


print(sum(bin_sigma_tot*bin_E)*(bin_lims[1]-bin_lims[0]))

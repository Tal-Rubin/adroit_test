import numpy as np
import xraylib as xrl

num_bins = 30

bin_lims = np.linspace(0,300,num_bins +1) #kev
bins = (bin_lims[1:]+bin_lims[:-1]) * 0.5

bin_sigma_tot = np.zeros(bins.shape)

for b,bin in bin_sigma_tot, bins:
	b = xrl.CS_Total(92,bin)
print(bin_sigma_tot)
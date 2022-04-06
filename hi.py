import numpy as np
import xraylib as xrl

num_bins = 30
min_E = 0
max_E = 600



bin_lims = np.linspace(min_E, max_E ,num_bins +1) #kev
bin_E = (bin_lims[1:]+bin_lims[:-1]) * 0.5

bin_sigma_tot = np.zeros(num_bins)

for i in range(num_bins):
	bin_sigma_tot[i] = xrl.CS_Total(92,bin_E[i])

print(bin_sigma_tot)
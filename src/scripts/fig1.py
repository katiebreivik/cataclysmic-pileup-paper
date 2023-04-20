import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import paths

max_distance = 1000.0
dat = pd.read_csv(paths.data / "dat_maxDistance_1000.txt")
dat = dat.rename(columns={'# m1[Msun]':'m1',
                          ' m2[Msun]':'m2',
                          ' inclination[rad]':'inc',
                          ' f_gw[Hz]':'f_gw',
                          ' x_gal[kpc]':'x', 
                          ' y_gal[kpc]':'y', 
                          ' z_gal[kpc]':'z',
                          ' Pala_reassigned' : 'pala'})

PT = pd.read_hdf(paths.data / 'Pala_2020_dat_combo.h5')


plt.figure(figsize=(6, 4))
plt.hist(1/(dat.f_gw/2) / 3600, bins=200, density=True, histtype='step', label=f'{max_distance} pc sample')
plt.hist(PT.porb/60, density=True, histtype='step', bins=50, label='Pala+2020 histogram')
plt.xlabel('porb [hr]', size=12)
plt.ylabel('count density', size=12)
plt.legend(prop={'size':12})
plt.tick_params('both', labelsize=10)
plt.xlim(1, 5.5)
plt.tight_layout()
plt.savefig(paths.figures / "fig1.pdf")



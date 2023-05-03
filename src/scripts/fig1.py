import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import paths

dark_green = '#264653'
teal = '#2a9d8f'
yellow =  '#e9c46a'
orange = '#f4a261'
salmon = '#e76f51'


# set the default font and fontsize
plt.rc('font', family='serif')
plt.rcParams['text.usetex'] = False
fs = 12

# update various fontsizes to match
params = {'figure.figsize': (6, 4),
          'legend.fontsize': fs,
          'axes.labelsize': fs,
          'xtick.labelsize': 0.7 * fs,
          'ytick.labelsize': 0.7 * fs}
plt.rcParams.update(params)

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
plt.hist(PT.porb/60, density=True, histtype='step', bins=50, label='Pala+2020', lw=1.5, color=orange)
plt.hist(1/(dat.f_gw/2) / 3600, bins=200, density=True, histtype='step', label=f'{int(max_distance/1000)} kpc', color=teal, lw=1.5)
plt.xlabel('porb [hr]', size=12)
plt.ylabel('count density', size=12)
plt.legend(prop={'size':12})
plt.tick_params('both', labelsize=10)
plt.xlim(1, 5.5)
plt.tight_layout()
plt.savefig(paths.figures / "fig1.pdf")



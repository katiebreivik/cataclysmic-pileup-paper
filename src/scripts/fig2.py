import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import paths
import astropy.units as u
from legwork import strain, psd, utils
from astropy.coordinates import SkyCoord
from CV_pop_create import sample_kpc_population
import seaborn as sns

# set the default font and fontsize
plt.rc('font', family='serif')
plt.rcParams['text.usetex'] = False
fs = 12
dark_green = '#264653'
teal = '#2a9d8f'
yellow =  '#e9c46a'
orange = '#f4a261'
salmon = '#e76f51'

# update various fontsizes to match
params = {'figure.figsize': (6, 4),
          'legend.fontsize': fs,
          'axes.labelsize': fs,
          'xtick.labelsize': 0.7 * fs,
          'ytick.labelsize': 0.7 * fs}
plt.rcParams.update(params)

n_realizations = 1

## THESE ARE FIXED FOR THIS STUDY
max_distance = 1000.0 # u.kpc
mu_m1 = 0.7
sigma_m1 = 0.001
sigma_m2 = 0.001

frequency_range=np.logspace(-5, 0, 1000) * u.Hz
LISA = psd.lisa_psd(frequency_range, approximate_R=False, confusion_noise=None)

h_f = []
m_c = []
f_gw = []

h_f_local = []
m_c_local = []
f_gw_local = []
for ii in range(n_realizations):
    dat = sample_kpc_population(max_distance, mu_m1, sigma_m1, sigma_m2)
    local_mask = dat[:,7] == 2

    mc = mc = utils.chirp_mass(dat[:,0] * u.Msun, dat[:,1] * u.Msun)
    c = SkyCoord(dat[:,4], dat[:,5], dat[:,6], unit=u.kpc, frame='galactocentric')
    dist = c.icrs.distance

    ASD = strain.h_0_n(m_c=mc, f_orb=dat[:,2]/2 * u.Hz,  
                       ecc=np.zeros(len(mc)), dist=dist, 
                       n=2, position=None, polarisation=None, 
                       inclination=None, interpolated_g=None) * np.sqrt(4 * 3.155e7)

    h_f.extend(ASD)
    m_c.extend(mc.value)
    f_gw.extend(dat[:,2])
    
    h_f_local.extend(ASD[local_mask])
    m_c_local.extend(mc.value[local_mask])
    f_gw_local.extend(dat[local_mask,2])


fig = plt.figure(figsize=(6, 4))

#sns.kdeplot(x=np.log10(f_gw), y=np.log10(np.array(h_f).flatten()), levels=5, bw_adjust=1.5, color="blue", fill=False, alpha=1, label='1 kpc, 50 realizations')
#sns.kdeplot(x=np.log10(f_gw_local), y=np.log10(np.array(h_f_local).flatten()), levels=5, bw_adjust=1.5, color="grey", fill=False, alpha=1, label='150 pc, 50 realizations')
local_mask = dat[:,7] == 2
Pala_mask = dat[:,7] == 1

plt.scatter(np.log10(dat[Pala_mask, 2]), np.log10(ASD[Pala_mask]), label='Pala+2020', s=40, edgecolors='black', linewidths=0.75, c=teal)
plt.scatter(np.log10(f_gw), np.log10(h_f), label='1 kpc sample', s=10, alpha=0.5, zorder=0, c=orange)
plt.plot(np.log10(frequency_range.value), np.log10(LISA.value**0.5), lw=2, c='black', ls='--', label='instrument noise')   
plt.xlim(-4.5, -2.5)
plt.ylim(-20, -16)
plt.xlabel('GW frequency [Hz]', size=12)
plt.ylabel('ASD [Hz$^{-1/2}$]', size=12)
plt.tick_params('both', labelsize=10)
plt.legend(prop={'size':12})
plt.tight_layout()

plt.savefig(paths.figures / "fig2.pdf")

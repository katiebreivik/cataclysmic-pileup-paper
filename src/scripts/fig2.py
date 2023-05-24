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
dat = pd.read_csv(paths.data / "dat_maxDistance_1000_final.txt")
dat = dat.rename(columns={'# m1[Msun]':'m1',
                          ' m2[Msun]':'m2',
                          ' inclination[rad]':'inc',
                          ' f_gw[Hz]':'f_gw',
                          ' x_gal[kpc]':'x', 
                          ' y_gal[kpc]':'y', 
                          ' z_gal[kpc]':'z',
                          ' Pala_reassigned' : 'pala'})

mc = utils.chirp_mass(dat.m1.values * u.Msun, dat.m2.values * u.Msun)
c = SkyCoord(dat.x.values, dat.y.values, dat.z.values, unit=u.kpc, frame='galactocentric')
dist = c.icrs.distance
ASD = strain.h_0_n(m_c=mc, f_orb=dat.f_gw.values/2 * u.Hz,  
                   ecc=np.zeros(len(mc)), dist=dist, 
                   n=2, position=None, polarisation=None, 
                   inclination=None, interpolated_g=None) * np.sqrt(4 * 3.155e7)

fig = plt.figure(figsize=(6, 4))

dat_Pala = dat.loc[dat.pala == 1]
dat_150 = dat.loc[dat.pala == 2]
mc_Pala = utils.chirp_mass(dat_Pala.m1.values * u.Msun, dat_Pala.m2.values * u.Msun)
c_Pala = SkyCoord(dat_Pala.x.values, dat_Pala.y.values, dat_Pala.z.values, unit=u.kpc, frame='galactocentric')
dist_Pala = c_Pala.icrs.distance
ASD_Pala = strain.h_0_n(m_c=mc_Pala, f_orb=dat_Pala.f_gw.values/2 * u.Hz,  
                        ecc=np.zeros(len(mc_Pala)), dist=dist_Pala, 
                        n=2, position=None, polarisation=None, 
                        inclination=None, interpolated_g=None) * np.sqrt(4 * 3.155e7)

mc_150 = utils.chirp_mass(dat_150.m1.values * u.Msun, dat_150.m2.values * u.Msun)
c_150 = SkyCoord(dat_150.x.values, dat_150.y.values, dat_150.z.values, unit=u.kpc, frame='galactocentric')
dist_150 = c_150.icrs.distance
ASD_150 = strain.h_0_n(m_c=mc_150, f_orb=dat_150.f_gw.values/2 * u.Hz,
                        ecc=np.zeros(len(mc_150)), dist=dist_150,
                        n=2, position=None, polarisation=None,
                        inclination=None, interpolated_g=None) * np.sqrt(4 * 3.155e7)


plt.scatter(dat_Pala.f_gw, ASD_Pala.flatten(), label='Pala+2020', s=40, edgecolors='black', linewidths=0.75, c=orange, zorder=3, rasterized=True)
plt.scatter(dat_150.f_gw, ASD_150.flatten(), label='150 pc sample', s=30, c=orange, zorder=1, rasterized=True)
plt.scatter(dat.f_gw, ASD.flatten(), label='1 kpc sample', s=10, alpha=0.5, zorder=0, c=teal, rasterized=True)
plt.plot(frequency_range.value, LISA.value**0.5, lw=2, c='black', ls='--', label='instrument noise')   
plt.xlim(10**(-4.5), 10**(-2.5))
plt.ylim(10**(-20), 10**(-16))
plt.xlabel('GW frequency [Hz]', size=12)
plt.ylabel('ASD [Hz$^{-1/2}$]', size=12)
plt.tick_params('both', labelsize=10)
plt.xscale('log')
plt.yscale('log')
plt.legend(prop={'size':12})
plt.tight_layout()

plt.savefig(paths.figures / "fig2.pdf")

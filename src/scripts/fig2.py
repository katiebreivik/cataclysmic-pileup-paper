import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import paths
import scipy
import astropy.units as u
from legwork import strain, psd, utils
from astropy.coordinates import SkyCoord


dat = pd.read_csv(paths.data / "dat_maxDistance_1000.txt")
dat = dat.rename(columns={'# m1[Msun]':'m1',
                          ' m2[Msun]':'m2',
                          ' inclination[rad]':'inc',
                          ' f_gw[Hz]':'f_gw',
                          ' x_gal[kpc]':'x', 
                          ' y_gal[kpc]':'y', 
                          ' z_gal[kpc]':'z',
                          ' Pala_reassigned' : 'pala'})

mc = utils.chirp_mass(dat.m1.values * u.Msun, dat.m2.values * u.Msun)
c = SkyCoord(dat.x, dat.y, dat.z, unit=u.kpc, frame='galactocentric')
dist = c.icrs.distance

ASD = strain.h_0_n(m_c=mc, f_orb=dat['f_gw'].values/2 * u.Hz,  
                   ecc=np.zeros(len(mc)), dist=dist, 
                   n=2, position=None, polarisation=None, 
                   inclination=None, interpolated_g=None) * np.sqrt(4 * 3.155e7)
frequency_range=np.logspace(-5, 0, 1000) * u.Hz
LISA = psd.lisa_psd(frequency_range, approximate_R=False, confusion_noise='robson19')
ind, = np.where(dat.pala == 1)
indD, = np.where(dat.pala==2)


fig = plt.figure(figsize=(8, 6))

# Compute the 2D histogram
hist, x_edges, y_edges = np.histogram2d(np.squeeze(dat.f_gw), np.squeeze(ASD), bins=50)

# Define the contour levels
levels = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
          
# Smooth the data using a Gaussian filter
sigma = 2  # Adjust this parameter to control the smoothing strength
smoothed_hist = scipy.ndimage.gaussian_filter(hist, sigma)

# Plot the contour lines
plt.contour(x_edges[:-1], y_edges[:-1], smoothed_hist.T, levels=levels)


# Plot the rest
plt.plot(frequency_range, LISA**0.5, lw=1, c='black')   
plt.scatter(dat.f_gw, ASD, s=5, alpha= 0.5, label='1000 pc sample')
plt.scatter(dat.loc[dat.pala == 2].f_gw, ASD[indD], s=20, label='150pc')
plt.scatter(dat.loc[dat.pala == 1].f_gw, ASD[ind], s=20, label='Pala+2020')
plt.legend(prop={'size':12})
plt.yscale('log')
plt.xscale('log')
plt.xlim(9e-5, 1e-3)
plt.ylim(1e-18, 3e-17)
plt.xlabel('GW frequency [Hz]', size=12)
plt.ylabel('ASD [Hz$^{-1/2}$]', size=12)
plt.tick_params('both', labelsize=10)
plt.tight_layout()

plt.savefig(paths.figures / "fig2.pdf")

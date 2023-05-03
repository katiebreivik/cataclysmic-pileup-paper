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


# Generate CDF of the orbital period distribution

porbGrid = np.linspace(60/60, 11, 100000) # hours
x = np.log10(porbGrid*60) # log10 of porb in minutes
C1 = 1.885
C2 = 2.269
x1 = x - C1
x2 = x - C2

y1 = 0.8344
cut1 = 76.78
cut2 = 2.10*60
cut3 = 3.09*60
cut4 = 10.58*60


cdf1 = x*0
cdf2 = 1.98 * ( (x1 + x1**(3/4)) / (1 + x1 + x1**2) )
cdf3 = y1
cdf4 = 0.26 * (x2 + x2**(3/4)) / (1 + x2 + x2**2) + y1
cdf5 = x*0 + 1

# Generate the model CDF
modelCDF = np.zeros(len(porbGrid))
ix = porbGrid < (cut1/60)
modelCDF[ix] = cdf1[ix]

ix = (porbGrid >= (cut1/60)) & (porbGrid <= (cut2/60))
modelCDF[ix] = cdf2[ix]

ix = (porbGrid > (cut2/60)) & (porbGrid < (cut3/60))
modelCDF[ix] = cdf3

ix = porbGrid >= (cut3/60)
modelCDF[ix] = cdf4[ix]

ix = porbGrid >= (cut4/60)
modelCDF[ix] = cdf5[ix]



fig = plt.figure(figsize=(6, 4))
ax1 = fig.add_subplot(111)
ax1.hist(PT.porb/60, density=True, histtype='step', bins=50, label='Pala+2020', lw=1.5, color=orange)
ax1.hist(1/(dat.f_gw/2) / 3600, bins=200, density=True, histtype='step', label=f'{int(max_distance/1000)} kpc', color=teal, lw=1.5)
ax1.set_xlabel('porb [hr]', size=12)
ax1.set_ylabel('normalized counts', size=12)
ax1.tick_params('both', labelsize=10)
ax1.set_xlim(1, 11)

# generate CDF of PT.porb
palaCDF = np.cumsum(1.0 / len(PT.porb) * np.arange(len(PT.porb)))

#palaCDF = np.cumsum(PT.porb/60) / np.sum(PT.porb/60) 
# Normalize the CDF
palaCDF = palaCDF / np.max(palaCDF)
#plot PT.porb vs palaCDF on right y-axis using dotted line
ax2 = plt.gca().twinx()
PT = PT.sort_values(by='porb')
ax2.step(PT.porb/60, palaCDF, color=orange, linestyle='dotted', label='Pala+2020 CDF', lw=1.5)
# make staircase plot


# plot modelCDF on right y-axis using solid line
ax2.plot(porbGrid, modelCDF, color='black', linestyle='dotted', label='Eq.5 CDF', lw=1.5)


ax2.set_ylabel('CDF', size=12)
# set ax2 ylim between 0 and 1
ax2.set_ylim(0, 1.01)
ax2.tick_params('both', labelsize=10)
#plt.legend(loc='center right', prop={'size':10})

lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, prop={'size':12}, loc='center right')

#ax2.legend(prop={'size':12}, loc='center right')


plt.tight_layout()
plt.savefig(paths.figures / "fig1.pdf")



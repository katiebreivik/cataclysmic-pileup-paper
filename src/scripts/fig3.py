import numpy as np
import matplotlib.pyplot as plt
import paths 

fig = plt.figure(figsize=(8, 6))
                   
data=np.loadtxt(paths.data / 'lisa_data_w_cvs.dat')
residual=np.loadtxt(paths.data / 'lisa_residual_w_cvs.dat')
cv_confusion=np.loadtxt(paths.data / 'psd_w_cvs.dat')
wdwd_confusion=np.loadtxt(paths.data / 'psd_no_cvs.dat')
noise=np.loadtxt(paths.data / 'noise.dat')

plt.plot(data[:,0],np.sqrt(data[:,1]*data[:,1]+data[:,2]*data[:,2]),label='data',color="gray")   
plt.plot(residual[:,0],np.sqrt(residual[:,1]*residual[:,1]+residual[:,2]*residual[:,2]),label='residual',color="lightgray")   
plt.plot(cv_confusion[:,0],np.sqrt(cv_confusion[:,1]),label='cv confusion',color="tab:blue", linewidth=3)   
plt.plot(wdwd_confusion[:,0],np.sqrt(wdwd_confusion[:,1]),label='wdwd confusion',color="tab:orange", linewidth=3)   
plt.plot(noise[:,0],np.sqrt((noise[:,3]-noise[:,4])/2),label='instrument noise',color="black", linestyle='dashed')   

plt.legend(prop={'size':12}, loc='upper left')
plt.yscale('log')
plt.xscale('log')
plt.xlim(1e-4, 4e-3)
plt.ylim(1e-24, 1e-18)

plt.xlabel('GW frequency [Hz]', size=12)
plt.ylabel('ASD [Hz$^{-1/2}$]', size=12)
plt.tick_params('both', labelsize=10)
plt.tight_layout()

plt.savefig(paths.figures / "fig3.pdf")

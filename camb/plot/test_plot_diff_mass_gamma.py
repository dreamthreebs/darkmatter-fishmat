import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc # latex in matplotlib
rc('text',usetex=True)

def log_Cubic_interpolate(xx,yy):
    from scipy.interpolate import CubicSpline
    logx=np.log10(xx)
    logy=np.log10(yy)
    Cubic_interp=CubicSpline(logx,logy)
    log_Cubic_interp=lambda zz: np.power(10, Cubic_interp(np.log10(zz)))
    return log_Cubic_interp

mass_len=50
mass_start=1.01e-5
mass_end=5e3
dm_mass_set=np.geomspace(mass_start,mass_end,mass_len)

sig_gamma_cvl=np.load('../data/sig_gamma_data/sig_gamma_diff_mass.npy')

xs=np.geomspace(mass_start,mass_end,1000)

cs=log_Cubic_interpolate(dm_mass_set, sig_gamma_cvl[:,3])
plt.loglog(xs,cs(xs),color='pink',linestyle='-',linewidth=1)
cs=log_Cubic_interpolate(dm_mass_set, sig_gamma_cvl[:,2])
plt.loglog(xs,cs(xs),color='blue',linestyle='-',linewidth=1)
cs=log_Cubic_interpolate(dm_mass_set, sig_gamma_cvl[:,1])
plt.loglog(xs,cs(xs),color='yellow',linestyle='-',linewidth=1)
cs=log_Cubic_interpolate(dm_mass_set, sig_gamma_cvl[:,0])
plt.loglog(xs,cs(xs),color='purple',linestyle='-',linewidth=1)



# plt.legend(['test'],loc='upper right')

# plt.title('ali l from 30 to 620, pico l from 10 to 1000')
plt.xlabel(r'$m_\chi\left[GeV\right]$')
plt.ylabel(r'$\left\langle \sigma v\right\rangle/ m_\chi \left[cm^3s^{-1}GeV^{-1}\right] $')
# plt.ylim((1e-28,1e-26))
plt.xlim((1e-5,5e3))
plt.title(r'constraints derived by fisher matrix')
# plt.savefig('./figure/constraints_diff_mass_without_planck/1.png',dpi=300)
plt.show()



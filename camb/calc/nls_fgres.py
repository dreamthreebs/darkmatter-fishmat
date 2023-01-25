import numpy as np
import os
from consts import *

def noise_dls2cls(dls,ells,noise_lmin,noise_lmax):
    cls=np.zeros((ells,2))
    for i in np.arange(2):
        for l in np.arange(noise_lmin,noise_lmax):
            cls[l,i]=dls[l,i]*(2*np.pi)/(l*(l+1))
    return cls

def fgres_dls2cls(dls,ells,fgres_lmin,fgres_lmax):
    cls=np.zeros((ells,3))
    for i in np.arange(3):
        for l in np.arange(fgres_lmin,fgres_lmax):
            cls[l,i]=dls[l,i]*(2*np.pi)/(l*(l+1))
    return cls

def ali_noise_level(ells): # remember to set noise to (ells,2) matrix
    ali_noise_dls=np.loadtxt('./input/ali/Ali_noise.dat',usecols=(1,2))
    insert_first_thirty_rows=np.zeros((30,2))
    half_Noise=np.insert(ali_noise_dls,0,insert_first_thirty_rows,axis=0)
    insert_last_many_rows=np.zeros((ells-len(half_Noise),2))
    noisedls=np.insert(half_Noise,len(half_Noise),insert_last_many_rows,axis=0)
    ali_noise=noise_dls2cls(noisedls,ells,30,621)
    return ali_noise

def ali_fg_res(ells):
    # from matplotlib import pyplot as plt
    TTFgres=np.load('./input/ali/TT_fg.npy')
    EEFgres=np.load('./input/ali/EE_fg.npy')
    TEFgres=np.load('./input/ali/TE_fg.npy')
    half_res=np.stack((TTFgres,EEFgres,TEFgres),axis=1)
    insert_last_many_rows=np.zeros((ells-len(half_res),3))
    fgresdls=np.insert(half_res,len(half_res),insert_last_many_rows,axis=0)
    ali_fgres=fgres_dls2cls(fgresdls,ells,1,len(TTFgres))
    return ali_fgres

def pico_noise_level(ells):
    pico_noise_TT=np.load('./input/pico/Nl_TT_results.npy')
    pico_noise_EE=np.load('./input/pico/Nl_EE_results.npy')
    half_noise=np.stack((pico_noise_TT,pico_noise_EE),axis=1)
    insert_last_many_rows=np.zeros((ells-len(half_noise),2))
    noisedls=np.insert(half_noise,len(half_noise),insert_last_many_rows,axis=0)
    pico_noise=noise_dls2cls(noisedls, ells, 1, len(pico_noise_TT))
    return pico_noise

def pico_fg_res(ells):
    TTFgres=np.load('./input/pico/fg_TT_results.npy')
    EEFgres=np.load('./input/pico/fg_EE_results.npy')
    TEFgres=np.load('./input/pico/fg_TE_results.npy')
    half_res=np.stack((TTFgres,EEFgres,TEFgres),axis=1)
    insert_last_many_rows=np.zeros((ells-len(half_res),3))
    fgresdls=np.insert(half_res,len(half_res),insert_last_many_rows,axis=0)
    pico_fgres=fgres_dls2cls(fgresdls,ells,1,len(TTFgres))
    return pico_fgres

def check_nls_and_fgres():
    from matplotlib import pyplot as plt
    ls=np.arange(ells)
    plt.figure(1)
    plt.loglog(ls,ls*ls*ali_noise_level(ells)[:,0])
    plt.loglog(ls,ls*ls*pico_noise_level(ells)[:,0])

    plt.figure(2)
    plt.loglog(ls,ls*ls*ali_noise_level(ells)[:,1])
    plt.loglog(ls,ls*ls*pico_noise_level(ells)[:,1])

    plt.figure(3)
    plt.semilogx(ls,ls*ls*ali_fg_res(ells)[:,0])
    plt.semilogx(ls,ls*ls*pico_fg_res(ells)[:,0])

    plt.figure(4)
    plt.semilogx(ls,ls*ls*ali_fg_res(ells)[:,1])
    plt.semilogx(ls,ls*ls*pico_fg_res(ells)[:,1])

    plt.figure(5)
    plt.semilogx(ls,ls*ls*ali_fg_res(ells)[:,2])
    plt.semilogx(ls,ls*ls*pico_fg_res(ells)[:,2])
    plt.show()

def log_Cubic_interpolate(xx,yy):
    from scipy.interpolate import CubicSpline
    logx=np.log10(xx)
    logy=np.log10(yy)
    Cubic_interp=CubicSpline(logx,logy)
    log_Cubic_interp=lambda zz: np.power(10, Cubic_interp(np.log10(zz)))
    return log_Cubic_interp

def zero_noise(ells):
    return np.zeros((ells,2))

def zero_fgres(ells):
    return np.zeros((ells,3))

if __name__=="__main__":
    os.chdir("../")
    print(os.getcwd())

    # check_nls_and_fgres()
    nls=zero_noise(ells)
    fgres=zero_fgres(ells)
    print(f"{nls}\n {nls.shape}")
    print(f"{fgres}\n {fgres.shape}")



import numpy as np
from matplotlib import pyplot as plt

# ls=np.arange(30,621)
# a=np.load('./TT_fg.npy')
# print(a.shape)

ls=np.arange(2501)

# # ali_TT_fgres=np.load('./TT_fg.npy')
# pico_TT_fgres=np.load('./fg_TT_results.npy')
# plt.figure(1)
# # plt.loglog(ls,ali_TT_fgres,color='g')
# plt.loglog(ls,pico_TT_fgres,color='k')
# # ali_EE_fgres=np.load('EE_fg.npy')
# pico_EE_fgres=np.load('./fg_EE_results.npy')
# plt.figure(2)
# # plt.loglog(ls,ali_EE_fgres,color='g')
# plt.loglog(ls,pico_EE_fgres,color='k')
# # ali_TE_fgres=np.load('TE_fg.npy')
# pico_TE_fgres=np.load('fg_TE_results.npy')
# plt.figure(3)
# # plt.loglog(ls,ali_TE_fgres,color='g')
# plt.loglog(ls,pico_TE_fgres,color='k')
# plt.show()

def pico_fg_res():
    # from matplotlib import pyplot as plt
    TTFgres=np.load('fg_TT_results.npy')
    TEFgres=np.load('fg_TE_results.npy')
    EEFgres=np.load('fg_EE_results.npy')
    half_res=np.stack((TTFgres,EEFgres,TEFgres),axis=1)
    insert_last_many_rows=np.zeros((2501-len(half_res),3))
    fg_resdls=np.insert(half_res,len(half_res),insert_last_many_rows,axis=0)
    fg_res=fg_res_dls2cls(fg_resdls)
    return fg_res

def fg_res_dls2cls(dls):
    cls=np.zeros((2501,3))
    for i in np.arange(3):
        for l in np.arange(1,1024):
            cls[l,i]=dls[l,i]*(2*np.pi)/(l*(l+1))
    return cls

def Ali_Noise_level(): # remember to set noise to (2501,2) matrix
    ali_noise_dls=np.loadtxt('Ali_noise.dat',usecols=(1,2))
    print(ali_noise_dls)
    insert_first_thirty_rows=np.zeros((30,2))
    half_Noise=np.insert(ali_noise_dls,0,insert_first_thirty_rows,axis=0)
    insert_last_many_rows=np.zeros((2501-len(half_Noise),2))
    Noised=np.insert(half_Noise,len(half_Noise),insert_last_many_rows,axis=0)
    Noise=noise_dls2cls(Noised)
    return Noise

def noise_dls2cls(dls):
    cls=np.zeros((2501,2))
    for i in np.arange(2):
        for l in np.arange(30,620):
            cls[l,i]=dls[l,i]*(2*np.pi)/(l*(l+1))
    return cls

def pico_noise_level():
    pico_noise_TT=np.load('./Nl_TT_results.npy')
    pico_noise_EE=np.load('./Nl_EE_results.npy')
    half_noise=np.stack((pico_noise_TT,pico_noise_EE),axis=1)
    insert_last_many_rows=np.zeros((2501-len(half_noise),2))
    noise_dls=np.insert(half_noise,len(half_noise),insert_last_many_rows,axis=0)
    pico_noise=noise_dls2cls(noise_dls)
    return pico_noise

Noisecls=Ali_Noise_level()
Noisecls1=pico_noise_level()
# fgrescls=pico_fg_res()

# plt.figure(1)
# plt.loglog(ls,Noisecls[:,0],color='b')
# plt.loglog(ls,Noisecls1[:,0],color='k')
# plt.figure(2)
# plt.loglog(ls,Noisecls[:,1],color='b')
# plt.loglog(ls,Noisecls1[:,1],color='k')
# plt.show()






# Test on foreground residual 
import numpy as np
from matplotlib import pyplot as plt

TTFgres=np.load('TT_fg.npy')
ls=np.arange(30,621)
plt.figure(1)
plt.loglog(ls,TTFgres[30:621])
TEFgres=np.load('TE_fg.npy')
plt.figure(2)
plt.loglog(ls,TEFgres[30:621])
EEFgres=np.load('EE_fg.npy')
plt.figure(3)
plt.loglog(ls,EEFgres[30:621])
plt.show()

# ali_noise_dls=np.loadtxt('Ali_noise.dat',usecols=(1,2))
# print(ali_noise_dls[0:30,0])
# print(ali_noise_dls[0:30,1])

def ali_fg_res():
    from matplotlib import pyplot as plt
    TTFgres=np.load('TT_fg.npy')
    TEFgres=np.load('TE_fg.npy')
    EEFgres=np.load('EE_fg.npy')
    # half_res=np.stack(TTFgres,TEFgres,axis=1)
    half_res=np.stack((TTFgres,TEFgres,EEFgres),axis=1)
    insert_last_many_rows=np.zeros((2501-len(half_res),3))
    fg_resdls=np.insert(half_res,len(half_res),insert_last_many_rows,axis=0)
    fg_res=fg_res_dls2cls(fg_resdls)
    print(fg_res.shape)
    return fg_res
    ls=np.arange(30,621)

    # plt.figure(1)
    # plt.loglog(ls,fg_res[30:621,0])
    # plt.figure(2)
    # plt.loglog(ls,fg_res[30:621,1])
    # plt.figure(3)
    # plt.loglog(ls,fg_res[30:621,2])
    # plt.show()

    # ls=np.arange(0,len(fg_res[:,0]))
    # plt.figure(1)
    # plt.loglog(ls,fg_res[:,0])
    # plt.figure(2)
    # plt.loglog(ls,fg_res[:,1])
    # plt.figure(3)
    # plt.loglog(ls,fg_res[:,2])
    # plt.show()

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

def fg_res_dls2cls(dls):
    cls=np.zeros((2501,3))
    for i in np.arange(3):
        for l in np.arange(1,1024):
            cls[l,i]=dls[l,i]*(2*np.pi)/(l*(l+1))
    return cls

def main():
    fg_res=ali_fg_res()

if __name__=='__main__':
    main()

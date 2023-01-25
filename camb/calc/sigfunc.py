import numpy as np
import os
from consts import *
from pdfunc import *
from nls_fgres import *
import datetime

def load_all_pd():
    DM_Pann_CLprime=np.load('./data/pd_data/DM_Pann_CLprime.npy')
    DM_Gamma_CLprime=np.load('./data/pd_data/DM_Gamma_CLprime.npy')
    thetastarmc_CLprime=np.load('./data/pd_data/thetastarmc_CLprime.npy')
    ombh2_CLprime=np.load('./data/pd_data/ombh2_CLprime.npy')
    omch2_CLprime=np.load('./data/pd_data/omch2_CLprime.npy')
    As_CLprime=np.load('./data/pd_data/As_CLprime.npy')
    ns_CLprime=np.load('./data/pd_data/ns_CLprime.npy')
    optical_depth_CLprime=np.load('./data/pd_data/optical_depth_CLprime.npy')
    return DM_Pann_CLprime,DM_Gamma_CLprime,thetastarmc_CLprime,ombh2_CLprime,omch2_CLprime,As_CLprime,ns_CLprime,optical_depth_CLprime


def load_basic_6params_pd():
    thetastarmc_CLprime=np.load('./data/pd_data/thetastarmc_CLprime.npy')
    ombh2_CLprime=np.load('./data/pd_data/ombh2_CLprime.npy')
    omch2_CLprime=np.load('./data/pd_data/omch2_CLprime.npy')
    As_CLprime=np.load('./data/pd_data/As_CLprime.npy')
    ns_CLprime=np.load('./data/pd_data/ns_CLprime.npy')
    optical_depth_CLprime=np.load('./data/pd_data/optical_depth_CLprime.npy')
    return thetastarmc_CLprime,ombh2_CLprime,omch2_CLprime,As_CLprime,ns_CLprime,optical_depth_CLprime

def check_all_pd():
    from matplotlib import pyplot as plt
    DM_Pann_CLprime,DM_Gamma_CLprime,thetastarmc_CLprime,ombh2_CLprime,omch2_CLprime,As_CLprime,ns_CLprime,optical_depth_CLprime=load_all_pd()
    ls=np.arange(ells)
    list=[DM_Pann_CLprime,DM_Gamma_CLprime,thetastarmc_CLprime,ombh2_CLprime,omch2_CLprime,As_CLprime,ns_CLprime,optical_depth_CLprime]
    fig, axs = plt.subplots(8,3)
    for i in np.arange(3):
        for index,pd in enumerate(list):
            axs[index,i].semilogx(ls,ls*ls*pd[:,i])
            if i==0:
                axs[index,0].set_title(['pann','gamma','thetastarmc','ombh2','omch2','As','ns','tau'][index],{'fontsize': 10})
    plt.show()

def get_TT_fisher_matrix(pd,cls,nls,fgres,lmin,lmax,fsky):
    FM_TT=np.zeros((params_num,params_num))
    Cell=np.zeros((1,1))
    Cellprime_i=np.zeros((1,1))
    Cellprime_j=np.zeros((1,1))
    for i in np.arange(params_num):
        for j in np.arange(params_num):
            for l in np.arange(lmin,lmax):
                Cell[0,0]=cls[l,0]+nls[l,0]+fgres[l,0]
                Cell_inv=np.linalg.inv(Cell)
                Cellprime_i[0,0]=pd[l,0,i]
                Cellprime_j[0,0]=pd[l,0,j]
                Mul=np.matmul(Cell_inv,Cellprime_i)
                Mult=np.matmul(Mul,Cell_inv)
                Multi=np.matmul(Mult,Cellprime_j)
                FM_TT[i,j]+=(2*l+1)*np.trace(Multi)/2
   # print(FM_TT)
    FM_TT*=fsky
    FI=np.linalg.inv(FM_TT)
    sigma=np.zeros((params_num,1))
    #get covariance
    for i in np.arange(params_num):
        sigma[i]=np.sqrt(FI[i,i])
        # print(sigma[i])
    return sigma

def get_EE_fisher_matrix(pd,cls,nls,fgres,lmin,lmax,fsky):
    FM_EE=np.zeros((params_num,params_num))
    Cell=np.zeros((1,1))
    Cellprime_i=np.zeros((1,1))
    Cellprime_j=np.zeros((1,1))
    for i in np.arange(params_num):
        for j in np.arange(params_num):
            for l in np.arange(lmin,lmax):
                Cell[0,0]=cls[l,1]+nls[l,1]+fgres[l,1]
                Cell_inv=np.linalg.inv(Cell)
                Cellprime_i[0,0]=pd[l,1,i]
                Cellprime_j[0,0]=pd[l,1,j]
                Mul=np.matmul(Cell_inv,Cellprime_i)
                Mult=np.matmul(Mul,Cell_inv)
                Multi=np.matmul(Mult,Cellprime_j)
                FM_EE[i,j]+=(2*l+1)*np.trace(Multi)/2
            # if i==j==3:
                # FM_EE+=1/(0.013**2)
   # print(FM_EE)
    FM_EE*=fsky
    FI=np.linalg.inv(FM_EE)
    sigma=np.zeros((params_num,1))
    #get covariance
    for i in np.arange(params_num):
        sigma[i]=np.sqrt(FI[i,i])
        # print(sigma[i])
    return sigma

def get_TE_fisher_matrix(pd,cls,nls,fgres,lmin,lmax,fsky):
    FM_TE=np.zeros((params_num,params_num))
    Cell=np.zeros((1,1))
    Cellprime_i=np.zeros((1,1))
    Cellprime_j=np.zeros((1,1))
    for i in np.arange(params_num):
        for j in np.arange(params_num):
            for l in np.arange(lmin,lmax):
                Cell[0,0]=np.sqrt((cls[l,2]+fgres[l,2])**2+(cls[l,0]+nls[l,0]+fgres[l,0])*(cls[l,1]+nls[l,1]+fgres[l,1]))
                Cell_inv=np.linalg.inv(Cell)
                Cellprime_i[0,0]=pd[l,2,i]
                Cellprime_j[0,0]=pd[l,2,j]
                Mul=np.matmul(Cell_inv,Cellprime_i)
                Mult=np.matmul(Mul,Cell_inv)
                Multi=np.matmul(Mult,Cellprime_j)
                FM_TE[i,j]+=(2*l+1)*np.trace(Multi)
            # if i==j==3:
                # FM_TE[i,j]+=1/(0.013**2)
   # print(FM_TE)
    FM_TE*=fsky
    FI=np.linalg.inv(FM_TE)
    sigma=np.zeros((params_num,1))
    #get covariance
    for i in np.arange(params_num):
        sigma[i]=np.sqrt(FI[i,i])
        # print(sigma[i])
    return sigma

def get_combined_fisher_matrix(pd,cls,nls,fgres,lmin,lmax,fsky):
    FM=np.zeros((params_num,params_num))
    Cell=np.zeros((2,2))
    Cellprime_i=np.zeros((2,2))
    Cellprime_j=np.zeros((2,2))
    for i in np.arange(params_num):
        for j in np.arange(params_num):
            for l in np.arange(lmin,lmax):
                Cell[0,0]=cls[l,0]+nls[l,0]+fgres[l,0]
                Cell[1,0]=cls[l,2]+fgres[l,2]
                Cell[0,1]=cls[l,2]+fgres[l,2]
                Cell[1,1]=cls[l,1]+nls[l,1]+fgres[l,1]
    #            Cell[2,2]=totCL[l,2]
                Cell_inv=np.linalg.inv(Cell)
                Cellprime_i[0,0]=pd[l,0,i]
                Cellprime_i[1,0]=pd[l,2,i]
                Cellprime_i[0,1]=pd[l,2,i]
                Cellprime_i[1,1]=pd[l,1,i]
    #            Cellprime_i[2,2]=pd[l,2,i]
                Cellprime_j[0,0]=pd[l,0,j]
                Cellprime_j[1,0]=pd[l,2,j]
                Cellprime_j[0,1]=pd[l,2,j]
                Cellprime_j[1,1]=pd[l,1,j]
    #            Cellprime_j[2,2]=pd[l,2,j]
                Mul=np.matmul(Cell_inv,Cellprime_i)
                Mult=np.matmul(Mul,Cell_inv)
                Multi=np.matmul(Mult,Cellprime_j)
                FM[i,j]+=(2*l+1)*np.trace(Multi)/2
   # print(FM)
    FM*=fsky
    FI=np.linalg.inv(FM)
    sigma=np.zeros((params_num,1))
    #get covariance
    for i in np.arange(params_num):
        sigma[i]=np.sqrt(FI[i,i])
        # print(sigma[i])
    return sigma

def check_fisher_pann():
    DM_Pann_CLprime,DM_Gamma_CLprime,thetastarmc_CLprime,ombh2_CLprime,omch2_CLprime,As_CLprime,ns_CLprime,optical_depth_CLprime=load_all_pd()
    pd_Pann=np.zeros((ells,3,params_num))
    pd_Pann[:,:,0]=ombh2_CLprime
    pd_Pann[:,:,1]=omch2_CLprime
    pd_Pann[:,:,2]=thetastarmc_CLprime
    pd_Pann[:,:,3]=optical_depth_CLprime
    pd_Pann[:,:,4]=ns_CLprime
    pd_Pann[:,:,5]=As_CLprime
    pd_Pann[:,:,6]=DM_Pann_CLprime # choose one of clprime
    cls=initial_totCL()
    zero_nls=np.zeros((ells,2))
    zero_fgres=np.zeros((ells,3))
    lmin=10
    lmax=ells
    fsky=0.37
    sig_pann_TT=get_TT_fisher_matrix(pd_Pann, cls, zero_nls, zero_fgres, lmin, lmax, fsky)
    sig_pann_EE=get_EE_fisher_matrix(pd_Pann, cls, zero_nls, zero_fgres, lmin, lmax, fsky)
    sig_pann_TE=get_TE_fisher_matrix(pd_Pann, cls, zero_nls, zero_fgres, lmin, lmax, fsky)
    sig_pann_combined=get_combined_fisher_matrix(pd_Pann, cls, zero_nls, zero_fgres, lmin, lmax, fsky)
    print(f"sig_pann_TT is :{sig_pann_TT}")
    print(f"sig_pann_EE is :{sig_pann_EE}")
    print(f"sig_pann_TE is :{sig_pann_TE}")
    print(f"sig_pann_combined is :{sig_pann_combined}")

def check_fisher_gamma():
    DM_Pann_CLprime,DM_Gamma_CLprime,thetastarmc_CLprime,ombh2_CLprime,omch2_CLprime,As_CLprime,ns_CLprime,optical_depth_CLprime=load_all_pd()
    pd_Gamma=np.zeros((ells,3,params_num))
    pd_Gamma[:,:,0]=ombh2_CLprime
    pd_Gamma[:,:,1]=omch2_CLprime
    pd_Gamma[:,:,2]=thetastarmc_CLprime
    pd_Gamma[:,:,3]=optical_depth_CLprime
    pd_Gamma[:,:,4]=ns_CLprime
    pd_Gamma[:,:,5]=As_CLprime
    pd_Gamma[:,:,6]=DM_Gamma_CLprime
    cls=initial_totCL()
    zero_nls=np.zeros((ells,2))
    zero_fgres=np.zeros((ells,3))
    lmin=10
    lmax=ells
    fsky=0.37
    sig_gamma_TT=get_TT_fisher_matrix(pd_Gamma, cls, zero_nls, zero_fgres, lmin, lmax, fsky)
    sig_gamma_EE=get_EE_fisher_matrix(pd_Gamma, cls, zero_nls, zero_fgres, lmin, lmax, fsky)
    sig_gamma_TE=get_TE_fisher_matrix(pd_Gamma, cls, zero_nls, zero_fgres, lmin, lmax, fsky)
    sig_gamma_combined=get_combined_fisher_matrix(pd_Gamma, cls, zero_nls, zero_fgres, lmin, lmax, fsky)
    print(f"sig_gamma_TT is :{sig_gamma_TT}")
    print(f"sig_gamma_EE is :{sig_gamma_EE}")
    print(f"sig_gamma_TE is :{sig_gamma_TE}")
    print(f"sig_gamma_combined is :{sig_gamma_combined}")

def calc_pann_two_sigma_error(cls,nls,fgres,ells,lmin,lmax,fsky,filename):
    # from matplotlib import pyplot as plt
    thetastarmc_CLprime,ombh2_CLprime,omch2_CLprime,As_CLprime,ns_CLprime,optical_depth_CLprime=load_basic_6params_pd()
    pd_Pann=np.zeros((ells,3,params_num))
    pd_Pann[:,:,0]=ombh2_CLprime
    pd_Pann[:,:,1]=omch2_CLprime
    pd_Pann[:,:,2]=thetastarmc_CLprime
    pd_Pann[:,:,3]=optical_depth_CLprime
    pd_Pann[:,:,4]=ns_CLprime
    pd_Pann[:,:,5]=As_CLprime

    mass_len=50
    mass_start=1e-5
    mass_end=5e3
    dm_mass_set=np.geomspace(mass_start,mass_end,mass_len)

    sig_Pann=np.zeros((mass_len,4))
    for i,DM_mass in enumerate(dm_mass_set):
        print('this is loop:', i)
        xcordinate,dp,ls,min_step_mat,DM_Pann_CLprime=DM_Pann_prime(DM_mass, ells, length=30, start_footstep=1e-26, end_footstep=1e-31)
        pd_Pann[:,:,6]=DM_Pann_CLprime
        sig_Pann[i,0]=2*get_TT_fisher_matrix(pd_Pann, cls, nls, fgres, lmin, lmax, fsky)[params_num-1]
        sig_Pann[i,1]=2*get_EE_fisher_matrix(pd_Pann, cls, nls, fgres, lmin, lmax, fsky)[params_num-1]
        sig_Pann[i,2]=2*get_TE_fisher_matrix(pd_Pann, cls, nls, fgres, lmin, lmax, fsky)[params_num-1]
        sig_Pann[i,3]=2*get_combined_fisher_matrix(pd_Pann, cls, nls, fgres, lmin, lmax, fsky)[params_num-1]
    np.save(filename,sig_Pann)

def test_plot():
    from matplotlib import pyplot as plt
    xs=np.geomspace(1e-5,5e3,1000)
    cs=log_Cubic_interpolate(dm_mass_set, sig_Pann[:,0])
    plt.loglog(xs,cs(xs),color='b')
    plt.xlabel('dark matter mass')
    plt.ylabel('energy injection')
    cs=log_Cubic_interpolate(dm_mass_set, sig_Pann[:,1])
    plt.loglog(xs,cs(xs),color='g')
    cs=log_Cubic_interpolate(dm_mass_set, sig_Pann[:,2])
    plt.loglog(xs,cs(xs),color='y')
    cs=log_Cubic_interpolate(dm_mass_set, sig_Pann[:,3])
    plt.loglog(xs,cs(xs),color='r')
    plt.legend(['TT','EE','TE','TT+EE+TE'],loc='upper right')
    plt.ylim((1e-28,1e-25))
    plt.xlim((1e-5,5e3))
    plt.title('CVL condition lmax=4000')
    plt.show()

def calc_pann_two_sigma_error_with_prior_pd_pann(cls,nls,fgres,ells,lmin,lmax,fsky,filename,pd_pann_diff_mass):
    # from matplotlib import pyplot as plt
    thetastarmc_CLprime,ombh2_CLprime,omch2_CLprime,As_CLprime,ns_CLprime,optical_depth_CLprime=load_basic_6params_pd()
    pd_Pann=np.zeros((ells,3,params_num))
    pd_Pann[:,:,0]=ombh2_CLprime
    pd_Pann[:,:,1]=omch2_CLprime
    pd_Pann[:,:,2]=thetastarmc_CLprime
    pd_Pann[:,:,3]=optical_depth_CLprime
    pd_Pann[:,:,4]=ns_CLprime
    pd_Pann[:,:,5]=As_CLprime

    mass_len=50
    mass_start=1e-5
    mass_end=5e3
    dm_mass_set=np.geomspace(mass_start,mass_end,mass_len)

    sig_Pann=np.zeros((mass_len,4))
    for i,DM_mass in enumerate(dm_mass_set):
        print('this is loop:', i)
        pd_Pann[:,:,6]=pd_pann_diff_mass[:,:,i]
        sig_Pann[i,0]=2*get_TT_fisher_matrix(pd_Pann, cls, nls, fgres, lmin, lmax, fsky)[params_num-1]
        sig_Pann[i,1]=2*get_EE_fisher_matrix(pd_Pann, cls, nls, fgres, lmin, lmax, fsky)[params_num-1]
        sig_Pann[i,2]=2*get_TE_fisher_matrix(pd_Pann, cls, nls, fgres, lmin, lmax, fsky)[params_num-1]
        sig_Pann[i,3]=2*get_combined_fisher_matrix(pd_Pann, cls, nls, fgres, lmin, lmax, fsky)[params_num-1]
    np.save(filename,sig_Pann)



if __name__=="__main__":
    print(os.getcwd())
    os.chdir("../")
    print(os.getcwd())

    # set_all_params()
    # check_all_pd()

    # check_fisher_pann()
    # check_fisher_gamma()

    # calc_pann_two_sigma_error(cls=initial_totCL(), nls=zero_noise(ells), fgres=zero_fgres(ells), ells=ells, lmin=10, lmax=ells, fsky=1,filename='../data/sigma_data/CVL_fsky1_lmin10_lmaxells.npy')

    # pd_pann_diff_mass=np.load('./data/pd_pann_50_data/pd_pann_50_diff_mass.npy')
    # calc_pann_two_sigma_error_with_prior_pd_pann(cls=initial_totCL(), nls=zero_noise(ells), fgres=zero_fgres(ells), ells=ells, lmin=10, lmax=4000, fsky=1, filename='./calc/test_sig.npy', pd_pann_diff_mass=pd_pann_diff_mass)
    # print(datetime.datetime.now())

    mass_len=50
    mass_start=1e-5
    mass_end=5e3
    DM_mass_set=np.geomspace(mass_start,mass_end,mass_len)


    sig_Pann=np.load('./calc/test_sig.npy')
    from matplotlib import pyplot as plt
    xs=np.geomspace(1e-5,5e3,1000)
    cs=log_Cubic_interpolate(DM_mass_set, sig_Pann[:,0])
    plt.loglog(xs,cs(xs),color='b')
    plt.xlabel('dark matter mass')
    plt.ylabel('energy injection')
    cs=log_Cubic_interpolate(DM_mass_set, sig_Pann[:,1])
    plt.loglog(xs,cs(xs),color='g')
    cs=log_Cubic_interpolate(DM_mass_set, sig_Pann[:,2])
    plt.loglog(xs,cs(xs),color='y')
    cs=log_Cubic_interpolate(DM_mass_set, sig_Pann[:,3])
    plt.loglog(xs,cs(xs),color='r')
    plt.legend(['TT','EE','TE','TT+EE+TE'],loc='upper right')
    plt.ylim((1e-28,1e-25))
    plt.xlim((1e-5,5e3))
    plt.title('CVL condition lmax=4000')
    plt.show()






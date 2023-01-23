import numpy as np
import os
from consts import *
from pdfunc import *

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


def check_all_pd():
    from matplotlib import pyplot as plt
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


if __name__=="__main__":
    print(os.getcwd())
    os.chdir("../")
    print(os.getcwd())


    DM_Pann_CLprime,DM_Gamma_CLprime,thetastarmc_CLprime,ombh2_CLprime,omch2_CLprime,As_CLprime,ns_CLprime,optical_depth_CLprime=load_all_pd()

    # check_all_pd()

    pd_Pann=np.zeros((ells,3,params_num))
    pd_Pann[:,:,0]=ombh2_CLprime
    pd_Pann[:,:,1]=omch2_CLprime
    pd_Pann[:,:,2]=thetastarmc_CLprime
    pd_Pann[:,:,3]=optical_depth_CLprime
    pd_Pann[:,:,4]=ns_CLprime
    pd_Pann[:,:,5]=As_CLprime
    pd_Pann[:,:,6]=DM_Pann_CLprime # choose one of clprime

    # pd_Gamma=np.zeros((ells,3,params_num))
    # pd_Gamma[:,:,0]=ombh2_CLprime
    # pd_Gamma[:,:,1]=omch2_CLprime
    # pd_Gamma[:,:,2]=thetastarmc_CLprime
    # pd_Gamma[:,:,3]=optical_depth_CLprime
    # pd_Gamma[:,:,4]=ns_CLprime
    # pd_Gamma[:,:,5]=As_CLprime
    # pd_Gamma[:,:,6]=DM_Gamma_CLprime

    cls=initial_totCL()
    zero_nls=np.zeros((ells,2))
    zero_fgres=np.zeros((ells,3))
    lmin=10
    lmax=620
    fsky=0.37
    sigma_TT=get_TT_fisher_matrix(pd_Pann, cls, zero_nls, zero_fgres, lmin, lmax, fsky)
    sigma_EE=get_EE_fisher_matrix(pd_Pann, cls, zero_nls, zero_fgres, lmin, lmax, fsky)
    sigma_TE=get_TE_fisher_matrix(pd_Pann, cls, zero_nls, zero_fgres, lmin, lmax, fsky)
    sigma_combined=get_combined_fisher_matrix(pd_Pann, cls, zero_nls, zero_fgres, lmin, lmax, fsky)
    print(f"sigma_TT is :{sigma_TT}")
    print(f"sigma_EE is :{sigma_EE}")
    print(f"sigma_TE is :{sigma_TE}")
    print(f"sigma_combined is :{sigma_combined}")








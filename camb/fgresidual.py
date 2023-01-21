import fileinput
import numpy as np
import re
import os
import sys
import time
from matplotlib import pyplot as plt
from matplotlib import rc # latex in matplotlib
rc('text',usetex=True)

#initial value of six fundamental parameters
ombh2_0 = 0.02242
omch2_0 = 0.11933
As_0 = 2.105209331337507e-09 #scalar_amp(1)
ns_0 = 0.9665 #scalar_spectral_index(1)
optical_depth_0 = 0.0561 # re_optical_depth
hubble_0 = 67.66
thetastarmc_0=1.040997
DM_mass_0=1e2
DM_Pann_0=0
DM_Gamma_0=0
params_num=7
ells=4001
ls=np.arange(ells)

def initial_totCL():
    os.system('./camb test.ini')
    test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
    insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
    dlsn=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
    clsn=dls2cls(dlsn)
    return clsn


def set_l_max_scalar(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('l_max_scalar'):
            print('l_max_scalar = '+str(new_value))
        else:
            print(line.strip())

def get_thetastarmc_from_files():
    filename="output.txt"
    with fileinput.FileInput(filename) as f:
        for line in f :
            if line.startswith('100 theta (CosmoMC)'):
                thetastarmc=re.findall(r"\d+\.?\d*",line[line.find("="):])[0]
                return float(thetastarmc)

def dls2cls(dls):
    cls=np.zeros((ells,3))
    for i in np.arange(3):
        for l in np.arange(2,ells):
            cls[l,i]=dls[l,i]*(2*np.pi)/(l*(l+1))
    return cls

def set_hubble(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('hubble'):
            print('hubble = '+str(new_value))
        else:
            print(line.strip())

def set_ombh2(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('ombh2'):
            print('ombh2 = '+str(new_value))
        else:
            print(line.strip())

def set_omch2(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('omch2'):
            print('omch2 = '+str(new_value))
        else:
            print(line.strip())

def set_As(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('scalar_amp(1)'):
            print('scalar_amp(1) = '+str(new_value))
        else:
            print(line.strip())

def set_ns(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('scalar_spectral_index(1)'):
            print('scalar_spectral_index(1) = '+str(new_value))
        else:
            print(line.strip())

def set_optical_depth(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('re_optical_depth'):
            print('re_optical_depth = '+str(new_value))
        else:
            print(line.strip())

def set_DM_Pann(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('DM_Pann'):
            print('DM_Pann = '+str(new_value))
        else:
            print(line.strip())

def set_DM_Gamma(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('DM_Gamma'):
            print('DM_Gamma = '+str(new_value))
        else:
            print(line.strip())

def set_DM_mass(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('DM_mass'):
            print('DM_mass = '+str(new_value))
        else:
            print(line.strip())

print('thetastarmc_0= ',thetastarmc_0)
#set initial value of fundamental Model
set_l_max_scalar(ells-1)
set_ns(ns_0)
set_ombh2(ombh2_0)
set_omch2(omch2_0)
set_optical_depth(optical_depth_0)
set_As(As_0)
set_hubble(hubble_0)
set_DM_mass(DM_mass_0)
set_DM_Pann(DM_Pann_0)
set_DM_Gamma(DM_Gamma_0)


clsn=initial_totCL()

def DM_Pann_prime(params_value):
    set_DM_mass(params_value)
    set_DM_Gamma(0)
    set_DM_Pann(0)
    DM_Pann_0=0
    length=30
    start_footstep=1e-26
    end_footstep=1e-31
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    index=0
    dp=np.zeros((ells,3,length))
    difference=np.ones((ells,3,length))
    DM_Pann_CLprime=np.zeros((ells,3))
    ls=np.arange(ells)
    x1=DM_Pann_0
    set_DM_Pann(x1)
    os.system('./camb test.ini')
    test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
    insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
    dlsn=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
    clsn=dls2cls(dlsn)
    
    for h in xcordinate:
        x2=DM_Pann_0+h
        set_DM_Pann(x2)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        dlsp=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)# the output test_scalCls is Dl(i.e. l(l+1)C_l)
        clsp=dls2cls(dlsp)
        dp[:,:,index]=(clsp-clsn)/h
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,ells):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        print('loop=',index)
        index+=1
    '''
        show dp(derivative against params) at difference footstep to find stable point
    '''
    # take footstep at smallest difference as DM_Pann_CLprime
    sum=0
    min_step_mat=np.zeros((ells,3))
    delta_dls=np.zeros((ells,3))
    for i in np.arange(3):
        for l in np.arange(2,ells):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l,i]=minimum_step[0][0]
            DM_Pann_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/((ells-2)*3))
    print(best_h)
    set_DM_Pann(xcordinate[best_h])
    os.system('./camb test.ini')
    test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))
    insert_first_two_rows=np.array([[0,0,0],[0,0,0]])
    dlsd=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
    delta_dls=dlsd-dlsn
    DM_Pann_CLprime=dp[:,:,best_h-1]
    set_DM_Pann(0)
    return xcordinate,dp,ls,min_step_mat,DM_Pann_CLprime

# xcordinate,dp,ls,min_step_mat,DM_Pann_CLprime=DM_Pann_prime(DM_mass_0)

def DM_Gamma_prime(params_value):
    set_DM_Pann(0)
    set_DM_mass(params_value)
    set_DM_Gamma(0)
    DM_Gamma_0=0
    length=30
    start_footstep=1e-24
    end_footstep=1e-31
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    index=0
    dp=np.zeros((ells,3,length))
    difference=np.ones((ells,3,length))
    DM_Gamma_CLprime=np.zeros((ells,3))
    ls=np.arange(ells)
    x1=DM_Gamma_0
    set_DM_Gamma(x1)
    os.system('./camb test.ini')
    test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
    insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
    dlsn=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
    clsn=dls2cls(dlsn)
    
    for h in xcordinate:
        x2=DM_Gamma_0+h
        set_DM_Gamma(x2)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        dlsp=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        clsp=dls2cls(dlsp)
        dp[:,:,index]=(clsp-clsn)/h
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,ells):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        print('loop=',index)
        index+=1
    '''
        show dp(derivative against params) at different footstep to find stable point
    '''
    # take footstep at smallest difference as DM_Gamma_CLprime
    sum=0
    min_step_mat=np.zeros((ells,3))
    delta_dls=np.zeros((ells,3))
    for i in np.arange(3):
        for l in np.arange(2,ells):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l][i]=minimum_step[0][0]
            DM_Gamma_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/((ells-2)*3))
    print(best_h)
    set_DM_Gamma(xcordinate[best_h])
    os.system('./camb test.ini')
    test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))
    insert_first_two_rows=np.array([[0,0,0],[0,0,0]])
    dlsd=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
    delta_dls=dlsd-dlsn
    DM_Gamma_CLprime=dp[:,:,best_h]
    set_DM_Gamma(0)
    return xcordinate,dp,ls,min_step_mat,DM_Gamma_CLprime

# xcordinate,dp,ls,min_step_mat,DM_Gamma_CLprime=DM_Gamma_prime(DM_mass_0)

def thetastarmc_prime():
    ls=np.arange(ells)
    length=30
    start_footstep=2e-1
    end_footstep=1e-4
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    thetastar_tmp=np.zeros((length*2,1))
    hubble_tmp=np.zeros((length*2,1))
    dp=np.zeros((ells,3,length))
    difference=np.ones((ells,3,length))
    thetastarmc_CLprime=np.zeros((ells,3))
    index=0
    for h in xcordinate:
        x1=hubble_0-h*hubble_0
        x2=hubble_0+h*hubble_0
        hubble_tmp[index]=x1
        hubble_tmp[2*length-index-1]=x2
        set_hubble(x1)
        os.system('./camb test.ini | grep "100 theta" >output.txt')
        thetastar_tmp[index]=get_thetastarmc_from_files()
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_n=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_n=dls2cls(DL_n)
        set_hubble(x2)
        os.system('./camb test.ini | grep "100 theta" >output.txt')
        thetastar_tmp[2*length-index-1]=get_thetastarmc_from_files()
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_p=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_p=dls2cls(DL_p)
        dp[:,:,index]=(CL_p-CL_n)/(thetastar_tmp[2*length-index-1]-thetastar_tmp[index])
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,ells):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        index+=1
        print(index)
    sum=0
    min_step_mat=np.zeros((ells,3))
    delta_dls=np.zeros((ells,3))
    for i in np.arange(3):
        for l in np.arange(2,ells):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l][i]=minimum_step[0][0]
            thetastarmc_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/((ells-2)*3))
    print(best_h) # for this parameter, 10 as best_h is better, around 1.45e-2 
    thetastarmc_CLprime=dp[:,:,10]
    set_hubble(hubble_0)
    return thetastarmc_CLprime

# thetastarmc_CLprime=thetastarmc_prime()

#derivative against ombh2
def ombh2_prime():
    ls=np.arange(ells)
    length=30
    start_footstep=1e-1
    end_footstep=1e-7
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    dp=np.zeros((ells,3,length))
    difference=np.ones((ells,3,length))
    ombh2_CLprime=np.zeros((ells,3))
    index=0
    for h in xcordinate:
        x1=ombh2_0-h*ombh2_0
        x2=ombh2_0+h*ombh2_0
        set_ombh2(x1)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_n=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_n=dls2cls(DL_n)
        set_ombh2(x2)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_p=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_p=dls2cls(DL_p)
        dp[:,:,index]=(CL_p-CL_n)/(2*h*ombh2_0)
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,ells):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        index+=1
        print(index)
    sum=0
    min_step_mat=np.zeros((ells,3))
    delta_dls=np.zeros((ells,3))
    for i in np.arange(3):
        for l in np.arange(2,ells):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l][i]=minimum_step[0][0]
            ombh2_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/((ells-2)*3))
    print(best_h) # for this parameter, 10 as best_h is better, around 1.45e-2 
    ombh2_CLprime=dp[:,:,best_h]
    set_ombh2(ombh2_0)
    return ombh2_CLprime

# ombh2_CLprime=ombh2_prime()

def omch2_prime():
    ls=np.arange(ells)
    length=30
    start_footstep=1e-1
    end_footstep=1e-7
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    dp=np.zeros((ells,3,length))
    difference=np.ones((ells,3,length))
    omch2_CLprime=np.zeros((ells,3))
    index=0
    for h in xcordinate:
        x1=omch2_0-h*omch2_0
        x2=omch2_0+h*omch2_0
        set_omch2(x1)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_n=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_n=dls2cls(DL_n)
        set_omch2(x2)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_p=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_p=dls2cls(DL_p)
        dp[:,:,index]=(CL_p-CL_n)/(2*h*omch2_0)
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,ells):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        index+=1
        print(index)
    sum=0
    min_step_mat=np.zeros((ells,3))
    delta_dls=np.zeros((ells,3))
    for i in np.arange(3):
        for l in np.arange(2,ells):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l][i]=minimum_step[0][0]
            omch2_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/((ells-2)*3))
    print(best_h) # for this parameter, 10 as best_h is better, around 1.45e-2 
    omch2_CLprime=dp[:,:,best_h]
    set_omch2(omch2_0)
    return omch2_CLprime

# omch2_CLprime=omch2_prime()

def optical_depth_prime():
    ls=np.arange(ells)
    length=30
    start_footstep=1e-1
    end_footstep=1e-4
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    dp=np.zeros((ells,3,length))
    difference=np.ones((ells,3,length))
    optical_depth_CLprime=np.zeros((ells,3))
    index=0
    for h in xcordinate:
        x1=optical_depth_0-h*optical_depth_0
        x2=optical_depth_0+h*optical_depth_0
        set_optical_depth(x1)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_n=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_n=dls2cls(DL_n)
        set_optical_depth(x2)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_p=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_p=dls2cls(DL_p)
        dp[:,:,index]=(CL_p-CL_n)/(2*h*optical_depth_0)
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,ells):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        index+=1
        print(index)
    sum=0
    min_step_mat=np.zeros((ells,3))
    delta_dls=np.zeros((ells,3))
    for i in np.arange(3):
        for l in np.arange(2,ells):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l][i]=minimum_step[0][0]
            optical_depth_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/((ells-2)*3))
    print(best_h) # for this parameter, 10 as best_h is better, around 1.45e-2 
    optical_depth_CLprime=dp[:,:,best_h]
    set_optical_depth(optical_depth_0)
    return optical_depth_CLprime

# optical_depth_CLprime=optical_depth_prime()

def ns_prime():
    ls=np.arange(ells)
    length=30
    start_footstep=1e-1
    end_footstep=1e-7
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    dp=np.zeros((ells,3,length))
    difference=np.ones((ells,3,length))
    ns_CLprime=np.zeros((ells,3))
    index=0
    for h in xcordinate:
        x1=ns_0-h*ns_0
        x2=ns_0+h*ns_0
        set_ns(x1)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_n=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_n=dls2cls(DL_n)
        set_ns(x2)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_p=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_p=dls2cls(DL_p)
        dp[:,:,index]=(CL_p-CL_n)/(2*h*ns_0)
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,ells):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        index+=1
        print(index)
    sum=0
    min_step_mat=np.zeros((ells,3))
    delta_dls=np.zeros((ells,3))
    for i in np.arange(3):
        for l in np.arange(2,ells):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l][i]=minimum_step[0][0]
            ns_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/((ells-2)*3))
    print(best_h) # for this parameter, 10 as best_h is better, around 1.45e-2 
    ns_CLprime=dp[:,:,best_h]
    set_ns(ns_0)
    return ns_CLprime

# ns_CLprime=ns_prime()

def As_prime():
    ls=np.arange(ells)
    length=30
    start_footstep=1e-1
    end_footstep=1e-7
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    dp=np.zeros((ells,3,length))
    difference=np.ones((ells,3,length))
    As_CLprime=np.zeros((ells,3))
    index=0
    for h in xcordinate:
        x1=1e-10*np.exp((1-h)*np.log(1e10*As_0))
        x2=1e-10*np.exp((1+h)*np.log(1e10*As_0))
        set_As(x1)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_n=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_n=dls2cls(DL_n)
        set_As(x2)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_p=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_p=dls2cls(DL_p)
        dp[:,:,index]=(CL_p-CL_n)/(np.log(1e10*x2)-np.log(1e10*x1))
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,ells):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        index+=1
        print(index)
    sum=0
    min_step_mat=np.zeros((ells,3))
    delta_dls=np.zeros((ells,3))
    for i in np.arange(3):
        for l in np.arange(2,ells):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l][i]=minimum_step[0][0]
            As_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/((ells-2)*3))
    print(best_h) # for this parameter, 10 as best_h is better, around 1.45e-2 
    As_CLprime=dp[:,:,best_h]
    set_As(As_0)
    return As_CLprime

# As_CLprime=As_prime()

def DM_mass_prime():
    ls=np.arange(ells)
    length=30
    start_footstep=9e-1
    end_footstep=1e-2
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    dp=np.zeros((ells,3,length))
    difference=np.ones((ells,3,length))
    DM_mass_CLprime=np.zeros((ells,3))
    index=0
    for h in xcordinate:
        x1=(1-h)*DM_mass_0
        x2=(1+h)*DM_mass_0
        set_DM_mass(x1)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_n=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_n=dls2cls(DL_n)
        set_DM_mass(x2)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_p=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_p=dls2cls(DL_p)
        dp[:,:,index]=(CL_p-CL_n)/(2*h*DM_mass_0)
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,ells):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        index+=1
        print(index)
    sum=0
    min_step_mat=np.zeros((ells,3))
    delta_dls=np.zeros((ells,3))
    for i in np.arange(3):
        for l in np.arange(2,ells):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l][i]=minimum_step[0][0]
            DM_mass_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/((ells-2)*3))
    print(best_h) # for this parameter, 10 as best_h is better, around 1.45e-2 
    DM_mass_CLprime=dp[:,:,best_h]
    set_DM_mass(DM_mass_0)
    return dp,ls,xcordinate,min_step_mat,DM_mass_CLprime

# dp,ls,xcordinate,min_step_mat,DM_mass_CLprime=DM_mass_prime()

def get_TT_fisher_matrix(Derivative,lmin,lmax):
    FM_TT=np.zeros((params_num,params_num))
    Cell=np.zeros((1,1))
    Cellprime_i=np.zeros((1,1))
    Cellprime_j=np.zeros((1,1))
    for i in np.arange(params_num):
        for j in np.arange(params_num):
            for l in np.arange(lmin,lmax):
                Cell[0,0]=clsn[l,0]+Noisecls[l,0]+fgrescls[l,0]
                Cell_inv=np.linalg.inv(Cell)
                Cellprime_i[0,0]=Derivative[l,0,i]
                Cellprime_j[0,0]=Derivative[l,0,j]
                Mul=np.matmul(Cell_inv,Cellprime_i)
                Mult=np.matmul(Mul,Cell_inv)
                Multi=np.matmul(Mult,Cellprime_j)
                FM_TT[i,j]+=(2*l+1)*np.trace(Multi)/2
   #             if l==2016:
   #                 print(i,j)
            # if i==j==3:
                # FM_TT[i,j]+=1/(0.013**2)
   # print(FM_TT)
    FM_TT*=fsky
    FI=np.linalg.inv(FM_TT)
    sigma=np.zeros((params_num,1))
    #get covariance
    for i in np.arange(params_num):
        sigma[i]=np.sqrt(FI[i,i])
        # print(sigma[i])
    return sigma

def get_EE_fisher_matrix(Derivative,lmin,lmax):
    FM_EE=np.zeros((params_num,params_num))
    Cell=np.zeros((1,1))
    Cellprime_i=np.zeros((1,1))
    Cellprime_j=np.zeros((1,1))
    for i in np.arange(params_num):
        for j in np.arange(params_num):
            for l in np.arange(lmin,lmax):
                Cell[0,0]=clsn[l,1]+Noisecls[l,1]+fgrescls[l,1]
                Cell_inv=np.linalg.inv(Cell)
                Cellprime_i[0,0]=Derivative[l,1,i]
                Cellprime_j[0,0]=Derivative[l,1,j]
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

def get_TE_fisher_matrix(Derivative,lmin,lmax):
    FM_TE=np.zeros((params_num,params_num))
    Cell=np.zeros((1,1))
    Cellprime_i=np.zeros((1,1))
    Cellprime_j=np.zeros((1,1))
    for i in np.arange(params_num):
        for j in np.arange(params_num):
            for l in np.arange(lmin,lmax):
                Cell[0,0]=np.sqrt((clsn[l,2]+fgrescls[l,2])**2+(clsn[l,0]+Noisecls[l,0]+fgrescls[l,0])*(clsn[l,1]+Noisecls[l,1]+fgrescls[l,1]))
                Cell_inv=np.linalg.inv(Cell)
                Cellprime_i[0,0]=Derivative[l,2,i]
                Cellprime_j[0,0]=Derivative[l,2,j]
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

def get_combined_fisher_matrix(Derivative,lmin,lmax):
    FM=np.zeros((params_num,params_num))
    Cell=np.zeros((2,2))
    Cellprime_i=np.zeros((2,2))
    Cellprime_j=np.zeros((2,2))
    for i in np.arange(params_num):
        for j in np.arange(params_num):
            for l in np.arange(lmin,lmax):
                Cell[0,0]=clsn[l,0]+Noisecls[l,0]+fgrescls[l,0]
                Cell[1,0]=clsn[l,2]+fgrescls[l,2]
                Cell[0,1]=clsn[l,2]+fgrescls[l,2]
                Cell[1,1]=clsn[l,1]+Noisecls[l,1]+fgrescls[l,1]
    #            Cell[2,2]=totCL[l,2]
                Cell_inv=np.linalg.inv(Cell)
                Cellprime_i[0,0]=Derivative[l,0,i]
                Cellprime_i[1,0]=Derivative[l,2,i]
                Cellprime_i[0,1]=Derivative[l,2,i]
                Cellprime_i[1,1]=Derivative[l,1,i]
    #            Cellprime_i[2,2]=Derivative[l,2,i]
                Cellprime_j[0,0]=Derivative[l,0,j]
                Cellprime_j[1,0]=Derivative[l,2,j]
                Cellprime_j[0,1]=Derivative[l,2,j]
                Cellprime_j[1,1]=Derivative[l,1,j]
    #            Cellprime_j[2,2]=Derivative[l,2,j]
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

def cangjs_Noise_level():
    ls=np.arange(ells)
    flag='all'
    Noise=np.zeros((ells,3))
    CMB_temperature=2.725#K
    sensitivity_tem_95=1.45#muK  arcmin
    sensitivity_tem_150=1.45#muK arcmin
    sensitivity_pol_95=2.06#muK  arcmin
    sensitivity_pol_150=2.06#muK arcmin
    thetaFWHM_95=15.37*0.000291#rad
    thetaFWHM_150=9.73*0.000291#rad
    omega_minus_one_95_tem=(sensitivity_tem_95*0.000291)**2
    omega_minus_one_95_pol=(sensitivity_pol_95*0.000291)**2
    omega_minus_one_150_tem=(sensitivity_tem_150*0.000291)**2
    omega_minus_one_150_pol=(sensitivity_pol_150*0.000291)**2
    for i in np.arange(4):
        for l in np.arange(ells):
            if flag=='95Hz':
                if i==0:
                    Noise[l,i]=omega_minus_one_95_tem*np.exp((l*(l+1)*(thetaFWHM_95)**2)/(8*np.log(2)))
                if i==1:
                    Noise[l,i]=omega_minus_one_95_pol*np.exp((l*(l+1)*(thetaFWHM_95)**2)/(8*np.log(2)))
            if flag=='150Hz':
                if i==0:
                    Noise[l,i]=omega_minus_one_150_tem*np.exp((l*(l+1)*(thetaFWHM_150)**2)/(8*np.log(2)))
                if i==1:
                    Noise[l,i]=omega_minus_one_150_pol*np.exp((l*(l+1)*(thetaFWHM_150)**2)/(8*np.log(2)))
            if flag=='all':
                if i==0:
                    Noise[l,i]=1/(1/(omega_minus_one_95_tem*np.exp((l*(l+1)*(thetaFWHM_95)**2)/(8*np.log(2))))+1/(omega_minus_one_150_tem*np.exp((l*(l+1)*(thetaFWHM_150)**2)/(8*np.log(2)))))
                if i==1:
                    Noise[l,i]=1/(1/(omega_minus_one_95_pol*np.exp((l*(l+1)*(thetaFWHM_95)**2)/(8*np.log(2))))+1/(omega_minus_one_150_pol*np.exp((l*(l+1)*(thetaFWHM_150)**2)/(8*np.log(2)))))
    return Noise

def Ali_Noise_level(): # remember to set noise to (ells,2) matrix
    ali_noise_dls=np.loadtxt('Ali_noise.dat',usecols=(1,2))
    print(ali_noise_dls)
    insert_first_thirty_rows=np.zeros((30,2))
    half_Noise=np.insert(ali_noise_dls,0,insert_first_thirty_rows,axis=0)
    insert_last_many_rows=np.zeros((ells-len(half_Noise),2))
    Noised=np.insert(half_Noise,len(half_Noise),insert_last_many_rows,axis=0)
    Noise=noise_dls2cls(Noised)
    return Noise

def noise_dls2cls(dls):
    cls=np.zeros((ells,2))
    for i in np.arange(2):
        for l in np.arange(30,620):
            cls[l,i]=dls[l,i]*(2*np.pi)/(l*(l+1))
    return cls

def pico_noise_dls2cls(dls):
    cls=np.zeros((ells,2))
    for i in np.arange(2):
        for l in np.arange(10,1001):
            cls[l,i]=dls[l,i]*(2*np.pi)/(l*(l+1))
    return cls


def log_Cubic_interpolate(xx,yy):
    from scipy.interpolate import CubicSpline
    logx=np.log10(xx)
    logy=np.log10(yy)
    Cubic_interp=CubicSpline(logx,logy)
    log_Cubic_interp=lambda zz: np.power(10, Cubic_interp(np.log10(zz)))
    return log_Cubic_interp

def fg_res_dls2cls(dls):
    cls=np.zeros((ells,3))
    for i in np.arange(3):
        for l in np.arange(1,1024):
            cls[l,i]=dls[l,i]*(2*np.pi)/(l*(l+1))
    return cls



def ali_fg_res():
    # from matplotlib import pyplot as plt
    TTFgres=np.load('TT_fg.npy')
    TEFgres=np.load('TE_fg.npy')
    EEFgres=np.load('EE_fg.npy')
    half_res=np.stack((TTFgres,EEFgres,TEFgres),axis=1)
    insert_last_many_rows=np.zeros((ells-len(half_res),3))
    fg_resdls=np.insert(half_res,len(half_res),insert_last_many_rows,axis=0)
    fg_res=fg_res_dls2cls(fg_resdls)
    return fg_res

def pico_fg_res():
    # from matplotlib import pyplot as plt
    TTFgres=np.load('fg_TT_results.npy')
    TEFgres=np.load('fg_TE_results.npy')
    EEFgres=np.load('fg_EE_results.npy')
    half_res=np.stack((TTFgres,EEFgres,TEFgres),axis=1)
    insert_last_many_rows=np.zeros((ells-len(half_res),3))
    fg_resdls=np.insert(half_res,len(half_res),insert_last_many_rows,axis=0)
    fg_res=fg_res_dls2cls(fg_resdls)
    return fg_res

def pico_noise_level():
    pico_noise_TT=np.load('./Nl_TT_results.npy')
    pico_noise_EE=np.load('./Nl_EE_results.npy')
    half_noise=np.stack((pico_noise_TT,pico_noise_EE),axis=1)
    insert_last_many_rows=np.zeros((ells-len(half_noise),2))
    noise_dls=np.insert(half_noise,len(half_noise),insert_last_many_rows,axis=0)
    pico_noise=pico_noise_dls2cls(noise_dls)
    return pico_noise


"""
    main function
"""

##initial value of six fundamental parameters
#ombh2_0 = 0.02242
#omch2_0 = 0.11933
#As_0 = 2.105209331337507e-09 #scalar_amp(1)
#ns_0 = 0.9665 #scalar_spectral_index(1)
#optical_depth_0 = 0.0561 # re_optical_depth
#hubble_0 = 67.66
#thetastarmc_0=1.040997
#DM_mass_0=1e2
#DM_Pann_0=0
#DM_Gamma_0=0
#params_num=7
#ells=4001

#print('thetastarmc_0= ',thetastarmc_0)
##set initial value of fundamental Model
#set_l_max_scalar(ells-1)
#set_ns(ns_0)
#set_ombh2(ombh2_0)
#set_omch2(omch2_0)
#set_optical_depth(optical_depth_0)
#set_As(As_0)
#set_hubble(hubble_0)
#set_DM_mass(DM_mass_0)
#set_DM_Pann(DM_Pann_0)
#set_DM_Gamma(DM_Gamma_0)


# clsn=initial_totCL()

# DM_Pann_CLprime=DM_Pann_prime(DM_mass_0)
# DM_Gamma_CLprime=DM_Gamma_prime(DM_mass_0)
# ombh2_CLprime=ombh2_prime()
# omch2_CLprime=omch2_prime()
# thetastarmc_CLprime=thetastarmc_prime()
# optical_depth_CLprime=optical_depth_prime()
# ns_CLprime=ns_prime()
# As_CLprime=As_prime()

# np.savez('derivative_data',DM_Pann_CLprime=DM_Pann_CLprime,DM_Gamma_CLprime=DM_Gamma_CLprime,ombh2_CLprime=ombh2_CLprime,omch2_CLprime=omch2_CLprime,thetastarmc_CLprime=thetastarmc_CLprime,optical_depth_CLprime=optical_depth_CLprime,ns_CLprime=ns_CLprime,As_CLprime=As_CLprime)

'''
    load data from npz file !!!! when you have to calculate error at different fsky or noise level, please use this after derivative are calculated, be careful about the dark matter mass here, default value is 1e-4 GeV
'''

# npzdata=np.load('derivative_data.npz')
# DM_Pann_CLprime=npzdata['DM_Pann_CLprime']
# DM_Gamma_CLprime=npzdata['DM_Gamma_CLprime']
# ombh2_CLprime=npzdata['ombh2_CLprime']
# omch2_CLprime=npzdata['omch2_CLprime']
# thetastarmc_CLprime=npzdata['thetastarmc_CLprime']
# optical_depth_CLprime=npzdata['optical_depth_CLprime']
# ns_CLprime=npzdata['ns_CLprime']
# As_CLprime=npzdata['As_CLprime']


'''
    fisher matrix -> one sigma error
'''

# Derivative_Pann=np.zeros((ells,3,params_num))
# Derivative_Pann[:,:,0]=ombh2_CLprime
# Derivative_Pann[:,:,1]=omch2_CLprime
# Derivative_Pann[:,:,2]=thetastarmc_CLprime
# Derivative_Pann[:,:,3]=optical_depth_CLprime
# Derivative_Pann[:,:,4]=ns_CLprime
# Derivative_Pann[:,:,5]=As_CLprime
# Derivative_Pann[:,:,6]=DM_Pann_CLprime # choose one of clprime

# Derivative_Gamma=np.zeros((ells,3,params_num))
# Derivative_Gamma[:,:,0]=ombh2_CLprime
# Derivative_Gamma[:,:,1]=omch2_CLprime
# Derivative_Gamma[:,:,2]=thetastarmc_CLprime
# Derivative_Gamma[:,:,3]=optical_depth_CLprime
# Derivative_Gamma[:,:,4]=ns_CLprime
# Derivative_Gamma[:,:,5]=As_CLprime
# Derivative_Gamma[:,:,6]=DM_Gamma_CLprime


'''
    setup for Noise level fsky and foreground residual
'''

# Noisecls=np.zeros((ells,2)) # zero noise level
# fgrescls=np.zeros((ells,3))
# fsky=1.0

# Noisecls=Ali_Noise_level()
# fsky=0.37

# fgrescls=ali_fg_res()
# Noisecls=Ali_Noise_level()
# fsky=0.37

# fgrescls=pico_fg_res()
# fgrescls=np.zeros((ells,3))
# Noisecls=pico_noise_level()
# fsky=1.0



##################################################

# # one sigma error for Pann
# sigma_TT=2*get_TT_fisher_matrix(Derivative_Pann,10,4001)
# sigma_EE=2*get_EE_fisher_matrix(Derivative_Pann,10,4001)
# sigma_TE=2*get_TE_fisher_matrix(Derivative_Pann,10,4001)
# sigma_combined=2*get_combined_fisher_matrix(Derivative_Pann,10,4001)
# print('Pann_sigma_TT= ',sigma_TT)
# print('Pann_sigma_EE= ',sigma_EE)
# print('Pann_sigma_TE= ',sigma_TE)
# print('Pann_sigma_combined= ',sigma_combined)

# # one sigma error for gamma
# sigma_TT=get_TT_fisher_matrix(Derivative_Gamma,30,621)
# sigma_EE=get_EE_fisher_matrix(Derivative_Gamma,30,621)
# sigma_TE=get_TE_fisher_matrix(Derivative_Gamma,30,621)
# sigma_combined=get_combined_fisher_matrix(Derivative_Gamma,30,621)
# print('Gamma_sigma_TT= ',sigma_TT)
# print('Gamma_sigma_EE= ',sigma_EE)
# print('Gamma_sigma_TE= ',sigma_TE)
# print('Gamma_sigma_combined= ',sigma_combined)

'''
    fisher matrix at different dark matter mass. Take care of what error you are calculating: one sigma or two sigma error
'''

# time_0=time.time()
# # fgrescls=ali_fg_res()
# fgrescls=np.zeros((ells,3))
# # Noisecls=Ali_Noise_level()
# # fsky=0.5
# Noisecls=np.zeros((ells,2))
# # fgrescls=pico_fg_res()
# # Noisecls=pico_noise_level()
# fsky=1.0

# mass_len=50
# dm_mass_set=np.geomspace(1e-5,5e3,mass_len)

# npzdata=np.load('derivative_data.npz')
# ombh2_CLprime=npzdata['ombh2_CLprime']
# omch2_CLprime=npzdata['omch2_CLprime']
# thetastarmc_CLprime=npzdata['thetastarmc_CLprime']
# optical_depth_CLprime=npzdata['optical_depth_CLprime']
# ns_CLprime=npzdata['ns_CLprime']
# As_CLprime=npzdata['As_CLprime']

# Derivative_Pann=np.zeros((ells,3,params_num))
# Derivative_Pann[:,:,0]=ombh2_CLprime
# Derivative_Pann[:,:,1]=omch2_CLprime
# Derivative_Pann[:,:,2]=thetastarmc_CLprime
# Derivative_Pann[:,:,3]=optical_depth_CLprime
# Derivative_Pann[:,:,4]=ns_CLprime
# Derivative_Pann[:,:,5]=As_CLprime
# sigma_of_Pann=np.zeros((mass_len,4))
# for i in range(len(dm_mass_set)):
#     print('this is loop:', i)
#     DM_mass_0=dm_mass_set[i]
#     xcordinate,dp,ls,min_step_mat,DM_Pann_CLprime=DM_Pann_prime(DM_mass_0)
#     Derivative_Pann[:,:,6]=DM_Pann_CLprime
#     sigma_of_Pann[i,0]=2*get_TT_fisher_matrix(Derivative_Pann,10,4001)[params_num-1]
#     sigma_of_Pann[i,1]=2*get_EE_fisher_matrix(Derivative_Pann,10,4001)[params_num-1]
#     sigma_of_Pann[i,2]=2*get_TE_fisher_matrix(Derivative_Pann,10,4001)[params_num-1]
#     sigma_of_Pann[i,3]=2*get_combined_fisher_matrix(Derivative_Pann,10,4001)[params_num-1]
# print('calculate time is: ',time.time()-time_0)

# np.save('sigma1.npy',sigma_of_Pann)

######Remember dm_mass_set has to have the same length, so you need to check the length of dm_mass_seta and the corresponding sigma_of_Pann.npy

# sigma_of_Pann=np.load('sigma1.npy')
# xs=np.geomspace(1e-5,5e3,1000)
# cs=log_Cubic_interpolate(dm_mass_set, sigma_of_Pann[:,0])
# plt.loglog(xs,cs(xs),color='b')
# plt.xlabel('dark matter mass')
# plt.ylabel('energy injection')
# cs=log_Cubic_interpolate(dm_mass_set, sigma_of_Pann[:,1])
# plt.loglog(xs,cs(xs),color='g')
# cs=log_Cubic_interpolate(dm_mass_set, sigma_of_Pann[:,2])
# plt.loglog(xs,cs(xs),color='y')
# cs=log_Cubic_interpolate(dm_mass_set, sigma_of_Pann[:,3])
# plt.loglog(xs,cs(xs),color='r')
# plt.legend(['TT','EE','TE','combined'],loc='upper right')
# plt.ylim((1e-28,1e-25))
# plt.xlim((1e-5,5e3))
# plt.title('CVL condition lmax=4000')
# # plt.show()
# plt.savefig('constraint.png',dpi=300)

####plot for paper
#sigma_of_Pann_nofg=np.load('sigma7.npy')
#sigma_of_Pann_res=np.load('sigma6.npy')
#sigma_of_Pann_nofg_pico=np.load('sigma11.npy')
#sigma_of_Pann_res_pico=np.load('sigma10.npy')
#sigma_of_Pann_cv=np.load('sigma15.npy')

# xs=np.geomspace(1e-5,5e3,1000)
# cs=log_Cubic_interpolate(dm_mass_set, sigma_of_Pann_res[:,3])
# plt.loglog(xs,cs(xs),color='y',linestyle='-',linewidth=1)
# cs=log_Cubic_interpolate(dm_mass_set, sigma_of_Pann_nofg[:,3])
# plt.loglog(xs,cs(xs),color='y',linestyle='--',linewidth=1)
# cs=log_Cubic_interpolate(dm_mass_set, sigma_of_Pann_res_pico[:,3])
# plt.loglog(xs,cs(xs),color='b',linestyle='-',linewidth=1)
# cs=log_Cubic_interpolate(dm_mass_set, sigma_of_Pann_nofg_pico[:,3])
# plt.loglog(xs,cs(xs),color='b',linestyle='--',linewidth=1)
# cs=log_Cubic_interpolate(dm_mass_set, sigma_of_Pann_cv[:,3])
# plt.loglog(xs,cs(xs),color='k',linestyle='-',linewidth=1)

# plt.legend(['Ali with fg residual fsky=0.5','ali without fg residual fsky=0.5','pico with fg residual fsky=1','pico without fg residual fsky=1','CV fsky=1'],loc='upper right')
# plt.title('ali l from 30 to 620, pico l from 10 to 1000')
# plt.xlabel(r'$m_\chi\left[GeV\right]$')
# plt.ylabel(r'$\left\langle \sigma v\right\rangle/ m_\chi \left[cm^3s^{-1}GeV^{-1}\right] $')
# plt.ylim((1e-28,1e-26))
# plt.xlim((1e-5,5e3))
# plt.title(r'constraints derived by fisher matrix')
# plt.savefig('constraint6.png',dpi=300)

#sigma_of_Pann_nofg=np.load('sigma9.npy')
#sigma_of_Pann_res=np.load('sigma8.npy')
#xs=np.geomspace(1e-5,5e3,1000)
#cs=log_Cubic_interpolate(dm_mass_set, sigma_of_Pann_res[:,3])
#plt.loglog(xs,cs(xs),color='y',linestyle='-')
#cs=log_Cubic_interpolate(dm_mass_set, sigma_of_Pann_nofg[:,3])
#plt.loglog(xs,cs(xs),color='y',linestyle='--')
#plt.legend(['Ali with fg residual','Ali without fg residual'],loc='upper right')
#plt.xlabel(r'$m_\chi\left[GeV\right]$')
#plt.ylabel(r'$\left\langle \sigma v\right\rangle/ m_\chi \left[cm^3s^{-1}GeV^{-1}\right] $')
#plt.ylim((1e-28,1e-26))
#plt.xlim((1e-5,5e3))
#plt.title(r'constraint derived by fisher matrix, fsky=0.5, ali noise')
#plt.savefig('constraint2.png',dpi=300)



##########gamma

# npzdata=np.load('derivative_data.npz')
# ombh2_CLprime=npzdata['ombh2_CLprime']
# omch2_CLprime=npzdata['omch2_CLprime']
# thetastarmc_CLprime=npzdata['thetastarmc_CLprime']
# optical_depth_CLprime=npzdata['optical_depth_CLprime']
# ns_CLprime=npzdata['ns_CLprime']
# As_CLprime=npzdata['As_CLprime']

# Derivative_Gamma=np.zeros((2501,3,params_num))
# Derivative_Gamma[:,:,0]=ombh2_CLprime
# Derivative_Gamma[:,:,1]=omch2_CLprime
# Derivative_Gamma[:,:,2]=thetastarmc_CLprime
# Derivative_Gamma[:,:,3]=optical_depth_CLprime
# Derivative_Gamma[:,:,4]=ns_CLprime
# Derivative_Gamma[:,:,5]=As_CLprime
# sigma_of_Gamma=np.zeros((mass_len,4))
# for i in range(len(dm_mass_set)):
#     print('this is loop:', i)
#     DM_mass_0=dm_mass_set[i]
#     xcordinate,dp,ls,min_step_mat,DM_Gamma_CLprime=DM_Gamma_prime(DM_mass_0)
#     Derivative_Gamma[:,:,6]=DM_Gamma_CLprime
#     sigma_of_Gamma[i,0]=get_TT_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
#     sigma_of_Gamma[i,1]=get_EE_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
#     sigma_of_Gamma[i,2]=get_TE_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
#     sigma_of_Gamma[i,3]=get_combined_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
# print('calculate time is: ',time.time()-time_0)


'''
    plot derivative against powerspectrum
'''

# xcordinate,dp,ls,min_step_mat,DM_Pann_CLprime=DM_Pann_prime(1e-4)
# plt.figure(1)
# plt.semilogx(ls[10:],ls[10:]*(ls[10:]+1)*DM_Pann_CLprime[10:,0])
# plt.xlabel('ls')
# plt.ylabel('derivativeTT')
# plt.figure(2)
# plt.semilogx(ls[10:],ls[10:]*(ls[10:]+1)*DM_Pann_CLprime[10:,1])
# plt.xlabel('ls')
# plt.ylabel('derivativeEE')
# plt.figure(3)
# plt.semilogx(ls[10:],ls[10:]*(ls[10:]+1)*DM_Pann_CLprime[10:,2])
# plt.xlabel('ls')
# plt.ylabel('derivativeTE')

# xcordinate,dp,ls,min_step_mat,DM_Pann_CLprime=DM_Pann_prime(1e-4)
# plt.figure(1)
# plt.loglog(ls[2:],ls[2:]*(ls[2:]+1)*np.absolute(DM_Pann_CLprime[2:,0])/(2*np.pi))
# plt.xlabel('ls')
# plt.ylabel('DM_Pann_CLprimeTT')
# plt.figure(2)
# plt.loglog(ls[2:],ls[2:]*(ls[2:]+1)*np.absolute(DM_Pann_CLprime[2:,1])/(2*np.pi))
# plt.xlabel('ls')
# plt.ylabel('DM_Pann_CLprimeEE')
# plt.figure(3)
# plt.loglog(ls[2:],ls[2:]*(ls[2:]+1)*np.absolute(DM_Pann_CLprime[2:,2])/(2*np.pi))
# plt.xlabel('ls')
# plt.ylabel('DM_Pann_CLprimeTE')

'''
    plot partialCL/partial parameter against CL
'''

# plt.figure(1)
# plt.loglog(ls[10:],ls[10:]*(ls[10:]+1)*masked_pos[10:,0]/(clsn[10:,0]*2*np.pi),'-')
# plt.loglog(ls[10:],ls[10:]*(ls[10:]+1)*np.abs(masked_neg[10:,0])/(clsn[10:,0]*2*np.pi),'--')
# plt.xlabel('ls')
# plt.ylabel('DM_Pann_CLprimeTT')

# plt.figure(2)
# plt.loglog(ls[10:],ls[10:]*(ls[10:]+1)*masked_pos[10:,1]/(clsn[10:,1]*2*np.pi),'-')
# plt.loglog(ls[10:],ls[10:]*(ls[10:]+1)*np.abs(masked_neg[10:,1])/(clsn[10:,1]*2*np.pi),'--')
# plt.xlabel('ls')
# plt.ylabel('DM_Pann_CLprimeEE')

# plt.figure(3)
# plt.loglog(ls[10:],ls[10:]*(ls[10:]+1)*maskedTE_pos[10:]/(2*np.pi),'-')
# plt.loglog(ls[10:],np.abs(ls[10:]*(ls[10:]+1)*maskedTE_neg[10:]/(2*np.pi)),'--')
# # plt.ylim((1e25,1e34))
# plt.xlabel('ls')
# plt.ylabel('DM_Pann_CLprimeTE')

'''
   apply np.where to seperate positive powerspectrum and negative powerspectrum. in this process, negative portion in positive powerspectrum are set to zero, and the same action is done on negative powerspectrum. Then np.ma.masked_where function are applied to make zero part of the previous action into nan because plt.plot will not plot nan. In this way, tendency to zero when switching from positive to negative and vice versa can be avoided.
'''

#xcordinate,dp,ls,min_step_mat,DM_Pann_CLprime=DM_Pann_prime(1e-4)
#pos_Pann_CLprime=np.where(DM_Pann_CLprime<0,0,DM_Pann_CLprime)
#neg_Pann_CLprime=np.where(DM_Pann_CLprime>0,0,DM_Pann_CLprime)
#TE_pos=np.where(DM_Pann_CLprime[:,2]/clsn[:,2]<0,0,DM_Pann_CLprime[:,2]/clsn[:,2])
#TE_neg=np.where(DM_Pann_CLprime[:,2]/clsn[:,2]>0,0,DM_Pann_CLprime[:,2]/clsn[:,2])
#maskedTE_pos=np.ma.masked_where(TE_pos==0, TE_pos)
#maskedTE_neg=np.ma.masked_where(TE_neg==0, TE_neg)
#masked_pos=np.ma.masked_where(pos_Pann_CLprime==0, pos_Pann_CLprime)
#masked_neg=np.abs(np.ma.masked_where(neg_Pann_CLprime==0, neg_Pann_CLprime))
##**********#
#from matplotlib import rc # latex in matplotlib
#rc('text',usetex=True)
##**********#
#fig, axes = plt.subplots(nrows=3,ncols=1,sharex=True,figsize=(10,10))
## fig.suptitle('Title of this figure')
#yticks1=[1e25,1e27,1e29,1e31]
#ax1=axes[0]
#ax1.loglog(ls[10:],ls[10:]*(ls[10:]+1)*masked_pos[10:,0]/(clsn[10:,0]*2*np.pi),'-')
#ax1.loglog(ls[10:],ls[10:]*(ls[10:]+1)*np.abs(masked_neg[10:,0])/(clsn[10:,0]*2*np.pi),'--')
#ax1.set_ylabel(r'$\ell\left(\ell+1\right)/2\pi \cdot \frac{\partial C_\ell^{TT}}{\partial\theta_i}/C_\ell$',fontsize=10)
#ax1.set_yscale('log')
#ax1.set_yticks(yticks1)
##**********#
#yticks2=[1e27,1e28,1e29,1e30,1e31]
#ax2=axes[1]
#ax2.loglog(ls[10:],ls[10:]*(ls[10:]+1)*masked_pos[10:,1]/(clsn[10:,1]*2*np.pi),'-')
#ax2.loglog(ls[10:],ls[10:]*(ls[10:]+1)*np.abs(masked_neg[10:,1])/(clsn[10:,1]*2*np.pi),'--')
## ax2.set_title("Second subplot")
#ax2.set_yscale('log')
#ax2.set_yticks(yticks2)
#ax2.set_ylabel(r'$\ell\left(\ell+1\right)/2\pi \cdot \frac{\partial C_\ell^{EE}}{\partial\theta_i}/C_\ell$',fontsize=10)
##**********#
#yticks3=[1e26,1e28,1e30,1e32,1e34]
#ax3=axes[2]
#ax3.loglog(ls[10:],ls[10:]*(ls[10:]+1)*maskedTE_pos[10:]/(2*np.pi),'-')
#ax3.loglog(ls[10:],np.abs(ls[10:]*(ls[10:]+1)*maskedTE_neg[10:]/(2*np.pi)),'--')
#ax3.set_xscale('log')
#ax3.set_xticks([1e1,1e2,1e3])
#ax3.set_yscale('log')
#ax3.set_yticks(yticks3)
#ax3.set_ylabel(r'$\ell\left(\ell+1\right)/2\pi \cdot \frac{\partial C_\ell^{TE}}{\partial\theta_i}/C_\ell$',fontsize=10)
#ax3.set_xlabel(r'$\ell$')
##**********#
#plt.savefig('derivative.png',dpi=300)


'''
    fisher matrix at different noise level
'''

# Noisecls_0=Ali_Noise_level()
# # Noisecls_0=cangjs_Noise_level()
# fsky=0.37
# noise_low_fold=0
# noise_high_fold=100000
# noise_length=10
# nsfold=np.linspace(noise_low_fold,noise_high_fold,noise_length)
# sigma_of_Gamma=np.zeros((noise_length,4))
# sigma_of_Pann=np.zeros((noise_length,4))
# print("calculation on different noise level")
# for i in range(len(nsfold)):
#     print(i)
#     Noisecls=Noisecls_0*nsfold[i]
#     sigma_of_Gamma[i,0]=get_TT_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
#     sigma_of_Gamma[i,1]=get_EE_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
#     sigma_of_Gamma[i,2]=get_TE_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
#     sigma_of_Gamma[i,3]=get_combined_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
#     sigma_of_Pann[i,0]=get_TT_fisher_matrix(Derivative_Pann,10,1001)[params_num-1]
#     sigma_of_Pann[i,1]=get_EE_fisher_matrix(Derivative_Pann,10,1001)[params_num-1]
#     sigma_of_Pann[i,2]=get_TE_fisher_matrix(Derivative_Pann,10,1001)[params_num-1]
#     sigma_of_Pann[i,3]=get_combined_fisher_matrix(Derivative_Pann,10,1001)[params_num-1]
# plt.figure(1)
# plt.semilogy(nsfold,sigma_of_Pann[:,0],color='b')
# plt.semilogy(nsfold,sigma_of_Pann[:,1],color='g')
# plt.semilogy(nsfold,sigma_of_Pann[:,2],color='y')
# plt.semilogy(nsfold,sigma_of_Pann[:,3],color='r')
# plt.legend(['TT','EE','TE','combined'])
# plt.xlabel('noise fold compare to standard noise level')
# plt.ylabel('one sigma error for dm_pann at different noise level')
# plt.savefig('error_at_dif_noise_Pann.png',dpi=150)
# plt.close()
# plt.figure(2)
# plt.semilogy(nsfold,sigma_of_Gamma[:,0],color='b')
# plt.semilogy(nsfold,sigma_of_Gamma[:,1],color='g')
# plt.semilogy(nsfold,sigma_of_Gamma[:,2],color='y')
# plt.semilogy(nsfold,sigma_of_Gamma[:,3],color='r')
# plt.legend(['TT','EE','TE','combined'])
# plt.xlabel('noise fold compare to standard noise level')
# plt.ylabel('one sigma error for dm_gamma at different noise level')
# plt.savefig('error_at_dif_noise_Gamma.png',dpi=150)
# plt.close()

'''
    fisher matrix at different fsky
'''

# Noisecls=Ali_Noise_level()
# fsky_start=0.1
# fsky_end=1
# fsky_length=20
# sigma_of_Gamma=np.zeros((fsky_length,4))
# sigma_of_Pann=np.zeros((fsky_length,4))
# fs=np.linspace(fsky_start, fsky_end, fsky_length)
# print("calculation on different fsky:")
# for i in range(len(fs)):
#     fsky=fs[i]
#     print(i)
#     sigma_of_Gamma[i,0]=get_TT_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
#     sigma_of_Gamma[i,1]=get_EE_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
#     sigma_of_Gamma[i,2]=get_TE_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
#     sigma_of_Gamma[i,3]=get_combined_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
#     sigma_of_Pann[i,0]=get_TT_fisher_matrix(Derivative_Pann,10,1001)[params_num-1]
#     sigma_of_Pann[i,1]=get_EE_fisher_matrix(Derivative_Pann,10,1001)[params_num-1]
#     sigma_of_Pann[i,2]=get_TE_fisher_matrix(Derivative_Pann,10,1001)[params_num-1]
#     sigma_of_Pann[i,3]=get_combined_fisher_matrix(Derivative_Pann,10,1001)[params_num-1]
# plt.figure(1)
# plt.semilogy(fs,sigma_of_Pann[:,0],color='b')
# plt.semilogy(fs,sigma_of_Pann[:,1],color='g')
# plt.semilogy(fs,sigma_of_Pann[:,2],color='y')
# plt.semilogy(fs,sigma_of_Pann[:,3],color='r')
# plt.legend(['TT','EE','TE','combined'])
# plt.xlabel('fsky')
# plt.ylabel('one sigma error for dm_pann')
# plt.savefig('error_at_dif_fsky_pann.png',dpi=150)
# plt.close()
# plt.figure(2)
# plt.semilogy(fs,sigma_of_Gamma[:,0],color='b')
# plt.semilogy(fs,sigma_of_Gamma[:,1],color='g')
# plt.semilogy(fs,sigma_of_Gamma[:,2],color='y')
# plt.semilogy(fs,sigma_of_Gamma[:,3],color='r')
# plt.legend(['TT','EE','TE','combined'])
# plt.xlabel('fsky')
# plt.ylabel('one sigma error for dm_gamma')
# plt.savefig('error_at_dif_fsky_gamma.png',dpi=150)
# plt.close()

'''
    how sigma change with different lmax or lmin
'''

npzdata=np.load('derivative_data.npz')
DM_Pann_CLprime=npzdata['DM_Pann_CLprime']
DM_Gamma_CLprime=npzdata['DM_Gamma_CLprime']
ombh2_CLprime=npzdata['ombh2_CLprime']
omch2_CLprime=npzdata['omch2_CLprime']
thetastarmc_CLprime=npzdata['thetastarmc_CLprime']
optical_depth_CLprime=npzdata['optical_depth_CLprime']
ns_CLprime=npzdata['ns_CLprime']
As_CLprime=npzdata['As_CLprime']

Derivative_Pann=np.zeros((ells,3,params_num))
Derivative_Pann[:,:,0]=ombh2_CLprime
Derivative_Pann[:,:,1]=omch2_CLprime
Derivative_Pann[:,:,2]=thetastarmc_CLprime
Derivative_Pann[:,:,3]=optical_depth_CLprime
Derivative_Pann[:,:,4]=ns_CLprime
Derivative_Pann[:,:,5]=As_CLprime
Derivative_Pann[:,:,6]=DM_Pann_CLprime # choose one of clprime

Noisecls=np.zeros((ells,2)) # zero noise level
fgrescls=np.zeros((ells,3))
fsky=1.0

lmax_start=600
lmax_end=4000
lmax_length=100
lmax_set=np.arange(lmax_start,lmax_end,lmax_length)
# sigma_of_Gamma=np.zeros((lmax_length,4))
sigma_of_Pann=np.zeros((34,4))
lmin=10
print("calculation on different lmax:")
for i,lmax in enumerate(lmax_set):
    print(f'this is loop:{i},where lmax={lmax}')
    lmax=lmax
    # sigma_of_Gamma[i,0]=2*get_TT_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
    # sigma_of_Gamma[i,1]=2*get_EE_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
    # sigma_of_Gamma[i,2]=2*get_TE_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
    # sigma_of_Gamma[i,3]=2*get_combined_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
    sigma_of_Pann[i,0]=2*get_TT_fisher_matrix(Derivative_Pann,lmin,lmax)[params_num-1]
    sigma_of_Pann[i,1]=2*get_EE_fisher_matrix(Derivative_Pann,lmin,lmax)[params_num-1]
    sigma_of_Pann[i,2]=2*get_TE_fisher_matrix(Derivative_Pann,lmin,lmax)[params_num-1]
    sigma_of_Pann[i,3]=2*get_combined_fisher_matrix(Derivative_Pann,lmin,lmax)[params_num-1]

plt.figure(1)
plt.semilogy(lmax_set,sigma_of_Pann[:,0],color='b')
plt.semilogy(lmax_set,sigma_of_Pann[:,1],color='g')
plt.semilogy(lmax_set,sigma_of_Pann[:,2],color='y')
plt.semilogy(lmax_set,sigma_of_Pann[:,3],color='r')
plt.legend(['TT','EE','TE','combined'])
plt.xlabel('lmax')
plt.ylabel('sigma of dm pann')
# plt.ylim((1e-28,1e-26))
# plt.show()
plt.savefig('diff_lmax.png',dpi=300)

# lmin_start=10
# lmin_end=1000
# lmin_length=50
# lmin_set=np.arange(lmin_start,lmin_end,lmin_length)
# # sigma_of_Gamma=np.zeros((lmin_length,4))
# sigma_of_Pann=np.zeros((20,4))
# lmax=2500
# print("calculation on different lmin:")
# for i,lmin in enumerate(lmin_set):
#     print(f'this is loop:{i},where lmin={lmin}')
#     lmin=lmin
#     # sigma_of_Gamma[i,0]=2*get_TT_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
#     # sigma_of_Gamma[i,1]=2*get_EE_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
#     # sigma_of_Gamma[i,2]=2*get_TE_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
#     # sigma_of_Gamma[i,3]=2*get_combined_fisher_matrix(Derivative_Gamma,10,1001)[params_num-1]
#     sigma_of_Pann[i,0]=2*get_TT_fisher_matrix(Derivative_Pann,lmin,lmax)[params_num-1]
#     sigma_of_Pann[i,1]=2*get_EE_fisher_matrix(Derivative_Pann,lmin,lmax)[params_num-1]
#     sigma_of_Pann[i,2]=2*get_TE_fisher_matrix(Derivative_Pann,lmin,lmax)[params_num-1]
#     sigma_of_Pann[i,3]=2*get_combined_fisher_matrix(Derivative_Pann,lmin,lmax)[params_num-1]
# plt.figure(1)
# plt.semilogy(lmin_set,sigma_of_Pann[:,0],color='b')
# plt.semilogy(lmin_set,sigma_of_Pann[:,1],color='g')
# plt.semilogy(lmin_set,sigma_of_Pann[:,2],color='y')
# plt.semilogy(lmin_set,sigma_of_Pann[:,3],color='r')
# plt.legend(['TT','EE','TE','combined'])
# plt.xlabel(f'lmin')
# plt.ylabel(f'one sigma error for dm_pann')
# # plt.ylim((1e-28,1e-26))
# # plt.show()
# plt.savefig('diff_lmin.png',dpi=300)
# # plt.close()

### check the correctness of derivative

# xcordinate,dp,ls,min_step_mat,DM_Pann_CLprime=DM_Pann_prime(1e-4)

# plt.figure(1)
# plt.loglog(xcordinate,np.absolute(dp[200,2,:]))
# plt.xlabel('DM_Pann')
# plt.ylabel('dp')
# plt.figure(2)
# plt.semilogx(ls[2:],ls[2:]*ls[2:]*DM_Pann_CLprime[2:,1])
# plt.xlabel('ls')
# plt.ylabel('DM_Pann_CLprime')
# plt.figure(3)
# plt.plot(ls[2:],min_step_mat[2:,2])
# plt.xlabel('ls')
# plt.ylabel('min_step_mat')
# plt.show()






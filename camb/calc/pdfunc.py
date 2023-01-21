import numpy as np
from setparams import *
import os

def initial_totCL():
    os.system('./camb test.ini')
    test_scalCls=np.loadtxt('./test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
    insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
    dlsn=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
    clsn=dls2cls(dlsn)
    return clsn

def dls2cls(dls,ells):
    cls=np.zeros((ells,3))
    for i in np.arange(3):
        for l in np.arange(2,ells):
            cls[l,i]=dls[l,i]*(2*np.pi)/(l*(l+1))
    return cls

def DM_Pann_prime(DM_mass,ells,length,start_footstep,end_footstep):
    set_DM_mass(DM_mass)
    set_DM_Gamma(0)
    set_DM_Pann(0)
    DM_Pann_0=0
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    dp=np.zeros((ells,3,length))
    difference=np.ones((ells,3,length))
    DM_Pann_CLprime=np.zeros((ells,3))
    ls=np.arange(ells)
    x1=DM_Pann_0
    set_DM_Pann(x1)
    os.system('./camb test.ini')
    test_scalCls=np.loadtxt('./test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
    insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
    dlsn=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
    clsn=dls2cls(dlsn,ells)
    for index,h in enumerate(xcordinate):
        x2=DM_Pann_0+h
        set_DM_Pann(x2)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('./test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        dlsp=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)# the output test_scalCls is Dl(i.e. l(l+1)C_l)
        clsp=dls2cls(dlsp,ells)
        dp[:,:,index]=(clsp-clsn)/h
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,ells):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        print('loop=',index)
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
    print("best_step is:",best_h,'where h=',xcordinate[best_h])
    set_DM_Pann(xcordinate[best_h])
    os.system('./camb test.ini')
    test_scalCls=np.loadtxt('./test_scalCls.dat',usecols=(1,2,3))
    insert_first_two_rows=np.array([[0,0,0],[0,0,0]])
    dlsd=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
    delta_dls=dlsd-dlsn
    DM_Pann_CLprime=dp[:,:,best_h]
    set_DM_Pann(0)
    return xcordinate,dp,ls,min_step_mat,DM_Pann_CLprime

# xcordinate,dp,ls,min_step_mat,DM_Pann_CLprime=DM_Pann_prime(DM_mass_0,ells,length=30,start_footstep=1e-26,end_footstep=1e-31)

def DM_Gamma_prime(params_value):
    set_DM_mass(params_value)
    set_DM_Pann(0)
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

if __name__=="__main__":
    from setparams import *
    from consts import *

    print(os.getcwd())
    os.chdir("../")
    print(os.getcwd())

    set_l_max_scalar(ells-1)
    set_hubble(hubble_0)
    set_ombh2(ombh2_0)
    set_omch2(omch2_0)
    set_optical_depth(optical_depth_0)
    set_As(As_0)
    set_ns(ns_0)
    set_DM_Pann(DM_Pann_0)
    set_DM_Gamma(DM_Gamma_0)
    set_DM_mass(DM_mass_0)


    # xcordinate,dp,ls,min_step_mat,DM_Pann_CLprime=DM_Pann_prime(DM_mass_0,ells,length=30,start_footstep=1e-26,end_footstep=1e-31)







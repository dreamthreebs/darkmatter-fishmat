import numpy as np
from setparams import *
from consts import *
from pdfunc import *

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
show_all_params()

xcordinate,dp,ls,min_step_mat,DM_Pann_CLprime=DM_Pann_prime(DM_mass_0, ells, length=30, start_footstep=1e-26,end_footstep=1e-31)
checksave_pd_stability(xcordinate, dp, ls, min_step_mat, DM_Pann_CLprime,'./plot/figure/pdstability/DM_Pann_CLprime.png')
np.save('./data/pd_data/DM_Pann_CLprime.npy',DM_Pann_CLprime)

xcordinate,dp,ls,min_step_mat,DM_Gamma_CLprime=DM_Gamma_prime(DM_mass_0, ells, length=30, start_footstep=1e-24, end_footstep=1e-31)
checksave_pd_stability(xcordinate, dp, ls, min_step_mat, DM_Gamma_CLprime,'./plot/figure/pdstability/DM_Gamma_CLprime.png')
np.save('./data/pd_data/DM_Gamma_CLprime.npy',DM_Gamma_CLprime)

xcordinate,dp,ls,min_step_mat,thetastarmc_CLprime=thetastarmc_prime(hubble_0, ells, length=30, start_footstep=2e-1, end_footstep=1e-4)
checksave_pd_stability(xcordinate, dp, ls, min_step_mat, thetastarmc_CLprime,'./plot/figure/pdstability/thetastarmc_CLprime.png')
np.save('./data/pd_data/thetastarmc_CLprime.npy',thetastarmc_CLprime)

xcordinate,dp,ls,min_step_mat,ombh2_CLprime=ombh2_prime(ombh2_0, ells, length=30, start_footstep=1e-1, end_footstep=1e-7)
checksave_pd_stability(xcordinate, dp, ls, min_step_mat, ombh2_CLprime,'./plot/figure/pdstability/ombh2_CLprime.png')
np.save('./data/pd_data/ombh2_CLprime.npy',ombh2_CLprime)

xcordinate,dp,ls,min_step_mat,omch2_CLprime=omch2_prime(omch2_0, ells, length=30, start_footstep=1e-1, end_footstep=1e-7)
checksave_pd_stability(xcordinate, dp, ls, min_step_mat, omch2_CLprime,'./plot/figure/pdstability/omch2_CLprime.png')
np.save('./data/pd_data/omch2_CLprime.npy',omch2_CLprime)

xcordinate,dp,ls,min_step_mat,optical_depth_CLprime=optical_depth_prime(optical_depth_0, ells, length=30, start_footstep=1e-1, end_footstep=1e-4)
checksave_pd_stability(xcordinate, dp, ls, min_step_mat, optical_depth_CLprime,'./plot/figure/pdstability/optical_depth_CLprime.png')
np.save('./data/pd_data/optical_depth_CLprime.npy',optical_depth_CLprime)

xcordinate,dp,ls,min_step_mat,ns_CLprime=ns_prime(ns_0, ells, length=30, start_footstep=1e-1, end_footstep=1e-7)
checksave_pd_stability(xcordinate, dp, ls, min_step_mat, ns_CLprime,'./plot/figure/pdstability/ns_CLprime.png')
np.save('./data/pd_data/ns_CLprime.npy',ns_CLprime)

xcordinate,dp,ls,min_step_mat,As_CLprime=As_prime(As_0, ells, length=30, start_footstep=1e-1, end_footstep=1e-7)
checksave_pd_stability(xcordinate, dp, ls, min_step_mat, As_CLprime,'./plot/figure/pdstability/As_CLprime.png')
np.save('./data/pd_data/ns_CLprime.npy',As_CLprime)


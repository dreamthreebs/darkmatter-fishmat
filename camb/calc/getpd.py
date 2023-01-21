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

xcordinate,dp,ls,min_step_mat,DM_Pann_CLprime=DM_Pann_prime(DM_mass_0,ells,length=30,start_footstep=1e-26,end_footstep=1e-31)




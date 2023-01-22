ombh2_0 = 0.02242
omch2_0 = 0.11933
As_0 = 2.105209331337507e-09 #scalar_amp(1)
ns_0 = 0.9665 #scalar_spectral_index(1)
optical_depth_0 = 0.0561 # re_optical_depth
hubble_0 = 67.66
DM_mass_0=1e2
DM_Pann_0=0
DM_Gamma_0=0

thetastarmc_0=1.040997
params_num=7
ells=4001

def show_all_params():
    print(f"fiducial values are:")
    print(f"ombh2_0={ombh2_0}")
    print(f"omch2_0={omch2_0}")
    print(f"As_0={As_0}")
    print(f"ns_0={ns_0}")
    print(f"optical_depth_0={optical_depth_0}")
    print(f"hubble_0={hubble_0}")
    print(f"thetastarmc_0={thetastarmc_0}")
    print(f"DM_Pann_0={DM_Pann_0}")
    print(f"DM_Gamma_0={DM_Gamma_0}")
    print(f"DM_mass_0={DM_mass_0}")
    print(f"ells={ells}")
    print(f"parameter_number={params_num}")

if __name__=="__main__":
    show_all_params()



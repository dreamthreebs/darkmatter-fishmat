import os
import numpy as np
from sigfunc import *


if __name__=="__main__":
    print(datetime.datetime.now())
    print("file position:",os.getcwd())
    os.chdir("../")
    print("work position:",os.getcwd())

    set_all_params()

    # check_all_pd()
    # check_fisher_pann()
    # check_fisher_gamma()

    # calc_pann_two_sigma_error(cls=initial_totCL(), nls=zero_noise(ells), fgres=zero_fgres(ells), ells=ells, lmin=10, lmax=ells, fsky=1,filename='./data/sigma_data/CVL_fsky1_lmin1_lmaxells.npy')

    # calc_pann_two_sigma_error(cls=initial_totCL(), nls=ali_noise_level(ells), fgres=ali_fg_res(ells), ells=ells, lmin=30, lmax=620, fsky=0.4, filename='./data/sigma_data/ali_noise_fgres_fsky37_lmin30_lmax620.npy')
    # calc_pann_two_sigma_error(cls=initial_totCL(), nls=ali_noise_level(ells), fgres=zero_fgres(ells), ells=ells, lmin=30, lmax=620, fsky=0.4, filename='./data/sigma_data/ali_noise_nofgres_fsky40_lmin30_lmax620.npy')

    # calc_pann_two_sigma_error(cls=initial_totCL(), nls=ali_noise_level(ells), fgres=zero_fgres(ells), ells=ells, lmin=30, lmax=3000, fsky=0.4, filename='./data/sigma_data/ali_noise_nofgres_fsky40_lmin30_lmax3000.npy')
    # calc_pann_two_sigma_error(cls=initial_totCL(), nls=pico_noise_level(ells), fgres=zero_fgres(ells), ells=ells, lmin=10, lmax=3000, fsky=0.7, filename='./data/sigma_data/pico_noise_nofgres_fsky37_lmin10_lmax3000.npy')

    # calc_pann_two_sigma_error(cls=initial_totCL(), nls=ali_noise_level(ells), fgres=zero_fgres(ells), ells=ells, lmin=30, lmax=4000, fsky=0.4, filename='./data/sigma_data/ali_noise_nofgres_fsky40_lmin30_lmax4000.npy')
    # calc_pann_two_sigma_error(cls=initial_totCL(), nls=pico_noise_level(ells), fgres=zero_fgres(ells), ells=ells, lmin=10, lmax=4000, fsky=0.7, filename='./data/sigma_data/pico_noise_nofgres_fsky37_lmin10_lmax4000.npy')

#     calc_pann_two_sigma_error(cls=initial_totCL(), nls=ali_noise_level(ells), fgres=zero_fgres(ells), ells=ells, lmin=30, lmax=1000, fsky=0.4, filename='./data/sigma_data/ali_noise_nofgres_fsky40_lmin30_lmax1000.npy')
#     calc_pann_two_sigma_error(cls=initial_totCL(), nls=pico_noise_level(ells), fgres=zero_fgres(ells), ells=ells, lmin=10, lmax=1000, fsky=0.7, filename='./data/sigma_data/pico_noise_nofgres_fsky37_lmin10_lmax1000.npy')

    pd_gamma_diff_mass=np.load('./data/pd_gamma_50_data/pd_gamma_50_diff_mass.npy')  # todo
    calc_gamma_two_sigma_error_with_prior_pd_gamma(cls=initial_totCL(), nls=zero_noise(ells), fgres=zero_fgres(ells), ells=ells, lmin=10, lmax=4000, fsky=1, filename='./data/sig_gamma_data/CVL_fsky1_lmin1_lmaxells.npy', pd_gamma_diff_mass=pd_gamma_diff_mass)

    calc_gamma_two_sigma_error_with_prior_pd_gamma(cls=initial_totCL(), nls=ali_noise_level(ells), fgres=ali_fg_res(ells), ells=ells, lmin=30, lmax=620, fsky=0.4, filename='./data/sig_gamma_data/ali_noise_fgres_fsky40_lmin30_lmax620.npy', pd_gamma_diff_mass=pd_gamma_diff_mass)
    calc_gamma_two_sigma_error_with_prior_pd_gamma(cls=initial_totCL(), nls=ali_noise_level(ells), fgres=zero_fgres(ells), ells=ells, lmin=30, lmax=620, fsky=0.4, filename='./data/sig_gamma_data/ali_noise_nofgres_fsky40_lmin30_lmax620.npy', pd_gamma_diff_mass=pd_gamma_diff_mass)
    calc_gamma_two_sigma_error_with_prior_pd_gamma(cls=initial_totCL(), nls=ali_noise_level(ells), fgres=zero_fgres(ells), ells=ells, lmin=30, lmax=1000, fsky=0.4, filename='./data/sig_gamma_data/ali_noise_nofgres_fsky40_lmin30_lmax1000', pd_gamma_diff_mass=pd_gamma_diff_mass)
    calc_gamma_two_sigma_error_with_prior_pd_gamma(cls=initial_totCL(), nls=ali_noise_level(ells), fgres=zero_fgres(ells), ells=ells, lmin=30, lmax=3000, fsky=0.4, filename='./data/sig_gamma_data/ali_noise_nofgres_fsky40_lmin30_lmax3000', pd_gamma_diff_mass=pd_gamma_diff_mass)
    calc_gamma_two_sigma_error_with_prior_pd_gamma(cls=initial_totCL(), nls=ali_noise_level(ells), fgres=zero_fgres(ells), ells=ells, lmin=30, lmax=4000, fsky=0.4, filename='./data/sig_gamma_data/ali_noise_nofgres_fsky40_lmin30_lmax4000', pd_gamma_diff_mass=pd_gamma_diff_mass)

    calc_gamma_two_sigma_error_with_prior_pd_gamma(cls=initial_totCL(), nls=pico_noise_level(ells), fgres=zero_fgres(ells), ells=ells, lmin=10, lmax=1000, fsky=0.7, filename='./data/sig_gamma_data/pico_noise_nofgres_fsky70_lmin10_lmax1000', pd_gamma_diff_mass=pd_gamma_diff_mass)
    calc_gamma_two_sigma_error_with_prior_pd_gamma(cls=initial_totCL(), nls=pico_noise_level(ells), fgres=zero_fgres(ells), ells=ells, lmin=10, lmax=3000, fsky=0.7, filename='./data/sig_gamma_data/pico_noise_nofgres_fsky70_lmin10_lmax3000', pd_gamma_diff_mass=pd_gamma_diff_mass)
    calc_gamma_two_sigma_error_with_prior_pd_gamma(cls=initial_totCL(), nls=pico_noise_level(ells), fgres=zero_fgres(ells), ells=ells, lmin=10, lmax=4000, fsky=0.7, filename='./data/sig_gamma_data/pico_noise_nofgres_fsky70_lmin10_lmax4000', pd_gamma_diff_mass=pd_gamma_diff_mass)

    print(datetime.datetime.now())

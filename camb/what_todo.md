# Calculate constrains that CMB powerspectrum can give, for dark matter annilation and decay parameters
only consider annilation situation
dark matter mass from 1e-5 -> 5e3 
use 7 parameters: ombh2, omch2, 100thetastarmc, optical depth, ln(10^10As), ns, dark matter annilation energy injection/dark matter decay energy injection
want plot a diagram: x is dark matter mass, y is energy injection, x footstep is 100, y is given by derivative-> one-error-sigma
need to do decay situation!!!!

lmin=30,lmax=620
sigma1.npy dark matter mass from 1e-5 -> 5e3 one sigma error footstep=100 alinoise nofgres fsky=0.37
sigma2.npy dark matter mass from 1e-5 -> 5e3 one sigma error footstep=50 alinoise nofgres fsky=0.37
sigma3.npy dark matter mass from 1e-5 -> 5e3 one sigma error footstep=10 alinoise alifgres fsky=0.37
sigma4.npy dark matter mass from 1e-5 -> 5e3 two sigma error footstep=50 alinoise alifgres fsky=0.37 constraint.png
sigma5.npy dark matter mass from 1e-5 -> 5e3 two sigma error footstep=50 alinoise alifgres fsky=1 constraint1.png
sigma6.npy dark matter mass from 1e-5 -> 5e3 two sigma error footstep=50 alinoise alifgres fsky=0.5 constraint2.png
sigma7.npy dark matter mass from 1e-5 -> 5e3 two sigma error footstep=50 alinoise nofgres fsky=0.5 constraint3.png
sigma8.npy dark matter mass from 1e-5 -> 5e3 two sigma error footstep=50 piconoise picofgres fsky=1.0
sigma9.npy dark matter mass from 1e-5 -> 5e3 two sigma error footstep=50 piconoise nofgres fsky=1.0

lmin=10,lmax=1000
sigma10.npy dark matter mass from 1e-5 -> 5e3 two sigma error footstep=50 piconoise picofgres fsky=1.0 
sigma11.npy dark matter mass from 1e-5 -> 5e3 two sigma error footstep=50 piconoise nofgres fsky=1.0 in dm_fisher
sigma15.npy dark matter mass from 1e-5 -> 5e3 two sigma error footstep=50 nonoise nofgres fsky=1.0  in random

lmin=10,lmax=2500
sigma13.npy dark matter mass from 1e-5 -> 5e3 two sigma error footstep=50 piconoise picofgres fsky=1.0 in dm_fisher1
sigma14.npy dark matter mass from 1e-5 -> 5e3 two sigma error footstep=50 piconoise nofgres fsky=1.0  in dm_fisher2

fgresidual.py plus the foreground residue

dm_fisher is different from random folder on lmax_scalar, where dm_fisher have lmax_scalar=4001, random have lmax_scalar=2501

derivative.npz also store derivative information on 8 parameters
sigma1.npy dark matter mass from 1e-5 -> 5e3 two sigma error footstep=50 nonoise nofgres fsky=1.0 in dm_fisher


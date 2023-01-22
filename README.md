# Fishmat
Fisher matrix method to calculate constraint on dark matter parameters
## Dependency
ifort(intel fortran), icc(intel c++), python(numpy, matplotlib, scipy)
## Structure
./camb/calc/ : all codes for calculation partial derivative(pd), and one error sigma(sig)  
./camb/data/ : pd data calculated by getpd.py and sigma data calculated by getsig.py  
./camb/input/ : noise level(nls) and foreground residual(fgres) for ali, pico, and planck  
./camb/plot/ : plot all figures and save in figure folder  
## How to use
make Hyrec; make camb
run ./camb/calc/getpd.py to get all pd data in ./camb/data/pd_data/*CLprime.npy  
run ./camb/calc/getsig.py to get all one sigma error in ./camb/data/sigma_data/sigma*.npy  
run ./camb/plot/plot*.py to get different figures of what you want  

'''
l_max_scalar, hubble parameter, ombh2, omch2, optical_depth, As, ns, DM_Pann, DM_Gamma, DM_mass are set in this file
'''

import fileinput
import re

def set_l_max_scalar(new_value):
    for line in fileinput.input('./test.ini',inplace=True):
        if line.startswith('l_max_scalar'):
            print('l_max_scalar = '+str(new_value))
        else:
            print(line.strip())

def set_hubble(new_value):
    for line in fileinput.input('./test.ini',inplace=True):
        if line.startswith('hubble'):
            print('hubble = '+str(new_value))
        else:
            print(line.strip())

def set_ombh2(new_value):
    for line in fileinput.input('./test.ini',inplace=True):
        if line.startswith('ombh2'):
            print('ombh2 = '+str(new_value))
        else:
            print(line.strip())

def set_omch2(new_value):
    for line in fileinput.input('./test.ini',inplace=True):
        if line.startswith('omch2'):
            print('omch2 = '+str(new_value))
        else:
            print(line.strip())

def set_As(new_value):
    for line in fileinput.input('./test.ini',inplace=True):
        if line.startswith('scalar_amp(1)'):
            print('scalar_amp(1) = '+str(new_value))
        else:
            print(line.strip())

def set_ns(new_value):
    for line in fileinput.input('./test.ini',inplace=True):
        if line.startswith('scalar_spectral_index(1)'):
            print('scalar_spectral_index(1) = '+str(new_value))
        else:
            print(line.strip())

def set_optical_depth(new_value):
    for line in fileinput.input('./test.ini',inplace=True):
        if line.startswith('re_optical_depth'):
            print('re_optical_depth = '+str(new_value))
        else:
            print(line.strip())

def set_DM_Pann(new_value):
    for line in fileinput.input('./test.ini',inplace=True):
        if line.startswith('DM_Pann'):
            print('DM_Pann = '+str(new_value))
        else:
            print(line.strip())

def set_DM_Gamma(new_value):
    for line in fileinput.input('./test.ini',inplace=True):
        if line.startswith('DM_Gamma'):
            print('DM_Gamma = '+str(new_value))
        else:
            print(line.strip())

def set_DM_mass(new_value):
    for line in fileinput.input('./test.ini',inplace=True):
        if line.startswith('DM_mass'):
            print('DM_mass = '+str(new_value))
        else:
            print(line.strip())

def get_thetastarmc_from_files():
    filename="./output.txt"
    with fileinput.FileInput(filename) as f:
        for line in f :
            if line.startswith('100 theta (CosmoMC)'):
                thetastarmc=re.findall(r"\d+\.?\d*",line[line.find("="):])[0]
                return float(thetastarmc)

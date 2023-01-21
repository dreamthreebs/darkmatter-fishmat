Files for PBH Hawking Radiation and Decay/Annihilation of Dark Matter

EFF.tar.gz contains files needed, untar EFF.tar.gz if EFF directory is corrupt
Input parameters documented in input_DM.dat and input_PBH.dat
Replace history.c and history.h with one of the following to switch between models

Contains following processes:

1 --- PBH
      For Black Holes in [2E13,1E18] g mass range with 5 different initial Kerr spins: [0,0.25,0.5,0.75,0.999]
      Files:
             history_PBH.c
             history_PBH_K1.h    ----> a=0
             history_PBH_K2.h    ----> a=0.25
             history_PBH_K3.h    ----> a=0.5
             history_PBH_K4.h    ----> a=0.75
             history_PBH_K5.h    ----> a=0.999
             history_PBH_K6.h    ----> a=0.9999
      The Following are EFF induced by Primary dNdEdt
             history_PBH_K1B.h
             history_PBH_K2B.h
             history_PBH_K3B.h
             history_PBH_K4B.h
             history_PBH_K5B.h
             history_PBH_K6B.h

2 --- PBH_Lite
      For Schwarzschild Black Holes in [1E15,1E17] g mass range
      Files:
             history_PBH_Lite.c
             history_PBH_Lite.h

3 --- DM_Elec
      For Dark Matter decay/annhilation into electron+positron pair, with 2 DM distribution models.
      HMG    : homogeneous distribution
      MODEL1 : DM destribution described in ArXiv.2002.03380
      Files:
             history_Elec.c
             history_Elec_HMG.h      --> HMG
             history_Elec_MODEL1.h   --> MODEL1

4 --- DM_Phot
      Same as process 3 but for gamma channel
      Files:
             history_Phot.c
             history_Phot_HMG.h      --> HMG 
             history_Phot_MODEL1.h   --> MODEL1

5 --- DM_SM
      For Dark Matter decay/annhilation into all EM Standard Model channels except for processes 3 and 4.
      Homogeneous distribution used as default
      Files:
             history_SM.c
             history_SM_HMG.h
      Usage:
             Default channel is Higgs, to change channel, replace 'Higgs' in h file by one of following strings:
             [ Muon   Tau   Q   Charm   Bottom   Top   W   Z   Gluon   Higgs ]

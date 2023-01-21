    !---------------------------------------------------------------------------------------------------
    ! Recombination module for CAMB, using HyRec
    ! Author: Antony Lewis
    !---------------------------------------------------------------------------------------------------
	!
	! P D Meerburg 2015, added P_ann (dark matter annihilation coefficient) [cm^3/s/GeV]
	!					 added fine structure constant variation
	!				     added electron mass variation
    module Recombination
    use constants
    use AMLUtils
    implicit none
    private

    type RecombinationParams
    !place to define spurious other recombination params
	real(dl) :: Pann  ! Dark matter annihilation rate
	real(dl) :: DM_Gamma  ! Dark matter decay parameter
	real(dl) :: DM_mass  ! Dark matter particle mass
	real(dl) :: BH_mass  ! Blackhole mass
	real(dl) :: DM_BH_frac  ! ratio of dark matter mass that is blackhole
    real(dl) :: FineSRat !Fine structure ratio (default 1)
    real(dl) :: EMassRat !ELectron mass ratio (default 1)
    integer  :: Rfback
    end type RecombinationParams

    character(LEN=*), parameter :: Recombination_Name = 'HyRec'

    public RecombinationParams, Recombination_xe, Recombination_tm,   &
    Recombination_ReadParams, Recombination_SetDefParams, Recombination_init, &
    Recombination_Validate, Recombination_Name


    contains

    subroutine Recombination_ReadParams(R, Ini)
    use IniFile
    Type(RecombinationParams) :: R
    Type(TIniFile) :: Ini
    !read other possible params
        R%Pann = Ini_Read_Double_File(Ini, 'DM_Pann', 0.0D0) !input in cm^3/s/GeV
        R%DM_Gamma = Ini_Read_Double_File(Ini, 'DM_Gamma', 0.0D0) !input in 1/s
        R%DM_mass = Ini_Read_Double_File(Ini, 'DM_mass', 1.0D0) !input in GeV
        R%BH_mass = Ini_Read_Double_File(Ini, 'BH_mass', 1.0D100) !input in grams
        R%DM_BH_frac = Ini_Read_Double_File(Ini, 'DM_BH_frac', 0.0D0) !input unitless
	R%FineSRat = Ini_Read_Double_File(Ini, 'FineS', 1.0D0)
	R%EMassRat = Ini_Read_Double_File(Ini, 'EMass', 1.0D0)
	
    end subroutine Recombination_ReadParams

    

    subroutine Recombination_SetDefParams(R)
    type (RecombinationParams) ::R


    end subroutine Recombination_SetDefParams


    subroutine Recombination_Validate(R, OK)
    Type(RecombinationParams), intent(in) :: R
    logical, intent(inout) :: OK


    end subroutine Recombination_Validate



    function Recombination_tm(a)
    real(dl), intent(in) :: a
    real(dl) Recombination_tm, hyrec_tm
    external hyrec_tm

    Recombination_tm =  hyrec_tm(a);

    end function Recombination_tm


    function Recombination_xe(a)
    real(dl), intent(in) :: a
    real(dl) Recombination_xe,hyrec_xe
    external hyrec_xe

    Recombination_xe = hyrec_xe(a);

    end function Recombination_xe


	!modded PDM 2015
    subroutine Recombination_init(Recomb, OmegaC, OmegaB, OmegaN, Omegav, h0inp, tcmb, yp, num_nu)
    use AMLUtils
    implicit none
    Type (RecombinationParams), intent(in) :: Recomb
    real(dl), intent(in) :: OmegaC, OmegaB, OmegaN, OmegaV, h0inp, tcmb, yp, num_nu
    external rec_build_history_camb
	
    call rec_build_history_camb(OmegaC, OmegaB, OmegaN, Omegav, h0inp, tcmb, yp, num_nu, & 
    			Recomb%Pann, Recomb%DM_Gamma, Recomb%DM_mass, Recomb%BH_mass, Recomb%DM_BH_frac, &
                        Recomb%FineSRat, Recomb%EMassRat, Recomb%Rfback)

    end subroutine Recombination_init




    end module Recombination


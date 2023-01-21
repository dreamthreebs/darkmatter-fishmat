/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         history.h: functions for recombination history                                        */
/*                                                                                               */
/*************************************************************************************************/


/**** Structure for cosmological parameters. 
      Include information on starting and ending redshit and timestep  
      Added May 2012: the period of Hydrogen recombination is evaluated right away so
      tables of radiation field can be reduced to only what is needed and memory is saved. 
****/

/*Efficiency Tables are comprised of Hydrogen ionization, Lyman-Alpha, and Heating*/

/*****Hydrogen Ionization Efficiency Tables*****/

#define ANNIHILATION_ENERGY_AXIS_SIZE			62
#define ANNIHILATION_REDSHIFT_AXIS_SIZE			63

#define ANNIHILATION_H_IONIZATION_ENERGY_AXIS_FILE	"./EFF/M_ANN_1.txt"
#define ANNIHILATION_H_IONIZATION_REDSHIFT_AXIS_FILE	"./EFF/Redshift_1.txt"
#define ANNIHILATION_H_IONIZATION_EFF_TABLE_FILE	"./EFF/Higgs_HIon_HMG.txt"

#define ANNIHILATION_ALPHA_ENERGY_AXIS_FILE		"./EFF/M_ANN_2.txt"
#define ANNIHILATION_ALPHA_REDSHIFT_AXIS_FILE		"./EFF/Redshift_2.txt"
#define ANNIHILATION_ALPHA_EFF_TABLE_FILE		"./EFF/Higgs_LyA_HMG.txt"

#define ANNIHILATION_HEATING_ENERGY_AXIS_FILE		"./EFF/M_ANN_3.txt"
#define ANNIHILATION_HEATING_REDSHIFT_AXIS_FILE		"./EFF/Redshift_3.txt"
#define ANNIHILATION_HEATING_EFF_TABLE_FILE		"./EFF/Higgs_Heat_HMG.txt"

/*****Decay Efficiency Tables*****/

#define DECAY_ENERGY_AXIS_SIZE				62
#define DECAY_REDSHIFT_AXIS_SIZE			63

#define DECAY_H_IONIZATION_ENERGY_AXIS_FILE		"./EFF/M_DEC_1.txt"
#define DECAY_H_IONIZATION_REDSHIFT_AXIS_FILE		"./EFF/Redshift_4.txt"
#define DECAY_H_IONIZATION_EFF_TABLE_FILE		"./EFF/Higgs_HIon_DEC.txt"

#define DECAY_ALPHA_ENERGY_AXIS_FILE			"./EFF/M_DEC_2.txt"
#define DECAY_ALPHA_REDSHIFT_AXIS_FILE			"./EFF/Redshift_5.txt"
#define DECAY_ALPHA_EFF_TABLE_FILE			"./EFF/Higgs_LyA_DEC.txt"

#define DECAY_HEATING_ENERGY_AXIS_FILE			"./EFF/M_DEC_3.txt"
#define DECAY_HEATING_REDSHIFT_AXIS_FILE		"./EFF/Redshift_6.txt"
#define DECAY_HEATING_EFF_TABLE_FILE			"./EFF/Higgs_Heat_DEC.txt"

/*****Blackhole Efficiency Tables*****/

#define BH_ENERGY_AXIS_SIZE                             200
#define BH_REDSHIFT_AXIS_SIZE                           63

#define BH_H_IONIZATION_ENERGY_AXIS_FILE                "./EFF/PBH_Temperature_1.txt"
#define BH_H_IONIZATION_REDSHIFT_AXIS_FILE              "./EFF/Redshift_7.txt"
#define BH_H_IONIZATION_EFF_TABLE_FILE                  "./EFF/PBH_HIon_K1.txt"

#define BH_ALPHA_ENERGY_AXIS_FILE                       "./EFF/PBH_Temperature_2.txt"
#define BH_ALPHA_REDSHIFT_AXIS_FILE                     "./EFF/Redshift_8.txt"
#define BH_ALPHA_EFF_TABLE_FILE                         "./EFF/PBH_LyA_K1.txt"

#define BH_HEATING_ENERGY_AXIS_FILE                     "./EFF/PBH_Temperature_3.txt"
#define BH_HEATING_REDSHIFT_AXIS_FILE                   "./EFF/Redshift_9.txt"
#define BH_HEATING_EFF_TABLE_FILE                       "./EFF/PBH_Heat_K1.txt"

/**********/

typedef struct {
   double T0;                          /* CMB temperature today in K*/
   double obh2, omh2, odeh2, okh2;     /* density parameters */
   double Y;                           /* primordial helium abundance */
   double Nnueff;                      /* effective number of neutrinos */

   double fsR, meR;              /* alpha_fs/alpha_fs(today) and me/me(today) (Added April 2012)*/
     
   double pann;                 /* Dark matter annihilation parameter (added January 2015 ) */
   double dm_gamma;             /* Dark matter decay parameter (added February 2016 ) */

/* Added Sept 2016 */
   double dm_mass;		/* Dark matter mass in GeV */
   double bh_mass;		/* Black hole mass in grams */
   double dm_bh_frac;		/* Fraction of dark matter mass that is black hole */
/****/

   /* Secondary parameters, to avoid recalculating every time */
   double nH0;                  /* density of hydrogen today in m^{-3} */  
   double fHe;                  /* Helium fraction by number */

   long nz;                     /* total number of redshift steps */
   long izH0;                   /* index when H recombination starts to be considered */  
   double zH0;                  /* Redshift at which H recombination starts (zH0 = z[izH0]) */
   long nzrt;                   /* number of redshift steps while radiative transfer is computed */
} REC_COSMOPARAMS;

/**** Structure for switches for additional output ****/

typedef struct {
  int print_spec, print_21cm;
  char file_spec[100];
  char file_21cm[100];
  int Nz_spec, Nz_21cm;
  double zmin_spec, zmax_spec, zmin_21cm, zmax_21cm;
} REC_EXTRAS;

/**** Structure for radiation field outputs. Added January 2015 ****/

typedef struct {
  double **Dfnu;
  double *Dfminus_Ly[3];
} RAD_OUTPUTS;

/**** Structure for Effective Efficiency Tables ****/

typedef struct {
   double *energy_axis;
   double *redshift_axis;
   double **efficiency;
} EFFECTIVE_EFF;

void allocate_rad_output(RAD_OUTPUTS *rad, long int nzrt);
void free_rad_output(RAD_OUTPUTS *rad);





void rec_get_cosmoparam(FILE *fin, REC_COSMOPARAMS *param, REC_EXTRAS *extras);
double rec_HubbleConstant(REC_COSMOPARAMS *param, double z);
double rec_Tmss(double xe, double Tr, double H, double fHe, 
                double fsR, double meR, double dE_dtdV, double nh);
double rec_dTmdlna(double xe, double Tm, double Tr, double H, double fHe, 
                   double fsR, double meR, double dE_dtdV, double nh);
void rec_get_xe_next1_He(REC_COSMOPARAMS *param, double z_in, double *xHeII, 
                         double *dxHeIIdlna_prev, double *dxHeIIdlna_prev2, int *post_saha); 
double rec_xH1s_postSaha(REC_COSMOPARAMS *param, unsigned iz_out, double z_out, double xHeII_out, 
                         HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params,
		         double **Dfminus_hist, double *Dfminus_Ly_hist[], double **Dfnu_hist, int *post_saha, 
                         EFFECTIVE_EFF *ionization_annihilation, EFFECTIVE_EFF *ionization_decay, EFFECTIVE_EFF *ionization_blackhole, 
                         EFFECTIVE_EFF *alpha_annihilation, EFFECTIVE_EFF *alpha_decay, EFFECTIVE_EFF *alpha_blackhole);
void get_rec_next2_HHe(REC_COSMOPARAMS *param, unsigned iz_in, double z_in, double Tm_in, double *xH1s, double *xHeII,
                       HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params, double **Dfminus_hist, double *Dfminus_Ly_hist[], double **Dfnu_hist,
                       double *dxHIIdlna_prev,  double *dxHeIIdlna_prev, double *dxHIIdlna_prev2, double *dxHeIIdlna_prev2, int *post_saha, 
                       EFFECTIVE_EFF *ionization_annihilation, EFFECTIVE_EFF *ionization_decay, EFFECTIVE_EFF *ionization_blackhole, 
                       EFFECTIVE_EFF *alpha_annihilation, EFFECTIVE_EFF *alpha_decay, EFFECTIVE_EFF *alpha_blackhole);
void rec_get_xe_next1_H(REC_COSMOPARAMS *param, double z_in, double xe_in, double Tm_in, double *xe_out,
                        HRATEEFF *rate_table, unsigned iz, TWO_PHOTON_PARAMS *twog_params, double **Dfminus_hist, double *Dfminus_Ly_hist[], 
                        double **Dfnu_hist, double *dxedlna_prev, double *dxedlna_prev2, int *post_saha, 
                        EFFECTIVE_EFF *ionization_annihilation, EFFECTIVE_EFF *ionization_decay, EFFECTIVE_EFF *ionization_blackhole, 
                        EFFECTIVE_EFF *alpha_annihilation, EFFECTIVE_EFF *alpha_decay, EFFECTIVE_EFF *alpha_blackhole);
void rec_get_xe_next2_HTm(int func_select, REC_COSMOPARAMS *param, double z_in, double xe_in, double Tm_in, double *xe_out, double *Tm_out,
                          HRATEEFF *rate_table, unsigned iz, TWO_PHOTON_PARAMS *twog_params, double **Dfminus_hist, double *Dfminus_Ly_hist[], 
                          double **Dfnu_hist, double *dxedlna_prev, double *dTmdlna_prev, double *dxedlna_prev2, double *dTmdlna_prev2, 
                          EFFECTIVE_EFF *ionization_annihilation, EFFECTIVE_EFF *ionization_decay, EFFECTIVE_EFF *ionization_blackhole, 
                          EFFECTIVE_EFF *alpha_annihilation, EFFECTIVE_EFF *alpha_decay, EFFECTIVE_EFF *alpha_blackhole, 
                          EFFECTIVE_EFF *thermal_annihilation, EFFECTIVE_EFF *thermal_decay, EFFECTIVE_EFF *thermal_blackhole);
void rec_build_history(REC_COSMOPARAMS *param, HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params,
                        EFFECTIVE_EFF *eff_ionization_annihilation, EFFECTIVE_EFF *eff_ionization_decay, EFFECTIVE_EFF *eff_ionization_blackhole, 
                        EFFECTIVE_EFF *eff_alpha_annihilation, EFFECTIVE_EFF *eff_alpha_decay, EFFECTIVE_EFF *eff_alpha_blackhole,
                        EFFECTIVE_EFF *eff_thermal_annihilation, EFFECTIVE_EFF *eff_thermal_decay, EFFECTIVE_EFF *eff_thermal_blackhole, 
                        double *xe_output, double *Tm_output, double **Dfnu_hist, double *Dfminus_Ly_hist[3]);
void read_eff_eff(EFFECTIVE_EFF *effeff, char *ENERGY_AXIS_FILE, char *REDSHIFT_AXIS_FILE, char *EFFICIENCY_FILE, unsigned ENERGY_AXIS_SIZE, unsigned REDSHIFT_AXIS_SIZE);
double interp_eff_eff(EFFECTIVE_EFF *effeff, double energy, double z, unsigned ENERGY_AXIS_SIZE, unsigned REDSHIFT_AXIS_SIZE);
double dE_dtdV_func(double z, REC_COSMOPARAMS *param, EFFECTIVE_EFF *annihilation, EFFECTIVE_EFF *decay, EFFECTIVE_EFF *blackhole);
void free_eff_eff(EFFECTIVE_EFF *effeff, unsigned REDSHIFT_AXIS_SIZE);

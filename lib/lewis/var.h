/**********************************************************************/
/*                                                                    */
/*              ----- Basic Structure Header File -----               */
/*                                                                    */
/*     Author Information:                                            */
/*                                                                    */
/*          Gary A. Allen, Jr.                                        */
/*          NASA Ames Research Center                                 */
/*          Mail Stop 230-2                                           */
/*          Moffett Field CA 94035-1000                               */
/*                                                                    */
/*          Telephone:  (650) 604-4228                                */
/*          FAX:        (650) 604-0350                                */
/*                                                                    */
/*          E-mail:  gallen@mail.arc.nasa.gov                         */
/*                                                                    */
/*                Version:  Mk 1.8, 22 June 2005                      */
/*                 Written by Gary A. Allen, Jr.                      */
/*                                                                    */
/* Notices                                                            */
/* -------                                                            */
/*                                                                    */
/* Copyright 1999-2002 U. S. Government as represented by the         */
/* Administrator of the National Aeronautics and Space Administration */
/*                                                                    */
/* This software may be used, copied, and provided to others only as  */
/* permitted under the terms of the contract or other agreement under */
/* which it was acquired from the U.S. Government. If you have        */
/* received this software for use under a Government Contract,        */
/* Cooperative Agreement or Grant the software is for federal         */
/* government research purposes only, and may not be further          */
/* distributed to third parties. Neither title to nor ownership of    */
/* the software is hereby transferred. This notice shall remain on    */
/* all copies of the software.                                        */
/*                                                                    */
/* Disclaimers                                                        */
/* -----------                                                        */
/*                                                                    */
/* This software is provided "as is" without any warranty of any      */
/* kind, either expressed, implied, or statutory, including, but not  */
/* limited to, any warranty that the software will conform to         */
/* specifications, any implied warranties of merchantability, fitness */
/* for a particular purpose, or freedom from infringement, any        */
/* warranty that the software will be error free, or any warranty     */
/* that documentation, if provided, will conform to the software. In  */
/* no event shall the U.S. Government, or the U.S. Government's       */
/* contractors or subcontractors, be liable for any damages,          */
/* including, but not limited to, direct, indirect, special or        */
/* consequential damages, arising out of, resulting from, or in any   */
/* way connected with this software, whether or not based upon        */
/* warranty, contract, tort, or otherwise, whether or not injury was  */
/* sustained by persons or property or otherwise, and whether or not  */
/* loss was sustained from, or arose out of the results of, or use    */
/* of, the software or services provided hereunder.                   */
/*                                                                    */
/* Recipient agrees to waive any and all claims against the U.S.      */
/* Government, and the U.S. Government's contractors and              */
/* subcontractors, and shall indemnify and hold harmless the U.S.     */
/* Government, and the U.S. Government's contractors and              */
/* subcontractors, for any damage that recipient may incur from       */
/* recipient's prior or future use of the provided software,          */
/* including any damages from products based on, or resulting from,   */
/* the use thereof.                                                   */
/*                                                                    */
/* If further release or distribution of this software or technical   */
/* data derived from this software is permitted, Recipient agrees to  */
/* obtain this identical disclaimer of warranty agreement and waiver  */
/* of claims, indemnification and hold harmless agreement with any    */
/* entities that receive the software or technical data derived from  */
/* the software.                                                      */
/*                                                                    */
/**********************************************************************/

#include <math.h>
#include <stdio.h>

#define FALSE   0
#define TRUE    1

#define Error	-1

/* Macro for stripping off garbage from a data line */
#define Strip(fp) for (is = 0; is < 132; is++) if (0xa == fgetc(fp)) break

#ifdef Six_DoF
#define NMAX    23      /* Total number of first order ODEs             */
#define NSTG    15      /* Number of ODEs including stag. point heating */
#define NDYN    13      /* Number of ODEs driving trajectory dynamics   */
#else
#define NMAX    17      /* Total number of first order ODEs             */
#define NSTG     9      /* Number of ODEs including stag. point heating */
#define NDYN     7      /* Number of ODEs driving trajectory dynamics   */
#endif
#define NCNT     6      /* Number of vectors in continuation matrix     */

#define M_grav_max  55
#define N_grav_max  55

/* Math constants */
#define Pi      	3.1415926535897932384626433832795028841972
#define Pi2             6.2831853071795864769252867665590057683944
#define Pi_2            1.5707963267948966192313216916397514420986
#define Rad70		1.2217304763960307038465835379420288994100
#define Sqrt2		1.4142135623730950488016887242096980785696

/* Unit conversions */
#define Deg2rad         1.7453292519943295769236907684886127134428e-2
#define Rad2deg         57.295779513082320876798154814105170332405
#define Meter2inch	39.370078740157480314960629921259842519685
#define Meter2foot	3.2808398950131233595800524934383202099738
#define Foot2meter	0.3048
#define Kelvin2Rankine	1.8
#define Rankine2Kelvin  0.5555555555555555555555555555555555555555
#define Km2meter	1000.0
#define Meter2km	0.001
#define Second2day	1.1574074074074074074074074074074074074075e-5
#define Second2minute   1.6666666666666666666666666666666666666667e-2
#define Day2second	86400.0
#define Au2meter	1.49597870691e11
#define Meter2au        6.68458712267e-12
#define Au2km  		1.49597870691e8
#define Km2au           6.68458712267e-9
#define Inch2cm		2.54

#define Small           1.0e-6
#define Very_small      1.0e-12

/* Univeral gas constant */
#define R_univ          8314.32  /* Joules/([kg mol] deg.K)   */
#define R_gram          8314.472 /* Joules/([kg mol] deg.K)   */
#define R_gurv          8.31441  /* Joules/([g mol] deg.K), Gurvich */
#define R_lewis         8.31451  /* Joules/([g mol] deg.K), Pg. 3, RP-1311 */

/* Solar flux at 1 au */
#define Solar_flux      0.1418          /* Watts/cm^2 */

/* Stefan-Boltzmann Constant                                         */
/*     previous value used was 5.67032e-12  Watts/(cm^2 deg.K^4)     */
/* http://physics.nist.gov/cgi-bin/cuu/Value?sigma|search_for=stefan */
#define Stef_Boltz      5.670400e-12 /* Watts/(cm^2 deg.K^4) */

/* Thermodynamic control parameters */
#define Rho_T           0
#define P_H             1
#define P_T             2
#define Enthalpy_only   3

/* File management return values */
#define Normal_Return   0
#define Abort		1
#define Failed_to_open  176
#define Truncated       177
#define Bad_format	178
#define Bad_check_sum   179
#define Write_protected 180

/* Turbulent transition flag values (sequence is important) */
#define Laminar         1
#define Transitional    2
#define Turbulent       3

/* Array pointers */
#define Lmnr		0
#define Turb  		1

/* Spherecone location flag values */
#define Nose            0
#define Cone            1
#define Ring            2

/* Run type values               */ 
/* Valid operands for "RunType": */
#define Run_normal      0
#define No_gravity      1
#define No_atmos        2
#define Julian_test	3
#define Top_test	4
#define Spin_test	5
#define N_run_type_list	6

#define Run_type_list {				\
	{"Run_normal",    Run_normal},		\
	{"No_gravity",    No_gravity},		\
	{"No_atmos",      No_atmos},		\
	{"Julian_test",   Julian_test}, 	\
	{"Top_test",      Top_test}, 		\
	{"Spin_test",     Spin_test}  		\
}

#define N_max_carlo_pnts	1000000
#define N_bin_max               10000
#define N_buffer		100
#define N_size			100

/* Completion codes -- Check "Run_to" structure consistency before changing */
#define Normal			1000
#define Help			1001
#define Help_tpl		1002
#define Help_hard_copy		1003
#define End_quietly		1004
#define Fiat_failed		14
#define Underflow		13
#define Converged		12
#define Post_processed          11
#define Time_out                10
#define Skip_hyperbolic         9
#define Skip_atmos              8
#define Entry_alt               7
#define Apoapsis                6
#define Dyn_700                 5
#define Mach_delay              4
#define Mach_end                3
#define Subsonic                2
#define Impact                  1
#define Failed			0	
#define Run_to_impact           -1
#define Run_to_subsonic         -2
#define Run_to_mach_end         -3
#define Run_to_mach_delay       -4
#define Run_to_dyn_700          -5
#define Run_to_apoapsis         -6
#define Run_to_entry_alt        -7

/* Beginning of warning messages for Complain */
#define Warning		800

/* File names */
#define	Traj_c		0
#define Maneuver_c	1
#define	Force_c		2
#define	Argv_c 		3
#define	Heat_c 		4
#define	Rkf45_c		5
#define	Aero_c 		6
#define	Spline_c 	7
#define	Subs_c 		8
#define	Output_c	9
#define	Bridge_c	10
#define Cmd_arg_c	11
#define Subs_fiat_c	12
#define Material_c      13
#define Fiat_c     	14
#define Io_fiat_c	15
#define Argv_fiat_c	16
#define Carlo_c         17
#define Therm_c		18
#define Grav_c		19
#define Dstr_c		20
#define Qbtnturb_cone_c 21
#define Qbtnturb_nose_c 22
#define Qdt_cone_c      23
#define Qdt_nose_c	24
#define Qdt_ring_c	25
#define Qreradia_c	26
#define Qturb_cone_c	27
#define Atmos_c		28
#define Unspecified     29
#define N_file_list     30

#define File_list {        \
	"traj.c         ", \
	"maneuver.c     ", \
	"force.c        ", \
	"argv.c         ", \
	"heat.c         ", \
	"rkf45.c        ", \
	"aero.c         ", \
	"spline.c       ", \
	"subs.c         ", \
	"output.c       ", \
	"bridge.c       ", \
	"cmd_arg.c      ", \
	"subs_fiat.c    ", \
	"material.c     ", \
	"fiat.c         ", \
	"io_fiat.c      ", \
	"argv_fiat.c    ", \
	"carlo.c        ", \
	"therm.c        ", \
	"grav.c         ", \
	"dstr.c         ", \
	"qbtnturb_cone.c", \
	"qbtnturb_nose.c", \
	"qdt_cone.c     ", \
	"qdt_nose.c     ", \
	"qdt_ring.c     ", \
	"qreradia.c     ", \
	"qturb_cone.c   ", \
	"atmos.c        ", \
	"NOT specified!!"} \

/* Input file Sufx op code list                          */
/* Tps must be 0, Bnd must be 1 and Str must be 2.       */
/* These defines must match those in fiat.h              */
#define Tps             0
#define Bnd             1
#define Str             2
#define Trj             3
#define Aer             4
#define Dsq             5
#define Dsp             6
#define Alp             7
#define Bnk             8
#define Ppc             9
#define Scr             10
#define Fat             11
#define Bpm             12
#define Ctl             13
#define Atm             14
#define Lst             15
#define Fre             16
#define Out             17
#define N_sufx_list     18

#define Sufx_list {     \
        {"tps", Tps},   \
        {"bnd", Bnd},   \
        {"str", Str},   \
        {"trj", Trj},   \
        {"aer", Aer},   \
        {"dsq", Dsq},   \
        {"dsp", Dsp},   \
        {"alp", Alp},   \
        {"bnk", Bnk},   \
        {"ppc", Ppc},   \
        {"scr", Scr},   \
        {"inp", Fat},   \
        {"bpm", Bpm},   \
        {"ctl", Ctl},   \
        {"atm", Atm},   \
        {"lst", Lst},	\
        {"fre", Fre},	\
        {"out", Out} 	\
}

/* Mollier table dimensions */
#define Ml_i    360
#define Ml_j    50

/* Planets                      */
/* Valid operands for "PLANET": */
/* The structure Planet_data must have consistency with this */
/* Planet_data is checked by Setup for consistency errors    */
#define Sun  		0
#define Mercury		1
#define Venus           2
#define Earth           3
#define Mars            4
#define Jupiter         5
#define Saturn          6
#define Uranus          7
#define Neptune         8
#define Pluto           9
#define Moon            10
#define Titan           11
#define Triton          12
#ifndef GrndTst
#define GrndTst         13
#endif
#define N_planet_list   14

#define Planet_list {                   \
	{"Sun",   	Sun},   	\
	{"Mercury", 	Mercury},	\
	{"Venus",       Venus},         \
	{"Earth",       Earth},         \
	{"Mars",        Mars},          \
	{"Jupiter",     Jupiter},       \
	{"Saturn",      Saturn},        \
	{"Uranus",      Uranus},        \
	{"Neptune",     Neptune},       \
	{"Pluto",       Pluto},		\
	{"Moon",  	Moon}, 		\
	{"Titan",       Titan},         \
	{"Triton",      Triton},	\
	{"GrndTst",     GrndTst}}

#define N_thrm_list	3
#define Hidden_dir	2
#define Thrm_list {         \
	"/opt/traj_data/",  \
	"/tmp/",            \
	"/.traj_data/"}

/* Command parsing structure */
typedef struct {
        char *mnemonic ;  /* Up to six letter mnemonic code  */
        int op_code ;     /* Corresponding code              */
} Parse ;

#define Uninit  -1
#define Air     3

/* Chemical formula structure                                            */
typedef struct {
        int species ;           /* Gas species identifier                */
        double mole_fract ;     /* mole fraction of the gas species      */
} Formula ;
#define N_species       54	/* Number of chemical species Traj supports */

/* Atmospheric gas thermodynamic data structure, equilibrium values */
typedef struct {
        int gas ;      		/* atmosphere gas type, e.g. Earth, Mars */
        int phase ;    		/* gas phase flag value                  */
	double t ;		/* temperature in deg.K                  */
	double p ;		/* pressure in Pascals                   */
	double rho ;		/* density in kg/meter^3                 */
	double h_cnv ;		/* H^0(T)-H^0(0) @ 298.15 in Gurvich     */
	double h ;		/* absolute enthalpy in Joules/kg        */
	double h_298_15 ;	/* enthalpy ref.@ 298.15 deg., Joules/kg */
        double h_d ;            /* heat of formation @ 298.15 deg., J/kg */
	double s ;     		/* entropy, Joules/(kg-mole)             */
	double cp ; 		/* specific heat Cp in Joules/(kg deg.K) */
	double gamma ;		/* ratio of specific heats               */
	double mw ;    		/* molecular weight in kg/(kg-moles)     */
	double mw_1 ;  		/* molecular weight in kg/(kg-moles)     */
	double c ; 		/* speed of sound in meters/sec          */
	double mu ;		/* viscosity in kg/(m sec)               */
	double cnd_equ ;	/* equil. conductivity in W/(m deg.K)    */
	double cnd_frz ;	/* frozen conductivity in W/(m deg.K)    */
	double prandtl ;	/* Frozen Prandtl number                 */
	double lewis_nmbr ;     /* Frozen Lewis number                   */
	Formula *frm ;		/* chemical model as mole fraction       */
} Thermo ;

/* Structure of aerodynamic coefficients [all dimensionless]                */
typedef struct {
	double L ;       /* Lift coefficient                                */
	double D ;       /* Drag coefficient                                */
	double N ;       /* Normal force coefficient                        */
	double A ;       /* Axial force coefficient                         */
	double l ;       /* Moment coefficient (x axis)                     */
	double m ;       /* Moment coefficient (y axis)                     */
	double n ;       /* Moment coefficient (z axis)                     */
	double L_alpha ; /* Lift coefficient alpha derivative coefficient   */
	double l0 ;      /* Roll roll driving moment coefficient            */
	double l_p ;     /* Roll damping derivative coefficient             */
	double m_alpha ; /* Static pitching moment derivative coefficient   */
	double m_q ;     /* Moment coefficient Q derivative coefficient     */
	double n_beta ;  /* Moment coefficient beta derivative coefficient  */
	double n_r ;     /* Moment coefficient R derivative coefficient     */
} Coeff ;

/* Maximum value of parameter alias values */
#define Max_cnv         pnt[0]	/* convective heat flux including blockage */
#define Max_rad         pnt[1]	/* radiative heat flux including blockage  */
#define Max_tot         pnt[2]	/* total heat flux including blockage      */
#define Max_cnv0        pnt[3]	/* convective heat flux without blockage   */
#define Max_rad0        pnt[4]	/* radiative heat flux without blockage    */
#define Max_tot0        pnt[5]	/* total heat flux without blockage        */
#define Max_t_wall      pnt[6]	/* adiabatic outer wall temperature        */
#define Max_t_wall_tps  pnt[7]  /* outer wall with conduction & ablation   */
#define Max_accel       pnt[8]	/* acceleration                            */
#define Max_p_dyn       pnt[9]  /* dynamic pressure                        */
#define Max_p_stg       pnt[10] /* stagnation pressure                     */
#define Min_rho         pnt[11] /* minimum freestream density              */
#define Max_rho         pnt[12] /* maximum freestream density              */
#define Min_alt         pnt[13]
#define I_max           14

/* Maximum value of parameter Structure */
typedef struct {
	double q ;      /* heat flux in Watts/cm^2                          */
	double value ;  /* actual value                                     */
	double time ;   /* time in seconds when the maximum occurred        */
	double alt ;    /* altitude in meters when the maximum occurred     */
	double v ;      /* velocity in meters/sec when the maximum occurred */
	double mach ;	/* free stream Mach number                          */
	double temp ;   /* wall temperature in deg.K                        */
	double B  ;  	/* ballistic coefficient, kg/m^2                    */
	double p ;      /* stagnation (wall) pressure in Pascals            */
	double ch ;     /* heat coefficient                                 */
} Max ;

#define Gam   		f[0]
#define V  		f[1]
#define Azi  		f[2]
#define Alt    		f[3]
#define Lat   		f[4]
#define Long  		f[5]
#define Arc   		f[6]
#define dGam		dfdt[0]
#define dV  		dfdt[1]
#define dAzi		dfdt[2]
#define dAlt		dfdt[3]
#define dLat		dfdt[4]
#define dLong		dfdt[5]
#define dArc		dfdt[6]

#ifdef Six_DoF
#define Pang		f[7]
#define Qang  		f[8]
#define Rang   		f[9]
#define Theta  		f[10]
#define Phi   		f[11]
#define Psi   		f[12]
#define Qcnv		f[13] 
#define Qrad		f[14]
#define Qtot		f[15]
#define Q_nose		f[16]
#define Q_cone		f[17]
#define Q_ring		f[18]
#define Qtot_base	f[19]
#define Mass_spall	f[20]
#define R_nose_spall	f[21]
#define R_base_spall	f[22]
#define dPang		dfdt[7]
#define dQang  		dfdt[8]
#define dRang  		dfdt[9]
#define dTheta 		dfdt[10]
#define dPhi   		dfdt[11]
#define dPsi   		dfdt[12]
#define qCnv		dfdt[13] 
#define qRad		dfdt[14]
#define qTot		dfdt[15]
#define q_Nose		dfdt[16]
#define q_Cone		dfdt[17]
#define q_Ring		dfdt[18]
#define qTot_base	dfdt[19]
#define dMass_spall	dfdt[20]
#define dR_nose_spall	dfdt[21]
#define dR_base_spall	dfdt[22]
#else
#define Qcnv		f[7] 
#define Qrad		f[8]
#define Qtot		f[9]
#define Q_nose		f[10]
#define Q_cone		f[11]
#define Q_ring		f[12]
#define Qtot_base	f[13]
#define Mass_spall	f[14]
#define R_nose_spall	f[15]
#define R_base_spall	f[16]
#define qCnv		dfdt[7] 
#define qRad		dfdt[8]
#define qTot		dfdt[9]
#define q_Nose		dfdt[10]
#define q_Cone		dfdt[11]
#define q_Ring		dfdt[12]
#define qTot_base	dfdt[13]
#define dMass_spall	dfdt[14]
#define dR_nose_spall	dfdt[15]
#define dR_base_spall	dfdt[16]
#endif

/* For use with distributed heat flux only */
#define dr_dt_spall 	dfdt[1]
#define dr_dt_ablate	dfdt[2]

/* Instant-in-time data structure */
typedef struct {
	/* Recessing geometry size parameters */
	int n_sphere ;          /* number of points modeling the nose      */
	int n_cone ;            /* number of points modeling the cone      */
	int n_sphere_cone ;     /* number of points modeling sphere-cone   */

	char desc[5] ;		/* trajectory description code letters     */
			 
	/* Time integration data */
	double time ;           /* Time from entry in seconds              */
	double julian_date ;    /* Time as a julian day number in days     */

	double f[NMAX];         /* Dependent variables, e.g. velocity      */
	double dfdt[NMAX] ;     /* First derivative (LHS)                  */
	double x[3] ;     	/* Dependent Cartesian inertial variables  */
	double dxdt[3] ; 	/* First derivative Cartesian variables    */
	double accel ;		/* sqrt(Lift^2 + Drag^2) / mass            */

	/* Aerodynamic data */
	double alpha_tot ;     	/* total angle-of-attack in radians        */
	double alpha ;		/* angle-of-attack in radians              */
	double alpha_trim ;	/* trim angle-of-attack in radians         */
	double beta ;		/* side-slip-angle in radians              */
	double bank ;		/* bank angle in radians                   */
	double roll ;           /* aerodynamic roll-angle in radians       */
	double mach ;           /* Mach number                             */
	double Re ;             /* Reynolds number                         */
	double v_atm ; 		/* atmosphere surface velocity, m/sec      */
	double v_err ;		/* velocity error from predefined values   */
	double p_dyn ; 		/* Dynamic pressure in Pascals             */
	double L_D ;            /* lift over drag ratio                    */
	double B  ;             /* ballistic coefficient, kg/m^2           */
	Coeff c ;		/* Aerodynamic coefficients                */

	/* Turbulent transition parameters */
	double theta_trns  ;
	double theta_turb ;
	double r_trns  ;
	double r_turb ;
	double dmass_dt_nose_spall ;
	double dmass_dt_cone_spall  ;
	double psi_rad_nose_trns ;

	/* Vehicle parameters that change due to ablation and spallation    */
	double mass  ;          /* vehicle total mass, kg.                  */
	double sphere_angle ;   /* dynamic nose sphere angle in radians     */
	double length ;		/* characteristic length in meters          */
	double x_center ;
	double r_nose ;		/* vehicle nose radius in meters            */
	double r_base ;		/* vehicle base radius in meters            */
	double r_base_ablate ;	/* base radius lost to ablation, meters     */
	double area ;		/* base area in meter^2                     */

        double qCnv0 ;          /* unblocked convective heat flux           */
	double qRad0 ;          /* unblocked radiative heat flux            */

	/* Thermodynamic data */
	Thermo strm  ;		/* Free stream thermodynamic state          */
	Thermo shock  ;		/* shock layer thermodynamic state          */
	Thermo edge  ;		/* Boundary layer edge thermodynamic state  */
	Thermo wall  ;		/* Stag. point wall thermodynamic state     */

	Formula frm_strm[N_species] ;
	Formula frm_wall[N_species] ;
	Formula frm_edge[N_species] ;
	Formula frm_shock[N_species] ;

	/* Orbital data */
	double r ;              /* radius from planet's center, m.          */
	double grav_r ;         /* gravity radial direction, m/sec^2        */
	double grav_phi ;	/* gravity phi direction, m/sec^2           */
	double p ;              /* semi-latus rectum in meters              */
	double a ;              /* semi-major axis in meters                */
	double e ;              /* eccentricity                             */
	double inc ;            /* inclination in radians                   */
	double asc ; 		/* ascending Node in radians                */
	double omg ; 		/* argument of Perigee in radians           */
	double nu ; 		/* true anomoly in radians                  */
	double mm ; 		/* mean anomoly in radians                  */
	double ra ;		/* right ascension, radians                 */
	double t_peri ;		/* time of periapsis in seconds             */

/* Structures common to both Fiat and Traj */
#ifdef Fiat_traj
	Var_fiat_def
#endif
} Var ;

/* Boolean switches (TRUE or FALSE) */
typedef struct {
	int help ;		/* command line control switch            */
	int null ;		/* command line control switch            */
	int number ;		/* command line control switch            */
	int last_cmd ;		/* command line control switch            */
	int tpl_first ;		/* command line control switch            */

	int above_atmos   ;     /* Is vehicle above sensible atmosphere?  */
	int above_entry ;       /* Is current altitude above entry?       */
	int aerocapture ;       /* Is Traj to solve for aerocapture?      */
	int aft_hlf_set ;       /* biconic aft half angle set             */
	int alpha_ascend ;	/* Angle-of_attack in aero. model ascends */
	int alpha_set ;		/* Is angle-of-attack set?                */
	int alpha_fixed ;	/* Is angle-of-attack fixed to one value? */
	int alpha_total_set ;	/* Is alpha set as TOTAL angle-of-attack? */
	int alpha_trim_set ;	/* Is trim angle-of-attack set?           */
	int alt_set ;		/* Is the altitude set?                   */
	int append_set ;	/* Append identifier text to file names   */
	int area_set ;		/* Is aerodynamic area set?               */
	int atm_gas_species ;	/* Does atmos. species vary with alt?     */
	int atm_rho_fct_set ;   /* Is free stream density multiplier set? */
	int atm_rho_fct_sigma_set ;  /* Is stand. deviation of above set? */
	int atmos_set ;		/* Is the Atmospheric model set?          */
	int atmos_default ;	/* Is the default Atmospheric model used? */
	int aux_ascend ;	/* Aux. parameter in aero. model ascends  */
	int B_set ;		/* True if Ballistic coef. is fixed       */ 
	int bank_set ;		/* Is Bank angle set?                     */
	int bank_rate_set ;	/* Is Bank angle rate set?                */
	int bank_accel_set ;	/* Is Bank angle acceleration set?        */
	int base_seiff ;	/* base radius based upon Al Seiff model? */
	int base_set ;		/* base radius set in command argument    */
	int beta_set ;		/* Is the side-slip angle set?            */
	int blockage_allowed ;	/* Is TPS blockage allowed?               */
	int blockag_set ;	/* Has Blockag command been set?          */
	int Cd_fct_set ;  	/* Is drag coeff., Cd multiplier set?     */
	int Cd_fct_sigma_set ;  /* Is the Cd mult. stand. deviation set?  */
	int Cm_alpha_set  ;     /* Is Cm_alpha override set?              */
	int cmd_seq ;		/* Does command line arg. have sequences? */
	int cea_direct ;        /* Use CEA thermo. routine directly?      */
	int cnv_heat_model_set;	/* stag. point conv. heating model read?  */
	int cnv_fct_set ;  	/* Is conv. heat flux multiplier set?     */
	int kepler_entry ;	/* Keplerian orbital element entry vector */
	int corner_set ;	/* corner radius set in command argument  */
	int correlation_set ;	/* Was Monte Carlo 2-D correlation set?   */
	int decel_trig_set ;	/* Was deceleration trigger value set?    */
	int dstr_set ;		/* Was distributed heating model read?    */
	int a_set ;             /* Was semi-major axis set?               */
	int e_set ;             /* Was eccentricity set?                  */
	int inc_set ;           /* Was inclination set?                   */
	int asc_set ; 		/* Was ascending Node set?                */
	int omg_set ; 		/* Was argument of Perigee set?           */
	int nu_set ; 		/* Was true anomoly set?                  */
	int half_angle_set ;	/* (forward) cone half angle set          */
	int emis_set ;		/* Emissivity selected in command arg.    */
	int event_triggered ;	/* A trajectory event was triggered.      */
	int fixed_solar_const ;	/* Assume APODES/LIESE F10.7 and Ap       */
	int post_file_atmos ;	/* Atmospheric model Is POST input file?  */
	int file_name_set[N_sufx_list] ;	/* external file name set */
	int floor_set ;		/* impact floor set?                      */
	int frame_set ;		/* Coordinate frame set                   */
	int gamma_set ;		/* Is velocity angle set?                 */
	int gamma_sigma_set ;	/* Is vel. angle standard deviation set?  */
	int geoid_set ;		/* planet model set                       */
	int plot_dstr ;		/* generate a distr. data graphics file   */
	int plot_3d ;		/* generate a Plot-3D format file         */
	int gravity_set ;	/* gravity model set                      */
	int h_max_set ;		/* Was max. allowed integ. step size set? */
	int h_set ;		/* enthalpy set in arg list?    [bridge]  */
	int h_out_set ;  	/* Output step size set                   */
	int heading_set ;	/* Is heading angle set?                  */
	int heading_sigma_set ;	/* heading angle standard deviation set?  */
	int hard_copy ;		/* log file page spacing for hard copy    */
	int heat_dstr ;		/* calculate distributed heating?         */
	int heat_method_set ;	/* surface heating methodology read?      */
	int heat_stag ;		/* calculate stagnation point heating?    */
	int hyperbolic ;	/* Is the trajectory hyperbolic?          */
	int landing_ellipse_set ;  /* Is the landing ellipse type set?    */
        int I_x_set ;          	/* Is I_x moment-of-inertia set?          */
        int I_y_set ;          	/* Is I_y moment-of-inertia set?          */
        int I_z_set ;          	/* Is I_z moment-of-inertia set?          */
        int I_xy_set ;         	/* Is I_xy product-of-inertia set?        */
        int I_xz_set ;         	/* Is I_xz product-of-inertia set?        */
        int I_yz_set ;         	/* Is I_yz product-of-inertia set?        */
        int I_prdct_set ;    	/* Is I_xy && I_xz && I_yz all set?       */ 
	int iter_max_set ;	/* Maximum number of iterestions set      */
	int jd_entry_set ;	/* Is entry time as Julian Day set?       */
	int jd_epoch_set ;	/* Is the epoch date set?                 */
	int L_D_set ;		/* True if lift/drag is fixed             */ 
	int lat_cntr_set ;	/* geocentric latitude set                */
	int lat_grph_set ;	/* geographic latitude set                */
	int length_set ;	/* characteristic length set in com. arg. */
	int long_set ;		/* longitude set in command argument      */
	int low_mach_model ;	/* low Mach # Cd aero. model available    */
	int mach_ascend ;	/* Mach number in aero. model ascends     */
	int mach_set ;		/* Mach number set                        */
	int mach_triggered ;	/* Was Mach number condition triggered?   */
	int mass_seiff ;	/* mass based upon Al Seiff model?        */
	int mass_set ;		/* Is mass set?                           */
	int mass_sigma_set ;	/* Is mass standard deviation set?        */
	int max_steps_set ;	/* maximum number of output steps set     */
	int mk_fiat_inp ;	/* Make an fiat.inp file?                 */
	int model_set ;		/* was aerodynamic model set?             */
	int monte_carlo ;	/* will this be a Monte Carlo simulation? */
	int msec ;		/* Is output time in milliseconds?        */
	int nose_set ;		/* nose radius set in command argument    */
	int number_bins_set ;	/* Is Monte Carlo number of bins set?     */
	int number_pnts_set ;	/* Is Monte Carlo number of points set?   */
	int output ;
	int output_q_load ;     /* Output to screen total heat loads      */
	int p_ang_set ;		/* Is X-axis angular acceleration set?    */
	int p_set ;		/* pressure set in arg list?    [bridge]  */
	int page_start ;	/* Start listing page at the beginning?   */
	int peak_g_set ;	/* was peak-g (deceleration) set?         */
	int photon_pressure ;	/* Include light pressure perturbations?  */
	int planet_set ;	/* Destination planet set in arg list     */
	int plot_aero_set ;	/* aerodynamic plots selected in com.arg. */
	int free_type_set ;	/* was the "free.plt" file type set?      */
	int traj_type_set ;	/* was the "traj.txt" file type set?      */
	int plot_type_set ;	/* was the "traj.plt" file type set?      */
	int post_process ;	/* Is Traj run as a postprocessor only?   */
	int post_prcss_p ;	/* Does postprocessing read-in pressure?  */
	int primed ;
	int prnt_tpl ;		/* State of Print to console TPL command? */
	int prnt_tpl_set ;	/* Was prnt_tpl set?                      */
	int prs_fct_set ;  	/* Is wall pressure multiplier set?       */
	int q_ang_set ;		/* Is Y-axis angular acceleration set?    */
	int r_ang_set ;		/* Is Z-axis angular acceleration set?    */
	int rad_fct_set ;  	/* Is rad. heat flux multiplier set?      */
	int rad_heat_model_set;	/* stag. point rad. heating model read?   */
	int Re_ascend ;		/* Reynolds number in aero. model ascends */
	int roll_ascend ;	/* Reynolds number in aero. model ascends */
	int reset_walnut ;	/* Reset initial shock strength on Walnut */
	int restart ;		/* Restart the integrator?                */
	int rho_set ;		/* density set in arg list?     [bridge]  */
	int rint_set ;          /* biconic interface radius set           */
	int roll_set ;		/* Is the vehicle roll angle set?         */
	int roll_rate_set ;	/* Is the vehicle roll rate set?          */
	int roll_accel_set ;	/* Is vehicle roll acceleration set?      */
	int rotating ;          /* Is coordinate frame rotating?          */
	int rotation_set ;	/* Planet rotation set in arg list        */
	int run_to_set ;	/* has the completion code been set?      */
	int run_type_set ;	/* has the run mode been manually set?    */
	int run_zoby ;		/* does the nose radius need correction?  */
	int s_set ;		/* entropy set in arg list?     [bridge]  */
	int skip_out_allowed ;	/* Is atmospheric skip-out allowed?       */
	int skip_out_set ;	/* Has SkipOut command been set?          */
	int tpl_start ;		/* TPL script properly started?           */
	int state_set ;		/* thermodynamic state set?     [bridge]  */
	int step_accurate ;     /* Is output time step to be accurate?    */
	int summary_set ;       /* Append final data to "summary.out"?    */
	int suppress_heading ;	/* Suppress printing heading in traj.out? */
	int t_set ;		/* temperature set in arg list? [bridge]  */
	int time_set ;		/* Is a time value set?                   */
	int t_max_set ;		/* Is the maximum allowed wall temp. set? */
	int one_heading_only ;	/* Print only top heading in traj.out?    */
	int scale_fct_set ;	/* Is the geometry factor set?            */
	int six_dof ;		/* Is a 6-Degree-of-Freedom (DoF) run?    */
	int tps_area_set ;	/* TPS area was calculated?               */
	int tps_mass_set ;	/* TPS mass was calculated?               */
	int transition_nose ;	/* Is there transition on the nose?       */
	int turbulent_nose ;	/* Is there turbulence on the nose?       */
	int transition_cone ;	/* Is there transition on the cone?       */
	int turbulent_cone ;	/* Is there turbulence on the cone?       */
	int v_inert ;		/* inertial(1) / hyp. excess(0) velocity  */
	int v_inf_set ;         /* velocity at infinity set in com. arg.? */
	int v_set ;		/* Is velocity or vel. mean value set?    */
	int v_sigma_set ;	/* Is velocity standard deviation set?    */
	int wall_model_set ;	/* was the outer wall heating model set?  */
	int decel_0_05g ;	/* has vehicle decelerated 0.05 G?        */      
        int x_cg_set ;		/* Is center-of-gravity X-coord. set?     */
        int y_cg_set ;		/* Is center-of-gravity Y-coord. set?     */
        int z_cg_set ;		/* Is center-of-gravity Z-coord. set?     */
	int x_cm_set  ;         /* Is aerodynamic center X-coord. set?    */
	int y_cm_set  ;         /* Is aerodynamic center Y-coord. set?    */
	int z_cm_set  ;         /* Is aerodynamic center Z-coord. set?    */

	int maneuver ; 		/* Are there maneuvers to do?             */

	/* Structures common to both Fiat and Traj */
#ifdef Fiat_traj
	Swt_fiat_def
#endif
} Swt ;

/* Maneuver script values and state switches */
typedef struct {
	int atmos ;		/* atmospheric model                   */
	int blockage_allowed ;	/* boolean state for TPS blockage      */
	int cnv_heat_model ;  	/* stag point convective heating model */
	int frame ;		/* coordinate frame type               */
	int dstr_type ;		/* type of distributed heating model   */
	int geoid_model ;	/* planet shape model                  */
	int gravity_model ;	/* gravity model                       */
	int heat_method ;	/* surface heating methodology         */
	int idcal ;		/* FIAT run mode                       */
	int iter_max ;		/* zero margin iteration maximum       */
	int landing_ellipse ;	/* type of landing ellipse calcuated   */
	int max_steps ;		/* maximum number of output steps      */
	int model ;		/* aerodynanic model                   */
        int number_bins ;       /* Number of bins for Monte Carlo      */
        int number_pnts ;       /* Number of points for Monte Carlo    */
	int planet ;		/* destination planet                  */
	int plt_file_type ;	/* "traj.plt" file type                */
	int fre_file_type ;	/* "free.txt" file type                */
	int trj_file_type ;	/* "traj.txt" file type                */
	int prnt_tpl ;		/* boolean state for sw->prnt_tpl      */
	int rad_heat_model ;    /* stag point convective heating model */
	int run_to ;		/* simulation completion condition     */
	int run_type ;		/* simulation run mode                 */
	int skip_out_allowed ;	/* boolean state for SkipOut           */
	int summary_type ;	/* type of summary being accumulated   */
	int type ;		/* maneuver type                       */
	int wall_model ;	/* outer wall heating model            */

	char alpha_intrp ;	/* Angle of attack interpolation meth. */
	char bank_intrp ;	/* Bank angle interpolation method     */
	char file_name[N_sufx_list][N_size] ;	/* input file name(s)  */

	double aft_hlf_angle ;	/* aft biconic half angle in degrees   */
	double alt ;		/* altitude in meters                  */
	double alpha ;     	/* "generic" angle-of-attack in deg.   */
	double alpha_trim ;    	/* trim angle-of-attack in degrees     */
	double apoapsis ;	/* aerocapture apoapsis altitude, m.   */
	double area ;		/* base area in meter^2                */
	double atm_rho_fct ;	/* free stream density multiplier      */
	double atm_rho_fct_sigma ; /* standard deviation for above     */
	double B ;		/* ballistic coefficient kg/m^2        */
	double bank ;		/* Bank angle, degrees                 */
	double bank_rate ;	/* Initial bank angle rate, deg/sec    */
	double bank_accel   ;	/* Bank angle acceleration, deg/sec^2  */
	double beta ;     	/* side-slip angle in degrees          */
	double blow ;		/* TPS ablation blowing parameter      */
	double Cd_fct ;         /* drag coefficient multiplier         */
	double Cd_fct_sigma ;   /* standard deviation for above        */
	double Cm_alpha ;	/* Cm alpha derivative overide value   */
	double cnv_fct ;        /* convective heat flux multiplier     */
	double correlation ;	/* Bivariate Normal correlation        */
	double decel_trig ;	/* Deceleration trigger value, G       */
	double a ;              /* semi-major axis in meters           */
	double e ;              /* eccentricity                        */
	double inc ;            /* inclination in radians              */
	double asc ; 		/* ascending Node in radians           */
	double omg ; 		/* argument of Perigee in radians      */
	double nu ; 		/* true anomoly in radians             */
	double floor ;         	/* impact altitude relative to geoid   */
	double gamma ;         	/* velocity angle to horixon in deg.   */
	double gamma_sigma ;   	/* standard deviation for above, deg.  */
	double h_max ;		/* max. allowed integ. stepsize, sec.  */
	double h_out ;		/* output stepsize in seconds          */
	double half_angle ;	/* cone half angle in degrees          */
	double heading ;      	/* heading angle in degrees            */
	double heading_sigma ; 	/* standard deviation for above, deg.  */
        double I_x ;            /* moment-of-inertia, kg-m^2           */
        double I_y ;            /* moment-of-inertia, kg-m^2           */
        double I_z ;            /* moment-of-inertia, kg-m^2           */
        double I_xy ;           /* product-of-inertia, kg-m^2          */
        double I_xz ;           /* product-of-inertia, kg-m^2          */
        double I_yz ;           /* product-of-inertia, kg-m^2          */
	double jd_entry ;	/* entry time as a Julian Day Number   */
	double jd_epoch ;	/* epoch as a Julian Day Number        */
	double L_D ;		/* lift over drag ratio                */
	double latitude ;	/* entry latitude in degrees           */
	double length ;		/* characteristic length in meters     */
	double longitude ;	/* entry longitude in degrees          */
	double mach ;		/* Mach number                         */
	double mass ;		/* vehicle mass in kilograms           */
	double mass_sigma ;	/* standard deviation for above, kg.   */
	double peak_g ;		/* maximum allowed deceleration in Gs  */
	double dp_ang ;		/* angular acceleration about X-axis   */
	double p_ang ;		/* angular velocity about X-axis       */
	double prs_fct ;        /* wall pressure multiplier            */
	double q_ang ;		/* angular velocity about Y-axis       */
	double r_ang ;		/* angular velocity about Z-axis       */
	double r_base ;		/* vehicle base radius in meters       */
	double r_corner ;	/* vehicle corner radius in meters     */
	double r_int ;		/* biconic cone interface radius, m.   */
	double r_nose ;		/* vehicle nose radius in meters       */
	double rad_fct ;        /* radiative heat flux multiplier      */
	double roll ;		/* vehicle roll angle in degrees       */
	double t_zero_mrgn ;	/* Bond Line fail. temperature, deg.K  */
	double t_zero_mrgn_sigma; /* stand. deviation for above, deg.K */
	double t_init ;		/* TPS entry temperature, deg.K        */
	double t_init_sigma ;	/* standard deviation for above, deg.K */
	double t_max ;		/* maximum allowed wall temperature    */
	double time ;		/* time in seconds                     */
	double scale_fct ;	/* vehicle dimensional scale factor    */
        double sigma_out ;      /* MonteCarlo output delta times sigma */
        double sigma_in ;       /* MonteCarlo input delta times sigma  */
	double thk_tps ;	/* initial TPS thickness in meters     */
	double thk_tps_sigma ;	/* standard deviation for above, m.    */
	double thk_bnd ;	/* bondlayer thickness in meters       */
	double thk_str ;	/* support structure thickness, meters */
	double tps_cp_fct ;	/* TPS material heat capacity factor   */
	double tps_cp_fct_sigma; /* standard deviation for above       */
	double tps_cnd_fct ;	/* TPS material conductivity factor    */
	double tps_cnd_fct_sigma;  /* standard deviation for above     */
	double tps_rho_fct ;	/* TPS material density factor         */
	double tps_rho_fct_sigma; /* standard deviation for above      */
	double v ;          	/* velocity in meters/sec              */
	double v_sigma ;       	/* standard deviation for above, m/sec */
	double x_cg ;		/* center-of-gravity X coordinate, m   */
	double y_cg ;		/* center-of-gravity Y coordinate, m   */
	double z_cg ;		/* center-of-gravity Z coordinate, m   */
	double x_cm ;		/* aerodynamic center X coordinate, m  */
	double y_cm ;		/* aerodynamic center Y coordinate, m  */
	double z_cm ;		/* aerodynamic center Z coordinate, m  */
	double y_thr_cpl ;	/* thermocouple location in meters     */

	Swt rd ;	        /* data read switch structure          */
} Maneuver ;

/* Aerocapture and Monte-Carlo result structure */
typedef struct {
	int completion_code ;
	int error_type ;
	int iterations ;        /* number of iterations                  */
	int sw_increasing ;     /* Is gamma increasing?                  */
	
	double accel_max ;	/* maximum deceleration                  */
	double atm_rho_fct ;	/* free stream density multiplier        */
	double bank ;		/* bank angle, degrees                   */
	double Cd_fct ;		/* drag coefficient, Cd multiplier       */
	double cos_bank ;	/* cos(bank angle)                       */
	double f[NMAX];         /* Dependent variables, e.g. velocity    */
	double mass ;           /* vehicle mass in kilograms             */
	double t_bnd_max ;	/* Maximum Bondline temperature, deg.R   */
	double t_init ;         /* TPS entry temperature, deg.R          */
	double time ;           /* time in seconds                       */
	double thk_tps ;        /* zero margin TPS thickness, meters     */
	double tps_cp_fct ;	/* TPS material heat capacity multiplier */
	double tps_cnd_fct ;	/* TPS material conductivity multiplier  */
	double tps_rho_fct ;	/* TPS material density multiplier       */
} Rslt ;

/* Jupiter atmospheric models */
#define Galileo			6

/* Structure for using free stream velocity as independent variable */
/* and finding:                                                     */
/*             Pressure on wall behind oblique shockwave            */
/*             Angle from stagnation point for turbulent transition */
typedef struct {
	double u_inf ;	/* Velocity in meters/sec */
	double f ;
} Tb_u_inf ; 

/* Structure defining radiation intensity */
typedef struct {
	double u_inf ;		/* Velocity in meters/sec */
	double aintsy[4] ;	/* was [23][4] */
	double expn[3]  ; 	/* was [23][3] */
} Rad_ints ;

/* Data for Atmosphere table */
typedef struct {
        double alt ;            /* altitude from surface in meters */
        double p ;              /* pressure in Pascals             */
        double rho ;            /* density in kg/meter^3           */
        double t ;              /* temperature in degrees K.       */
        double c ;              /* speed of sound in meters/sec    */
        double mu ;             /* viscosity in kg/(m sec)         */
} Atmos ;

/* Data for Atmospheric reference constants */
typedef struct {
        int model ;             /* Model number                        */
        int planet ;            /* Planet number                       */
        int numbr_of_entries ;  /* Number of altitude entries in model */
        Atmos *tbl ;            /* Atmospheric model table pointer     */
        double h ;              /* isentropic scale height, meters     */
        double t ;              /* reference temperature, deg.K        */
        double p ;              /* reference pressure in Pascals       */
        double rho ;            /* reference density in kg/meter^3     */
        double c ;              /* reference speed-of-sound in m/sec   */
        double mol ;            /* ref. molecular weight, kg/kg-mole   */
} Atmos_ref ;

/* Atmospheric models                                       */
/* Atmos_ref_list in "atmos.h" is checked for consistency   */
/* Valid operands for "Atmosph":                            */
#define US_1976         0
/*      External        1 */	
/*      No_atmos        2 */
#define Pathfinder      3
/*      Default         4 */
#define Viking_1        5
#define Cospar_NS       6
#define Pluat2          7
#define Venus_cospar    8
#define MER_A		9
#define MER_B		10
#define Mars_05		11
#define Titan_min	12
#ifndef GrndTst
#define GrndTst         13
#endif
#define Titan_rec	14
#define Nept_atmos	15
#define Fire            16
#define Jup_galileo	17
#define Jup_jae		18
#ifndef POST
#define POST            19
#endif
#define Venus_low_t	20
#define Venus_high_t	21
#define Venus_GRAM	22
#define Isothermal	23
#define Saturn_orton	24
#define Uranus_saic	25
#define Valles		26
#define N_atmos_list    27

#define Atmos_list {                      \
	{"US_1976",       US_1976},       \
	{"External",      External},      \
	{"No_atmos",      No_atmos},      \
	{"Pathfinder",    Pathfinder},    \
	{"Default",       Default},       \
	{"Viking_1",      Viking_1},      \
	{"Cospar_NS",     Cospar_NS},     \
	{"Pluat2",        Pluat2},        \
	{"Venus_cospar",  Venus_cospar},  \
	{"MER_A",         MER_A},         \
	{"MER_B",         MER_B},         \
	{"Mars_05",       Mars_05},       \
	{"Titan_min",     Titan_min},     \
	{"GrndTst",       GrndTst},       \
	{"Titan_rec",     Titan_rec},     \
	{"Nept_atmos",    Nept_atmos},    \
	{"Fire",  	  Fire},          \
	{"Jup_galileo",   Jup_galileo},   \
	{"Jup_jae",       Jup_jae},       \
	{"POST",          POST},          \
	{"Venus_low_t",   Venus_low_t},   \
	{"Venus_high_t",  Venus_high_t},  \
	{"Venus_GRAM",    Venus_GRAM},    \
	{"Isothermal",    Isothermal},    \
	{"Saturn_orton",  Saturn_orton},  \
	{"Uranus_saic",   Uranus_saic},   \
	{"Valles",        Valles}         \
}

/* Surface heating calculation methodology */
/* Valid operands for "Heating":           */
#define No_heating        0     /* Trajectory only without heat flux         */
#define Heat_flux_only 	  1	/* Calculate the heat flux only.             */
#define Heat_load      	  2	/* Integrate heat flux to yield heat load.   */
#define Fiat_zero_margin  3	/* Iterate Fiat to calculate zero TPS margin */
#define Fiat_post_proc    4	/* Run Fiat internally with fiat decoupled   */
#define Fiat_file_gen     5	/* Radiative equilibrium makes a Fiat file   */
#define Fiat_coupled   	  6	/* Fiat computes mat. response on inner loop */
#define N_heat_list   	  7

#define Heat_list {               			\
	{"No_heating", 		No_heating},		\
	{"Heat_flux_only", 	Heat_flux_only},	\
	{"Heat_load", 		Heat_load}, 		\
	{"Fiat_zero_margin",   	Fiat_zero_margin}, 	\
	{"Fiat_post_proc",   	Fiat_post_proc}, 	\
	{"Fiat_file_gen",   	Fiat_file_gen}, 	\
	{"Fiat_coupled", 	Fiat_coupled}}

/* Outer wall temperature/heat flux model */
/* Valid operands for "Wall":             */
#ifndef No_model
#define No_model          0
#endif
#define Cold	          1	/* Cold outer wall (wall enthalpy is zero) */
#define Fixed_temp        2	/* Outer wall temperature is a fixed value */
#define Rad_equilibrium	  3	/* Assume radiative equilibrium            */
#define Heat_flux_indep	  4	/* Heat flux independent of wall temp.     */
#define N_wall_model_list 5

#define Wall_model_list {                         \
	{"No_model", 		No_model},        \
	{"Cold", 		Cold},            \
	{"Fixed_temp",		Fixed_temp},      \
	{"Rad_equilibrium", 	Rad_equilibrium}, \
	{"Heat_flux_indep",	Heat_flux_indep}}

/* Stagnation point convective heat flux model */
/* Valid operands for "CnvHeat":               */
/*      No_model          0        Surface heating not modelled              */
#define External          1	/* Model read-in as an external files        */
#define Tauber            2	/* Mike Tauber models                        */
#define Wright            3	/* Mike Wright models                        */
#define Catalytic      	  4	/* Generic catalytic Fay-Riddel model used   */
#define Noncatalytic   	  5	/* Generic noncatalytic Fay-Riddel model     */
#define N_cnv_model_list  6

#define Cnv_model_list {              		\
	{"No_model", 		No_model},	\
	{"External", 		External},	\
	{"Tauber",		Tauber}, 	\
	{"Wright",		Wright}, 	\
	{"Catalytic", 		Catalytic}, 	\
	{"Noncatalytic",	Noncatalytic}}

/* Stagnation point radiative heat flux model */
/* Valid operands for "RadHeat":              */
/*      No_model          0        Surface heating not modelled              */
/*      External          1	   Model read-in as an external files        */
/*	Tauber            2	   Mike Tauber models                        */
/*	Wright            3	   Mike Wright models                        */
#define Default        	  4	/* Generic noncatalytic Fay-Riddel model     */
#define Obsolete   	  5	/* Obsolete model retained for legacy        */
#define N_rad_model_list  6

#define Rad_model_list {              		\
	{"No_model", 		No_model},	\
	{"External", 		External},	\
	{"Tauber",		Tauber}, 	\
	{"Wright",		Wright},	\
	{"Default",		Default},	\
	{"Obsolete",		Obsolete}}

/* Distributed heat flux model types */
/* Valid operands for "Dstrbt":      */
/*      No_model        0          Distributed heat flux NOT modelled.     */
#define Internal        1	/* Model is compiled inside Traj.          */
/*      Tauber          2	   Mike Tauber models                      */
#define Dinesh          3  	/* Dinesh Prabhu's model                   */
#define Axisym          4	/* External axisymmetric model.            */
#define Two_D      	5	/* External two dimensional model.         */
#define Three_D	        6	/* External three dimensional model.       */
#define N_dstr_list   	7

#define Dstr_list {                \
	{"No_model",	No_model}, \
	{"Internal", 	Internal}, \
	{"Tauber",	Tauber},   \
	{"Dinesh",	Dinesh},   \
	{"Axisym", 	Axisym},   \
	{"Two_D",   	Two_D},    \
	{"Three_D", 	Three_D}}

/* Gravitational models           */
/* Valid operands for "Gravity":  */
#define Dipole          0
/*      No_gravity      1 */
#define J2              2
#define J4              3
#define Oblate          4
#define N_gravity_list  5

#define Gravity_list {			\
	{"Dipole",      Dipole},        \
	{"No_gravity",  No_gravity},    \
        {"J2",      	J2},            \
        {"J4",      	J4},            \
        {"Oblate",      Oblate}}

/* Landing ellipse Monte-Carlo description */
/* Valid operands for "LandElp":           */
/*      Not_used          0    	   Do not perform a landing ellipse study    */
#define Equatorial        1	/* Vehicle trajectory only along equator     */
#define Polar             2	/* Vehicle trajectory is polar orbit based   */
#define Altitude_only     3	/* Only look at altitude variation           */
/*      Default        	  4        Do a generic/default landing ellipse      */
#define N_LandElp_list 	  5

#define LandElp_list {              		\
	{"Not_used", 		Not_used},	\
	{"Equatorial", 		Equatorial},	\
	{"Polar",		Polar}, 	\
	{"Altitude_only",	Altitude_only},	\
	{"Default",		Default}}

/* Planet "geoid" models        */
/* Valid operands for "Geoid":  */
#define No_shape        0
#define Sphere          1
#define Ellipsoid       2
#define Special         3
#define Spheroid        4
#define N_geoid_list 	5

#define Geoid_list {			\
	{"No_shape",  No_shape}, 	\
	{"Sphere",    Sphere},  	\
	{"Ellipsoid", Ellipsoid},	\
	{"Special",   Special},		\
        {"Spheroid",  Spheroid}}

/* Coordinate frames            */
/* Valid operands for "Frame":  */
#define Rotating 	0
#define Inertial	1
#define N_frame_list	2

#define Frame_list {			\
	{"Rotating",    Rotating},	\
	{"Inertial",    Inertial}}

/* Simulation completion state   */
/* Valid operands for "Run_to":  */
/* define Failed	0 */
/* define Impact 	1 */
/* define Subsonic	2 */
/* define Mach_end	3 */
/* define Mach_delay	4 */
/* define Dyn_700 	5 */
/* define Apoapsis	6 */
/* define Entry_alt	7 */
#define N_run_to_list	8

#define Run_to_list {			\
	{"Failed",     Failed},		\
	{"Impact",     Impact},		\
	{"Subsonic",   Subsonic},	\
	{"Mach_end",   Mach_end},	\
	{"Mach_delay", Mach_delay},	\
	{"Dyn_700",    Dyn_700},	\
	{"Apoapsis",   Apoapsis},	\
	{"Entry_alt",  Entry_alt}}


/* Plot file type codes           */
/* Valid operands for "PltType":  */
#define No_Plot		0
#define SharpL1		1
#define PltHeat		2
#define Show_Cm		3
#define Default		4
#define PltTauber	5
#define PltJulian	6
#define PltTop		7
#define Plt6DoF		8
#define PltStrm		9
#define PltFiat		10
#define PltCarlo	11
#define PltPost		12
#define PltXyz		13
#define PltTec		14
#define N_Plt_type	15

#define Plt_type {			\
	{"No_Plot",   No_Plot},		\
	{"SharpL1",   SharpL1},		\
	{"PltHeat",   PltHeat},		\
	{"Show_Cm",   Show_Cm},		\
	{"Default",   Default},		\
	{"PltTauber", PltTauber},	\
	{"PltJulian", PltJulian},	\
	{"PltTop",    PltTop},		\
	{"Plt6DoF",   Plt6DoF},		\
	{"PltStrm",   PltStrm},	        \
	{"PltFiat",   PltFiat},         \
	{"PltCarlo",  PltCarlo},	\
	{"PltPost",   PltPost},		\
	{"PltXyz",    PltXyz},		\
	{"PltTec",    PltTec}}

/* Operands to the "GuidLaw" guidance law op code.  */
#define Unguided       	0
#define RCS_active      1  
/*      Apollo          2 */
#define P63             3   /* do not renumber */
#define P64             4   /* do not renumber */
#define P67             5   /* do not renumber */
#define Bank_flip       6   /* do not renumber */
#define Bank_sign_chg   7   /* do not renumber */
#define Bank_cnst_alt   8   /* do not renumber */
#define Alpha_bank_mono	9
#define Alpha_bank_linr	10
#define Alpha_bank_besl	11
#define Alpha_mono	12
#define Alpha_linr	13
#define Alpha_besl	14
#define Bank_mono	15
#define Bank_linr	16
#define Bank_besl	17
#define Under_shoot	18
#define Over_shoot	19
#define Mid_gamma	20
#define N_guidlaw_list  21

#define GuidLaw_list {                        \
	{"Unguided",        Unguided},        \
	{"RCS_active",      RCS_active},      \
	{"Apollo",          Apollo},          \
	{"P63",             P63},             \
	{"P64",             P64},             \
	{"P67",             P67},             \
	{"Bank_flip",       Bank_flip},       \
	{"Bank_sign_chg",   Bank_sign_chg},   \
	{"Bank_cnst_alt",   Bank_cnst_alt},   \
	{"Alpha_bank_mono", Alpha_bank_mono}, \
	{"Alpha_bank_linr", Alpha_bank_linr}, \
	{"Alpha_bank_besl", Alpha_bank_besl}, \
	{"Alpha_mono",      Alpha_mono},      \
	{"Alpha_linr",      Alpha_linr},      \
	{"Alpha_besl",      Alpha_besl},      \
	{"Bank_mono",       Bank_mono},       \
	{"Bank_linr",       Bank_linr},       \
	{"Under_shoot",     Under_shoot},     \
	{"Over_shoot",      Over_shoot},      \
	{"Mid_gamma",       Mid_gamma}}



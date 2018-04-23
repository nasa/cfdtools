/****************************************************************/
/*                                                              */
/*                   ----- Traj Header File -----               */
/*                                                              */
/*  Author Information:                                         */
/*                                                              */
/*       Gary A. Allen, Jr.                                     */
/*       NASA Ames Research Center                              */
/*       Mail Stop 230-2                                        */
/*       Moffett Field CA 94035-1000                            */
/*                                                              */
/*       Telephone:  (650) 604-4228                             */
/*       FAX:        (650) 604-0350                             */
/*                                                              */
/*       E-mail:  gallen@mail.arc.nasa.gov                      */
/*                                                              */
/*      National Aeronautics and Space Administration           */
/*      No Copyright claimed in USA under Title 17-USC          */
/*                 All other rights reserved                    */
/*                                                              */
/*            Version:  Mk 1.10, 10 July 2006                   */
/*              Written by Gary A. Allen, Jr.                   */
/*                                                              */
/****************************************************************/
#define Rot2Inert	0
#define Inert2Rot	1

#define Con_tau		0.5655506
#define Con_del_t	1.2136524

#define EPSILON         2.2e-16
#define Gamma_limit     1.5707		/* 89.994 degrees */

/* Minimum velocity for a spacecraft out of the atmosphere in m/sec */
#define V_min		10.0

/* Physical constants */
/* Earth surface acceleration of gravity      */
/* 9.80665 m/sec is exactly 1 G by definition */
#define Grav_accel      9.80665 	/*   meters/sec^2  */

/* Listing type codes */
#define No_listing	0
#define Vel_alt		1
#define Orbital_elem	2
#define Sharp_listing	3
#define Angle_listing	4
#define Blockage	5
#define Heading_incld	6
#define Xyz		7

/* Control flags for "veh_reset" */
#define Cold_start	0
#define Reset_swt	1
#define Reset_mnv	2

/* Definitions for command line argument sequences */
#define N_seq    	7

/* "summary.out" data types */
#define No_summary_data 	0
/*      Aerocapture		1   */
/*      Peak_heat_flux		2   */
/*      Fiat_zero_margin	3   */
#define Traj_cmd_gen		4
#define Peak_conditions 	5
#define Heat_thick	  	6
#define Bal_thick	   	7
#define Bal_ACC_thick	  	8
#define V_thick	 		9
#define Vary_mass_rnose		10
#define Temp_thick		11
#define Load_thick		12
#define Lat_long  		13
#define Excel_todd		14
#define N_summary_list   	15

#define Summary_list {                             \
	{"No_summary_data",	No_summary_data},  \
	{"Aerocapture",         Aerocapture},      \
	{"Peak_heat_flux",     	Peak_heat_flux},   \
	{"Fiat_zero_margin", 	Fiat_zero_margin}, \
	{"Traj_cmd_gen",     	Traj_cmd_gen},     \
	{"Peak_conditions",    	Peak_conditions},  \
	{"Heat_thick",  	Heat_thick},       \
	{"Bal_thick",     	Bal_thick},        \
	{"Bal_ACC_thick",	Bal_ACC_thick},    \
	{"V_thick",  		V_thick},          \
	{"Vary_mass_rnose",	Vary_mass_rnose},  \
	{"Temp_thick",   	Temp_thick},       \
	{"Load_thick",   	Load_thick},       \
	{"Lat_long",     	Lat_long},         \
	{"Excel_todd",   	Excel_todd}        \
}

/* Special data types */
#define No_special_data	0
#define Aerocapture	1
#define Peak_heat_flux	2
#define Entry_Alt_Rtrn	3
#define N_special_data	4

#ifndef Not_used
#define Not_used    0
#endif 

/* Maneuver Description values   */
#define Nose        0	/* nose radius                             */
#define Veloc       1	/* velocity or mean value                  */
#define Base        2   /* base radius                             */
#define Gamma       3	/* velocity vector angle to local horizon  */
#define Mass        4	/* vehicle mass                            */
#define B_coef      5	/* Ballistic coefficient                   */
#define ScalFct     6	/* Vehicle dimensions scale factor         */
#define Area        7	/* aerodynamic area                        */
#define MassSg      8	/* vehicle mass standard deviation         */
#define Aerodyn     9	/* aerodynamic model                       */
#define Bank       10	/* bank angle                              */
#define Altitud    11	/* altitude                                */
#define Action     12	/* Action in following parameter           */
#define Alpha      13	/* angle of attack                         */
#define VelocSg    14	/* velocity standard deviation (sigma)     */
#define GammaSg    15	/* velocity angle standard deviation       */
#define Corner     16   /* corner radius                           */
#define Length     17   /* Reynolds number reference length        */
#define NumPnts    18	/* Number of points (for Monte Carlo)      */
#define HlfAngl    19	/* Cone half angle                         */
#define Afthlf     20	/* Biconic aft half angle                  */
#define Rint       21	/* Biconic junction to cones radius        */       
#define Lft_Drg    22	/* Lift over Drag ratio                    */
#define OptiStp    23   /* Output time step for use by optimizer   */
#define H_max      24   /* Maximum allowed integ. step size        */
#define PLANET     25   /* destination planet, e.g. Mars           */
#define Atmosph    26	/* Atmospheric model                       */
#define TimeStp    27   /* Output time step                        */
#define Step_mx    28   /* Maximum number of output time steps     */
#define Frame      29   /* Frame of reference, e.g. inertial       */
#define Geoid      30   /* Planet gravity equal potential surface  */
#define Gravity    31   /* Gravity model, e.g. dipole              */
#define LatGrph    32   /* Geographic latitude                     */
#define Longtud    33   /* Longitude                               */
#define Heading    34	/* Heading angle                           */
#define HedngSg    35	/* Heading angle standard deviation        */
#define File       36	/* Maneuver script file name               */
#define Heating    37	/* Select heating model.                   */
#define PltType    38	/* Plot file type made after completion?   */
#define PrntTPL    39	/* Print TPL commands to console?          */
#define Floor      40	/* Impact altitude relative to geoid.      */
#define Run_to     41	/* Simulation termination condition        */
#define Dstrbt     42   /* Use distributed heat flux model         */
#define Init_t     43   /* Initial vehicle cold soak temperature   */
#define InitTSg    44   /* Cold soak temp. standard deviation      */
#define Itr_max    45   /* Maximum allowed iterations              */
#define BndFl_t    46   /* Bondline failure temperature            */
#define BndFlSg    47   /* Bondline fail. temp., stand. deviation  */
#define Thr_cpl    48   /* Thermal couple location                 */
#define Thk_tps    49   /* Initial TPS thickness                   */
#define TkTpsSg    50   /* Init. TPS thick., standard deviation    */
#define Thk_bnd    51   /* Bond layer thickness                    */
#define Thk_str    52   /* Support structure thickness             */
#define Tps_Cnd    53   /* TPS material conductivity multiplier    */
#define TpCndSg    54   /* TPS material conductivity multiplier    */
#define Tps_Cp     55   /* TPS material heat capacity multiplier   */
#define TpsCpSg    56   /* TPS material heat capacity multiplier   */
#define TpsMode    57   /* FIAT solution mode - pyrolysis or not?  */
#define Tps_Rho    58   /* TPS material density multiplier         */
#define TpRhoSg    59   /* TPS material density multiplier         */
#define Mach       60   /* Trajectory termination Mach number      */
#define LatCntr    61   /* Geocentric latitude                     */
#define Aerocpt    62   /* Select aerocapture apoapsis altitude    */
#define AlphTrm    63	/* Trim angle of attack                    */
#define Beta       64	/* Side Slip angle                         */
#define RollAng    65   /* Roll angle                              */
#define CnvHeat    66   /* Convective stag. point heating model    */
#define RadHeat    67   /* Radiative stag. point heating model     */
#define RunType    68   /* Manually set run mode for code testing  */
#define RollRat    69   /* Set roll rate about X-axis              */
#define RollAcl    70   /* Set roll rate acceleration about X-axis */
#define T_max      71   /* Maximum allowed wall temperature        */
#define Schedul    72	/* Read control function schedule file     */
#define NumBins    73	/* Number of data bins (for Monte Carlo)   */
#define Corltn     74	/* Bivariate normal correlation            */
#define Decelrt    75	/* Deceleration trigger value              */
#define Wall       76	/* Wall temp./heating model                */
#define AtmRho     77	/* Free stream density multiplier          */
#define AtRhoSg    78	/* Free stream density mult. standard dev. */
#define Ix         79	/* Moment-of-inertia for 6-DoF             */
#define Iy         80	/* Moment-of-inertia for 6-DoF             */
#define Iz         81	/* Moment-of-inertia for 6-DoF             */
#define Ixy        82	/* Product-of-inertia for 6-DoF            */
#define Ixz        83	/* Product-of-inertia for 6-DoF            */
#define Iyz        84	/* Product-of-inertia for 6-DoF            */
#define CdFct      85	/* Drag coefficient, Cd multiplier         */
#define CdFctSg    86	/* Cd multiplier standard deviation        */
#define LandElp    87	/* Calculate landing Ellipse (Monte Carlo) */
#define P_ang      88	/* Angular accelleration about X-axis      */
#define Q_ang      89	/* Angular accelleration about Y-axis      */
#define R_ang      90	/* Angular accelleration about Z-axis      */
#define X_cm       91	/* Aerodynamic center on X-axis            */
#define Y_cm       92	/* Aerodynamic center on Y-axis            */
#define Z_cm       93	/* Aerodynamic center on Z-axis            */
#define X_Cg       94	/* Center-of-gravity on X-axis             */
#define Y_Cg       95	/* Center-of-gravity on Y-axis             */
#define Z_Cg       96	/* Center-of-gravity on Z-axis             */
#define Blowing    97	/* TPS ablation blowing parameter          */
#define Cm_alph    98	/* Angle-of-attack derivative of Cm        */
#define CnvFct     99	/* Convective Heat flux multiplier         */
#define RadFct     100	/* Radiative Heat flux multiplier          */
#define PrsFct     101  /* wall pressure multiplier                */
#define Eccent     102  /* Eccentricity                            */
#define AscNode    103  /* Ascending node                          */
#define SemiMjr    104  /* Semi-major axis                         */
#define ArgPeri    105  /* Argument of periapsis                   */
#define Inclin     106  /* Inclination                             */
#define TruAnom    107  /* True anomaly                            */
#define JdEntry    108  /* Entry time as a Julian day number       */
#define Peak_g     109  /* Maximum allowed deceleration            */
#define SkipOut    110  /* Atmospheric skip out allowed            */
#define Time       111	/* time from entry                         */
#define Summary    112  /* Maximum allowed deceleration            */
#define N_mnv_dscp 113

#define Mnv_dscp {              \
	{"Nose",    Nose},      \
	{"Veloc",   Veloc},     \
	{"Base",    Base},      \
	{"Gamma",   Gamma},     \
	{"Mass",    Mass},      \
	{"B_coef",  B_coef},    \
	{"ScalFct", ScalFct},   \
	{"Area",    Area},      \
	{"MassSg",  MassSg},    \
	{"Aerodyn", Aerodyn},   \
	{"Bank",    Bank},      \
	{"Altitud", Altitud},   \
	{"Action",  Action},    \
	{"Alpha",   Alpha},     \
	{"VelocSg", VelocSg},   \
	{"GammaSg", GammaSg},   \
	{"Corner",  Corner},    \
	{"Length",  Length},    \
	{"NumPnts", NumPnts},   \
	{"HlfAngl", HlfAngl},   \
	{"Afthlf",  Afthlf},    \
	{"Rint",    Rint},      \
	{"Lft_Drg", Lft_Drg},   \
	{"OptiStp", OptiStp},	\
	{"H_max",   H_max},	\
	{"PLANET",  PLANET},	\
	{"Atmosph", Atmosph},	\
	{"TimeStp", TimeStp},	\
	{"Step_mx", Step_mx},	\
	{"Frame",   Frame},	\
	{"Geoid",   Geoid},	\
	{"Gravity", Gravity},	\
	{"LatGrph", LatGrph},	\
	{"Longtud", Longtud},	\
	{"Heading", Heading},   \
	{"HedngSg", HedngSg},   \
	{"File",    File},      \
	{"PltType", PltType}, 	\
	{"Heating", Heating},	\
	{"PrntTPL", PrntTPL},   \
	{"Floor",   Floor},     \
	{"Run_to",  Run_to},    \
	{"Dstrbt",  Dstrbt},    \
	{"Init_t",  Init_t},    \
	{"InitTSg", InitTSg},   \
	{"Itr_max", Itr_max},   \
	{"BndFl_t", BndFl_t},   \
	{"BndFlSg", BndFlSg},   \
	{"Thr_cpl", Thr_cpl},   \
	{"Thk_tps", Thk_tps},   \
	{"TkTpsSg", TkTpsSg},   \
	{"Thk_bnd", Thk_bnd},   \
	{"Thk_str", Thk_str},   \
	{"Tps_Cnd", Tps_Cnd},   \
	{"TpCndSg", TpCndSg},   \
	{"Tps_Cp",  Tps_Cp},    \
	{"TpsCpSg", TpsCpSg},   \
	{"TpsMode", TpsMode},   \
	{"Tps_Rho", Tps_Rho},   \
	{"TpRhoSg", TpRhoSg},   \
	{"Mach",    Mach},      \
	{"LatCntr", LatCntr},	\
	{"Aerocpt", Aerocpt},	\
        {"AlphTrm", AlphTrm},   \
        {"Beta",    Beta},	\
        {"RollAng", RollAng},	\
	{"CnvHeat", CnvHeat},	\
	{"RadHeat", RadHeat},	\
	{"RunType", RunType},	\
	{"RollRat", RollRat},	\
	{"RollAcl", RollAcl},	\
	{"T_max",   T_max},	\
	{"Schedul", Schedul},   \
	{"NumBins", NumBins},   \
	{"Corltn",  Corltn},	\
	{"Decelrt", Decelrt},	\
	{"Wall",    Wall},	\
	{"AtmRho",  AtmRho},	\
	{"AtRhoSg", AtRhoSg},	\
	{"I_x",     Ix},	\
	{"I_y",     Iy},	\
	{"I_z",     Iz},	\
	{"I_xy",    Ixy},	\
	{"I_xz",    Ixz},	\
	{"I_yz",    Iyz},	\
	{"CdFct",   CdFct},	\
	{"CdFctSg", CdFctSg},	\
	{"LandElp", LandElp},	\
	{"P_ang",   P_ang},	\
	{"Q_ang",   Q_ang},	\
	{"R_ang",   R_ang},	\
	{"X_cm",    X_cm},	\
	{"Y_cm",    Y_cm},	\
	{"Z_cm",    Z_cm},	\
	{"X_Cg",    X_Cg},	\
	{"Y_Cg",    Y_Cg},	\
	{"Z_Cg",    Z_Cg},	\
	{"Blowing", Blowing},	\
	{"Cm_alph", Cm_alph},	\
	{"CnvFct",  CnvFct},	\
	{"RadFct",  RadFct},	\
	{"PrsFct",  PrsFct},	\
	{"Eccent",  Eccent},	\
	{"AscNode", AscNode},	\
	{"SemiMjr", SemiMjr},	\
	{"ArgPeri", ArgPeri},	\
	{"Inclin",  Inclin},	\
	{"TruAnom", TruAnom},	\
	{"JdEntry", JdEntry},	\
	{"Peak_g",  Peak_g},	\
	{"SkipOut", SkipOut},	\
	{"Time",    Time},      \
	{"Summary", Summary}}

/* Vehicle models -- Be sure to update                            */
/* Chck_data_base and Dflt_data_base are checked for consistency  */                                  
#ifndef No_model
#define No_model	 0
#endif 
#define Sphere		 1
#define AFE		 2
#define Pathfinder	 3
#define Sphere_sect	 4
#define Sphere_cone	 5
#define Biconic		 6
#define Viking_seiff	 7
#define Galileo_seiff	 8
#define Galileo_jae	 9
#define Pioneer_seiff	10
#define X38		11
#define Viking		12
#ifndef GrndTst
#define GrndTst		13
#endif 
#define Viking_air	14
#define AMaRV		15
#define Fire		16
#define Apollo		17
#define CEV   		18
#ifndef POST
#define POST            19
#endif 
#define PII		20
#define Julian		21
#define External_aero	22  /* File type marker, duplicating number below */
#define Ma_Alpha_Cl_Cd	22
#define Ma_Alpha_Cn_Ca	23
#define Ma_Alpha_Ca_Cn	24
#define Ma_Alpha_logRe	25
#define Ma_Alpha_Re	26
#define Ma_logRe_Cd	27
#define Ma_Re_Cd	28
#define Ma_Cd		29
#define N_model_list    30

#define Model_list {                        \
	{"No_model",       No_model},       \
	{"Sphere",         Sphere},         \
	{"AFE",            AFE},            \
	{"Pathfinder",     Pathfinder},     \
	{"Sphere_sect",    Sphere_sect},    \
	{"Sphere_cone",    Sphere_cone},    \
	{"Biconic",        Biconic},        \
	{"Viking_seiff",   Viking_seiff},   \
	{"Galileo_seiff",  Galileo_seiff},  \
	{"Galileo_jae",    Galileo_jae},    \
	{"Pioneer_seiff",  Pioneer_seiff},  \
	{"X38",            X38},            \
	{"Viking",         Viking},         \
	{"GrndTst",        GrndTst},        \
	{"Viking_air",     Viking_air},     \
	{"AMaRV",          AMaRV},          \
	{"Fire",           Fire},           \
	{"Apollo",         Apollo},         \
	{"CEV",            CEV},            \
	{"POST",           POST},           \
	{"PII",            PII},            \
	{"Julian",         Julian},         \
	{"Ma_Alpha_Cl_Cd", Ma_Alpha_Cl_Cd}, \
	{"Ma_Alpha_Cn_Ca", Ma_Alpha_Cn_Ca}, \
	{"Ma_Alpha_Ca_Cn", Ma_Alpha_Ca_Cn}, \
	{"Ma_Alpha_logRe", Ma_Alpha_logRe}, \
	{"Ma_Alpha_Re",    Ma_Alpha_Re},    \
	{"Ma_logRe_Cd",    Ma_logRe_Cd},    \
	{"Ma_Re_Cd",       Ma_Re_Cd},       \
	{"Ma_Cd",          Ma_Cd}           \
}

/*      TRUE     1                                              */   
#define Reqrd	 1	/* Required parameter                   */
/*      FALSE    0                                              */
#define Unusd    0 	/* Parameter unused by this geometry    */
#define Deriv	-1	/* Parameter can be derived from others */ 
#define AerDy	-2	/* Aerodynamic parameter                */
#define Specl   -3 	/* Special parameter                    */
#define Fixed   -4 	/* Fixed default value                  */
#define Chgbl   -5 	/* Changeable default value             */
#define Combn   -6 	/* Combined influence parameter         */
#define Vacnt   -7 	/* Vacant parameter                     */

/* Vehicle model description */
#define Chck_data_base {                                                \
{"No", /* no aerodynamic model */                                       \
 Unusd, Unusd, Reqrd, Specl, Unusd, Specl, Unusd, Unusd, Unusd, Unusd,  \
 Reqrd, Reqrd, Reqrd, Unusd, FALSE, Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 No_model},                                                             \
{"Sphere based upon experimental data",                                 \
 Unusd, Unusd, Deriv, Specl, Unusd, Deriv, Unusd, Unusd, Unusd, Unusd,  \
 Deriv, Reqrd, Deriv, Unusd, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Sphere},                                                               \
{"Aeroassisted Flight Experiment (AFE)",                                \
 Unusd, AerDy, Chgbl, Specl, AerDy, Chgbl, Unusd, Unusd, Unusd, Specl,  \
 Chgbl, Chgbl, Chgbl, Unusd, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 AFE},                                                                  \
{"Pathfinder Aeroshell [Martian Atmosphere type]",                      \
 Unusd, AerDy, Chgbl, Specl, AerDy, Chgbl, Chgbl, Fixed, Reqrd, Specl,  \
 Chgbl, Chgbl, Chgbl, Unusd, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Pathfinder},                                                           \
{"Modified Newtonian spherical section",                                \
 Unusd, AerDy, Deriv, Specl, AerDy, Deriv, Unusd, Unusd, Unusd, Specl,  \
 Deriv, Reqrd, Reqrd, Unusd, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Sphere_sect},                                                          \
{"Modified Newtonian Sphere-Cone",                                      \
 Unusd, AerDy, Deriv, Specl, AerDy, Deriv, Unusd, Deriv, Specl, Specl,  \
 Deriv, Reqrd, Reqrd, Unusd, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Sphere_cone},                                                          \
{"Modified Newtonian Axisymmetric Biconic",                             \
 Deriv, AerDy, Deriv, Specl, AerDy, Deriv, Unusd, Deriv, Unusd, Specl,  \
 Specl, Reqrd, Reqrd, Deriv, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Biconic},                                                              \
{"Viking Lander Aeroshell (Al Seiff version)",                          \
 Unusd, AerDy, Chgbl, Specl, AerDy, Chgbl, Chgbl, Fixed, Unusd, Specl,  \
 Chgbl, Chgbl, Chgbl, Unusd, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Viking_seiff},                                                         \
{"Galileo Aeroshell (Al Seiff version)",                                \
 Unusd, AerDy, Chgbl, Specl, AerDy, Chgbl, Chgbl, Fixed, Unusd, Specl,  \
 Chgbl, Chgbl, Chgbl, Unusd, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Galileo_seiff},                                                        \
{"Galileo Aeroshell (Mike Tauber / JAE version)",                       \
 Unusd, AerDy, Chgbl, Specl, AerDy, Chgbl, Chgbl, Fixed, Unusd, Specl,  \
 Chgbl, Chgbl, Chgbl, Unusd, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Galileo_jae},                                                          \
{"Pioneer Venus Aeroshell (Al Seiff version)",                          \
 Unusd, AerDy, Chgbl, Specl, AerDy, Chgbl, Chgbl, Fixed, Unusd, Specl,  \
 Chgbl, Chgbl, Chgbl, Unusd, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Pioneer_seiff},                                                        \
{"X-38 hypersonic full stall",                                          \
 Unusd, AerDy, Chgbl, Specl, AerDy, Chgbl, Unusd, Unusd, Unusd, Specl,  \
 Reqrd, Chgbl, Chgbl, Unusd, FALSE, Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 X38},                                                                  \
{"Viking Lander Aeroshell [Martian atmosphere type]",                   \
 Unusd, AerDy, Chgbl, Specl, AerDy, Chgbl, Chgbl, Fixed, Unusd, Specl,  \
 Deriv, Chgbl, Chgbl, Unusd, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Viking},                                                               \
{"Ground Test Facility Model",                                          \
 Fixed, AerDy, Chgbl, Specl, AerDy, Chgbl, Unusd, Fixed, Unusd, Specl,  \
 Specl, Chgbl, Chgbl, Chgbl, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 GrndTst},                                                              \
{"Viking Lander Aeroshell [Terrestrial atmoshere type]",                \
 Unusd, AerDy, Chgbl, Specl, AerDy, Chgbl, Chgbl, Fixed, Unusd, Specl,  \
 Deriv, Chgbl, Chgbl, Unusd, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Viking_air},                                                           \
{"Advanced Maneuverable Reentry Vehicle (AMaRV)",                       \
 Fixed, AerDy, Chgbl, Specl, AerDy, Chgbl, Unusd, Fixed, Unusd, Specl,  \
 Specl, Chgbl, Chgbl, Chgbl, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 AMaRV},                                                                \
{"Fire Flight Aeroshell",                                               \
 Fixed, AerDy, Chgbl, Specl, AerDy, Chgbl, Unusd, Fixed, Unusd, Specl,  \
 Specl, Chgbl, Chgbl, Chgbl, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Fire},                                                                 \
{"Apollo Command Module",                                               \
 Unusd, AerDy, Chgbl, Specl, AerDy, Chgbl, Chgbl, Fixed, Unusd, Specl,  \
 Deriv, Chgbl, Chgbl, Unusd, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Apollo},                                                               \
{"CEV Command Module",                                                  \
 Unusd, AerDy, Chgbl, Specl, AerDy, Chgbl, Chgbl, Fixed, Unusd, Specl,  \
 Deriv, Chgbl, Chgbl, Unusd, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 CEV},                                                                  \
{"POST format external model            ",                              \
 Unusd, Reqrd, Reqrd, Specl, Reqrd, Specl, Unusd, Unusd, Unusd, Specl,  \
 Reqrd, Reqrd, Reqrd, Unusd, FALSE, Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 POST},                                                                 \
{"Pershing-II Reentry Vehicle",                                         \
 Fixed, AerDy, Chgbl, Specl, AerDy, Chgbl, Unusd, Fixed, Unusd, Specl,  \
 Specl, Chgbl, Chgbl, Chgbl, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 PII},                                                                  \
{"Julian Allen Model for code validation",                              \
 Fixed, AerDy, Chgbl, Specl, AerDy, Chgbl, Unusd, Fixed, Unusd, Specl,  \
 Specl, Chgbl, Chgbl, Chgbl, TRUE,  Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Julian},                                                               \
{"Cd,Cl = f(Mach Number, Alpha)         ",                              \
 Unusd, Reqrd, Reqrd, Specl, Reqrd, Specl, Unusd, Unusd, Unusd, Specl,  \
 Reqrd, Reqrd, Reqrd, Unusd, FALSE, Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Ma_Alpha_Cl_Cd},                                                       \
{"Cn,Ca = f(Mach Number, Alpha)         ",                              \
 Unusd, Reqrd, Reqrd, Specl, Reqrd, Specl, Unusd, Unusd, Unusd, Specl,  \
 Reqrd, Reqrd, Reqrd, Unusd, FALSE, Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Ma_Alpha_Cn_Ca},                                                       \
{"Ca,Cn = f(Mach Number, Alpha)         ",                              \
 Unusd, Reqrd, Reqrd, Specl, Reqrd, Specl, Unusd, Unusd, Unusd, Specl,  \
 Reqrd, Reqrd, Reqrd, Unusd, FALSE, Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Ma_Alpha_Ca_Cn},                                                       \
{"Cd,Cl = f(Mach Number, Alpha, log(Re))",                              \
 Unusd, Reqrd, Reqrd, Specl, Reqrd, Specl, Unusd, Unusd, Unusd, Specl,  \
 Reqrd, Reqrd, Reqrd, Unusd, FALSE, Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Ma_Alpha_logRe},                                                       \
{"Cd,Cl = f(Mach Number, Alpha, Re)     ",                              \
 Unusd, Reqrd, Reqrd, Specl, Reqrd, Specl, Unusd, Unusd, Unusd, Specl,  \
 Reqrd, Reqrd, Reqrd, Unusd, FALSE, Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Ma_Alpha_Re},                                                          \
{"Cd = f(Mach Number, log(Re))          ",                              \
 Unusd, Unusd, Reqrd, Specl, Unusd, Specl, Unusd, Unusd, Unusd, Unusd,  \
 Reqrd, Reqrd, Reqrd, Unusd, FALSE, Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Ma_logRe_Cd},                                                          \
{"Cd = f(Mach Number, Re)               ",                              \
 Unusd, Unusd, Reqrd, Specl, Unusd, Specl, Unusd, Unusd, Unusd, Unusd,  \
 Reqrd, Reqrd, Reqrd, Unusd, FALSE, Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Ma_Re_Cd},                                                             \
{"Cd = f(Mach Number)                   ",                              \
 Unusd, Unusd, Reqrd, Specl, Unusd, Specl, Unusd, Unusd, Unusd, Unusd,  \
 Reqrd, Reqrd, Reqrd, Unusd, FALSE, Vacnt, Vacnt, Vacnt, Vacnt, Vacnt,  \
 Ma_Cd}}

/* Vehicle parameters check list structure */
typedef struct {
	char *name ;	/* Vehicle name / aerodynamic description   0 */
	int aft_hlf ; 	/* biconic aft half angle                   1 */
	int alpha ;	/* angle-of-attack useable                  2 */
	int area ;	/* aerodynamic characteristic area          3 */
	int B ;		/* ballistic coefficient                    4 */
	int bank ;	/* Bank angle                               5 */
	int base ;	/* base radius                              6 */
	int corner ;	/* corner radius                            7 */
	int half_angle;	/* (forward) cone half angle                8 */
	int heat_dstr ;	/* Distributed heating                      9 */
	int L_D ;	/* Lift-over-drag ratio                    10 */
	int length ;	/* characteristic length (Reynolds number) 11 */
	int mass ;	/* mass set                                12 */
	int nose ;	/* nose radius                             13 */
	int rint ;      /* biconic interface radius                14 */
	int round ;     /* is vehicle base area round?             15 */
	int unused_1  ; /* place holder                            16 */
	int unused_2  ; /* place holder                            17 */
	int unused_3  ; /* place holder                            18 */
	int unused_4  ; /* place holder                            19 */
	int unused_5  ; /* place holder                            20 */
	int model ;	/* Aerodynamic model name for consistency check */
} Chck ;

/* Default vehicle parameter database */
#define Dflt_data_base {                                                \
/* No Model                                         */ {No_model,       \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,      \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000},	\
/* Sphere based upon experimental data              */ {Sphere,         \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,      \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000},     \
/* AFE - Aeroassisted Flight Experiment             */ {AFE,            \
 1057.8449320,   0.0000000,   0.0000000,  36.7000000, 398000.0000,      \
   18.3500000,   0.0000000,   0.0000000,  20.0000000,   0.0000000},     \
/* Pathfinder - Mars 1997 lander Aeroshell          */ {Pathfinder,     \
    5.5154586,   0.0000000,  70.0000000,   2.6500000, 585.3000000,      \
    1.3250000,   0.0662500,   0.0000000,   0.6638000,   0.0000000},     \
/* Newtonian spherical section                      */ {Sphere_sect,    \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,      \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000},     \
/* Newtonian sphere-cone                            */ {Sphere_cone,    \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,      \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000},     \
/* Newtonian biconic                                */ {Biconic,        \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,      \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000},     \
/* Viking based upon Al Seiff's aerodynamics        */ {Viking_seiff,   \
    9.8408720,   0.0000000,  70.0000000,   3.5397440, 980.0000000,      \
    1.7698720,   0.0254000,   0.0000000,   0.8763000,   0.0000000},     \
/* Galileo based upon Al Seiff's aerodynamics       */ {Galileo_seiff,  \
    1.2566151,   0.0000000,  44.8600000,   1.2649000, 338.9300000,      \
    0.6324500,   0.0000001,   0.0000000,   0.2220000,   0.0000000},     \
/* Galileo based on Mike Tauber's JAE program       */ {Galileo_jae,    \
    1.2566151,   0.0000000,  44.8600000,   1.2649000, 338.9300000,      \
    0.6324500,   0.0000001,   0.0000000,   0.2220000,   0.0000000},     \
/* Pioneer Venus Large Probe, Seiff's aerodynamics  */ {Pioneer_seiff,  \
    1.58903,     0.0000000,  45.0000000,   1.4223900, 316.4830000,      \
    0.7111950,   0.0000001,   0.0000000,   0.3632200,   0.0000000},     \
/* X-38, a winged entry vehicle                     */ {X38,            \
   21.6720000,   0.0000000,   0.0000000,   8.4125000, 9344.000000,      \
    0.0000000,   0.0000000,   0.0000000,   0.3048000,   0.0000000},     \
/* Viking - Mars lander Aeroshell, Mars atmosphere  */ {Viking,         \
    9.8408720,   0.0000000,  70.0000000,   3.5000000, 980.0000000,      \
    1.7698720,   0.0254000,   0.0000000,   0.8763000,   0.0000000},     \
/* Ground Test Facility Model                       */ {GrndTst,        \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,      \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000},     \
/* Viking - Mars 1976 lander Aeroshell based on air */ {Viking_air,     \
    9.8408720,   0.0000000,  70.0000000,   3.5000000, 980.0000000,      \
    1.7698720,   0.0254000,   0.0000000,   0.8763000,   0.0000000},     \
/* AMaRV - Advanced Maneuverable Reentry Vehicle    */ {AMaRV,          \
    0.2680483,   6.0000000,  10.4000000,   2.0792694, 470.0000000,      \
    0.2921000,   0.0000000,   0.1460500,   0.0233680,   0.0000000},     \
/* Fire Flight aeroshell in air                     */ {Fire,           \
    0.3543,      0.0000000,   0.0,         0.6716,     86.568,          \
    0.3358,      0.0102,      0.0,         0.9347,      0.0000000},     \
/* Apollo Command Module in air (AS-501)            */ {Apollo,         \
   12.07958,     0.0,         0.0,         3.92176,    5424.51,         \
    1.96088,     0.1778,      0.0,         4.6609,      0.0000000},     \
/* CEV Command Module in air                        */ {CEV,            \
   79.45984,     0.0000000,   0.0,         5.0292,     7418.0,          \
    2.5146,      0.2750,      0.0,         6.03504,     0.0000000},     \
/* POST External File                               */ {POST,           \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,      \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000},     \
/* Pershing-II Reentry Vehicle                      */ {PII,            \
    0.2680483,   6.0000000,  10.4000000,   2.0792694, 470.0000000,      \
    0.2921000,   0.0000000,   0.1460500,   0.0233680,   0.0000000},     \
/* Julian - Dummy code validation aerodynamic model */ {Julian,         \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,      \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000},     \
/* External file, Cd,Cl = f(Mach #, Alpha)          */ {Ma_Alpha_Cl_Cd, \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,      \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000},     \
/* External file, Cn,Ca = f(Mach #, Alpha)          */ {Ma_Alpha_Cn_Ca, \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,      \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000},     \
/* External file, Ca,Cn = f(Mach #, Alpha)          */ {Ma_Alpha_Ca_Cn, \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,      \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000},     \
/* External file, Cd,Cl = f(Mach #, Alpha, log(Re)) */ {Ma_Alpha_logRe, \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,      \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000},     \
/* External file, Cd,Cl = f(Mach #, Alpha, Re)      */ {Ma_Alpha_Re,    \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,      \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000},     \
/* External file, Cd = f(Mach #, log(Re))           */ {Ma_logRe_Cd,    \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,      \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000},     \
/* External file, Cd = f(Mach #, Re)                */ {Ma_Re_Cd,       \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,      \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000},     \
/* External file, Cd = f(Mach #)                    */ {Ma_Cd,          \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000,      \
    0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.0000000}}

/* Definitions for POST format input file */
#define N_comma    	50
#define Not_used        0
#define Atemt        	1
#define Prest         	2
#define Cst            	3
#define Denst          	4
#define Cdt		5
#define Clt		6
#define N_post_cmd_list	7

#define Post_cmd_list {            \
	{"Not_used",    Not_used}, \
	{"atemt",       Atemt},    \
	{"prest",       Prest},    \
	{"cst",         Cst},      \
	{"denst",       Denst},    \
	{"cdt",		Cdt},      \
	{"clt",         Clt}}

/* Command line argument sequence struction */
typedef struct {
	int n_seq ;		/* number of elements in seq[]           */
	int type ;		/* command type, e.g. Gamma = 3          */
	double vct[N_comma] ;	/* command sequence vector               */ 
} Seq ;

/* Default value structure */
typedef struct {
	int model ;		/* Aerodynamic model name                */
	double area ;		/* base area in meter^2                  */
	double aft_hlf_angle ;	/* aft biconic half angle in degrees     */
	double half_angle ;	/* cone half angle in degrees            */
	double length ;		/* characteristic length in meters       */
	double mass ;		/* mass in kilogram                      */
	double r_base ;		/* vehicle base radius in meters         */
	double r_corner ;	/* vehicle corner radius in meters       */
	double r_int ;		/* biconic cone interface radius, m.     */
	double r_nose ;		/* vehicle nose radius in meters         */
	double unused    ;	/* place holder                          */
} Dflt ;

/* Integrator and Maneuver control */
#define N_scripts	10
typedef struct {
        int number_bins ;         /* Number of bins for Monte Carlo      */
        int number_pnts ;         /* Number of points for Monte Carlo    */
	int max_steps ;           /* Maximum number of output steps      */
	int pnt_script ;          /* Script pointer < N_scripts          */

	int mnv_max[N_scripts] ;  /* Number of maneuvers for script      */
	int mnv_pnt[N_scripts] ;  /* Pointer to line in script           */

        double correlation ;      /* MonteCarlo correlation 2-D          */
	double atm_rho_fct_sigma; /* free stream dens. mult. stand. dev. */
	double Cd_fct_sigma ;     /* Cd coefficient mult. standard dev.  */
        double gamma_sigma ;      /* vel. angle standard deviation, rad. */
        double heading_sigma ;    /* heading angle stand. dev., rad.     */
        double mass_sigma ;       /* vehicle mass stand. dev., kg        */
	double h ;                /* Integrator step size in seconds     */
	double h_min ;            /* minimum integ. step size, sec.      */
	double h_max ;            /* maximum integ. step size, sec.      */
	double h_out ;            /* output step size in seconds         */
	double h_trig ;           /* step size to trigger event, sec.    */
        double v_sigma ;          /* vel. standard deviation, meters/sec */

	Maneuver *mnv_addr[N_scripts]; /* Maneuver script base address    */

	/* Structures common to both Fiat and Traj */
#ifdef Fiat_traj
	Control_fiat_def
#endif
} Control ;

/* Command line constants */
typedef struct {
	int atmos ;		/* atmospheric model                    */
	int geoid_model ;	/* planet shape model                   */
	int gravity_model ;	/* gravitational model                  */
	int planet ;		/* destination planet, 3=Earth, 4=Mars  */
	int frame ;		/* coordinate frame                     */
	int iter_max ;		/* maximum allowed iterations           */
	int state ;		/* thermodynamic state                  */
	int summary_type ;	/* "summary.out" file data type         */
	int third_body ;        /* most significant third-body grav.    */
	double half_angle ;	/* cone half angle at entry, degrees    */
	double atm_rho_fct ;	/* Atmospheric density variation factor */
	double alt ;		/* altitude in kilometers               */
	double v ;		/* inertial entry velocity, meters/sec  */
	double v_inf ;		/* hyperbolic velocity @ inf., m./sec   */
	double alpha ;		/* angle-of-attack in degrees           */
	double bank ;		/* bank angle in degrees                */
	double dbank_dt ;	/* coning angle rate in deg/sec         */
	double floor ;		/* impact altitude (m.) rel. to geoid   */
        double gamma ;          /* velocity angle to horixon in deg.    */
        double heading ;        /* heading angle in degrees             */
        double jd_entry ;       /* entry time, Julian Day Number, days  */
	double jd_epoch ;       /* epoch date as a Julian Day Number    */
        double latitude ;       /* areographic latitude in degrees      */
        double longitude ;      /* areographic longitude in degrees     */
	double B  ;             /* ballistic coefficient, kg/m^2        */
	double L_D ;            /* lift over drag ratio                 */
	double Cd_fct ;    	/* Cd coefficient multiplier            */
	double p ;              /* semi-latus rectum in meters          */
        double a ;              /* semi-major axis in meters            */
	double e ;              /* eccentricity                         */
	double inc ;            /* inclination in radians               */
	double asc ;            /* ascending Node in radians            */
	double omg ;            /* argument of Perigee in radians       */
	double mm ; 		/* mean anomaly in radians              */
	double nu ;             /* true anomoly in radians              */
	double dp_ang ;		/* X-axis angular accel. in rads/sec^2  */
	double p_ang ;		/* X-axis angular velocity in rads/sec  */
	double q_ang ;		/* Y-axis angular velocity in rads/sec  */
	double r_ang ;		/* Z-axis angular velocity in rads/sec  */
	double r_peri ;		/* periapsis in meters                  */
	double t_peri ;         /* time of periapsis in seconds         */
	Thermo *strm  ;          /* Free stream thermodynamic state     */

	double mass ;		/* Vehicle mass at entry, kg.           */
	double length ;		/* characteristic length in meters      */
	double r_nose ;		/* Nose radius at entry, meters         */
	double r_base ;		/* Base radius at entry, meters         */
	double area ;		/* base area in meter^2                 */
	double scale_fct ;	/* vehicle geometry scale factor        */
 
	/* Structures common to both Fiat and Traj */
#ifdef Fiat_traj
	Entry_fiat_def
#endif
} Entry ;

/* External file aerodynamics Structure */
typedef struct {
	double mach ;		/* Mach number, indep. paramter            */
	double alpha ;		/* angle-of-attack, indep. parameter       */
	double roll ;		/* aero. roll angle, indep. parameter      */
	double logRe ;		/* log10(Reynold Number), indep. parameter */
	double Ca ;		/* Axial coefficient, dependent parameter  */
	double Cn ;		/* Normal coefficient, dependent parameter */
	double Cd ;		/* Drag coefficient, dependent parameter   */
	double Cl ;		/* Lift coefficient, dependent parameter   */
	double Cm ;		/* Moment coefficient, dependent parameter */
	double Ch ;		/* Hinge moment coef., dependent parameter */
	double dCa_mach ;	/* local derivative dCa/dmach              */
	double dCn_mach ;	/* local derivative dCn/dmach              */
	double dCd_mach ;	/* local derivative dCd/dmach              */
	double dCl_mach ;	/* local derivative dCl/dmach              */
	double dCm_mach ;	/* local derivative dCm/dmach              */
	double dCh_mach ;	/* local derivative dCh/dmach              */
	double dCa_alpha ;	/* local derivative dCa/dalpha             */
	double dCn_alpha ;	/* local derivative dCn/dalpha             */
	double dCd_alpha ;	/* local derivative dCd/dalpha             */
	double dCl_alpha ;	/* local derivative dCl/dalpha             */
	double dCm_alpha ;	/* local derivative dCm/dalpha             */
	double dCh_alpha ;	/* local derivative dCm/dalpha             */
	double dCa_logRe ;	/* local derivative dCa/dlogRe             */
	double dCn_logRe ;	/* local derivative dCn/dlogRe             */
	double dCd_logRe ;	/* local derivative dCd/dlogRe             */
	double dCl_logRe ;	/* local derivative dCl/dlogRe             */
	double dCm_logRe ;	/* local derivative dCm/dlogRe             */
	double dCh_logRe ;	/* local derivative dCh/dlogRe             */
	double dCa_roll ;	/* local derivative dCa/droll              */
	double dCn_roll ;	/* local derivative dCn/droll              */
	double dCd_roll ;	/* local derivative dCd/droll              */
	double dCl_roll ;	/* local derivative dCl/droll              */
	double dCm_roll ;	/* local derivative dCm/droll              */
	double dCh_roll ;	/* local derivative dCh/droll              */
} Aero ;

/* External trajectory envelope structure */
typedef struct {
	double v ;      /* envelope velocity in meters/sec */
	double alt ;    /* envelope altitude in meters     */
} Envel ;

/* Lookup table structure for Seiff aerodynamic model */
typedef struct {
	double alt ;            /* altitude surface in meters        */
	double v ;              /* velocity, meters/sec              */
	double mass ;           /* vehicle mass in kilograms         */
	double r_base ;         /* base radius, meters               */
	double Cd ;             /* drag coefficient                  */
} Seiff ;

/* Structure for distributed heat flux tables */
typedef struct {
        double x ;              /* X coord. of surface point in meters  */
        double y ;              /* Y coord. of surface point in meters  */
        double z ;              /* Z coord. of surface point in meters  */
        double s ;              /* distance from Stag. point in meters  */
        double q_coef ;         /* heat flux, Leading Coefficient MKS   */
        double q_rho_exp ;      /* heat flux, density (kg/m^3) exponent */
        double q_v_exp ;        /* heat flux, velocity (m/sec) exponent */
        double p_coef ;         /* pressure, leading Coefficient MKS    */
        double p_rho_exp ;      /* pressure, density (kg/m^3) exponent  */
        double p_v_exp ;        /* pressure, velocity (m/sec) exponent  */
} Dstr ;

/* Distributed heat flux and wall pressure output */
typedef struct {
	int flow_state ;        /* was "kk" */ 

        double x ;              /* X coord. of surface point in meters  */
        double y ;              /* Y coord. of surface point in meters  */
        double z ;              /* Z coord. of surface point in meters  */
	double r_cyl ;		/* radius from vehicle axis of symmetry */
	double rb ;
        double s ;              /* distance from Stag. point in meters  */
	double theta ; 

        double thk_tps ; 	/* TPS thickness in meters              */
	double t_bnd_max ;	/* Maximum bond line temperature, deg.K */
	double emis ;		/* Heat shield emissivity               */
        double f[NMAX];         /* x[] = heat load                      */
        double dfdt[NMAX] ;     /* First derivative of x[], heat flux   */

	double qCnv0 ;          /* unblocked convective heat flux       */
	double qRad0 ;          /* unblocked radiative heat flux        */

	Thermo wall ;		/* Local wall thermodynamic state       */
} Dist ;

/* Constants for describing the vehicle */
typedef struct {
	int completion_code ;	/* code returned at simulation end       */
	int dstr_type ;		/* specifies distributed heat flux model */
	int error_source ;	/* error source file for epitaph         */
	int error_type ;	/* error type used for epitaph           */
	int error_int ;		/* error specific integer for epitaph    */
	int heat_method ;	/* surface heating analysis method       */
	int cnv_heat_model ;	/* stag. point convective heating model  */
	int rad_heat_model ;	/* stag. point radiative heating model   */
	int il ;		/* number of Angle of Attack splines     */
	int jl ;		/* number of Mach Number splines         */
	int i_seq ;		/* argument list sequence vector index   */
	int landing_ellipse ;	/* specifies landing ellipse solution    */
	int listing_type ;	/* specifies output listing type         */
	int model ; 		/* specifies vehicle model               */
	int n_alpha ;		/* number of angle-of-attack fit points  */
	int n_alpha_schd   ;	/* number of points for alpha schedule   */
	int n_bank_schd   ;	/* number of points for bank schedule    */
	int n_continuum_aero ;  /* number of entries in cust. aero.      */
	int n_free_mole_aero ;	/* number of entries in free mole. aero. */
	int n_alpha_free ;	/* # free-mole. angle-of-attack points   */
	int n_envel ;		/* number of trajectory envelope points  */
	int n_mach ;		/* number of Mach number fit points      */
	int n_roll ;		/* number of aero. roll fit points       */
	int n_seiff_tab ;	/* number of entries in Seiff table      */
	int i_pnts   ;		/* # points in x dir. for dstr heat flux */
	int j_pnts   ;		/* # points in y dir. for dstr heat flux */
	int k_pnts   ;		/* # points in z dir. for dstr heat flux */
	int n_total ;		/* total # = i_pnts * j_pnts * k_pnts    */
	int n_Re ;		/* number of Reynolds Number fit points  */
	int iter_max ;		/* maximum allowed iterations            */
	int plt_file_type ;	/* specifies plt_file_type type          */
	int run_type ; 		/* specifies run model approximation     */
	int special_data ;	/* specifies special data type           */
	int sphere_cone_type ;	/* precise sphere-cone (half angle) type */
	int thrm_dir_id ;	/* thermodynamics source identification  */
	int wall_model ;	/* temp./heat flux model at outer wall   */  

	char alpha_intrp ;	/* alpha schedule interpolation method   */
	char bank_intrp ;	/* bank schedule interpolation method    */
	char error_char ;	/* error specific char. shown in epitaph */
	char *error_text ;	/* error specific text shown in epitaph  */
	char file_name[N_sufx_list][N_size] ;	/* input file name(s)    */
	char append[N_size] ;	/* unique text to append to file names   */
	char tpl_lbl[N_size];	/* first line of TPL input deck          */
	char ppc_lbl[N_size];	/* first line of Ppc data deck           */
	char mat_lbl[N_size];	/* first line of TPS material data deck  */

	double cnv_fct ;	/* convective heat flux muliplier        */
	double rad_fct ;	/* radiative heat flux muliplier         */
	double prs_fct ;	/* wall pressure muliplier               */
	double emis ; 		/* static heat shield emissivity         */
	double error_float ;	/* error spec. float. point for epitaph  */
	double r_corner ;	/* vehicle corner radius in meters       */
	double r_int ;		/* biconic cone interface radius, m.     */
	double half_angle ;	/* cone half angle in radians            */
	double aft_hlf_angle ;	/* aft biconic half angle in radians     */
	double deceleration ;	/* deceleration trigger value, m/sec^2   */
	double mach_end ;	/* trajectory termination Mach number    */
	double peak_g ;		/* maximum allowed deceleration in Gs    */
	double t_wall ;		/* First guess for wall temp., deg.K.    */
	double t_max ;		/* maximum allowed wall temp., deg.K.    */
	double time_delay ;	/* time delay after trigger event, sec.  */
	double tps_area ;	/* TPS surface area in meters^2          */
	double tps_mass ;	/* TPS forebody mass in kg.              */
	double p_dyn_min ;	/* Aerodynamics minimum dynamic pressure */
	double atm_rho_fct ;	/* Atmospheric density variation factor  */

	double bank_target ; 	/* Target bank angle for RCS maneuver    */

	double I_x ;		/* moment-of-inertia, kg-m^2             */
	double I_y ;		/* moment-of-inertia, kg-m^2             */
	double I_z ;		/* moment-of-inertia, kg-m^2             */
	double I_xy ;		/* product-of-inertia, kg-m^2            */
	double I_xz ;		/* product-of-inertia, kg-m^2            */
	double I_yz ;		/* product-of-inertia, kg-m^2            */
	double I_pl ;		/* Traj inertia parameter                */
	double I_qm ;		/* Traj inertia parameter                */
	double I_rn ;		/* Traj inertia parameter                */
	double I_pq ;		/* Traj inertia parameter                */
	double I_pr ;		/* Traj inertia parameter                */
	double I_qr ;		/* Traj inertia parameter                */

	double x_cg ;		/* Center-of-gravity on X-axis           */
	double y_cg ;		/* Center-of-gravity on Y-axis           */
	double z_cg ;		/* Center-of-gravity on Z-axis           */
	double x_cm ;		/* Aerodynamic-center on X-axis          */
	double y_cm ;		/* Aerodynamic-center on Y-axis          */
	double z_cm ;		/* Aerodynamic-center on Z-axis          */

	/* Displacement between longitudinal accelerometer and origin    */
	double x_x ;		/* x coord. displacement, meters         */
	double y_x ;		/* y coord. displacement, meters         */
	double z_x ;		/* z coord. displacement, meters         */

	/* Displacement between transverse accelerometer and origin      */
	double x_y ;		/* x coord. displacement, meters         */
	double y_y ;		/* y coord. displacement, meters         */
	double z_y ;		/* z coord. displacement, meters         */

	/* Displacement between normal accelerometer and origin          */
	double x_z ;		/* x coord. displacement, meters         */
	double y_z ;		/* y coord. displacement, meters         */
	double z_z ;		/* z coord. displacement, meters         */

	double Re_trns ;	/* Transition Reynolds Number            */
	double Cd_fct ; 	/* drag coefficient, Cd multiplier       */
	double turb_dstnc ;	/* transition to turbulence distance     */
	double heat_of_vapor ;	/* TPS heat of vaporization, Joules/kg   */
	double tps_density ;	/* TPS density, kg/m^3                   */
	double ring_width ;	/* width of probe base ring, meters      */

	double test_a ;
	double test_b ;
	double test_c ;
	double test_e ;
	double test_k ;
	double test_p ;
	double test_u1 ;
	double test_u2 ;
	double test_u3 ;

	double x_intcpt_old ;
	double x_intcpt ;
	double x0_center ;

	double B  ;             /* ballistic coefficient, kg/m^2         */
	double L_D ;            /* lift over drag ratio                  */
	double apoapsis  ;      /* aerocapture apoapsis altitude, meters */

	double *mach_pnt ;	/* Mach number fit point vector          */ 
	double *alpha_pnt ;	/* angle-of-attack fit point vector      */
	double *roll_pnt ;	/* aerodynamic roll fit point vector     */
	double *alpha_free_pnt; /* angle-of-attack (free mole.) pointer  */
	double **logRe_pnt ;	/* log(Reynolds number) fit point vector */
	double *t_alpha_schd ;	/* time for alpha control schedule       */
	double *alpha_schd ;	/* angle-of-attack control schedule      */
	double *t_bank_schd ;	/* time for bank angle control schedule  */
	double *bank_schd ;	/* bank angle control schedule           */

	Aero *continuum_aero ;	/* Continuum aerodynamics table pointer  */
	Aero *free_mole_aero ;	/* Free mole. aerodynamics table pointer */
	Envel *envel ;		/* Sharp trajectory envelope pointer     */
	Seiff *seiff_tab ;	/* Seiff aerodynamic model table pointer */
	Seq seq[N_seq] ;	/* argument list sequence vector         */

	Dstr *dstr ;		/* distr. heat flux data table pointer   */
	Max pnt[I_max] ;	/* Maximum/minimum parameter values      */

	/* Structures common to both Fiat and Traj */
#ifdef Fiat_traj
	Vehicle_fiat_def
#endif
/* Kepler elements appended below -- to be removed */
        double wtmol ;
} Vehicle ;

/* Kepler definitions below -- to be removed */
#define Sphere_jnk      200
#define Plate           100
#define Request         1

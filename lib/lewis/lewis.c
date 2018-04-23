#ifdef Fiat_traj
#include "fiat.h"
#else
#include "var.h"
#endif
#include "lewis.h"

/* #define DEBUG */

int lewis(int i_state, int *p_iter, Thermo *atm)
/***************************************************************************/
/*                                                                         */
/*         --- Thermodynamic State from Lewis Code Algorithms ---          */
/*                                                                         */
/*   This subroutine calculates the equilibrium composition of a gas       */
/*   by Gibbs free energy minimization.  The technique is that of          */
/*   Gordon and McBride, NASA SP-273                                       */
/*                                                                         */
/*   Heading of the specific Gordon and McBride version that this code     */
/*   derives from:                                                         */
/*                                                                         */
/*   NASA-LEWIS CHEMICAL EQUILIBRIUM PROGRAM CEA,  DEC. 12, 1996           */
/*               BY  BONNIE MCBRIDE AND SANFORD GORDON                     */
/*   REFS: NASA RP-1311, PART I, 1994 AND NASA RP-1311, PART II, 1996      */
/*                                                                         */
/*   The last revision date of the Gordon and McBride code that this code  */
/*   was derived from was 23 April 2001.                                   */
/*                                                                         */
/*   This C code derivative uses the NASA Lewis thermodynamic data curve   */
/*   fits.  This code was inspired by a FORTRAN version by Grant Palmer.   */
/*                                                                         */
/*            C Language                        FORTRAN                    */
/*   ----------------------------          ------------------              */
/*   spc[j]->mole_kg_mix / 1000.0              En(j,Npt)                   */
/*                                                                         */
/*                                                                         */
/*   i_state = P_T   for constant pressure and temperature                 */
/*             P_H   for constant pressure and enthalpy                    */ 
/*             P_S   for constant pressure and entropy                     */
/*             Rho_T for constant density and temperature                  */
/*             Rho_H for constant density and enthalpy                     */
/*             Rho_S for constant density and entropy                      */
/*                                                                         */
/*              National Aeronautics and Space Administration              */
/*             No Copyright claimed in USA under Title 17-USC              */
/*                        All other rights reserved                        */
/*                                                                         */
/*                Version:  Mk 2.1,      12 November 2009                  */
/*                      Written by Gary A. Allen, Jr.                      */
/*                                                                         */
/***************************************************************************/
{
register int i, j, k, m, n ;

int bnd, n_elem, n_spc, elem, sw_phi_tabulated ;
int spc_map[N_species], sw_ionized, i_err ;
int elem_remap[N_elem], elm_map[N_elem], n_react ;
int react_list[N_species], sw_n_fnd, sw_m_fnd ;
int sw_C_set, sw_Cgr_set ;
int sw_h2o_set, sw_Wtr_set ;
int sw_el_set, sw_suppress_condensate ;

int sw_pegged = FALSE ;		/* suppress Gcc warning */
int sw_gibbs = FALSE ;		/* suppress Gcc warning */
int n_mat = 0 ;			/* suppress Gcc warning */
int iter = 0 ;			/* suppress Gcc warning */
int sw_finish = 0 ;		/* suppress Gcc warning */
int i_start_cnd = 0 ;		/* suppress Gcc warning */
int i_extra = 0 ;		/* suppress Gcc warning */
int i_for = 0 ;			/* suppress Gcc warning */
int i_against = 0 ;		/* suppress Gcc warning */
int iter_max = 0 ;		/* suppress Gcc warning */

#define Cnvrg_crt	1.0e-5
double dlnn_spc ; 	/* convergence parameter */

double rnln, dlv_dlt_p, dlv_dlp_t, stem ;
double tmpln, g[N_species][N_species] ;
double rl, log_t, cp_frzn, cp_reac, cp  ;
double dlnn, exp_dlnn, mij, arg ;
double bvtr[N_max], bzero[N_max], amtrx[N_max][N_species] ;
double b_save[N_max], a[N_species][N_elem] ;
double mole_fract[N_species], phi, x_phi, cond[N_species] ;
double kt_e, kt_e_2, kt_e_3, lambda, debye, ionic ;
double x_psi, psi, mole_ij, delr_h[N_species], cnd_rct_spc[N_species] ;
double visc[N_species][N_species], accuracy, acc_1, acc_2 ;
double stx[N_species][N_species][N_species], rt_pd, cnd_rct ;
double stc_m, stc_n, denom, totn, wtmol ;

double rmlfr, log_rmlfr, log_p ;
static double p_ref = 1.0e5 ;	/* pascals, reference condition */

#ifdef DEBUG
double sum ;
#endif /* DEBUG */

double one_mw = 0.0 ; 		/* suppress Gcc warning */
double mw = 0.0 ; 		/* suppress Gcc warning */

double t = 0.0 ; 		/* suppress Gcc warning */
double dlnnmax = 0.0 ;   	/* suppress Gcc warning */
double relax_rate = 0.0 ; 	/* suppress Gcc warning */
double ambda = 0.0 ;        	/* suppress Gcc warning */
double p = 0.0 ;		/* suppress Gcc warning */
double rho = 0.0 ;		/* suppress Gcc warning */
double omega_ion = 0.0 ;	/* suppress Gcc warning */
double omega_neut = 0.0 ;	/* suppress Gcc warning */
double omega_el = 0.0 ;		/* suppress Gcc warning */
double omega = 0.0 ;		/* suppress Gcc warning */

Formula *frm ;
Thermo orig, scr ;
Therm_data *spc[N_species], *s, *si, *sj, *sm, *sn, *elec ;
static Therm_data base[N_species] = {
	Thermo_data_base_1,
	Thermo_data_base_2,
	Thermo_data_base_3
} ;

static Visc_inter_data inter_act[N_interact] = Visc_interaction_data_base ;
Visc_inter_data *x ;

/* Function Prototypes */
int condensed(Therm_data*, Thermo*) ;
int solve_lewis(int, double[][N_species], double[], int) ;
int mixture(int,int,Therm_data*[],Thermo*) ;

/* supress a warning */
orig.p = 0.0 ;
orig.t = 0.0 ;
orig.h_298_15 = 0.0 ;
orig.s = 0.0 ;
orig.rho = 0.0 ;


/* pass a pointer */
frm = atm->frm ;

/* See if condensed spieces was selected and if so then locate it */
k = 0 ; /* get rid of a compiler warning */
for (n = 0; n < N_species ; ++n) spc_map[n] = Jnk ;
for (n = 0; n < N_species ; ++n) {
	if (frm[n].species == End) break ;
}

i_err = FALSE ;	/* intially assume success */
n_spc = n ;
n_elem = 0 ;
sw_C_set = FALSE ;
sw_Cgr_set = FALSE ;
sw_Wtr_set = FALSE ;
sw_h2o_set = FALSE ;
for (n = 0; n < n_spc ; ++n) {
	if (base[frm[n].species].sw_elem) ++n_elem ;
	if (frm[n].species == C) sw_C_set = TRUE ;
	if (frm[n].species == Cgr) {
		k = n ; 	/* remember where condensed species is */
		sw_Cgr_set = TRUE ;
	}
	if (frm[n].species == H2O) sw_h2o_set = TRUE ;
	if (frm[n].species == Wtr) {
		k = n ; 	/* remember where condensed species is */
		sw_Wtr_set = TRUE ;
	}
}

/* Set matrix size */
n_mat = n_elem + 1;
i_start_cnd = n_elem ;

/* Swap condensed species with whatever was next to the end */
if (sw_Wtr_set || sw_Cgr_set) {
	frm[k].species = frm[n_spc-1].species ;
	if (sw_Wtr_set) {
		if (!sw_h2o_set) return(Bad_phase_input) ;
		frm[n_spc-1].species = Wtr ;
	}
	else {
		if (!sw_C_set) return(Bad_phase_input) ;
		frm[n_spc-1].species = Cgr ;
	}
	frm[n_spc].mole_fract = frm[n_spc-1].mole_fract ;
	frm[n_spc-1].mole_fract = frm[k].mole_fract ;
	frm[k].mole_fract = frm[n_spc].mole_fract ;
	if ((i_state == P_T)||(i_state == Rho_T)) i_extra = n_elem + 1 ;
	else i_extra = n_elem + 2 ;
	++n_mat ; 
}
else i_extra = n_elem ;

/* condensate removal loop */
sw_suppress_condensate = FALSE ;
for (;;) {

if (sw_suppress_condensate) {
	i_extra = n_elem ;
	--n_mat ;
	--n_spc ;
	frm[n_spc].species = End ;
	sw_suppress_condensate = FALSE ;
}

/* Initialize the species pointers */
sw_el_set = FALSE ;
for (n = 0; n < n_spc ; ++n) {
	/* Deal with the Electron species singularity */
	if (frm[n].species == El) {
		sw_el_set = TRUE ;
		if (frm[n].mole_fract <= 0.0) frm[n].mole_fract = 1.0e-16 ;
	}
	mole_fract[frm[n].species] = frm[n].mole_fract ;
}

for (n = n_spc + 1; n < N_species ; ++n) frm[n].species = Jnk ;
for (n = 0; n < N_species ; ++n) {
        if (frm[n].species != Jnk) spc_map[frm[n].species] = n ;
}

/* Deal with bad parameters */
if (n_spc > N_species) return(Too_many_species) ;

/* Assign species */
for (n = 0; n < N_species ; ++n) {
        if (spc_map[n] != Jnk) spc[spc_map[n]] = &base[n] ;
}

/* Construct the element mapping table and scan for ions */
m = 0 ;
sw_ionized = FALSE ;
for (n = 0; n < n_spc; ++n) {
	s = spc[n] ;
	if (s->type == Ion) sw_ionized = TRUE ;
	for (k = 0; k < 3; ++k) {
		if ((elem = s->comp[k].species) < 0) break ;
		if (m > 0) {
			for (j = 0; j < m; j++) {
				if (elem_remap[j] == elem) break ;
			}
			if (j < m) continue ;
		}
		elem_remap[m++] = elem ;
	}
}
if (m != n_elem) {
	printf("m = %d, n_elem = %d\n", m, n_elem) ;
	return(Too_many_elements) ;
}

/* Deal with bad ionized species input */
if (((!sw_el_set) && sw_ionized)||(sw_el_set && (!sw_ionized))) {
	return(Bad_ion_input) ;
}

/* Assign elements */
for (n = 0; n < n_elem; ++n) elm_map[elem_remap[n]] = n ;

for (n = 0; n < n_spc; ++n) {
	for (j = 0; j < n_elem; ++j) a[n][j] = 0.0 ;
	for (k = 0; k < 3; ++k) {
		if ((j = spc[n]->comp[k].species) < 0) break ;
		a[n][elm_map[j]] = (double)(spc[n]->comp[k].stoic) ;
	}
}

/* Determine initial thermodynamic state */
switch(i_state) {
case P_T:	/* Constant Pressure and Temperature */
	sw_gibbs = TRUE ;
	iter_max = 150 ;
	relax_rate = 0.01 ; /* tuned for Make_mollier */
	if ((t = atm->t) <= 0.0) return(Bad_temperature) ;
	if ((p = atm->p) <= 0.0) return(Bad_pressure) ;
	orig.p = p ; 
	orig.t = t ; 
	break ;

case P_H:	/* Constant Pressure and Enthalpy */
	sw_gibbs = FALSE ;
	relax_rate = 0.008 ;
	iter_max = 100 ;
	n_mat++ ;
	t = 12000.0 ; /* Temperature first guess */
	if ((p = atm->p) <= 0.0) return(Bad_pressure) ;
	orig.p = p ; 
	orig.h_298_15 = atm->h_298_15 ; 
	break ;

case P_S:	/* Constant Pressure and Entropy */
	sw_gibbs = FALSE ;
	relax_rate = 0.01 ; /* tuned to Jupiter entry */
	iter_max = 60 ;
	n_mat++ ;
	t = 12000.0 ; /* Temperature first guess */
	if ((p = atm->p) <= 0.0) return(Bad_pressure) ;
	orig.p = p ; 
	orig.s = atm->s ; 
	break ;

case Rho_T:	/* Constant Density and Temperature */
	sw_gibbs = TRUE ;
	relax_rate = 0.007 ;
	iter_max = 150 ;
	if ((t = atm->t) <= 0.0) return(Bad_temperature) ;
	if ((rho = atm->rho) <= 0.0) return(Bad_density) ;
	p = 1.0e5 ; /* first guess */
	orig.rho = rho ; 
	orig.t = t ; 
	break ;

case Rho_H:	/* Constant Density and Enthalpy */
	sw_gibbs = FALSE ;
	relax_rate = 0.02 ;
	iter_max = 60 ;
	n_mat++ ;
	if ((rho = atm->rho) <= 0.0) return(Bad_density) ;
	t = 12000.0 ; /* Temperature first guess */
	p = 1.0e5 ; /* first guess */
	t = 12000.0 ; /* Temperature first guess */
	orig.rho = rho ; 
	orig.h_298_15 = atm->h_298_15 ; 
	break ;

case Rho_S:	/* Constant Density and Entropy */
	sw_gibbs = FALSE ;
	relax_rate = 0.01 ; /* tuned to Jupiter entry */
	iter_max = 60 ;
	n_mat++ ;
	if ((rho = atm->rho) <= 0.0) return(Bad_density) ;
	t = 12000.0 ; /* Temperature first guess */
	p = 1.0e5 ; /* first guess */
	orig.rho = rho ; 
	orig.s = atm->s ; 
	break ;

default:
	return(Unknown_state) ;
} /* end of initial state Switch block */

/*   calculate initial molar mass in kg/mole */
mw = 0.0 ;
for (n = 0; n < n_spc ; ++n) {
	j = frm[n].species ;
	spc[n]->mole_fr = mole_fract[j] ;
	mw += spc[n]->mole_fr * spc[n]->mole_wt;
}
one_mw = 1.0 / mw ;

for (n = 0; n < n_spc ; ++n) {
	spc[n]->mass_fr = spc[n]->mole_fr * spc[n]->mole_wt * one_mw ;
}

/* convert mole fractions to moles per kilogram mixture */
for (n = 0; n < n_spc; ++n) spc[n]->mole_kg_mix = spc[n]->mole_fr / mw;

/*   set elemental conservation constants for Gibbs */
/*   free energy minimization. */
for (m = 0; m < n_elem; ++m) {
	bzero[m] = 0.0 ;
	for (n = 0; n < n_spc; ++n) {
		bzero[m] += a[n][m] * spc[n]->mole_kg_mix ;
	}
}

if (n_mat > N_max) return(Too_many_elements) ;

sw_finish = FALSE ;
for (iter = 1; iter <= iter_max; ++iter) {

/*   Calculate the species specific heat, enthalpy, and gibbs free */
/*   energy for the current temperature.                           */
/*   Thermodynamic quantities are based on curve fits to Cp        */
scr.mw = mw * 1000.0 ;
scr.t = t ;
scr.p = p ;
sw_pegged = mixture(FALSE, n_spc, spc, &scr) ;

/*   load the matrices and solve for the thermodynamic updates */
/*   fill [a] matrix.  see NASA SP-273 page 228. */
for (n = 0; n < n_mat; ++n) for (m = 0; m < n_mat; ++m) amtrx[n][m] = 0.0 ;
amtrx[i_extra][i_extra] = -one_mw ;
bvtr[i_extra] = one_mw;
for (m = 0; m < n_elem; ++m) bvtr[m] = bzero[m] ;

#ifdef DEBUG
printf("bug-0: bvtr[0] = bzero[0] = %e\n", bzero[0]) ;
#endif /* DEBUG */

for (n = 0; n < n_spc; ++n) {
	s = spc[n] ;
	if ((s->type == Solid)||(s->type == Liquid)) {
#ifdef DEBUG
printf("Solid/Liquid mole [%d] fraction %e s->mole_kg_mix %e 1/mw %e\n",
		n, s->mole_kg_mix * mw, s->mole_kg_mix, 1/mw) ; 
#endif /* DEBUG */

/* temporary -- put gibbs free energy test here */
if (sw_Wtr_set) {
		rmlfr = (base + H2O)->mole_kg_mix * mw;
		if (rmlfr <= 1.0e-10) log_rmlfr = -23.02585093 ;
		else log_rmlfr = log(rmlfr) ;
		log_p = log(scr.p / p_ref) ;
#ifdef DEBUG
		printf("Water->gibbs_free %e\n", (base + Wtr)->gibbs_free) ;
		printf("Stream->gibbs_free %e\n",  (base + H2O)->gibbs_free
			+ log_p + log_rmlfr) ;
#endif /* DEBUG */
		if (((base + H2O)->gibbs_free
                        + log_p + log_rmlfr) < (base + Wtr)->gibbs_free) {
#ifdef DEBUG
/*
			sw_suppress_condensate = TRUE ;
			break ;
*/
			printf("suppress condensate\n") ; /* temporary */
#endif /* DEBUG */
		}
#ifdef DEBUG
		else printf("do NOT suppress condensate\n") ; /* temporary */
#endif /* DEBUG */
}

if (sw_Cgr_set) {
		rmlfr = (base + C)->mole_kg_mix * mw;
		if (rmlfr <= 1.0e-10) log_rmlfr = -23.02585093 ;
		else log_rmlfr = log(rmlfr) ;
		log_p = log(scr.p / p_ref) ;
#ifdef DEBUG
		printf("Cgr->gibbs_free %e\n", (base + Cgr)->gibbs_free) ;
		printf("C->gibbs_free %e\n",  (base + C)->gibbs_free
			+ log_p + log_rmlfr) ;
#endif /* DEBUG */
		if (((base + C)->gibbs_free
                        + log_p + log_rmlfr) < (base + Cgr)->gibbs_free) {
#ifdef DEBUG
/*
			sw_suppress_condensate = TRUE ;
			break ;
*/
			printf("suppress condensate\n") ; /* temporary */
#endif /* DEBUG */
		}
#ifdef DEBUG
		else printf("do NOT suppress condensate\n") ; /* temporary */
#endif /* DEBUG */
}
		for (m = 0; m < n_elem; ++m) {
    			amtrx[i_start_cnd][m] = a[n][m] * 1000.0 ; 
    			amtrx[m][i_start_cnd] = a[n][m] * 1000.0 ; 
			bvtr[m] -= a[n][m] * s->mole_kg_mix * 1000.0 ; /* fix */
		}
		amtrx[i_start_cnd][i_start_cnd] = 0.0 ;
		a[i_start_cnd][i_extra] = 0.0 ;
		a[i_extra][i_start_cnd] = 0.0 ; 
		bvtr[i_start_cnd] = s->gibbs_free * 1000.0 ;
		continue ;
	}
#ifdef DEBUG
	else {
printf("Gas mole [%d] fraction %e s->mole_kg_mix %e 1/mw %e\n",
		n, s->mole_kg_mix * mw, s->mole_kg_mix, 1/mw) ; 
	}
#endif /* DEBUG */

	for (m = 0; m < n_elem; ++m) {
    		amtrx[i_extra][m] += a[n][m] * s->mole_kg_mix ;
    		amtrx[m][i_extra] += a[n][m] * s->mole_kg_mix ;
		for (j = 0; j < n_elem; ++j) {
			amtrx[j][m] += a[n][j] * a[n][m] * s->mole_kg_mix;
		}

		/* Fill [b] vector */
		bvtr[m] += a[n][m] * s->mole_kg_mix * (s->gibbs_free - 1.0) ;
	}
	amtrx[i_extra][i_extra] += s->mole_kg_mix ;
	bvtr[i_extra] += s->mole_kg_mix * (s->gibbs_free - 1.0) ;
}

if (sw_suppress_condensate) break ;

if (sw_finish) {
	if (iter >= iter_max) i_err = Not_converged ;
	else *p_iter = iter ;
	break;
}

/* Additional matrix elements */
switch(i_state) {
case P_H:	/* Constant Pressure and Enthalpy */
case Rho_H:	/* Constant Density and Enthalpy */
	bvtr[i_extra + 1] = atm->h_298_15 / (R_lewis * t) ;
	for (n = 0; n < n_spc; ++n) {
		s = spc[n] ;
		for (m = 0; m < n_elem; ++m) {
			amtrx[i_extra + 1][m] += a[n][m] * s->mole_kg_mix
				* s->enthalpy;
			amtrx[m][i_extra + 1] += a[n][m] * s->mole_kg_mix
				* s->enthalpy;
		}
		amtrx[i_extra + 1][i_extra] += s->mole_kg_mix * s->enthalpy;
		amtrx[i_extra][i_extra + 1] += s->mole_kg_mix * s->enthalpy;
		amtrx[i_extra + 1][i_extra + 1] += s->mole_kg_mix * (s->cp
			+ s->enthalpy * s->enthalpy);
		bvtr[i_extra + 1] += s->mole_kg_mix * s->enthalpy
			* (s->gibbs_free - 1.0) ;
	}
	break ;

case P_S:	/* Constant Pressure and Entropy */
case Rho_S:	/* Constant Density and Entropy */
	bvtr[i_extra + 1] = (atm->s / R_lewis) + one_mw ;
	for (n = 0; n < n_spc; ++n) {
		s = spc[n] ;
		for (m = 0; m < n_elem; ++m) {
			amtrx[i_extra + 1][m] += a[n][m] * s->mole_kg_mix
				* s->entropy;
			amtrx[m][i_extra + 1] += a[n][m] * s->mole_kg_mix
				* s->enthalpy;
		}
		amtrx[i_extra + 1][i_extra] += s->mole_kg_mix * s->entropy;
		amtrx[i_extra][i_extra + 1] += s->mole_kg_mix * s->enthalpy;
		amtrx[i_extra + 1][i_extra + 1] += s->mole_kg_mix * (s->cp
			+ s->enthalpy * s->entropy);

		bvtr[i_extra + 1] += s->mole_kg_mix * (s->entropy
			* s->gibbs_free - s->entropy - 1.0);
	}
	break ;

case P_T:	/* Constant Pressure and Temperature */
case Rho_T:	/* Constant Density and Temperature */
	/* Row/column n_elem + 1 not used */
	break ;

/* Unknown thermodynamic state requested so abort */
default:
	return(Unknown_state) ;
} /* end of additional element Switch block */ 

#ifdef DEBUG
printf("G(i,k) A:  \n") ;
for (m = 0; m < n_mat; ++m) {
	for (j = 0; j < n_mat; ++j) {
		printf(" %13.6e", amtrx[j][m]/1000.0) ;
	}
	/* last column of g(i,k) */
	printf(" %13.6e\n", bvtr[m]/1000.0) ;
}
#endif /* DEBUG */

/*  solve the matrix for species and temperature updates. */
if (solve_lewis(TRUE, amtrx, bvtr, n_mat)) return(Singular) ;

#ifdef DEBUG
for (m = 0; m < n_mat; ++m) {
	printf("soltn[%d] = %e\n", m, bvtr[m]) ;
}
#endif /* DEBUG */

/*   Select the weighting factor.  this is used to prevent the program */
/*   from overshooting too far on the first few iterations. */
if ((rl = iter * relax_rate) >= 1.0) ambda = 1.0 ;
else ambda = rl ;

#ifdef DEBUG
ambda = 1.0 ; 
  sum = fabs(bvtr[i_extra]) * 5.0 ;
  if (sum > 2.0) ambda = 2.0 / sum ; 

printf("abmda %f\n", ambda) ;
#endif /* DEBUG */

dlnnmax = 0.0 ;
i_for = 0 ;
i_against = 0 ;

/*   apply update corrections and check for convergence    */
/*   update species moles per gram mixture                 */
/*   The "bvtr[i_extra + 1]" error appeared in this routine */  
for (n = 0; n < n_spc; ++n) {
	s = spc[n] ;
	if ((s->type == Solid)||(s->type == Liquid)) {
		s->mole_kg_mix += ambda * bvtr[i_start_cnd] ;
	}
	else {
		dlnn = -s->gibbs_free + bvtr[i_extra] ;
		if (!sw_gibbs) dlnn += s->enthalpy * bvtr[i_extra + 1] ;
		for (m = 0; m < n_elem; ++m) dlnn += a[n][m] * bvtr[m] ;
	
		if (dlnn < -50.0) exp_dlnn = 0.0 ;
		else {
			if (dlnn > 50.0) exp_dlnn = 5.18e21 ;
			else exp_dlnn = exp(dlnn) ;
		}
		dlnn_spc = fabs(1.0 - exp_dlnn);
		if (dlnn_spc < Cnvrg_crt) ++i_for ;
		else ++i_against ;

		if (dlnnmax < dlnn_spc) dlnnmax = dlnn_spc ;
		if ((rnln = ambda * dlnn) <= -50.0) rnln = -50.0 ;
		if (rnln >= 50.0) rnln = 50.0 ;
		if ((s->mole_kg_mix *= exp(rnln)) <= 1.0e-10) {
			s->mole_kg_mix = 1.0e-10 ;
		}
	}
}

/* Update molar mass */
if ((arg = ambda * bvtr[i_extra]) <= -50.0) arg = -50.0 ;
if (arg >= 50.0) arg = 50.0 ;
if ((one_mw *= exp(arg)) <= 1.0e-10) one_mw = 1.0e-10 ;
mw = 1.0 / one_mw ;

/* Update nonconstant parameters */
switch(i_state) {
case P_T:	/* Constant Pressure and Temperature */
	break ;

case P_H:	/* Constant Pressure and Enthalpy */
case P_S:	/* Constant Pressure and Entropy */
	if ((tmpln = ambda * bvtr[i_extra + 1]) <= -50.0) tmpln = -50.0 ;
	if (tmpln >= 50.0) tmpln = 50.0 ;
	t *= exp(tmpln);
	if (t < 100.0) t = 100.0 ;
        else if (t > 20000.0) t = 20000.0 ;
	break ;

case Rho_T:	/* Constant Density and Temperature */
	p = rho * R_lewis * t / mw ;
	break ;

case Rho_H:	/* Constant Density and Enthalpy */
case Rho_S:	/* Constant Density and Entropy */
	tmpln = ambda * bvtr[i_extra + 1];
	t *= exp(tmpln);
	if (t < 100.0) t = 100.0 ;
        else if (t > 20000.0) t = 20000.0 ;
	p = rho * R_lewis * t / mw ;
	break ;

default:
	return(Unknown_state) ;
} /* end of parameter update Switch block */

if ((iter == iter_max - 1)||(dlnnmax < Cnvrg_crt)||(i_for > 2 * i_against)) {
	sw_finish = TRUE ;
}
} /* end of iteration FOR loop */
if (sw_finish) break ;
printf("Condensed species must be removed\n") ;
} /* end of condensate removal loop */

/* Solve for isentropic exponent */
amtrx[i_extra][i_extra] = 0.0 ;
n_mat = i_extra + 1;
for (m = 0; m < n_mat; ++m) bvtr[m] = 0.0 ;
totn = 0.0 ;
for (n = 0; n < n_spc; ++n) {
	s = spc[n] ;

	/* Do not add non-gas phase terms */
	if ((s->type == Solid)||(s->type == Liquid)) {
		bvtr[i_start_cnd] = - s->enthalpy * 1000.0 ;
		totn += s->mole_kg_mix * 1000.0 ;
		continue ;
	}
	else totn += s->mole_kg_mix ;

	for (m = 0; m < i_extra; ++m) {
		bvtr[m] -= a[n][m] * s->mole_kg_mix * s->enthalpy;
	}
	bvtr[i_extra] -= s->mole_kg_mix * s->enthalpy;
}
for (m = 0; m < n_mat; ++m) b_save[m] = bvtr[m] ;

#ifdef DEBUG
printf("G(i,k) B:  \n") ;
for (m = 0; m < n_mat; ++m) {
	for (j = 0; j < n_mat; ++j) {
		printf(" %13.6e", amtrx[j][m]/1000.0) ;
	}
	/* last column of g(i,k) */
	printf(" %13.6e\n", bvtr[m]/1000.0) ;
}
#endif /* DEBUG */

if (solve_lewis(TRUE, amtrx, bvtr, n_mat)) return(Singular) ;

#ifdef DEBUG
for (m = 0; m < n_mat; ++m) {
	printf("soltn[%d] = %e\n", m, bvtr[m]) ;
}
#endif /* DEBUG */

/*  compute final mixture enthalpy and entropy */
atm->mw = mw * 1000.0 ;
if (fabs(totn) > Small) atm->mw_1 = 1000.0 / totn ;
else atm->mw_1 = 0.0 ;
atm->t = t ;
log_t = log(t) ; /* for viscosity calculations */
atm->p = p ;
sw_pegged = mixture(TRUE, n_spc, spc, atm) ;

cp_reac = 0.0 ;
cp_frzn = 0.0 ;
for (m = 0; m < n_mat; ++m) cp_reac -= b_save[m] * bvtr[m] ;
for (n = 0; n < n_spc; ++n) {
	s = spc[n] ;
	/* condensed species does not contribute to reaction Cp */
	if ((s->type == Solid)||(s->type == Liquid)) {
		cp_frzn += s->mole_kg_mix * 1000.0 * s->cp ;
	}
	else {
		cp_reac += s->mole_kg_mix * s->enthalpy * s->enthalpy ;
		cp_frzn += s->mole_kg_mix * s->cp ;
	}
}
 
cp = cp_frzn + cp_reac ;
atm->cp = cp * R_lewis ;
dlv_dlt_p = 1.0 + bvtr[i_extra] ;

for (m = 0; m < n_mat; ++m) bvtr[m] = 0.0 ;
for (n = 0; n < n_spc; ++n) {
	s = spc[n] ;

	if ((s->type == Solid)||(s->type == Liquid)) continue ;
	for (m = 0; m < n_elem; ++m) {
    		bvtr[m] += a[n][m] * s->mole_kg_mix ;
	}
	bvtr[i_extra] += s->mole_kg_mix ;
}

#ifdef DEBUG
printf("G(i,k) C:  \n") ;
for (m = 0; m < n_mat; ++m) {
	for (j = 0; j < n_mat; ++j) {
		printf(" %13.6e", amtrx[j][m]/1000.0) ;
	}
	/* last column of g(i,k) */
	printf(" %13.6e\n", bvtr[m]/1000.0) ;
}
#endif /* DEBUG */

solve_lewis(FALSE, amtrx, bvtr, n_mat) ;

dlv_dlp_t = -1.0 + bvtr[i_extra] ;
atm->gamma = -1 / (dlv_dlp_t + (dlv_dlt_p * dlv_dlt_p / (cp * mw))) ;

#ifdef DEBUG
printf("\ndlv_dlt_p %e\n", dlv_dlt_p) ;
printf("dlv_dlp_t %e\n\n", dlv_dlp_t) ;
#endif /* DEBUG */

atm->h_cnv = 0.0 ;
atm->h_d = 0.0 ;

/* Calculate conversion heat of formation and total gas */
/* heat of formation as enthalpy.                       */
for (n = 0; n < n_spc; ++n) {
	s = spc[n] ;
	if ((s->type == Solid)||(s->type == Liquid)) {
		s->mole_fr = s->mole_kg_mix * 1000.0 / totn ;
		s->mass_fr = s->mole_kg_mix * 1000.0 * s->mole_wt;
		atm->h_cnv += s->mole_kg_mix * 1000.0 * s->h_cnv ;
		atm->h_d += s->mole_kg_mix * 1000.0 * s->del_h_298_15 ;
	}
	else {
		s->mole_fr = s->mole_kg_mix / totn ;
		s->mass_fr = s->mole_kg_mix * s->mole_wt;
		atm->h_cnv += s->mole_kg_mix * s->h_cnv ;
		atm->h_d += s->mole_kg_mix * s->del_h_298_15 ;
	}
	mole_fract[frm[n].species] = s->mole_fr ;
} /* end of mass and mole fraction FOR loop */

for (i = 0; i < n_spc; ++i) {
	for (j = 0; j < n_spc; ++j) visc[i][i] = 0.0 ;
}

/* If gas is ionized then calculate the global parameters */
if (sw_ionized) {
	elec = base + El ; /* special pointer to electron gas */
	kt_e = (Boltz_lewis * t * 1.0e10) / E_2 ;
	kt_e_2 = kt_e * kt_e ;
	kt_e_3 = kt_e * kt_e_2 ;
	debye = 9.0 * R_lewis * t * kt_e_3 / (40.0 * Pi * Avgdr_lewis
		* elec->mole_fr)  ;
	ionic = pow(9.0 * debye, 2.0/3.0) ;
	lambda = sqrt(debye + ionic) ;
	if (lambda < 2.718281828) lambda = 2.718281828 ;
	omega_ion = 1.36e-6 * log(lambda) / kt_e_2 ;
	omega_el = 1.29e-6 * log(lambda) / kt_e_2 ;
	omega_neut = 1.0e-28 * exp(6.776 - 0.4 * log_t) ;
}

/* Compute species viscosity, in micropoises and conductivity */
for (n = 0; n < n_spc; ++n) {
	s = spc[n] ;

	/* Ions require special treatment */
	if (s->type == Ion) {
		/* Calculate viscosity from cross section equation */
		visc[n][n] = 5.0 * sqrt(10.0 * Pi * s->mole_wt
			* Boltz_lewis * t / Avgdr_lewis) / (omega_ion * 16.0) ;

		/* Calculate conductivity based upon above viscosity */
		cond[n] = visc[n][n] * R_lewis * (0.00375 + 0.00132
			* (s->cp - 2.5)) / s->mole_wt ;
		continue ;
	}

	/* Use crude species viscosity and conductivity */
	/* approximation if untabulated.                */
	if (!s->n_bands) {
		/* CN and C2 viscosity are calculated here */ 
		omega  = 1.0e-26 * log(50.0 * pow(s->mole_wt * 1000.0, 4.6)
			/ pow(t, 1.4)) ;
		if (omega < 1.0e-26) omega = 1.0e-26 ;
		visc[n][n] = 5.0 * sqrt(1.0e5 * s->mole_wt * Boltz_lewis * t
			/ (Pi * Avgdr_lewis)) / (omega * 16.0) ;

		/* Calculate conductivity based upon above viscosity */
		cond[n] = visc[n][n] * R_lewis * (0.00375 + 0.00132
			* (s->cp - 2.5)) / s->mole_wt ;
		continue ;
	}

	/* Search for temperature band */
	for (;;) {
	/* Deal with temperature below the data base */
	if (t <= s->t_low[0]) {
		bnd = 0 ;
		break ;
	}

	/* Deal with temperature above the data base  */
	if (t >= s->t_hgh[s->n_bands - 1]) {
		bnd = s->n_bands - 1 ;
		break ;
	}

	/* Temperature is within the data base so search for the viscosity */
	/*    temperature band for this species.                           */
	for (bnd = 0; bnd < s->n_bands; ++bnd) {
		if ((s->t_low[bnd] <= t) && (t <= s->t_hgh[bnd])) break ;
	}
	break ;
	/*NOTREACHED*/
	} /* end of temperature band inf. FOR loop */

	visc[n][n] = exp(s->visc_poly[0][bnd] * log_t 
		+ (s->visc_poly[1][bnd] + s->visc_poly[2][bnd] / t) / t
		+ s->visc_poly[3][bnd]) ;

	cond[n] = exp(s->cond_poly[0][bnd] * log_t 
		+ (s->cond_poly[1][bnd] + s->cond_poly[2][bnd] / t) / t
		+ s->cond_poly[3][bnd]) ;
	
} /* end of FOR loop calculating individual species viscosity */

/* Calculate the viscosity species interaction in micropoises */
for (i = 0; i < n_spc; ++i) {
si = spc[i] ;

for (j = 0; j < n_spc; ++j) {
if (i == j) break ;
sj = spc[j] ;
sw_phi_tabulated = FALSE ;

/* Search table if viscous combinations tabulated */
if (si->n_combo_pnt) {

	/* Scan viscous combination species table */
	for (k = 0; k < si->n_combo_pnt; k++) {

		/* Does combination species match table entry? */
		if (frm[j].species == si->combo_species[k]) {
			sw_phi_tabulated = TRUE ;
			x = &inter_act[si->combo_pnt[k]] ;

			/* Determine the temperature band */
			for(;;) {
			/* Deal with temperature below the data base */
			if (t <= x->t_low[0]) {
				bnd = 0 ;
				break ;
			}

			/* Deal with temperature above the data base   */
			if (t >= x->t_hgh[x->n_bands - 1]) {
				bnd = x->n_bands - 1 ;
				break ;
			}

			/* Temperature is within the data base so search for */
			/*   the temperature band for this species.          */
			for (bnd = 0; bnd < x->n_bands; ++bnd) {
				if ((x->t_low[bnd] <= t)
					&& (t <= x->t_hgh[bnd])) {
					break ;
				}
			}
			break ;
			/*NOTREACHED*/
			} /* end of temperature band inf. FOR loop */
			visc[i][j] = exp(x->visc_poly[0][bnd] * log_t
				+ (x->visc_poly[1][bnd]
				+ x->visc_poly[2][bnd] / t) / t
				+ x->visc_poly[3][bnd]) ;
			break ;
		} /* end of table entry match routine */
	} /* end of visc. comb. table scan FOR loop */
} /* end of tabulated viscous combo. search */

if (!sw_phi_tabulated) {
if (sw_ionized) {
switch(si->type) {
case Ion:
	switch(sj->type) {
	case Ion:
		omega = omega_ion ;
		break ;

	case El:
		omega = omega_el ;
		break ;

	case Neutral:
	case Liquid:  /* temporary */
	case Solid:  /* temporary */
		omega = omega_neut ;
		break ;
	}
	break ;

case El:
	switch(sj->type) {
	case Ion:
	case El:
		omega = omega_el ;
		break ;

	case Neutral:
	case Liquid:  /* temporary */
	case Solid:  /* temporary */
		omega = omega_neut ;
		break ;
	}
	break ;

case Neutral:
case Liquid:  /* temporary */
case Solid:  /* temporary */
	switch(sj->type) {
	case Ion:
	case El:
		omega = omega_neut ;
		break ;

	case Neutral:
	case Liquid:  /* temporary */
	case Solid:  /* temporary */
		if ((fabs(visc[i][i]) < Very_small)
			||(fabs(visc[j][j]) < Very_small)) {

			visc[i][j] = 0.0 ;
			visc[j][i] = 0.0 ;
			continue ;
		}
		/* replace with rigid sphere routine */
		mij = sqrt(sqrt(si->mole_wt) / visc[i][i])
			+ sqrt(sqrt(sj->mole_wt) / visc[j][j]) ;
		omega = 5.0 * sqrt(20.0 * Pi * Boltz_lewis * t / Avgdr_lewis)
			* mij * mij / (16.0 * 5.656854) ;
		break ;
	}
	break ;

default:
	return(Logic_error) ;
} /* end of Not-tabulated Switch block */

visc[i][j] = 5.0 * sqrt(20.0 * Pi * si->mole_wt * sj->mole_wt * Boltz_lewis * t
	/ (Avgdr_lewis * (si->mole_wt + sj->mole_wt))) / (omega * 16.0) ;
} /* end of Ionized block */

/* If viscous combination species not tabulated then use Wilke's Rule */
/* based on equation (5.3), pg. 21 of RP-1311.  The Xi term in the    */
/* denominator of equ.(5.3) is a typographical error.  "Hypersonic    */
/* and High Temperature Gas Dynamics" by Anderson shows the correct   */
/* equation in equ.(16.20), pg. 597.                                  */
else {
	mij = sqrt(visc[i][i] / sqrt(si->mole_wt)) + sqrt(visc[j][j]
		/ sqrt(sj->mole_wt)) ;
	if (fabs(mij) < Very_small) visc[i][j] = 0.0 ;
	else {
		visc[i][j] = 4.0 * visc[i][i] * visc[j][j] * sqrt(2.0
			/ (si->mole_wt + sj->mole_wt)) / (mij * mij) ;
	}
} /* end of Wilke's Rule block */
} /* end of Not-tabulated block */

/* Visc[][] matrix is symmetrical */
visc[j][i] = visc[i][j] ;

} /* end of "j" FOR loop for viscosity species interaction */
} /* end of "i" FOR loop for viscosity species interaction */

/* Calculate the total mixed viscosity and frozen conductivity*/
atm->mu = 0.0 ;
atm->cnd_frz = 0.0 ;
for (i = 0; i < n_spc; ++i) {
	if (visc[i][i] < Visc_cutoff) continue ;
	si = spc[i] ;
	if ((si->type == Solid)||(si->type == Liquid)) continue ;
	x_phi = 0.0 ;
	x_psi = 0.0 ;
	for (j = 0; j < n_spc; ++j) {
		if (visc[i][j] < Visc_cutoff) continue ;
		sj = spc[j] ;
		if ((sj->type == Solid)||(sj->type == Liquid)) continue ;
		mole_ij = si->mole_wt + sj->mole_wt ;
		stem = sj->mole_wt / (visc[i][j] * mole_ij) ;
		if (i == j) {
			phi = 1.0 ;
			psi = 1.0 ;
		}
		else {
			phi = 2.0 * visc[i][i] * stem ;
			psi = phi * (1.0 + (2.41 * (si->mole_wt - sj->mole_wt)
				* (si->mole_wt - 0.142 * sj->mole_wt)
				/ (mole_ij * mole_ij))) ;
		}

		x_phi += sj->mole_fr * phi ;
		x_psi += sj->mole_fr * psi ;
	} /* end of "j" FOR loop for mixed viscosity */
	atm->mu  += si->mole_fr * visc[i][i] / x_phi ;
	atm->cnd_frz += si->mole_fr * cond[i] / x_psi ;
} /* end of "i" FOR loop for mixed viscosity */

/* Calculate heat of reaction for determining equilibrium conductivity */
j = 0 ;
for (i = 0; i < n_spc; i++) {
	si = spc[i] ;
	if (si->sw_elem) continue ;
	if ((si->type == Solid)||(si->type == Liquid)) continue ;
	
	react_list[j] = i ;
	delr_h[j] = -si->enthalpy ;
	for (k = 0; k < 3; k++) {
		if (si->comp[k].species == End) break ;
		delr_h[j] += base[si->comp[k].species].enthalpy
			* ((double)(si->comp[k].stoic)) ;
	}
	j++ ;
}
n_react = j ;

/* load the "ax" matrix */
for (i = 0; i < n_react; i++) {
si = spc[react_list[i]] ;
if ((si->type == Solid)||(si->type == Liquid)) continue ;
for (m = 0; m < n_spc - 1; ++m) {
sm = spc[m] ;
if ((sm->type == Solid)||(sm->type == Liquid)) continue ;
for (n = m + 1; n < n_spc; ++n) {
	sn = spc[n] ;
	if ((sn->type == Solid)||(sn->type == Liquid)) continue ;
	sw_n_fnd = FALSE ;
	sw_m_fnd = FALSE ;
	stc_m = 0.0 ;
	stc_n = 0.0 ;

	if (m == react_list[i]) {
		stc_m = -1.0 ; 
		sw_m_fnd = TRUE ;
		for (j = 0; j < 3; j++) {
			if (si->comp[j].species == End) break ;
			if (si->comp[j].species == frm[n].species) {
				stc_n = (double)(si->comp[j].stoic) ;
				sw_n_fnd = TRUE ;
				break ;
			}
		} /* end of "j" reaction component FOR loop */
	}
	else {
	if (n == react_list[i]) {
		stc_n = -1.0 ;
		sw_n_fnd = TRUE ;
		for (j = 0; j < 3; j++) {
			if (si->comp[j].species == End) break ;
			if (si->comp[j].species == frm[m].species) {
				stc_m = (double)(si->comp[j].stoic) ;
				sw_m_fnd = TRUE ;
				break ;
			}
		} /* end of "j" reaction component FOR loop */
	}
	}
	if (!sw_m_fnd) {
	for (j = 0; j < 3; j++) {
		if (si->comp[j].species == End) break ;
		if (si->comp[j].species == frm[m].species) {
			stc_m = (double) si->comp[j].stoic  ;
			sw_m_fnd = TRUE ;
			break ;
		}
	} /* end of "j" reaction component FOR loop */
	}
	if (!sw_n_fnd) {
	for (j = 0; j < 3; j++) {
		if (si->comp[j].species == End) break ;
		if (si->comp[j].species == frm[n].species) {
			stc_n = (double) si->comp[j].stoic ;
			sw_n_fnd = TRUE ;
			break ;
		}
	} /* end of "j" reaction component FOR loop */
	}
	stx[i][m][n] = (sn->mole_fr * stc_m) - (sm->mole_fr * stc_n) ; 

#ifdef DEBUG
	printf("frm[react_list[%d]].species %d\n",
		i, frm[react_list[i]].species) ;
	printf("   stx[%d][%d][%d] = %e\n", i, m, n, stx[i][m][n]) ;
	printf("   sn->mole_fr %e sm->mole_fr %e\n", sn->mole_fr, sm->mole_fr) ;
	printf("   stc_m %e stc_n %e\n", stc_m, stc_n) ;
#endif /* DEBUG */

} /* end of "n" element FOR loop */
} /* end of "m" element FOR loop */
} /* end of "i" reaction FOR loop */

/* load the "g" matrix */
for (i = 0; i < n_react; i++) {
si = spc[react_list[i]] ;
if ((si->type == Solid)||(si->type == Liquid)) continue ;
for (j = i; j < n_react; j++) {
	sj = spc[react_list[j]] ;
	if ((sj->type == Solid)||(sj->type == Liquid)) continue ;
	g[i][j] = 0.0 ;

	for (m = 0; m < n_spc - 1; ++m) {
		sm = spc[m] ;
		if ((sm->type == Solid)||(sm->type == Liquid)) continue ;
		for (n = m + 1; n < n_spc; ++n) {
			if (visc[m][n] < Visc_cutoff) continue ;
			sn = spc[n] ;
			if ((sn->type == Solid)
				||(sn->type == Liquid)) continue ;
			if ((denom = sm->mole_fr * sn->mole_fr) < 1.0e-20) {
				continue ;
			}
			/* Equ. 5.13, pg. 22, RP-1311 */
			rt_pd = 5.0 * sm->mole_wt * sn->mole_wt * 1000.0
				/ (3.0 * 1.1 * visc[m][n]
				* (sm->mole_wt + sn->mole_wt)) ;
			/* Binary diffusion coefficient, cm^2/sec          */
			/* Only 90% accurate compared to CRC manual values */
		/*	d[m][n] = t * R_lewis / (p * rt_pd) ; */ 
			g[i][j] += rt_pd * stx[i][m][n] * stx[j][m][n] 
				/ denom ;
		}
	} /* end of "m" element FOR loop */
	g[j][i] = g[i][j] ; /* matrix is symmetric */
} /* end of "j" species FOR loop */
} /* end of "i" species FOR loop */

for (i = 0; i < n_react; i++) cnd_rct_spc[i] = delr_h[i] ;

if (solve_lewis(TRUE, g, cnd_rct_spc, n_react)) return(Singular) ;

cnd_rct = 0.0 ;
for (i = 0; i < n_react; i++) cnd_rct += delr_h[i] * cnd_rct_spc[i] ;
cnd_rct *= R_lewis ;

atm->cnd_equ = atm->cnd_frz + cnd_rct ;

atm->mu  *= 1.0e-7 ;  /* convert micropoises to kg/(m sec) */
atm->cnd_equ *= 1.0e-4 ;  /* convert microwatts / (cm sec) to Watts / (m sec) */
atm->cnd_frz *= 1.0e-4 ;  /* convert microwatts / (cm sec) to Watts / (m sec) */

/* Recalculate frozen heat capacity and molecular weight for gases only */
cp_frzn = 0.0 ;
wtmol = 0.0 ;
for (n = 0; n < n_spc; ++n) {
	s = spc[n] ;
	if ((s->type == Solid)||(s->type == Liquid)) continue ;
	cp_frzn += s->mole_kg_mix * mw * s->cp ;
	wtmol += s->mole_kg_mix * mw * s->mole_wt;
}

/* Set final values and return */
/* Frozen Prandtl number */
atm->prandtl = cp_frzn * R_lewis * atm->mu / (wtmol * atm->cnd_frz) ;
/* Assume a bogus Schmidt number = 1.0 for later correction */
atm->lewis_nmbr = cp_frzn * R_lewis * atm->mu / (wtmol * atm->cnd_frz) ;
atm->rho = p * mw / (t * R_lewis);
atm->c = sqrt(atm->gamma * t * R_lewis / mw) ; 

for (n = 0; n < n_spc ; ++n) frm[n].mole_fract = mole_fract[frm[n].species] ;

/* Do the final sanity test */
/* temporary 
if ((atm->gamma < 1.0)
	||(atm->gamma > 1.66666666666)
	||(atm->mw < 5.4e-4)) return(Junk_result) ;
*/

if (sw_Wtr_set) {
	(base + Wtr)->mole_fr = (base + H2O)->mole_fr ;
	condensed(base + Wtr, atm) ;
	if ((base + Wtr)->gibbs_free < (base + H2O)->gibbs_free) {
		printf("Liquid water exists!\n") ;
	}
	else printf("Water vapor exists!\n") ;
}

if (sw_Cgr_set) {
	(base + Cgr)->mole_fr = (base + C)->mole_fr ;
	condensed(base + Cgr, atm) ;
	if ((base + Cgr)->gibbs_free < (base + C)->gibbs_free) {
		printf("Graphite exists!\n") ;
	}
	else printf("Carbon gas exists!\n") ;
}

/* Check accuracy and restore the constant input state */
switch(i_state) {
case P_T:	/* Constant Pressure and Temperature */
	acc_1 = fabs((orig.p - atm->p) / orig.p) ;
	acc_2 = fabs((orig.t - atm->t) / orig.t) ;
	atm->p = orig.p ;
	atm->t = orig.t ;
	break ;

case P_H:	/* Constant Pressure and Enthalpy */
	acc_1 = fabs((orig.p - atm->p) / orig.p) ;
	acc_2 = fabs((orig.h_298_15 - atm->h_298_15) / orig.h_298_15) ;
	atm->p = orig.p ;
	atm->h_298_15 = orig.h_298_15 ;
	break ;

case P_S:	/* Constant Pressure and Entropy */
	acc_1 = fabs((orig.p - atm->p) / orig.p) ;
	acc_2 = fabs((orig.s - atm->s) / orig.s) ;
	atm->p = orig.p ;
	atm->s = orig.s ;
	break ;

case Rho_T:	/* Constant Density and Temperature */
	acc_1 = fabs((orig.rho - atm->rho) / orig.rho) ;
	acc_2 = fabs((orig.t - atm->t) / orig.t) ;
	atm->rho = orig.rho ;
	atm->t = orig.t ;
	break ;

case Rho_H:	/* Constant Density and Enthalpy */
	acc_1 = fabs((orig.rho - atm->rho) / orig.rho) ;
	acc_2 = fabs((orig.h_298_15 - atm->h_298_15) / orig.h_298_15) ;
	atm->rho = orig.rho ;
	atm->h_298_15 = orig.h_298_15 ;
	break ;

case Rho_S:	/* Constant Density and Entropy */
	acc_1 = fabs((orig.rho - atm->rho) / orig.rho) ;
	acc_2 = fabs((orig.s - atm->s) / orig.s) ;
	atm->rho = orig.rho ;
	atm->s = orig.s ;
	break ;

default:
	return(Unknown_state) ;
} /* end of accuracy check Switch block */
if (sw_pegged) return(Out_of_limits) ;
if (acc_1 > acc_2) accuracy = acc_1 ;
else accuracy = acc_2 ;
if (accuracy > 1.0e-3) return(Bad_accuracy) ;

return(i_err) ; /* normal return */
} /* --- end of the "Lewis" function --- */

int solve_lewis(int sw_reduce, double a[][N_species], double b[], int n) 
/***************************************************************/
/*                                                             */
/*           --- Solve a Square Linear System ---              */
/*                                                             */
/* Don't improve this numerical method if you want to continue */
/* reproducing the results of the original Gordon and McBride  */
/* code.  Yes!  There are better linear system solution        */
/* methods out there but don't use them unless you have a      */
/* better accuracy standard (there isn't one) than the         */
/* original Gordon and McBride code.                           */
/*                                                             */
/*       National Aeronautics and Space Administration         */
/*       No Copyright claimed in USA under Title 17-USC        */
/*                All other rights reserved                    */
/*                                                             */
/*            Version:  Mk 1.0,  16 January 1999               */
/*               Written by Gary A. Allen, Jr.                 */
/*                                                             */
/***************************************************************/
{
register int i, j, k ;
int n_1 ;

n_1 = n - 1;

if (sw_reduce) {
	for (i = 0; i < n_1; ++i) {
		/* do not replace "==" with "fabs() < Very_small" */
		if (a[i][i] == 0.0) return(TRUE) ;
		a[i][i] = 1.0 / a[i][i] ;
		for (j = i + 1; j < n; ++j) a[j][i] *= a[i][i];
		for (j = i + 1; j < n; ++j) for (k = i + 1; k < n; ++k)
			a[k][j] -= a[k][i] * a[i][j] ;
	}
	/* do not replace "==" with "fabs() < Very_small" */
	if (a[n_1][n_1] == 0.0) return(TRUE) ;
	a[n_1][n_1] = 1.0 / a[n_1][n_1] ;
}

for (i = 0; i < n_1; ++i) {
	for (j = i + 1; j < n; ++j) {
		b[j] -= a[j][i] * b[i];
	}
}
for (i = n_1; i >= 0; --i) {
	b[i] *= a[i][i] ;
	for (j = 0; j < i; ++j) b[j] -= a[j][i] * b[i];
}

return(FALSE) ;
} /* --- end of the "solve_lewis" subroutine --- */

int err_lewis(int sw_suppress_warning, int ret_code, int i_state, Thermo *atm)
/***************************************************************************/
/*                                                                         */
/*              --- Lewis Code Error Message Routine ---                   */
/*                                                                         */
/*             National Aeronautics and Space Administration               */
/*             No Copyright claimed in USA under Title 17-USC              */
/*                      All other rights reserved                          */
/*                                                                         */
/*                  Version:  Mk 1.0,  28 January 1999                     */
/*                     Written by Gary A. Allen, Jr.                       */
/*                                                                         */
/***************************************************************************/
{
int sw_abort ;

switch(ret_code) {
case Bad_temperature:
	printf("!!!! Temperature of %f deg.K is spurious.\n", atm->t) ;
	sw_abort = TRUE ;
	break ;

case Bad_density:
	printf("!!!! Density of %f kg/m^3 is spurious.\n", atm->rho) ;
	sw_abort = TRUE ;
	break ;

case Bad_pressure:
	printf("!!!! Pressure of %e Pascals is spurious.\n", atm->p) ;
	sw_abort = TRUE ;
	break ;

case Logic_error:
	printf("!!!! Logic error.\n") ;
	sw_abort = TRUE ;
	break ;

case Singular:
	printf("!!!! Matrix in \"Lewis\" is singular.\n") ;
	sw_abort = TRUE ;
	break ;

case Out_of_limits:
	sw_abort = FALSE ;
	if (sw_suppress_warning) break ;
	printf("Warning:  Temperature is beyond model limits!\n\n") ;
	break ;

case Bad_accuracy:
	sw_abort = FALSE ;
	if (sw_suppress_warning) break ;
	printf("Warning:  The result is inaccurate!\n\n") ;
	break ;

case Junk_result:
	sw_abort = FALSE ;
	if (sw_suppress_warning) break ;
	printf("Warning:  The result is NOT correct!\n\n") ;
	break ;

case Not_converged:
	sw_abort = FALSE ;
	if (sw_suppress_warning) break ;
	printf("Warning:  Solution did NOT converge.\n\n") ;
	break ;

case Too_many_elements:
	printf("!!!! Too many elements.\n") ;
	sw_abort = TRUE ;
	break ;

case Unknown_state:
	printf("!!!! Unknown thermodynamic state of %d selected.\n", i_state) ;
	sw_abort = TRUE ;
	break ;

case Bad_phase_input:
	printf("!!!! Inconsistant gas/liquid/solid input.\n") ;
	sw_abort = TRUE ;
	break ;

case Bad_ion_input:
	printf("!!!! Inconsistant ionized species selected.\n") ;
	sw_abort = TRUE ;
	break ;

default:
	printf("!!!! Unknown error code: %d\n", ret_code) ;
	sw_abort = TRUE ;
} /* end of Switch block */ 

if (sw_abort) printf("!!!! Function \"Lewis\" aborted!\n") ;

return(sw_abort) ;
} /* --- end of the "err_lewis" function --- */

int mixture(int sw_sum, int n_spc, Therm_data *base[], Thermo *atm)
/******************************************************/
/*                                                    */
/*   --- Determine Mixture Enthalpy and Entropy ---   */
/*                                                    */
/*    National Aeronautics and Space Administration   */
/*    No Copyright claimed in USA under Title 17-USC  */
/*             All other rights reserved              */
/*                                                    */
/*          Version:  Mk 1.0,  12 June 1999           */
/*            Written by Gary A. Allen, Jr.           */
/*                                                    */
/******************************************************/
{
register int i, n ;

int sw_pegged, bnd, num_of_levels  ;

double entropy, rmlfr, log_rmlfr, enthalpy, mw, t  ;
double p, arg, q_trans, denom, q_int, log_p, log_t ;
double q_int_bar, q_bar_q, q_int_bar_2, h, del_h_298_15 ;

static double p_ref = 1.0e5 ;	/* pascals, reference condition */
static double t_ref = 298.15 ;	/* deg.K,   reference condition */

static Level n_level[37] = N_level ;	 /* Atomic nitrogen [neutral]      */
static Level el_level = {2.0, 0.0} ;	 /* Electron                       */
static Level npp_level[30] = Npp_level ; /* Double ionized atomic nitrogen */
static Level opp_level[30] = Opp_level ; /* Double ionized atomic oxygen   */
static Level cpp_level[30] = Cpp_level ; /* Double ionized atomic carbon   */

Therm_data *s ;

Level *level ;
static Level *level_table[6] = {
	n_level,
	n_level,
	&el_level,
	npp_level,
	opp_level,
	cpp_level
} ;

/* Function Prototypes */
int condensed(Therm_data*, Thermo*) ;

/* Initialize input parameters */
mw = atm->mw / 1000.0 ;
t = atm->t ;
p = atm->p ;
log_p = log(p / p_ref) ;
sw_pegged = FALSE ;

/* Initialize mixture entropy and enthalpy */
entropy = 0.0 ;
enthalpy = 0.0 ;

/* Sum all of the species enthalpies and entropies to find the mixture */
for (n = 0; n < n_spc; ++n) {

s = base[n] ;

rmlfr = s->mole_kg_mix * mw;
if (rmlfr <= 1.0e-10) log_rmlfr = -23.02585093 ;
else log_rmlfr = log(rmlfr) ;

if (s->part) {
	
	arg = 2.0 * Pi * s->mole_wt * 1000.0 *  Boltz / (Avgdr * Planck
		* Planck) ;

	q_trans = pow(arg * t_ref, 1.5) * Boltz * t_ref / p_ref ;
	denom = Boltz * t_ref / (Planck * Spd_of_lght * 100.0) ;

	/* Calculate reference temperature enthalpy using only the   */
	/* first term of the partition function (all that is needed) */
	level = level_table[s->part] ;
	q_int = level->g * exp(-level->nu / denom) ; 
	q_int_bar = level->g * (level->nu / denom) * exp(-level->nu / denom) ; 
	q_bar_q = q_int_bar / q_int ;
	del_h_298_15 = R_gurv * t_ref * ((5.0 / 2.0) + q_bar_q) ;

	q_trans = pow(arg * t, 1.5) * Boltz * t / p ;
	denom = Boltz * t / (Planck * Spd_of_lght * 100.0) ;

	q_int = 0.0 ;
	q_int_bar = 0.0 ;
	q_int_bar_2 = 0.0 ;
	num_of_levels = 27 ; /* temporary */
	num_of_levels = s->num_of_levels ;
	for (i = 0; i < num_of_levels; i++) {
		level = level_table[s->part] + i ;

		if ((arg = level->nu / denom) > 500.0) arg = 500.0 ;
		q_int += level->g * exp(-arg) ; 
		q_int_bar += level->g * arg * exp(-arg) ; 
		q_int_bar_2 += level->g * arg * arg * exp(-arg) ; 
	}
	q_bar_q = q_int_bar / q_int ;
	h = R_gurv * t * ((5.0 / 2.0) + q_bar_q) ;
	s->phi = R_gurv * log(q_int * q_trans) ;
	s->q_int = q_int ;
	s->cp = R_gurv * ((5.0 / 2.0) + (q_int_bar_2 / q_int)
		- q_bar_q * q_bar_q) / R_lewis ;
	s->entropy = ((s->phi + h / t) / R_lewis) - log_p - log_rmlfr ;
	s->enthalpy = (h + s->del_h_298_15 - del_h_298_15) / (t * R_lewis) ;

} /* end of partition function routine */

/* Partition function not tabulated so use the Gordon and McBride fits */
else { 

	if ((s->type == Solid)||(s->type == Liquid)) {
		sw_pegged = condensed(s, atm) ;
	}
	else {
	/* Select temperature band */
	for(;;) {
		/* Deal with temperature below the database */
		if (t <= 100.0) {
			/* t = 100.0 ; */
			t = fabs(t) ;
			bnd = 0 ;
			sw_pegged = TRUE ;
			break ;
		}
		/* Deal with temperature above the database */
		if (t >= 20000.0) {
			t = 20000.0 ;
			bnd = 2 ;
			sw_pegged = TRUE ;
			break ;
		}

		/* Determine temperature band and set reusable parameters */
		/* Band:  100.0 < t <= 1000.0 */
		if (t <= 1000.0) {
			bnd = 0 ;
			break ;
		}

		/* Band:  1000.0 < t <= 6000.0 */
		if (t <= 6000.0) {
			bnd = 1;
			break ;
		}

		/* Band:  6000.0 < t <= 20000.0 */
		bnd = 2;
		break ;
		/*NOTREACHED*/
	} /* end of temperature band selection For loop */

	log_t = log(t) ;

	s->enthalpy = (-s->a[0][bnd] / t + s->a[1][bnd] * log_t
		+ s->b[0][bnd]) / t + s->a[2][bnd]
		+ t * (s->a[3][bnd] / 2.0 + t * (s->a[4][bnd] / 3.0
		+ t * (s->a[5][bnd] / 4.0 + t * s->a[6][bnd] / 5.0))) ;

	s->entropy = ((-s->a[0][bnd] / (t * 2.0) - s->a[1][bnd]) / t
		+ s->a[2][bnd] * log_t + t * (s->a[3][bnd]
		+ t * (s->a[4][bnd] / 2.0 + t * (s->a[5][bnd] / 3.0
		+ t * s->a[6][bnd] / 4.0))) + s->b[1][bnd]) - log_p
		- log_rmlfr ;

	s->cp = (s->a[0][bnd] / t + s->a[1][bnd]) / t + s->a[2][bnd]
		+ t * (s->a[3][bnd] + t * (s->a[4][bnd] + t * (s->a[5][bnd]
		+ t * s->a[6][bnd]))) ;
	}
}
if (sw_sum) {

	if ((s->type == Solid)||(s->type == Liquid)) {
		enthalpy += s->mole_kg_mix * 1000.0
			* R_lewis * t * s->enthalpy ;
		entropy += s->mole_kg_mix * 1000.0 * R_lewis * s->entropy ;
	}
	else {
		enthalpy += s->mole_kg_mix * R_lewis * t * s->enthalpy ;
		entropy += s->mole_kg_mix * R_lewis * s->entropy ;
	}

#ifdef DEBUG
	printf("%f, n = %d\n", s->mole_wt * 1000.0, n) ;
	printf("              h %f\n", s->enthalpy) ;
	if ((s->type == Solid)||(s->type == Liquid)) {
		printf("              s %f\n", s->entropy) ;
	}
	else {
		printf("              s %f\n", s->entropy + log_p
			+ log_rmlfr) ;
	}
	printf("             Cp %f\n", R_lewis * s->cp) ;
	printf("----------------------\n") ;
#endif /* DEBUG */
}
else s->gibbs_free = s->enthalpy - s->entropy ;
} /* end species FOR loop */

/* Output results to data structure */
if (sw_sum) {
	atm->h_298_15 = enthalpy ;
	atm->s = entropy ;
}

return(sw_pegged) ;
} /* --- end of the "mixture" subroutine --- */

int make_mollier(int i_planet, Thermo h_p[][Ml_j],
	Thermo rho_t[][Ml_j], Thermo p_t[][Ml_j])
/************************************************************/
/*                                                          */
/*        National Aeronautics and Space Administration     */
/*        No Copyright claimed in USA under Title 17-USC    */
/*                 All other rights reserved                */
/*                                                          */
/*            Version:  Mk 1.1,  11 November 2009           */
/*               Written by Gary A. Allen, Jr.              */
/*                                                          */
/************************************************************/

{
register int n ;
register int i = 0, j = 0 ;	/* suppress Gcc bogus warnings */

int i_state = 0, i_err = 0 ;	/* suppress Gcc bogus warnings */
int iter = 0 ;			/* suppress Gcc bogus warnings */

double t ;
double schmidt_number = 0.7 ;
double enthalpy = 0.0 ;		/* suppress Gcc bogus warnings */
double h_step = 0.0 ;		/* suppress Gcc bogus warnings */
double del_p = 0.0 ;		/* suppress Gcc bogus warnings */
double rho_exp_start = 0.0 ;	/* suppress Gcc bogus warnings */
double h_ref = 0.0  ;		/* suppress Gcc bogus warnings */

Thermo ref ;
Thermo *atm = 0 ; 		/* suppress Gcc bogus warnings */

Formula frm[N_species] ;
static Formula formula[N_planet_list][N_species] = Planet_atmos_mole_fraction ;

/* Function Prototypes */
int lewis(int, int*, Thermo*) ;
int err_lewis(int,int,int,Thermo*) ;

switch(i_planet) {
case Titan:
	h_step = 5.0e5 ;
	rho_exp_start = -9.0 ;
	del_p = 0.15 ;
	schmidt_number = 0.7 ;
	break ;

case Triton:
	h_step = 5.0e5 ;
	rho_exp_start = -9.0 ;
	del_p = 0.2 ;
	schmidt_number = 0.7 ;
	break ;

case Venus:
	h_step = 5.0e5 ;
	rho_exp_start = -8.0 ;
	del_p = 0.2 ;
	schmidt_number = 0.7 ;
	break ;

case Earth:
	h_step = 5.0e5 ;
	rho_exp_start = -9.0 ;
	del_p = 0.2 ;
	schmidt_number = 0.6 ;
	break ;

case Mars:
	h_step = 5.0e5 ;
	rho_exp_start = -9.0 ;
	del_p = 0.18 ;
	schmidt_number = 0.7 ;
	break ;

case Jupiter:
	h_step = 5.0e6 ;
	rho_exp_start = -9.0 ;
	del_p = 0.2 ;
	schmidt_number = 0.4 ;
	break ;

case Saturn:
	h_step = 5.0e6 ;
	rho_exp_start = -9.0 ;
	del_p = 0.2 ;
	schmidt_number = 0.4 ;
	break ;

case Uranus:
	h_step = 5.0e6 ;
	rho_exp_start = -9.0 ;
	del_p = 0.2 ;
	schmidt_number = 0.4 ;
	break ;

case Neptune:
	h_step = 5.0e6 ;
	rho_exp_start = -9.0 ;
	del_p = 0.2 ;
	schmidt_number = 0.4 ;
	break ;

case Pluto:
	h_step = 5.0e5 ;
	rho_exp_start = -9.0 ;
	del_p = 0.2 ;
	schmidt_number = 0.7 ;
	break ;

default:
	printf("!!!! Logic error in \"Make_mollier\"\n") ;
	return(TRUE) ;
} /* end planet selection Switch block */

/* Load scratch structure */
for (n = 0; n < N_species ; ++n) {
	frm[n].species = formula[i_planet][n].species ;
	frm[n].mole_fract = formula[i_planet][n].mole_fract ;
}

i_state = P_T ;
ref.t = 298.15 ;
ref.p = 1.0e5 ;
ref.frm = frm ;
lewis(i_state, &iter, &ref) ;
h_ref = ref.h_cnv - ref.h_298_15 ;

/* Constant Pressure and Enthalpy */
i_state = P_H ;
for (i = 0; i < Ml_i; i++) {
	enthalpy = ((double) i) * h_step ;

	for (j = 0; j < Ml_j; j++) {
		atm = h_p[i] + j ;
		atm->h = enthalpy ;
		atm->h_298_15 = enthalpy - h_ref ;
		atm->p = pow(10.0, (((double)j) * del_p) - 2.0) ;
		for (n = 0; n < N_species ; ++n) {
			frm[n].species = formula[i_planet][n].species ;
			frm[n].mole_fract = formula[i_planet][n].mole_fract ;
		}
		atm->frm = frm ;
		
		if ((i_err = lewis(i_state, &iter, atm))) {
			if (err_lewis(TRUE, i_err, i_state, atm)) {
				printf("!!!! h = %8.1e, press = %8.1e pa\n",
					atm->h_298_15, atm->p) ;  
				return(TRUE) ;
			}
		}
		atm->frm = (Formula *)0 ;
		atm->lewis_nmbr /= schmidt_number ;
	}
} /* end of "i", P_H enthalpy FOR loop */

/* Constant Density and Temperature */
i_state = Rho_T ;
for (i = 0; i < Ml_i; i++) {
	t = 200.0 + ((double) i) * 55.0 ;

	for (j = 0; j < Ml_j; j++) {
		atm = rho_t[i] + j ;
		atm->t = t ;
		atm->rho = 1.2250 * pow(10.0, (((double)j) * 0.2) 
			+ rho_exp_start) ;
		for (n = 0; n < N_species ; ++n) {
			frm[n].species = formula[i_planet][n].species ;
			frm[n].mole_fract = formula[i_planet][n].mole_fract ;
		}
		atm->frm = frm ;

		if ((i_err = lewis(i_state, &iter, atm))) {
			if (err_lewis(TRUE, i_err, i_state, atm)) {
				printf("!!!! t = %7f, rho = %8.1e kg/m^3\n",
					atm->t, atm->rho) ;  
				return(TRUE) ;
			}
		}

		atm->frm = (Formula *)0 ;
		atm->h = atm->h_298_15 + h_ref ;
		atm->lewis_nmbr /= schmidt_number ;
	}
} /* end of "i", Rho_T temperature FOR loop */

/* Constant Pressure and Temperature */
i_state = P_T ;
for (i = 0; i < Ml_i; i++) {
	t = 200.0 + ((double) i) * 55.0 ;

	for (j = 0; j < Ml_j; j++) {
		atm = p_t[i] + j ;
		atm->p = pow(10.0, (((double)j) * del_p) - 2.0) ;
		atm->t = t ;
		for (n = 0; n < N_species ; ++n) {
			frm[n].species = formula[i_planet][n].species ;
			frm[n].mole_fract = formula[i_planet][n].mole_fract ;
		}
		atm->frm = frm ;
		if ((i_err = lewis(i_state, &iter, atm))) {
			if (err_lewis(TRUE, i_err, i_state, atm)) {
				printf("!!!! t = %7f deg.K, press = %8.1e pa\n",
					atm->t, atm->p) ;  
				return(TRUE) ;
			}
		}
		atm->frm = (Formula *)0 ;
		atm->h = atm->h_298_15 + h_ref ;
		atm->lewis_nmbr /= schmidt_number ;
	}
} /* end of "i", P_T temperature FOR loop */

return(FALSE) ;
} /* --- end of function "Make_mollier" --- */

int condensed(Therm_data *s, Thermo *atm)
/******************************************************/
/*                                                    */
/*   Determine Condensed Species Enthalpy & Entropy   */
/*                                                    */
/*    National Aeronautics and Space Administration   */
/*   No Copyright claimed in USA under Title 17-USC   */
/*              All other rights reserved             */
/*                                                    */
/*          Version:  Mk 1.0,  12 April 2000          */
/*            Written by Gary A. Allen, Jr.           */
/*                                                    */
/******************************************************/
{
int sw_pegged, bnd  ;

double t, log_t ;

/* Initialize input parameters */
t = atm->t ;
sw_pegged = FALSE ;

/* Partition function not tabulated so use the Gordon and McBride fits */
/* Select temperature band */
for(;;) {
	/* Deal with temperature below the database */
	if (t <= 200.0) {
		t = 200.0 ;
		bnd = 0 ;
		sw_pegged = TRUE ;
		break ;
	}
	/* Deal with temperature above the database */
	if (t >= 5000.0) {
		t = 5000.0 ;
		bnd = 2 ;
		sw_pegged = TRUE ;
		break ;
	}

	/* Determine temperature band and set reusable parameters */
	/* Band:  200.0 < t <= 600.0 */
	if (t <= 600.0) {
		bnd = 0 ;
		break ;
	}

	/* Band:  600.0 < t <= 2000.0 */
	if (t <= 2000.0) {
		bnd = 1;
		break ;
	}

	/* Band:  2000.0 < t <= 5000.0 */
	bnd = 2;
	break ;
	/*NOTREACHED*/
} /* end of temperature band selection For loop */

log_t = log(t) ;

s->enthalpy = (-s->a[0][bnd] / t + s->a[1][bnd] * log_t
	+ s->b[0][bnd]) / t + s->a[2][bnd]
	+ t * (s->a[3][bnd] / 2.0 + t * (s->a[4][bnd] / 3.0
	+ t * (s->a[5][bnd] / 4.0 + t * s->a[6][bnd] / 5.0))) ;

/* log of pressure and mole fraction not included because this is solid */ 
s->entropy = (-s->a[0][bnd] / (t * 2.0) - s->a[1][bnd]) / t
	+ s->a[2][bnd] * log_t + t * (s->a[3][bnd]
	+ t * (s->a[4][bnd] / 2.0 + t * (s->a[5][bnd] / 3.0
	+ t * s->a[6][bnd] / 4.0))) + s->b[1][bnd] ;

s->cp = (s->a[0][bnd] / t + s->a[1][bnd]) / t + s->a[2][bnd]
	+ t * (s->a[3][bnd] + t * (s->a[4][bnd] + t * (s->a[5][bnd]
	+ t * s->a[6][bnd]))) ;

#ifdef DEBUG
printf("Condensed  h %f\n", s->enthalpy) ;
printf("             %f\n", R_lewis * s->enthalpy) ;
printf("Condensed  s %f\n", s->entropy) ;
printf("             %f\n", R_lewis * s->entropy) ;
printf("Condensed Cp %f\n", R_lewis * s->cp) ;
#endif /* DEBUG */

s->gibbs_free = s->enthalpy - s->entropy ;

return(sw_pegged) ;
} /* --- end of the "condensed" subroutine --- */

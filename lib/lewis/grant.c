#include "var.h"
#include "lewis.h"

/* MAIN PROGRAM --- MAIN PROGRAM --- MAIN PROGRAM --- MAIN PROGRAM */
int main (int argc, char *argv[])
/***************************************************************************/
/*                                                                         */
/*          -----  Gibbs Free Energy Minimization Program  -----           */
/*                                                                         */
/*   This subroutine calculates the equilibrium composition of a gas       */
/*   by Gibbs free energy minimization.  The technique is that of          */
/*   Gordon and McBride, NASA SP-273                                       */
/*                                                                         */
/*                Version:  Mk 1.2,    12 November 2009                    */
/*                      (c) Copyright 2009 by NASA                         */
/*                     Written by Gary A. Allen, Jr.                       */
/*                                                                         */
/***************************************************************************/

{
register int i, n, is ;

int i_err, iter, ret_code ;
int i_state = 0 ; /* removes a warning */

double mole_fr, mole_fr_sum ;
double arg_1, arg_2, mass_fract ;

char buff[132] ;

FILE *fp_in = 0 ;

static Tagged_parse species[N_species] = Species_list ;
Tagged_parse state[N_state] = State_list ;
Formula frm[N_species + 1] ;

Therm_data spc[N_species] = {
	Thermo_data_base_1,
	Thermo_data_base_2,
	Thermo_data_base_3
} ;

Thermo atm ;

/* Function Prototypes */
int lewis(int, int*, Thermo*) ;
int err_lewis(int,int,int,Thermo*) ;

if (argc != 2) {
     printf("!!!! Wrong number of arguments on the input line. \n");
     printf("argc = %d argv[1] = %s\n", argc, argv[1]) ;
     return(0) ;
}
	
if ((fp_in = fopen(argv[1],"r")) == 0) {
     printf("!!!! Could not open input file:  \"%s\"\n", argv[1]);
     return(0) ;
}

ret_code = fscanf(fp_in, "%s %lf %lf", buff, &arg_1, &arg_2) ;
Strip(fp_in) ;

/* Find the species in the species list and lookup the op_code */
for (i = 0; i < N_state; i++) {
	if (!strncmp(buff, state[i].mnemonic, Max_size)) {
		i_state = state[i].op_code ;
		break ;
	}
}

/* Abort if the For loop hit the limit */
if (i >= N_state) {
	printf("!!!! Unknown thermodynamic state: %s\n", buff) ;
	printf("!!!! Program aborted!\n") ;
	return(0) ;
}

/* Print heading */
printf("Free energy minimization\n\n") ;
printf("%s held constant.\n\n", state[i_state].tag) ;
printf("  Initial Value\n") ;
printf("  -------------\n") ;

/* Print initial thermodynamic state */
switch(i_state) {
case Rho_T:
	atm.rho = arg_1 ;
	atm.t = arg_2 ;
	printf("density             = %e kg/m^3\n", atm.rho) ;
	printf("temperature         = %f deg.K\n\n", atm.t) ;
	break ;

case Rho_H:
	atm.rho = arg_1 ;
	atm.h_298_15 = arg_2 ;
	printf("density             = %e kg/m^3\n", atm.rho) ;
	printf("mixture enthalpy    = %e Joules/kg\n\n", atm.h_298_15) ;
	break ;

case Rho_S:
	atm.rho = arg_1 ;
	atm.s = arg_2 ;
	printf("density             = %e kg/m^3\n", atm.rho) ;
	printf("mixture entropy     = %e Joules/kg-mole\n\n", atm.s) ;
	break ;

case P_T:
	atm.p = arg_1 ;
	atm.t = arg_2 ;
	printf("pressure            = %e Pascals\n", atm.p) ;
	printf("temperature         = %f deg.K\n\n", atm.t) ;
	break ;

case P_H:
	atm.p = arg_1 ;
	atm.h_298_15 = arg_2 ;
	printf("pressure            = %e Pascals\n", atm.p) ;
	printf("mixture enthalpy    = %e Joules/kg\n\n", atm.h_298_15) ;
	break ;

case P_S:
	atm.p = arg_1 ;
	atm.s = arg_2 ;
	printf("pressure            = %e Pascals\n", atm.p) ;
	printf("mixture entropy     = %e Joules/kg-mole\n\n", atm.s) ;
	break ;

default:
	printf("!!!! Logic error!\n") ;
	printf("!!!! Program aborted!\n") ;
	return(0) ;
} /* end of initial thermodynamic state Switch block */

for (n = 0; n < N_species; ++n) frm[n].mole_fract = 0.0 ;
mole_fr_sum = 0.0 ;
atm.frm = frm ;
for (n = 0; n < N_species; ++n) {
	if (EOF == fscanf(fp_in, "%s %lf", buff, &mole_fr)) break  ;
	Strip(fp_in) ;

	mole_fr_sum += mole_fr ;

	/* Find the species in the species list and lookup the op_code */
	for (i = 0; i < N_species; i++) {
		if (!strncmp(buff, species[i].mnemonic, Max_size)) {
			frm[n].species = species[i].op_code ;
			frm[n].mole_fract = mole_fr ;
			break ;
		}
	}

	/* Abort if the For loop hit the limit */
	if (i >= N_species) {
		printf("!!!! Unknown species: %s\n", buff) ;
		printf("!!!! Program aborted!\n") ;
		return(0) ;
	}
} /* end of species For loop */

/* Mark the end of the formula */
frm[n].species = End ;

/* Check that input mole fractions are sane */
if (fabs(1.0 - mole_fr_sum) > 1.0e-15) {
	printf("!!!! Input mole fractions must add up to one.\n") ;
	printf("!!!! Difference: %e\n", fabs(1.0 - mole_fr_sum)) ;
	printf("!!!! Program aborted!\n") ;
	return(0) ;
}

if ((i_err = lewis(i_state, &iter, &atm))) {
	if (err_lewis(FALSE, i_err, i_state, &atm)) return(0) ;
}
else printf("Solution converged in %d steps.\n\n", iter) ;

/*   Print out the species */
for (n = 0; n < N_species; ++n) {
	if (frm[n].species == End) break ;
	mass_fract = frm[n].mole_fract * spc[frm[n].species].mole_wt * 1000.0
		/  atm.mw_1 ;

	printf("%23s mass fraction %e mole fraction %e\n",
		species[frm[n].species].tag, mass_fract, frm[n].mole_fract) ;
}
printf("\n") ;

printf("  Final Value\n") ;
printf("  -----------\n") ;
printf("pressure                 = %e Pascals\n", atm.p) ;
printf("density                  = %e kg/m^3\n", atm.rho) ;
printf("temperature              = %f deg.K\n", atm.t) ;
printf("molar mass               = %f\n", atm.mw) ;
printf("                           %f\n", atm.mw_1) ;
printf("speed of sound           = %f m/sec\n", atm.c) ;
printf("specific heat, Cp        = %f Joules/(kg deg.K)\n", atm.cp) ;
printf("isentropic exponent      = %f\n", atm.gamma) ;
printf("mixture enthalpy         = %e Joules/kg\n", atm.h_298_15) ;
printf("                           %e BTU/lbm\n", atm.h_298_15 * 4.299226e-4) ;
printf("mixture entropy          = %e Joules/kg-mole\n", atm.s) ;
printf("mixture viscosity        = %e kg/(m sec)\n", atm.mu) ;
printf("equilibrium conductivity = %e Watts/(m deg.K)\n", atm.cnd_equ) ;
printf("frozen conductivity      = %e Watts/(m deg.K)\n", atm.cnd_frz) ;
printf("frozen Prandtl number    = %f \n\n", atm.prandtl) ;

fclose(fp_in) ;
return(0) ;
} /* --- end of the Main routine --- */


#include "var.h"
#include "lewis.h"

/******************************************************************************/

int equilibrium_gas_
  (int*    iprint,            /* 0 => suppress printing; 1 => print details   */
   int*    istate,            /* Definition of which two state variables are  */
                              /* inputs:                                      */
                              /* 0 = Rho_T;  1 = P_H;    2 = P_T,             */
			      /* 3 = Rho_H;  4 = Rho_S;  5 = P_S              */
   double* state_1,           /* Corresponding input values                   */
   double* state_2,           /*                                              */
   int*    n_species,         /* Number of input gas species                  */
   int     ispecies[],        /* Corresponding "op. codes" in data structure  */
   double  mole_fractions[],  /* Input with starting guesses;                 */
                              /* output with equilibrium values               */
   double* pressure,          /* Output mixture properties in SI units        */
   double* density,           /*                                              */
   double* temperature,       /*                                              */
   double* enthalpy,          /*                                              */
   double* entropy,           /*                                              */
   double* speed_of_sound)    /* What else?                                   */

/* Returned value is 0 if no problem is detected, else the lewis status code. */

/* DESCRIPTION:

      This routine provides a means for transferring inputs & outputs between
   a Fortran 90 calling program and Gary Allen's implementation of CEA-type
   equilibrium gas calculations.  See also equilibrium_composition.f90.

   HISTORY:

   01/16/99- G.A.Allen, Jr. Function lewis.c and a driving program, used by the
   03/15/04                 nozzle throat conditions utilities via system calls.
   10/08/08- D.A.Saunders   Initial implementation of this C interface to avoid
   10/28/08                 system calls, with possible use by DPLR in mind too.
   10/29/08    "     "      Added iprint argument to control printing.
   09/06/11    "     "      Gary updated lewis.h, var.h, and lewis.c; lewis no
                            longer has Formula *f as an argument, and err_lewis
                            has a new suppress_warning argument.
   AUTHOR:

   David Saunders ELORET Corporation/NASA Ames Research Center, Mtn. View, CA */

/******************************************************************************/

{ /* Start "equilibrium_gas" function */

/* Local variables */

int	i, i_err, iter, n;
double	mass_fract, mole_fract, mole_fr_sum;

/* Local structures */

static Tagged_parse species[N_species] = Species_list;
Tagged_parse state[N_state] = State_list;

Formula frm[N_species + 1];

Therm_data spc[N_species] = {
	Thermo_data_base_1,
	Thermo_data_base_2,
	Thermo_data_base_3
			    };
Thermo	atm;

/* Function Prototypes */

int	lewis (int, int*, Thermo*);
int	err_lewis (int, int, int, Thermo*);

/* Execution */

/* printf ("iprint: %d\n", *iprint);       */
/* printf ("istate: %d\n", *istate);       */
/* printf ("n_species: %d\n", *n_species); */
/* printf ("state_1: %e\n", *state_1);     */
/* printf ("state_2: %e\n", *state_2);     */

   if (*iprint)
{
	/* printf ("\nFree energy minimization, case %d\n\n", *istate); */
	printf ("\n%s held constant.\n\n", state[*istate].tag);
	printf ("Initial Values\n");
	printf ("--------------\n");
}

/* Which two state variables are specified? */

switch (*istate)
{
case Rho_T:
	atm.rho = *state_1;
	atm.t   = *state_2;
	if (*iprint)
	{
		printf ("density          = %e kg/m^3\n",  atm.rho);
		printf ("temperature      = %f deg.K\n\n", atm.t);
	}
	break;

case P_H:
        atm.p        = *state_1;
        atm.h_298_15 = *state_2;
        if (*iprint)
        {
                printf ("pressure         = %e Pascals\n", atm.p);
                printf ("mixture enthalpy = %e Joules/kg\n\n", atm.h_298_15);
        }
        break;

case P_T:
        atm.p = *state_1;
        atm.t = *state_2;
        if (*iprint)
        {
                printf ("pressure         = %e Pascals\n", atm.p);
                printf ("temperature      = %f deg.K\n\n", atm.t);
        }
        break;

case Rho_H:
	atm.rho      = *state_1;
	atm.h_298_15 = *state_2;
	if (*iprint)
	{
		printf ("density          = %e kg/m^3\n",  atm.rho);
		printf ("mixture enthalpy = %e Joules/kg\n\n", atm.h_298_15);
	}
	break;

case Rho_S:
	atm.rho = *state_1;
	atm.s   = *state_2;
	if (*iprint)
	{
		printf ("density          = %e kg/m^3\n",  atm.rho);
		printf ("mixture entropy  = %e Joules/kg-mole\n\n", atm.s);
	}
	break;

case P_S:
	atm.p = *state_1;
	atm.s = *state_2;
	if (*iprint)
	{
		printf ("pressure         = %e Pascals\n", atm.p);
		printf ("mixture entropy  = %e Joules/kg-mole\n\n", atm.s);
	}
	break;

default:
	printf ("!!! Bad input value for istate: %d\n", *istate);
	printf ("!!! Aborting.\n");  fflush (stdout);
	return (0);

} /* end switch (*istate) */

/* Check that the input mole fractions are sane. */

for (n = 0; n < N_species; ++n) frm[n].mole_fract = 0.0;

atm.frm = frm;
mole_fr_sum = 0.0;
for (n = 0; n < *n_species; ++n)
{
	mole_fr_sum      += mole_fractions[n];
	frm[n].mole_fract = mole_fractions[n];
	frm[n].species    = ispecies[n];
/*      printf ("   n: %d, ispecies[n]: %d\n", n, ispecies[n]); */
}

/* printf ("   n: %d, End: %d\n", n, End); fflush (stdout); */

frm[n].species = End;     /* Mark the end of the formula */

if (fabs (1.0 - mole_fr_sum) > 1.0e-6)
{
	printf ("!!! Input mole fractions must add up to 1.\n");
       	printf ("!!! Sum found: %e\n", mole_fr_sum);
       	printf ("!!! Aborting.\n");  fflush (stdout);
       	return (0);
}

/* Calculate the equilibrium composition. */

if ((i_err = lewis (*istate, &iter, &atm)))
{
	if (err_lewis (FALSE, i_err, *istate, &atm)) return (i_err);
}

if (*iprint) printf ("Solution converged in %d steps.\n\n", iter);

for (n = 0; n < *n_species; ++n)
{
	mole_fract = frm[n].mole_fract;  i = frm[n].species;
	mole_fractions[n] = mole_fract;
	mass_fract        = mole_fract * spc[i].mole_wt * 1000./ atm.mw;

	if (*iprint) printf ("%23s mass fraction %e mole fraction %e\n",
		species[i].tag, mass_fract, mole_fract);
}

*pressure       = atm.p;
*density        = atm.rho;
*temperature    = atm.t;
*enthalpy       = atm.h_298_15;
*entropy        = atm.s;
*speed_of_sound = atm.c;

if (*iprint)
{
	printf ("\nFinal Values\n");
	printf ("------------\n");
	printf ("pressure                 = %e Pascals\n", atm.p);
	printf ("density                  = %e kg/m^3\n", atm.rho);
	printf ("temperature              = %f deg.K\n", atm.t);
	printf ("molar mass               = %f\n", atm.mw);
	printf ("speed of sound           = %f m/sec\n", atm.c);
	printf ("specific heat, Cp        = %f Joules/(kg deg.K)\n", atm.cp);
	printf ("isentropic exponent      = %f\n", atm.gamma);
	printf ("mixture enthalpy         = %e Joules/kg\n", atm.h_298_15);
	printf ("                           %e BTU/lbm\n", atm.h_298_15 *
		4.299226e-4);
	printf ("mixture entropy          = %e Joules/kg-mole\n", atm.s);
	printf ("mixture viscosity        = %e kg/(m sec)\n", atm.mu);
	printf ("equilibrium conductivity = %e Watts/(m deg.K)\n", atm.cnd_equ);
	printf ("frozen conductivity      = %e Watts/(m deg.K)\n", atm.cnd_frz);
	printf ("frozen Prandtl number    = %f \n\n", atm.prandtl);
	fflush (stdout);
}

return (0);

} /* End of equilibrium_gas function */

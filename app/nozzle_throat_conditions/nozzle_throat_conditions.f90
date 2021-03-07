!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program nozzle_throat_conditions
!
!  Description:
!  ============
!
!     Calculate boundary conditions at an arc-jet nozzle throat that correspond
!     to a choice of distributions of two real gas flow variables along a radius
!     of the (circular) throat, from center to wall.  Results are output in the
!     form of a pointwise boundary condition for the DPLR2D flow solver.
!
!  Approach:
!  =========
!
!     This is a more thorough implementation than the earlier THROAT_CONDITIONS,
!     which was limited to air, didn't provide species mass fractions, sought a
!     target equilibrium Mach number above 1 (to simulate a "frozen" Mach number
!     of 1), and treated a list of one or more target pairs of specific total
!     enthalpy and mass flow rate.
!
!     Here, equilibrium gas compositions are calculated with a C implementation
!     of the relevant portions of the CEA program by Gordon and McBride dated
!     May 20, 1998.  Gary Allen wrote the C package at NASA Ames Research Center
!     and continues to maintain it.  Initially, a C driving program was invoked
!     here via a system call whenever a gas composition was required.  This has
!     now been streamlined with a C interface routine (equilibrium_gas.c) even
!     though its inefficiency was not much of an issue.
!
!  Flow Specification Options:
!  ===========================
!
!     Gary/CEA provides 6 options for the given pair of flow specifications.
!     These are preserved as options here, but normal usage is likely to specify
!     other combinations, the first of which prompted this program.
!
!     Option   Flow variables specified          Units
!
!     Ht_MF    total stagnation enthalpy         J/kg
!              mass flux                         kg/(m^2.s)
!              (safeguarded Newton iteration)
!     Ht_Ru    total stagnation enthalpy         J/kg
!              mass flux                         kg/(m^2.s)
!              (2 nested 1-variable iterations, no longer recommended)
!     Rho_T    density                           kg/m^3
!              temperature                       K
!     Rho_H    density                           kg/m^3
!              mixture enthalpy                  J/kg
!     Rho_S    density                           kg/m^3
!              mixture entropy                   J/kg-mole
!     P_T      pressure                          Pa
!              temperature                       K
!     P_H      pressure                          Pa
!              mixture enthalpy                  J/kg
!     P_S      pressure                          Pa
!              mixture entropy                   J/kg-mole
!
!     Ht_MF and Ht_Ru Algorithms:
!     ---------------------------
!
!     These options determine throat conditions at a specified frozen Mach
!     number (not necessarily 1) for given [distributions of] total stagnation
!     (reservoir) enthalpy and mass flow rate per unit area.  If the species
!     mass fractions for starting guess pressure and temperature are treated as
!     known, the energy balance and mass balance equations are separable and may
!     be solved via inner and outer iterations (the Ht_Ru option).
!
!     Alternatively, the Ht_MF option treats the same case as two general
!     nonlinear equations, solved via a safeguarded 2-variable Newton iteration
!     with central difference derivatives as in the earlier THROAT_CONDITIONS.
!
!     Until both approaches were implemented, it was not clear which would cope
!     best with the limited precision available from the equilibrium composition
!     calculations.  The clear winner proves to be the Newton iteration, which
!     uses about one fourth of the function evaluations in spite of the central
!     differencing.
!
!     The equations at a point r are:
!
!        f(1)  =  ho(r)  -  t1 Cpbar T  -  evbar(T)  -  sum (ci hfi)  =  0
!
!        f(2)  =  rho u(r)  -  P Mf sqrt (gamf / (Rbar T))            =  0
!
!     where, at point r, the frozen flow variables are:
!
!        ho      =  mixture total stagnation enthalpy  =  0.5 u**2 +  h
!        h       =  mixture enthalpy per unit mass     =  P / rho  +  e
!        rho     =  mixture density
!        P       =  mixture pressure
!        T       =  mixture temperature (translational T = vibrational temp. Tv)
!        e       =  mixture internal energy per unit mass (H is per unit volume)
!        Mf      =  frozen Mach number (not necessarily 1)
!        af      =  frozen speed of sound
!        u       =  Mf af  =  Mf * sqrt (gamf Rbar T)
!        t1      =  1  +  0.5 (gamf - 1) Mf**2
!        ci      =  mass fraction of species i
!        hfi     =  reference enthalpy at O K for species i
!        Mi      =  molecular weight for species i
!        gamf    =  Cpbar / Cvbar
!        Cpbar   =  mixture specific heat at constant pr.  =  Rbar  +  Cvbar
!        Cvbar   =  mixture specific heat at constant vol. =  sum (ci (Cvi^/Mi))
!        Rbar    =  mixture gas constant                   =  sum (ci  R^ / Mi)
!        evbar   =  mixture vibrational energy             =  sum (ci (evi^/Mi))
!        evi^    =  R^ (thetavi / (exp (thetavi / Tv) - 1))
!        thetavi =  characteristic temperature for species i vibrational energy
!        Cvi^    =  1.5 R^ or 2.5 R^ for atomic or diatomic species resp.
!        R^      =  universal gas constant
!
!     The two variables to be solved for are T and P, with the mass fractions
!     essentially converging as well during the iterative solution.
!
!     Ht_Ru and Ht_MF Starting Guesses:
!     ---------------------------------
!
!     For a specified mass fraction of argon (say 10%), fractions for N2 and
!     O2 are readily derived from the values 0.767 and 0.233 for pure air, with
!     other mass fractions set to 0.
!
!     For the current point in space, the previous solution serves as a good
!     starting guess.  For the first point r in space, starting guesses for P
!     and T are implemented as follows:
!
!        P  =  C (2 / (gamf + 1)) ^ (gamf / (gamf - 1)) rho u (r) ho(r)^n
!        T  =    (2 / (gamf + 1)) ^ To
!
!     where
!
!        C     =  3.41831 (SI units)
!        n     =  0.397
!        gamf  =  1.15 - 1.18  (lower for higher ho)
!        To    =  6500 K  if  ho  <  15 MJ/kg  else
!              =  7000 K  if  ho  <  20 MJ/kg  else
!              =  7500 K  if  ho  >  20 MJ/kg
!
!  Sample control file (standard input):
!  =====================================
!
!     ========================================================================
!     Control file for NOZZLE_THROAT_CONDITIONS (2 flow variables prescribed)
!     ========================================================================
!     V1_V2 specification
!     ------------------------------------------------------------------------
!     Ht_MF         ! Ht_MF, Ht_Ru, Rho_T, Rho_H, Rho_S, P_T, P_H, or P_S
!     1.0           ! Frozen Mach number
!     ------------------------------------------------------------------------
!     Iteration controls
!     ------------------------------------------------------------------------
!     1             ! Iterate for target bulk enthalpy?   0 = no | 1 = yes
!     2.9E+07       ! Target bulk enthalpy, J/kg
!     1             ! Iterate for target bulk mass flow rate?  0 | 1
!     90.           ! Target bulk mass flow rate, kg/s
!     ------------------------------------------------------------------------
!     V_1 specifications
!     ------------------------------------------------------------------------
!     Linear        ! Uniform, Linear, Parabolic, Sinusoid, Gaussian, Dataset
!     0.   3.E+07   ! (r, V1) at nozzle center-line
!     1.   3.E+07   ! (r, V1) at nozzle wall
!     999.          ! "Width" of distribution (usage depends on profile type)
!     none          ! Dataset file name or "none" or Sigmoid steepness
!     ------------------------------------------------------------------------
!     V_2 specifications
!     ------------------------------------------------------------------------
!     Uniform       ! Profile shape as for V_1 (which now includes Lorentzian)
!     0.   0.1      ! (r, V2) at nozzle center-line
!     1.   0.1      ! (r, V2) at nozzle wall
!     999.          ! "Width"  (most profiles)
!     none          ! Dataset file name or "none" or Sigmoid steepness
!     ------------------------------------------------------------------------
!     Mixture specifications
!     ------------------------------------------------------------------------
!     13            ! # species
!     N2   0.6903   ! Species names (DPLR) and mass fraction starting guesses
!     O2   0.2097
!     NO   0.0
!     N    0.0
!     O    0.0
!     Ar   0.1000
!     Ar+  0.0
!     N2+  0.0
!     O2+  0.0
!     NO+  0.0
!     N+   0.0
!     O+   0.0
!     e    0.0
!     ------------------------------------------------------------------------
!     Output specifications
!     ------------------------------------------------------------------------
!     21             ! # uniform pts. for initial BC evals.; < 0 for diagnostics
!     grid.dat       ! File containing target radii for interpolating to | none
!     61             ! DPLR's input profile BC code:  60 | 61 | 62
!     conditions.f   ! File name for output throat conditions (PLOT3D function)
!     conditions.dat ! Plottable data (Tecplot ASCII)
!
!  Output format (DPLR *.pbca file):
!  =================================
!
!     Initially, we write just the PLOT3D-type function file portion.
!
!  Usage Notes:
!  ============
!
!     (0) Run the code like this:  nozzle_throat_conditions < xxx.inp > xxx.log
!
!     (1) As indicated above, the Newton iteration of method Ht_MF is preferable
!         to the nested 1-variable iterations of method Ht_Ru.  Excessive
!         "step-halving iteration failed" messages normally mean something
!         unreasonable is being requested.
!     (2) When the Sigmoid shape function option was installed, the extra
!         control variable it needs threatened inclusion of an extra control
!         for all shape functions.  However, this was avoided by making use of
!         the (character) file name control provided for shapes defined by
!         discretized points.
!     (3) The Lorentz profile tends toward zero less rapidly than the Gaussian,
!         and is more peaky at r = 0 for the same "width" parameter.
!
!  History:
!  ========
!
!     06/09/04  DKP/DAS  Earlier THROAT_CONDITIONS:   Dinesh Prabhu/D. Saunders.
!     12/19/05  TG /DAS  Initial NOZZE_THROAT_CONDITIONS (Ht_Ru option).
!     12/20/05  TG /DAS  Added the Ht_MF alternative (Newton iteration) to the
!                        Ht_Ru option (non-derivative zero finder). Newton wins.
!     01/19/05  TG /DAS  Tahir noticed some glitches in the Tec variable names.
!     11/28/06  DKP/DAS  One test on "diagnostics" treated it as an integer.
!                        Lowered the convergence tolerance from 1.e-5 to 1.e-6.
!     11/29/06   "   "   Put the tolerance back to 1.e-5.  It seems the "Gary"
!                        part of the scheme can give different results for
!                        different initial estimates of the mass fractions.
!                        This needs to be investigated.
!     03/16/07  D.Prabhu Dinesh suggested the sigmoid shape function, which
!               DAS      required another control dubbed "steepness", entered
!                        here via the "dataset" file name inputs.
!     03/28/07   "   "   Scaling the variables (P and T) and the residuals
!                        (f(1) and f(2)) allows use of 1.e-6 tolerance on ||f||
!                        and overcomes mysteriously poor convergence and ragged
!                        pressures observed in a low-pressure case from Dinesh.
!     04/11/07   "   "   Bulk mass flow rate and bulk enthalpy options added;
!                        frozen Mach number is not assumed to be 1 (new input).
!     05/02/07   "   "   Uniform profiles need to work with the wall value, not
!                        the center value, because of the bulk flow options.
!     06/14/07   "   "   Added Gaussian profile option.
!     07/05/07   "   "   Added Lorentz  profile option.
!     10/29/08     DAS   Replaced the system call in equilibrium_composition.f90
!                        to Gary's CEA driving program with a call to interface
!                        routine equilibrium_gas.c, without change here except
!                        to update the documentation and to control the extra
!                        printing to standard output.
!     09/07/11      "    Maria Pulsonetti needed CO2, CO, C species in the top
!     as in the 3D       level tables (already handled at the lower levels).
!     code update,       WARNING:  Polyatomic species have more than one
!     03/25/10           characteristic vibrational temperature, so if there is
!                        any significant CO2 fraction, enthalpy and temperature
!                        calculations are only approximate.
!                        Gary's package has also been updated, affecting
!                        equilibrium_composition.f90 and equilibrium.c.
!     09/08/11      "    Incorporated Tahir's more thorough handling of the
!                        CO2 and NCO triatomics with three vibrational modes.
!                        Note that lewis.c does not handle NCO yet, though.
!     09/12/11      "    Redirecting standard output (unit 6) was found to
!                        interact badly with singular matrix diagnostics from
!                        the lewis function.  Therefore, open a log file on
!                        unit 7 and write to that.  Any C code diagnostics
!                        still come to the screen.
!
!                        [Later:  Separating the two print streams is too
!                        inelegant.  Therefore, we revert to putting the log
!                        file name on the command line.]
!
!                        Continuing if the lewis (CEA) calculations fail was
!                        originally intentional (in case a poor Newton step
!                        can be halved enough to achieve success).  But now,
!                        a non-zero return from lewis/equilibrium_composition
!                        is treated as fatal.
!                        Also, Gary pointed out that Ar+ should have been among
!                        the subset of species here all along (at least for high
!                        temperature cases).
!
!  Authors:
!  ========
!
!     Analysis:          Tahir Gokcen,   ELORET/NASA Ames Research Center, CA
!     Implementation:    David Saunders, ELORET/NASA Ames Research Center, CA
!                                   now: ERC, Inc./NASA ARC
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Local constants:

   integer, parameter ::    &
      itmax       = 30,     &  ! outer iteration limit
      lunctl      = 5,      &  ! control file (standard input/command line)
      lunlog      = 6,      &  ! 6 here can cause standard output f90/C mix-up
      lunin       = 1,      &  ! for possible datasets specifying V1 and/or V2
      lunoutf     = 2,      &  ! for throat-condition results
      lunoutg     = 3,      &  ! for throat-condition uniform grid if not input
      lunplot     = 4,      &  ! for throat-condition uniform grid if not input
      lunprint    = lunlog, &  ! lunprint = -lunlog suppresses Newton itn. print
      max_fun     = 40,     &  ! limit on # fn. evals. for 1-D zero finder
      max_options = 8,      &  ! number of specified state cases handled
      max_species = 21,     &  ! length of species names dictionary
      nv          = 2          ! # nonlinear equations (that might be) solved

   real, parameter ::       &
      half   = 0.5,         &
      one    = 1.0,         &
      pt90   = 0.90,        &  ! For a tighter line search
      third  = 1.0 / 3.0,   &
      tol    = 1.E-5,       &  ! Rel. tolerance (1-D zero finder & Newton itn.)
      two    = 2.0,         &
      zero   = 0.0,         &
      R_hat  = 8.31441e+03     ! Universal gas constant, J/(kmol.K)

   logical, parameter ::    &
      false = .false.,      &
      true  = .true.

   character, parameter ::  &
      dataset * 7 = 'DATASET', &  ! Prescribed option for V1 or V2
      innerid * 5 = 'INNER',   &  ! Identifier for inner iteration of ZERORC
      outerid * 5 = 'OUTER'       !      "      "  outer     "     "    "

!  Local variables:

   integer :: &
      i, ibc, ier, ios, iter, ip(nv), j, k, l, n, ndim, neval, nf, nfplot,     &
      ngrid, numfun1, numfun2, n_species, n_v1_shape, n_v2_shape,              &
      status1, status2, status_CEA

   real :: &
      ae, c_over_M, Cp_bar, Cv_bar, e, ev_bar, gamf, h, ho, h_ref_sum,         &
      P, pi, r, rho, R_bar, s, T, t1,                                          &
      alpha, cube_root_eps, dr, fMsq, fnorm, fnorm0, frozen_Mach, toler,       &
      bulk_enthalpy, target_bulk_enthalpy, bulk_mass, target_bulk_mass,        &
      fscale1, fx1, tol1, x1, xa1, xb1, xscale1,                               &
      fscale2, fx2, tol2, x2, xa2, xb2, xscale2,                               &
      r_center, r_wall, square_root_eps, u,                                    &
      v1, v1_center, v1_wall, v1_wall_0, v1_width, v1_steepness,               &
      v2, v2_center, v2_wall, v2_wall_0, v2_width, v2_steepness

   real :: &
      hold1(13), hold2(13), &   ! Used by 1-D zero finder
      mol_weights(max_species), ref_enthalpies(max_species),                   &
      species_multipliers(max_species), theta_vibs(max_species),               &
      theta_CO2(3), theta_NCO(3)

   real :: &
      AJ(nv,nv), aitch(nv), dx(nv), f(nv), ffd(-1:1,nv,nv), x(nv), xfd(nv),    &
      xlast(nv)

   real, allocatable, dimension (:) :: &
      Cv_hat, ev_hat, hf, M, mass_fraction, mass_fraction_0, mole_fractions,   &
      reval, rgrid, theta_v, v1eval, v2eval

   real, pointer, dimension (:) :: &  ! Passed as unallocated arguments
      v1_shape_r, v1_shape_v, &
      v2_shape_r, v2_shape_v

   real, allocatable, dimension (:,:) :: &
      feval, fgrid, fplot

   logical :: &
      analytic, bulk_enthalpy_specified, bulk_mass_specified, diagnostics,     &
      interpolate_to_grid, reversed, show_details

   character :: &
      filename * 80, specified_state * 5, &
      v1_shape * 9, v2_shape * 9

   character :: &
      main_options(max_options) * 5, species_names(max_species) * 4

   character, allocatable, dimension (:) :: &
      name_input * 4, name_gary * 4

!  Procedures:

!! external &
!!    lusolve,                  & ! Square system solver (LU-factorization)
!!    equilibrium_composition,  & ! Gas properties routine (calls C package)
!!    lookup2,                  & ! Table look-up utility
!!    upcase                      ! Upper case utility

!  Storage:

!  The main options include the 6 handled by Gary's C package.
!  Adding new ones in front of those would be sensible.

   data main_options  /  &
      'Ht_MF', 'Ht_Ru', 'Rho_T', 'Rho_H', 'Rho_S', 'P_T  ', 'P_H  ', 'P_S  ' /

!  The following tables should be consistent with each other, but the species
!  may be in any order, and additions can be made at the ends.
!  Reference: DPLR's /cfdinput/chemprops.spec file.

   data species_multipliers /  &
      2.5,        2.5,        2.5,        1.5,        1.5,        1.5,        &
      2.5,        2.5,        2.5,        1.5,        1.5,        1.5,        &
      3.5,        2.5,        1.5,        2.5,        2.5,        1.5,        &
      2.5,        3.5,        1.5/

   data species_names       /  &
      'N2  ',     'O2  ',     'NO  ',     'N   ',     'O   ',     'Ar  ',     &
      'N2+ ',     'O2+ ',     'NO+ ',     'N+  ',     'O+  ',     'e   ',     &
      'CO2 ',     'CO  ',     'C   ',     'C2  ',     'CN  ',     'C+  ',     &
      'CO+ ',     'NCO ',     'Ar+ '/

   data mol_weights         /  &
      28.016,     32.0,       30.008,     14.008,     16.0,       39.944,     &
      28.01545,   31.99945,   30.00745,   14.00745,   15.99945,   0.00055,    &
      44.01100,   28.01100,   12.01100,   24.02200,   26.01900,   12.01045,   &
      28.01045,   42.01900,   39.94345/

   data ref_enthalpies      /  &
      0.,         0.,         2.996123e6, 3.362161e7, 1.542000e7, 0.,         &
      5.370000e7, 3.637000e7, 3.283480e7, 1.340000e8, 9.756000e7, 0.,         &
     -8.932880e6,-4.0630824e6,5.9211889e7,3.4234785e7,1.679542e7, 1.4967366e8,&
      4.4200904e7,4.212400e6, 3.806812e7/

   data theta_vibs          /  &
      3.395000e3, 2.239000e3, 2.817000e3, 0.,         0.,         0.,         &
      3.175800e3, 2.741200e3, 3.421000e3, 0.,         0.,         0.,         &
      -1.0,       3.122000e3, 0.,         2.66870e3,  2.97610e3,  0.,         &
      3.188000e3, -2.0,       0./

!     The negative theta_vibs for CO2 and NCO are used to indicate that there is
!     more than one theta for those triatomics, as follows:

   data theta_CO2 / &
      1.91870e+03, 9.59660e+02, 3.38210e+03/

   data theta_NCO / &
      1.83600e+03, 7.67100e+02, 2.76800e+03/
      

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                   Execution:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   pi = asin (one) * two

!! No - stay with redirecting standard output to any file name on command line.
!! open (lunlog, file='nozzle_throat_conditions.log', status='unknown')

   call control_inputs ()  ! Internal procedure below


!  Set up work-space associated with the indicated species:
!  --------------------------------------------------------

   n = n_species
   allocate (Cv_hat(n), ev_hat(n), hf(n), M(n), mole_fractions(n), theta_v(n))

   do i = 1, n_species

      call lookup2 (max_species, species_names, name_input(i), n)

      if (n <= 0) then
         write (lunlog, '(/, 2a)') ' Unhandled species: ', name_input(i)
         go to 99
      end if

      Cv_hat(i)  = R_hat * species_multipliers(n)
      M(i)       = mol_weights(n)
      hf(i)      = ref_enthalpies(n)
      theta_v(i) = theta_vibs(n)

   end do


!  Set up a uniform grid reval(:) for the initial flow evaluations:
!  ----------------------------------------------------------------

   if (v1_shape(1:1) == 'U' .and. v2_shape(1:1) == 'U') neval = 2 !Both uniform

   allocate (reval(neval), v1eval(neval), v2eval(neval))

   dr = (r_wall - r_center) / real (neval - 1)
   do i = 1, neval
      reval(i) = r_center + real (i - 1) * dr
   end do
   reval(neval) = r_wall


!  Discretize the indicated shapes for V1 and V2 on this uniform grid:
!  -------------------------------------------------------------------

   call discretize_shape (neval, reval, v1_shape, v1_center, v1_wall, v1_width,&
                         v1_steepness, v1eval, analytic)
   if (ios /= 0) go to 99  ! Bad shape spec.

   if (.not. analytic) call redistribute (n_v1_shape, v1_shape_r, v1_shape_v,  &
                                          neval, reval, v1eval)

   call discretize_shape (neval, reval, v2_shape, v2_center, v2_wall, v2_width,&
                         v2_steepness, v2eval, analytic)
   if (ios /= 0) go to 99  ! Bad shape spec.

   if (.not. analytic) call redistribute (n_v2_shape, v2_shape_r, v2_shape_v,  &
                                          neval, reval, v2eval)

   ier = 0

!  Option to seek specified bulk flow rate:
!  ----------------------------------------

   if (bulk_mass_specified) call bulk_iteration (1)

   if (ier < 0) go to 99


!  Option to seek specified bulk enthalpy:
!  ---------------------------------------

   if (bulk_enthalpy_specified) call bulk_iteration (2)

   if (ier < 0) go to 99


!  Calculate remaining flow conditions from the nozzle center to the wall:
!  -----------------------------------------------------------------------

!  Store all possible plottable data, but some of it may be suppressed.

   nfplot = n_species + 8

   allocate (fplot(neval,nfplot), feval(neval,nf))

   if (ndim > 0) then

      v1 = v1eval(1)
      v2 = v2eval(1)

      call starting_guess (0)  ! For the first point only

   end if

   do k = 1, neval

      v1 = v1eval(k)
      v2 = v2eval(k)

      call mass_mole_conversion (1)  ! Mass fractions to mole fractions

      select case (ndim)

      case (0) ! Straight call to Gary's package

         if (k == 1 .or. k == neval) then  ! Force printing the details
            status_CEA = 1
         else
            status_CEA = 0
         end if

         call control_outputs ()

         call equilibrium_composition (specified_state, v1, v2,                &
                                       n_species, name_gary, mole_fractions,   &
                                       P, rho, T, h, s, ae, status_CEA)
         if (status_CEA /= 0) then
            write (lunlog, '(a, i4, a, 1p, 2e13.6)') &
               ' Bad equilibrium gas composition status:', status_CEA,         &
               '  v1, v2:', v1, v2
            write (lunlog, '(3a, i4)') &
               ' Specified state case: ', specified_state, '    Point #:', k
            go to 99
         end if
     
         call mass_mole_conversion (2)       ! Mole fractions to mass fractions

      case (1) ! 2 separable eqns. solved via nested non-deriv. zero finders

         status2 = 2              ! Initialize the outer iteration
         numfun2 = max_fun        ! Max. no. of fn. evals. allowed

         call starting_guess (2)  ! Set [xa2, xb2] containing root x2 of eq. (2)

 20      continue  ! Outer loop, updating P from equation (2)

            call zerorc (xa2, xb2, x2, fx2, tol2, numfun2, outerid, lunlog, &
                         hold2, status2)

            status_CEA = 0        ! Suppress printout unless requested
            call control_outputs ()

            if (diagnostics) write (lunlog, '(/, a, i4)') &
               ' ZERORC status for eq. (2): ', status2

            if (status2 < -1) then      ! Fatal error

               write (lunlog, '(/, a, i4)') &
                  ' Bad return from ZERORC for eq. (2). Point k: ', k
               go to 99

            else if (status2 < 0) then  ! Iteration limit reached

               write (lunlog, '(/, a, i4, a)') &
                  ' Max. # fn. evals. for eq. (2), pt. k =', k, '. Proceeding.'

            else if (status2 > 0) then  ! Evaluate f(2) by [re]solving eq. (1)

               status1 = 2              ! Initialize the inner iteration
               numfun1 = max_fun        ! Max. no. of fn. evals. allowed

               call starting_guess (1)  ! Set [xa1, xb1] for root x1 of eq. (1)

 21            continue

               call zerorc (xa1, xb1, x1, fx1, tol1, numfun1, innerid, lunlog, &
                            hold1, status1)

               if (diagnostics) write (lunlog, '(/, a, i4)') &
                  ' ZERORC status for eq. (1): ', status1

               if (status1 < -1) then      ! Fatal error

                  write (lunlog, '(/, a, i4)') &
                     ' Bad return from ZERORC for eq. (1). Point k: ', k
                  go to 99

               else if (status1 < 0) then  ! Iteration limit reached

                  write (lunlog, '(/, a, i4, a)') &
                  ' Max. # fn. evals. for eq. (1), pt. k =', k, '. Proceeding.'

               else if (status1 > 0) then  ! Evaluate the function from eq. (1)

                  status_CEA = 0        ! Suppress printout unless requested
                  call control_outputs ()

                  call residual_eq_1 (x1, fx1)
                  go to 21
 
               ! else status1 = 0 - a zero has been found.
 
               end if

               call residual_eq_2 (x2, fx2)

               go to 20  ! Another outer iteration?

            ! else status2 = 0 - the outer iteration has converged

            end if

      case (2) ! 2 nonlinear equations; Newton iteration with central diffs.

!        The iteration involves solving  J dx = f, where J is the matrix
!        of partial derivatives df(i)/dx(j), updating x as x - dx, and
!        repeating until || f || becomes arbitrarily small.

!        Safeguarding involves halving each pure Newton step until || f ||
!        is reduced from its value for the previous iteration, meaning that
!        the update step is actually x <-- x - alpha dx, where alpha <= 1.

!        Tests for convergence are subject to the choice of tolerance, which
!        is dependent upon the data scaling and precision.  The tolerance
!        is used relative to the norm of the starting guess in SI units.

!        N.B. Don't confuse the "v"s here (P and T) with the specified v1 & v2.

         ier      = 0
         iter     = 0
!!       x(1)     = P              ! Already set via call starting_guess (0)
!!       x(2)     = T   
         dx(:)    = zero           ! For iteration 0 printout
         alpha    = zero           !  "      "     "     "
         fnorm0   = 1.E+30         ! Big to avoid iteration 0 test
         fscale1  = v1             ! Since eqn. 1 has units of enthalpy
         fscale2  = v2
         xscale1  = x(1);  x(1) = one
         xscale2  = x(2);  x(2) = one
         aitch(:) = (one + x(:)) * cube_root_eps
!!!      toler  = tol * max (abs (x(1)), abs (x(2))) ! Should be relative to f,
!!!                              !               not x, but this should suffice
!!!      toler  = max (toler, one)  ! It doesn't!  Use 1.0 since equation 1 has
!!!                                 ! units of total enthalpy (very large)
!!!      toler    = square_root_eps ! Allow for limited precision from Gary
         toler    = 1.e-6           ! Sqrt. (eps) is too small
         reversed = .false.         ! Probably never needed

         status_CEA = 0             ! Suppress printout unless requested
         call control_outputs ()

         if (diagnostics) then
            write (lunlog, '(/, a, 1p, 2e15.6)') &
               ' Initial variables:                      ', x(:)
            write (lunlog, '(a, 1p, 2e15.6)') &
               ' Central difference steps:               ', aitch
            write (lunlog, '(a, 1p, e15.6)') &
               ' Convergence tolerance on residual norm: ', toler
         end if

  30     continue

!           Evaluate the residuals for the current variables.
!           Actually, to calculate central derivatives while we're at it,
!           we evaluate the functions for each x at x - h, x, x + h.

            mass_fraction_0(:) = mass_fraction(:)  ! Start consistently

            do j = 1, nv

               xfd(:) = x(:)

               do l = -1, 1

                  xfd(j) = x(j) + real (l) * aitch(j)
                  mass_fraction(:) = mass_fraction_0(:)

                  status_CEA = 0        ! Suppress printout unless requested
                  call control_outputs ()

                  call residual_2d (xfd(1)*xscale1, xfd(2)*xscale2, &
                                    ffd(l,1,j), ffd(l,2,j))

                  ffd(l,1,j) = ffd(l,1,j) / fscale1
                  ffd(l,2,j) = ffd(l,2,j) / fscale2
               end do

!              Partial derivatives dfi / dxj:

               do i = 1, nv
                  AJ(i,j) = (ffd(1,i,j) - ffd(-1,i,j)) / (two * aitch(j))
               end do

            end do

!           Set up the residual vector f of function values as the RHS of
!           the system J dx = f.  The norm of f should converge to ~0.

            f(1) = ffd(0,1,1)
            f(2) = ffd(0,2,2)

            fnorm = max (abs (f(1)), abs (f(2)))
            mass_fraction(:) = mass_fraction_0(:)

  40        continue

            status_CEA = 0        ! Suppress printout unless requested
            call control_outputs ()

!           Halve the step until || f || is reduced (except first time through).

            if (fnorm > fnorm0) then
               if (alpha > 0.001) then
!!!               alpha = half * alpha
                  alpha = pt90 * alpha
                  do i = 1, nv
                     x(i) = xlast(i) - alpha * dx(i)
                  end do
                  mass_fraction(:) = mass_fraction_0(:)

                  call residual_2d (x(1)*xscale1, x(2)*xscale2, f(1), f(2))

                  f(1) = f(1) / fscale1
                  f(2) = f(2) / fscale2
                  fnorm = max (abs (f(1)), abs (f(2)))
                  if (diagnostics) then
                     write (lunprint, &
                        '(a, i3, a, 1p, e9.2, a, e9.2, 3(a, 2e15.7))') &
                        ' Itn.', iter, ' ||f||:', fnorm, ' step:', alpha, &
                        ' dx:', alpha*dx(1), alpha*dx(2), ' x:', x(:), &
                        ' f:', f(:)
                  end if
                  go to 40
               else if (.not. reversed) then ! Direction must have been uphill
                  reversed = .true.
                  alpha = one
                  dx(:) = -dx(:)
                  x(:)  = xlast(:) - dx(:)

                  call residual_2d (x(1)*xscale1, x(2)*xscale2, f(1), f(2))

                  f(1) = f(1) / fscale1
                  f(2) = f(2) / fscale2
                  fnorm = max (abs (f(1)), abs (f(2)))
                  if (diagnostics) then
                     write (lunprint, &
                        '(a, i3, a, 1p, e9.2, a, e9.2, 3(a, 2e15.7))') &
                        ' Itr.', iter, ' ||f||:', fnorm, ' step:', alpha, &
                        ' dx:', alpha*dx(1), alpha*dx(2), ' x:', x(:), &
                        ' f:', f(:)
                  end if
                  go to 40
               end if
               go to 83
            end if

            if (diagnostics) then
               dx(:) = alpha * dx(:)
               write (lunprint, &
                  '(a, i3, a, 1p, e9.2, a, e9.2, 2(a, 2e16.8))') &
                  ' Itn.', iter, '  ||f||:', fnorm, '  step:', alpha, &
                  '  dx:', dx(:), '  x:', x(:)
            end if

            if (fnorm < toler .and. iter > 0) go to 90  ! Success


            iter = iter + 1
            if (iter == itmax) go to 81

            fnorm0 = fnorm

!           Set up the LHS matrix for the next iteration.  AJ(i,j) = dfi/dxj.

!!!         AJ(1, 1) = ...  ! Already done above during the finite differencing
!!!         AJ(2, 1) = ...
!!!         AJ(1, 2) = ...
!!!         AJ(2, 2) = ...

!           Solve  J dx = f.

            call lusolve (nv, nv, AJ, f, ier)      ! Solution overwrites RHS

            if (ier /= 0) go to 84

            do i = 1, nv
               dx(i)    = f(i)
               xlast(i) = x(i)
               x(i)     = x(i) - dx(i)
            end do

            alpha = one
            reversed = false

         go to 30                           ! Do another iteration


!        Error handling.

  81     ier = 1
         write (lunlog, '(/, a)') ' Iteration limit reached.'
         if (fnorm < 10. * toler) go to 90

         ier = 2
!!!      go to 99   ! Bad data scaling can prevent the convergence test from
         go to 90   ! being satisfied even when the solution is still adequate.

  83     ier = 3
         write (lunlog, '(/, a)') ' Step-halving iteration failed.'
         go to 90

  84     ier = 4
         write (lunlog, '(/, a)') ' Singular matrix.'  ! Show solution anyway
!!!      go to 99

  90     continue    ! Wrap up, having converged or almost converged.

!        Ensure final values:

         mass_fraction(:) = mass_fraction_0(:)

         if (k == 1 .or. k == neval) then  ! Force printing the details
            status_CEA = 1
         else
            status_CEA = 0
         end if
         call control_outputs ()

         call residual_2d (x(1)*xscale1, x(2)*xscale2, f(1), f(2))

         f(1) = f(1) / fscale1
         f(2) = f(2) / fscale2

         if (diagnostics) write (lunlog, &
           '(a, 1p, 2e14.6, a, 2e14.6, a, 2e14.6)') &
           ' v1, v2:', v1, v2, '   x1, x2:', x(:), '   Residuals:', f(:)

         x(1) = x(1) * xscale1  ! For starting the next case
         x(2) = x(2) * xscale2

      end select  ! Over main options for specified 2 variables and soln. method

      mass_fraction_0(:) = mass_fraction(:)

!     Store the flow variables indicated by ibc for this point:

      call store_results (ibc, neval, nf, k, feval)   ! In feval(k,1:nf)

   end do ! Next location


!  Save the results on the uniform grid or redistribute them first if indicated:
!  -----------------------------------------------------------------------------

   if (interpolate_to_grid) then

      do i = 1, nf

         call redistribute (neval, reval, feval(1,i), ngrid, rgrid, fgrid(1,i))

      end do

      call save_results (lunoutf, lunplot, ngrid, nf, fgrid)

   else

      open  (lunoutg, file='uniform.g', status='unknown')
      write (lunoutg, '(2i3)') neval, 1
      write (lunoutg, '(1p, 6e18.10)') reval
      write (lunoutg, '(6f3.0)') (0., i = 1, neval)
      close (lunoutg)

      call save_results (lunoutf, lunplot, neval, nf, feval)

   end if

!  Done.


99 continue  ! Avoid possible system-dependent STOP behavior.

!  Internal procedures for program nozzle_throat_conditions:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine bulk_iteration (item)

!     Solve the 1-variable nonlinear equation that specifying the indicated bulk
!     flow quantity leads to.  The bulk enthalpy is dependent upon the bulk mass
!     flow (but not vice versa), so it should be determined second.
!     The one variable is a multiplier for the input wall level of the profile
!     for mass flux (item = 1) or total enthalpy (iterm = 2).  The center levels
!     of these quantities are fixed at the input values.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Argument:

      integer, intent (in) :: item  ! 1 = bulk flow rate, kg/s;
                                    ! 2 = bulk enthalpy, J/kg
!     Local constants:

      integer,   parameter :: maxfun   = 40
      real,      parameter :: tol      = 1.e-7
      character, parameter :: name * 8 = 'BULKITER'

!     Local variables:

      integer :: istat, lunout, numfun
      real    :: hold(13), toler, x, xa, xb

!     Execution:

      xa     = 0.5 ! Multiplier range to search
      xb     = 2.0
      toler  = tol * xa
      numfun = maxfun
      lunout = lunlog;  if (.not. diagnostics) lunout = -lunlog
      istat  = 2        ! Initialize the iteration

      do while (istat > 0)

         call zerorc (xa, xb, x, f, toler, numfun, name, lunout, hold, istat)

         if (istat > 0) then  ! Evaluate the function at this x

            select case (item)

               case (1) ! Mass flux profile is being adjusted

                  v2_wall = v2_wall_0 * x

               case (2) ! Total enthalpy profile is being adjusted

                  v1_wall = v1_wall_0 * x

            end select

!           Discretize the profiles for V1 = total enthalpy and V2 = mass flux:

            call discretize_shape (neval, reval, v1_shape, v1_center, v1_wall, &
                                   v1_width, v1_steepness, v1eval, analytic)
            if (ios /= 0) go to 99  ! Bad shape spec.

            if (.not. analytic) &
               call redistribute (n_v1_shape, v1_shape_r, v1_shape_v,          &
                                  neval, reval, v1eval)

            call discretize_shape (neval, reval, v2_shape, v2_center, v2_wall, &
                                   v2_width, v2_steepness, v2eval, analytic)
            if (ios /= 0) go to 99  ! Bad shape spec.

            if (.not. analytic) &
               call redistribute (n_v2_shape, v2_shape_r, v2_shape_v,          & 
                                  neval, reval, v2eval)

!           Integrate the bulk quantities (together, because ho depends on MF):

            call integrate_bulk_mass_and_enthalpy &
                        (neval, reval, v1eval, v2eval, bulk_mass, bulk_enthalpy)

            select case (item)

               case (1)

                  f = bulk_mass - target_bulk_mass

                  if (numfun == 0) write (lunlog, '(/, a, 1p, e15.7, a)')      &
                     ' Minimum bulk mass rate for this starting guess:    ',   &
                     bulk_mass, ' kg/s'
                  if (numfun == 1) write (lunlog, '(/, a, 1p, e15.7, a)')      &
                     ' Maximum bulk mass rate for this starting guess:    ',   &
                     bulk_mass, ' kg/s'

               case (2)

                  f = bulk_enthalpy - target_bulk_enthalpy

                  if (numfun == 0) write (lunlog, '(/, a, 1p, e15.7, a)')      &
                     ' Minimum bulk enthalpy for this starting guess:',        &
                     bulk_enthalpy, ' J/kg'
                  if (numfun == 1) write (lunlog, '(/, a, 1p, e15.7, a)')      &
                     ' Maximum bulk enthalpy for this starting guess:',        &
                     bulk_enthalpy, ' J/kg'

            end select

         end if

      end do

!     Wrap up:

      ier = istat

      if (istat >= -1) then ! 0 = converged); -1 = iteration limit - may be OK

!        Ensure the best x is active:

         select case (item)

            case (1) ! Mass flux profile is being adjusted

               v2_wall = v2_wall_0 * x

               write (lunlog, '(/, a, 1p, e15.7, a)')                          &
                  ' Adjusted mass flux profile wall value:     ', v2_wall,     &
                  ' kg/(m^2.s)'

            case (2) ! Total enthalpy profile is being adjusted

               v1_wall = v1_wall_0 * x

               write (lunlog, '(/, a, 1p, e15.7, a)')                          &
                 ' Adjusted total enthalpy profile wall value:', v1_wall,' J/kg'

         end select

         ier = 0

      else ! istat < -1) then      ! Fatal error

         write (lunlog, '(/, a, i3)') &
            ' BULK_ITERATION:  Bad return from ZERORC.  istat:', istat

         write (lunlog, '(2a)') ' Search is limited to [0.5, 2.0] for the',    &
            ' multiple of the initial profile at the wall.'

      end if

 99   return

      end subroutine bulk_iteration

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine control_inputs ()

!     Read the control file.  All variables are inherited from the host.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Echo control file (standard input) to the output log (standard output):

      do ! Until EOF
         read (lunctl, '(a)', iostat=ios) filename ! Handy buffer
         if (ios < 0) exit
         write (lunlog, '(1x, a)') trim (filename)
      end do

      rewind (lunctl)

!!    cube_root_eps   = epsilon (cube_root_eps) ** third ! For central diffs.
      cube_root_eps   =                   1.e-7 ** third ! In view of Gary pkg.
      square_root_eps = sqrt (epsilon (toler))

!     Read the control file:

      do i = 1, 5
         read (lunctl, *)
      end do

      read (lunctl, *) specified_state  ! Option descriptor

      call lookup2 (max_options, main_options, specified_state, n)

      if (n > 2) then  ! Just call Gary's package at each point r
         ndim = 0
      else if (n == 2) then
         ndim = 1      ! 2 separable equations (nested 1-variable iterations)
      else if (n == 1) then
         ndim = 2      ! 2 nonlinear equations (2-variable Newton iteration)
      else
         write (lunlog, '(/, 2a)') ' Unhandled main option: ', specified_state
         go to 99
      end if

      read (lunctl, *) frozen_Mach  ! Not necessarily 1

      fMsq = frozen_Mach**2

!     Iteration controls for bulk flow specs. (Ht_MF case only):

      do i = 1, 3
         read (lunctl, *)
      end do

      read (lunctl, *) i
      bulk_enthalpy_specified = i /= 0
      read (lunctl, *) target_bulk_enthalpy

      read (lunctl, *) i
      bulk_mass_specified = i /= 0
      read (lunctl, *) target_bulk_mass  ! N.B.: This should be performed first

!     Variable 1 specs:

      read (lunctl, *);  read (lunctl, *);  read (lunctl, *)
      read (lunctl, *) v1_shape
      read (lunctl, *) r_center, v1_center
      read (lunctl, *) r_wall,   v1_wall;   v1_wall_0 = v1_wall
      read (lunctl, *) v1_width
      read (lunctl, *) filename  ! or "none" or steepness control

      call upcase (v1_shape)

      v1_steepness = 1.
      if (trim (v1_shape) == dataset) then ! See internal procedure below

         call read_discretized_shape (n_v1_shape, v1_shape_r, v1_shape_v)
         if (ios /= 0) go to 99

      else ! Option for Sigmoid shape function that avoids another control

         if (ichar (filename(1:1)) < ichar ('A')) then  ! Assume it's numeric
            l = len_trim (filename)
            read (filename(1:l), *) v1_steepness
         end if

      end if

!     Variable 2 specs:

      read (lunctl, *);  read (lunctl, *);  read (lunctl, *)
      read (lunctl, *) v2_shape
      read (lunctl, *) r_center, v2_center
      read (lunctl, *) r_wall,   v2_wall;   v2_wall_0 = v2_wall
      read (lunctl, *) v2_width
      read (lunctl, *) filename

      call upcase (v2_shape)

      v2_steepness = 1.

      if (trim (v2_shape) == dataset) then ! See internal procedure below

         call read_discretized_shape (n_v2_shape, v2_shape_r, v2_shape_v)
         if (ios /= 0) go to 99

      else ! Option for Sigmoid shape function

         if (ichar (filename(1:1)) < ichar ('A')) then  ! Assume it's numeric
            l = len_trim (filename)
            read (filename(1:l), *) v2_steepness
         end if

      end if

!     Starting guesses for the species mass fractions:

      read (lunctl, *);  read (lunctl, *);  read (lunctl, *)
      read (lunctl, *) n_species

      nf = 4 + n_species  ! For 2-D and two temperatures (though Tv = T)

      n = n_species
      allocate (name_input(n), mass_fraction(n), mass_fraction_0(n))

      do i = 1, n_species
         read (lunctl, *, iostat=ios) name_input(i), mass_fraction(i)
         if (ios /= 0) then
            write (lunlog, '(/, a, i4)') &
               ' Trouble reading species name/mass_fraction; i =', i
            go to 99
         end if
      end do

      allocate (name_gary(n_species))

      l = len (name_input(1))

!     Convert from DPLR to C code convention (works for the likely species):

      do i = 1, n_species
         name_gary(i) = name_input(i)
         if (name_gary(i)(1:1) == 'e') name_gary(i)(1:2) = 'El'
         do j = 1, l
            if (name_gary(i)(j:j) == '+') name_gary(i)(j:j) = 'p'
            if (name_gary(i)(j:j) == '-') name_gary(i)(j:j) = 'n'
         end do
      end do

      read (lunctl, *);  read (lunctl, *);  read (lunctl, *)
      read (lunctl, *) neval  ! # uniform radii for initial calculations

      show_details = neval < 0    ! True means all possible printout;
      diagnostics  = show_details ! false means suppress most but not all;
      neval        = abs (neval)  ! the original "diagnostics" control variable
                                  ! is toggled accordingly.

!     Target radii for output results (1-D grid), if any:

      read (lunctl, *) filename

      interpolate_to_grid = filename(1:4) /= 'none'

      if (interpolate_to_grid) then

         l = len_trim (filename)
         open (lunin, file=filename(1:l), status='old', iostat=ios)
         if (ios /= 0) then
           write (lunlog, '(/, 2a)') ' Unable to open grid file ', filename(1:l)
           go to 99
         end if

!        A 2nd coordinate may be present or not.  I.e., 1-D and 2-D are both OK.

         read (lunin, *, iostat=ios) ngrid
         if (ios /= 0) then
            write (lunlog, '(/, 2a)') &
               ' Trouble reading # pts. in grid file ', filename(1:l)
            go to 99
         end if

         allocate (rgrid(ngrid), fgrid(ngrid,nf))

         read (lunin, *, iostat=ios) rgrid
         if (ios /= 0) then
            write (lunlog, '(/, 2a)') ' Trouble reading grid ', filename(1:l)
            go to 99
         end if
         close (lunin)

      end if

      read (lunctl, *) ibc       ! DPLR input profile BC code
      read (lunctl, *) filename  ! Output flow results (PLOT3D function file)

      open (lunoutf, file=filename, status='unknown')

      read (lunctl, *) filename  ! Plottable results (Tecplot ASCII)

      open (lunplot, file=filename, status='unknown')

      close (lunctl)

 99   return

      end subroutine control_inputs

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine control_outputs ()

!     Modularize the turning on and off of diagnostic printout from CEA.
!     If most diagnostics are suppressed, we still display the equilibrium gas
!     composition details for the solutions at the first and last spacial pts.,
!     but not during any finite differencing or Newton iteration evaluations.
!
!     Recall that full details are invoked by entering a negative number of
!     working grid points in the main control file, while the status argument
!     of equilibrium_composition is made use of as its printing control.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      diagnostics = show_details
      if (.not. show_details) diagnostics = status_CEA /= 0

      end subroutine control_outputs

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine discretize_shape (n, r, v_shape, v_center, v_wall, v_width,   &
                                   v_steepness, v, analytic)

!     Evaluate the indicated variable shape at pts. r(1:n) from center to wall.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,   intent (in)  :: n            ! # evaluations >= 2
      real,      intent (in)  :: r(n)         ! Abscissas for evaluations
      character, intent (in)  :: v_shape*(9)  ! Shape fn. descriptor, upper case
      real,      intent (in)  :: v_center     ! Variable values at the center ..
      real,      intent (in)  :: v_wall       ! .. and the wall
      real,      intent (in)  :: v_width      ! For sinusoid, nthdegree, sigmoid
      real,      intent (in)  :: v_steepness  ! For sigmoid
      real,      intent (out) :: v(n)         ! Discretized variable values
      logical,   intent (out) :: analytic     ! Signals whether to redistribute
                                              ! a discretized definition of V
!     Local variables:

      integer :: i
      real    :: a, b, rm, rscale, rshift, t

!     Execution:

      analytic = .true.

      select case (v_shape)

         case ('UNIFORM  ')

            v(:) = v_wall  ! Because this is what the bulk iteration adjusts

         case ('LINEAR   ')

            t = (v_wall - v_center) / (r(n) - r(1))

            do i = 1, n - 1
               v(i) = v_center + t * (r(i) - r(1))
            end do

         v(n) = v_wall  ! Just to be precise

         case ('PARABOLIC')

            t = (v_wall - v_center) / (r(n) - r(1))**2

            do i = 1, n - 1
               v(i) = v_center + t * (r(i) - r(1))**2
            end do

            v(n) = v_wall

         case ('MTHDEGREE')  ! Peak at mid-semi-axis for semi-ellipse case

            rm = (r(1) + r(n)) * half
            a  =  r(n) - rm
            b  =  v_center - v_wall

            v(1) = v_wall

            do i = 2, n - 1
               t = abs (r(i) - rm) / a
               v(i) = v_center - b * t**v_width
            end do

            v(n) = v_wall

         case ('NTHDEGREE')  ! Peak at center (r = 0)

            a = r(n) - r(1)
            b = v_center - v_wall

            do i = 1, n - 1
               t = (r(i) - r(1)) / a
               v(i) = v_wall + b * (one - t**v_width)
            end do

            v(n) = v_wall

         case ('SINUSOID ')

            rscale = (half * pi) / (r(n) - r(1))
            rshift = -rscale * r(1)

            do i = 1, n - 1
               t = rscale * r(i) + rshift
               v(i) = v_wall + (v_center - v_wall) * (cos (t))**v_width
            end do

            v(n) = v_wall

         case ('SIGMOID  ')

            b = v_wall;  a = v_center - v_wall
            v(1) = v_center

            do i = 2, n - 1
               t = v_steepness * (r(i) - v_width)
               v(i) = b + a / (one + exp (t))
            end do

            v(n) = b

         case ('GAUSSIAN ')  ! f(r) = a * e (-half r^2 / sigma^2)  +  b

!           Solve a pair of equations in a, b for given center and wall values:

            t      = exp (-half * (r(n) / v_width)**2)
            b      = one - t
            rscale = v_width * sqrt (two * pi)
            a      = rscale * (v_center - v_wall) / b
            b      = (v_wall - v_center * t) / b
            rscale = a / rscale

            v(1) = v_center

            do i = 2, n - 1
               v(i) = rscale * exp (-half * (r(i) / v_width)**2)  +  b
            end do

            v(n) = v_wall

         case ('LORENTZ  ')  ! f(r) = a * (gamma / (gamma^2 + r^2))  +  b

!           Solve a pair of equations in a, b for given center and wall values:

            t    = one + (v_width / r(n))**2
            a    = v_width * t * (v_center - v_wall)
            b    = t * v_wall  + (one - t) * v_center
            t    = a * v_width
            v(1) = v_center

            do i = 2, n - 1
               v(i) = t / (v_width**2 + r(i)**2)  +  b
            end do

            v(n) = v_wall

         case ('DATASET  ')

            analytic = .false.  ! Redistribute via another routine

         case default

            write (*, '(/, 2a)') &
               ' Unhandled flow variable shape spec.: ', v_shape
            ios = -1  ! Abort

      end select

      end subroutine discretize_shape

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine integrate_bulk_mass_and_enthalpy (n, r, v1, v2,               &
                                                   bulk_mass, bulk_enthalpy)

!     Integrate bulk mass flow and bulk enthalpy for an axisymmetric nozzle.
!     The functions v1 and v2 are assumed to be total stagnation enthalpy and
!     mass flow rate per unit area, respectively.  The bulk enthalpy is the
!     integral of rho u ho over the integral of rho u (w.r.t. area).
!     Trapezoidal integration would suffice, but the local spline quadrature
!     method allows (slight) extrapolation, so supply it with data at the cell
!     mid-points.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)  :: n             ! # data pts. along a nozzle radius
      real,    intent (in)  :: r(n)          ! Coordinates
      real,    intent (in)  :: v1(n), v2(n)  ! ho and mass flux data
      real,    intent (out) :: bulk_mass, bulk_enthalpy

!     Local constants:

      character, parameter :: method * 1 = 'B' ! Bessel method (loose fit)

!     Local variables:

      integer :: i
      real, allocatable :: x(:), f(:)

!     Execution:

!     Integrate over the entire nozzle cross section area, mass flux first:

!     dA(i) = pi (r(i+1)**2 - r(i)**2) = pi (r(i) + r(i+1)) dr(i)

      allocate (x(n-1), f(n-1))

      do i = 1, n - 1
         x(i) = (r(i) + r(i+1)) * half
         f(i) = (v2(i) + v2(i+1)) * x(i)
      end do

      call lcsquad (n - 1, x, f, r(1), r(n), method, bulk_mass)

      bulk_mass = pi * bulk_mass

!     Bulk enthalpy = integral of rho u ho / bulk mass

      do i = 1, n - 1
         f(i) = (v1(i) + v1(i+1)) * (v2(i) + v2(i+1)) * x(i)
      end do

      call lcsquad (n - 1, x, f, r(1), r(n), method, bulk_enthalpy)

      bulk_enthalpy = pi * bulk_enthalpy / (two * bulk_mass)

      write (lunlog, '(/, (a, 1p, e15.7, a))') &
         ' Approximate bulk enthalpy:      ', bulk_enthalpy, ' J/kg', &
         ' Approximate bulk mass flow rate:', bulk_mass,     ' kg/s'

      deallocate (x, f)

      end subroutine integrate_bulk_mass_and_enthalpy

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mass_mole_conversion (mode)

!     Switch between the mass fractions preferred here and the mole fractions
!     treated by Gary's C package.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, intent (in) :: mode  ! 1 for mass -> mole;  2 for mole -> mass

      integer :: i
      real    :: Mbar, sum

      sum = zero

      select case (mode)

         case (1) ! Mass fractions to mole fractions

            do i = 1, n_species
               mole_fractions(i) = mass_fraction(i) / M(i)
               sum = sum + mole_fractions(i)
            end do

            Mbar = one / sum  ! Mixture molecular weight

            mole_fractions(:) = Mbar * mole_fractions(:)

         case (2) ! Mole fractions to mass fractions

            do i = 1, n_species
               mass_fraction(i) = mole_fractions(i) * M(i)
               sum = sum + mass_fraction(i)
            end do

            Mbar = sum

            mass_fraction(:) = mass_fraction(:) / Mbar

         end select

      end subroutine mass_mole_conversion

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_discretized_shape (npts, r, v)

!     Read a 1-D distribution defining the indicated flow variable from nozzle
!     center-line to nozzle wall.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, intent (out)        :: npts    ! # points found
      real, pointer, dimension (:) :: r, v    ! 1-D shape distribution

      integer :: i, l

      l = len_trim (filename)

      open (lunin, file=filename(1:l), status='old', iostat=ios)

      if (ios /= 0) then
         write (lunlog, '(/, 2a)') ' Unable to open file ', filename(1:l)
         go to 99
      end if

!     Count the number of defining points (top to EOF):

      npts = 0
      do while (ios == 0)
         npts = npts + 1
         read (lunin, *, iostat=ios)
      end do
      npts = npts - 1
      rewind (lunin)

      allocate (r(npts), v(npts))

      do i = 1, npts
         read (lunin, *, iostat=ios) r(i), v(i)
         if (ios /= 0) then
            write (lunlog, '(/, 3a, i4)') &
               ' Trouble reading file ', filename(1:l), '; line #:', i
            go to 99
         end if
      end do

      close (lunin)

 99   return

      end subroutine read_discretized_shape

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine redistribute (n, r, v, neval, reval, veval)

!     Redistribute a 1-D dataset to specified abscissas.
!     The new abscissas are assumed to be contained within the data range.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, intent (in)  :: n            ! # points to be redistributed
      real,    intent (in)  :: r(n), v(n)   ! Distribution being interpolated
      integer, intent (in)  :: neval        ! Desired # points
      real,    intent (in)  :: reval(neval) ! Redistributed abscissas
      real,    intent (out) :: veval(neval) ! Reinterpolated ordinates

      real, allocatable :: unused_derivs (:)

      allocate (unused_derivs(neval))

      call lcsfit (n, r, v, .true., 'M', neval, reval, veval, unused_derivs)

      deallocate (unused_derivs)
      
      end subroutine redistribute

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine residual_eq_1 (x, fx)

!     Evaluate the residual for equation (1) of the specified main option.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      real, intent (in)  :: x
      real, intent (out) :: fx

!     Local variables:

      integer :: i
      real    :: P_specified, T_specified

!     Execution:

      select case (specified_state)

      case ('Ht_MF', 'Ht_Ru')  ! "x" = T in eq. (1)

!        f(1) = v1 - t1 Cpbar T - evbar - sum (ci hfi) = 0   where v1 = ho
!        f(2) = v2 - P sqrt (gamf / (Rbar T))          = 0     "   v2 = rho u

         P_specified = P  ! The pressure associated with the current T = x
         T_specified = x

!        Compute species densities, etc., for the current estimates of P, T:

         call mass_mole_conversion (1)  ! Mass fractions to mole fractions

         if (diagnostics) then
            write (lunlog, '(/, a, 1p, 2e14.6)') &
               ' Residual_eq_1:  P, T inputs to CEA: ', P_specified, T_specified
            if (k == 1) then
               write (lunlog, '(a, 1p, e25.17)') &            ! For direct input
               (name_gary(i), mole_fractions(i), i = 1, n_species)  ! to Gary's
            end if                                            ! driving program
            status_CEA = 1  ! Force CEA printout
         else
            status_CEA = 0
         end if

         call equilibrium_composition ('P_T  ', P_specified, T_specified,      &
                                       n_species, name_gary, mole_fractions,   &
                                       P, rho, T, h, s, ae, status_CEA)
         if (diagnostics) &
            write (lunlog, '(a, 1p, 3e14.6)') ' Output rho, h, ae: ', rho, h, ae

         if (status_CEA /= 0) then
            write (lunlog, '(2(a, i4), /, a, 1p, 2e13.6)')                     &
               ' Bad equilibrium composition status:', status_CEA,             &
               '  Point #:', k, ' P and T specified: ', P_specified, T_specified
!!!         go to 99
            stop   ! Hard to avoid, but better than keeping going
         end if

         call mass_mole_conversion (2)  ! Mole fractions to mass fractions

         if (diagnostics) &
            write (lunlog, '(a, /, 1p, 12e11.4)') ' Output mass fractions:',   &
            mass_fraction(:)

!        For this equilibrium composition, update the gas constants:

         call frozen_gas_constants (T_specified)

         t1 = one + half * (gamf - one) * fMsq
         fx = v1 - t1 * Cp_bar * T_specified - ev_bar - h_ref_sum

         if (diagnostics) write (lunlog, '(a, 1p, 3e14.6)') &
            ' Eq. (1) P, T, residual:   ', P, T_specified, fx

      end select

!! 99 return
      return

      end subroutine residual_eq_1

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine residual_eq_2 (x, fx)

!     Evaluate the residual for equation (2) of the specified main option.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, intent (in)  :: x
      real, intent (out) :: fx

      select case (specified_state)

      case ('Ht_MF', 'Ht_Ru')  ! "x" = P in eq. (2)

!        f(1) = v1 - t1 Cpbar T - evbar - sum (ci hfi) = 0   where v1 = ho
!        f(2) = v2 - P sqrt (gamf / (Rbar T))          = 0     "   v2 = rho u

         P    = x
         fx   = v2 - frozen_Mach * P * sqrt (gamf / (R_bar * T))

         if (diagnostics) write (lunlog, '(a, 1p, 3e14.6)') &
            ' Eq. (2) P, T, residual:   ', P, T, fx

      end select

      end subroutine residual_eq_2

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine residual_2d (x1, x2, f1, f2) 

!     Evaluate the relevant two residuals according to the main option. 
!     f1 is the residual for eq. (1);  f2 is the residual for eq. (2).

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, intent (in)  :: x1, x2
      real, intent (out) :: f1, f2

      select case (specified_state)

         case ('Ht_MF') ! x1 = P; x2 = T for this case

            P = x1
            call residual_eq_1 (x2, f1)  ! q.v.

            T = x2
            call residual_eq_2 (x1, f2)

      end select

      end subroutine residual_2d

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine frozen_gas_constants (T_specified)

!     This portion of residual_eq_1 has been modularized for reuse in
!     subroutine store_results.
!
!     For a given equilibrium composition, it calculates R_bar, Cp_bar, & gamf.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      real, intent (in) :: T_specified  ! Probably the output T from the
                                        ! equilibrium calculation would suffice
                                        ! for both calls, but this matches how
                                        ! it was before triatomics were treated
!     Local constants:

      real, parameter :: minus_one = -1.0, minus_two = -2.0, two = 2.0

!     Local variables:

      integer :: i

!     Execution:

       R_bar    = zero
      Cv_bar    = zero
      ev_bar    = zero
      h_ref_sum = zero

      do i = 1, n_species
         h_ref_sum = h_ref_sum + mass_fraction(i) * hf(i)
         c_over_M  = mass_fraction(i) / M(i)

         if (theta_v(i) > zero) then

            ev_hat(i) = R_hat * theta_v(i) / &
                          (exp (theta_v(i) / T_specified) - one)
            ev_bar    = ev_bar + c_over_M * ev_hat(i)

         else if (theta_v(i) == minus_one) then  ! CO2 triatomic

            ev_hat(i) = R_hat * ( &
                theta_CO2(1) / (exp (theta_CO2(1) / T_specified) - one) + &
            two*theta_CO2(2) / (exp (theta_CO2(2) / T_specified) - one) + &
                theta_CO2(3) / (exp (theta_CO2(3) / T_specified) - one))
            ev_bar    = ev_bar + c_over_M * ev_hat(i)

         else if (theta_v(i) == minus_two) then  ! NCO triatomic

            ev_hat(i) = R_hat * ( &
                theta_NCO(1) / (exp (theta_NCO(1) / T_specified) - one) + &
            two*theta_NCO(2) / (exp (theta_NCO(2) / T_specified) - one) + &
                theta_NCO(3) / (exp (theta_NCO(3) / T_specified) - one))
            ev_bar    = ev_bar + c_over_M * ev_hat(i)

         end if  ! Else theta_v(i) = zero - not a diatomic or triatomic

         Cv_bar = Cv_bar + c_over_M * Cv_hat(i)
          R_bar =  R_bar + c_over_M
      end do

       R_bar =  R_bar *  R_hat
      Cp_bar =  R_bar + Cv_bar
      gamf   = Cp_bar / Cv_bar

      if (diagnostics) then
         write (lunlog, '(a, 1p, e14.6)') '   gamf: ',  gamf
         write (lunlog, '(a, 1p, e14.6)') '  R_bar: ',  R_bar
         write (lunlog, '(a, 1p, e14.6)') ' Cp_bar: ', Cp_bar
         write (lunlog, '(a, 1p, e14.6)') ' Cv_bar: ', Cv_bar
         write (lunlog, '(a, 1p, e14.6)') ' ev_bar: ', ev_bar
      end if

      end subroutine frozen_gas_constants

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine save_results (lunf, lunp, n, nf, f)

!     Write a 2-D PLOT3D-type n x 1 function file for DPLR and a Tecplot file.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, intent (in) :: lunf, lunp, n, nf
      real,    intent (in) :: f(n,nf)

      character, parameter :: quotes * 1 = '"'

      integer :: i
      logical :: integrate

      write (lunf, '(3i3)') n, 1, nf
      do i = 1, nf
         write (lunf, '(1p, 6e18.10)') f(:,i)
      end do

      close (lunf)

      integrate = false

      write (lunp, '(a)') &
         'TITLE = "Nozzle throat conditions"', &
         'VARIABLES =', &
         '"Radius, m"'

      select case (specified_state) ! For the name of v1
      case ('Ht_MF', 'Ht_Ru')
         write (lunp, '(a)') '"Total enthalpy, J/kg"'
      case ('Rho_T', 'Rho_H', 'Rho_S')
         write (lunp, '(a)') '"Density, kg/m^3"'
      case ('P_T  ', 'P_H  ', 'P_S  ')
         write (lunp, '(a)') '"Pressure, Pa"'
      end select

      select case (specified_state) ! For the name of v2
      case ('Ht_MF', 'Ht_Ru')
         write (lunp, '(a)') '"Mass flux, kg/(m^2.s)"'
         integrate = true
      case ('Rho_T', 'P_T  ')
         write (lunp, '(a)') '"Temperature, K"'
      case ('Rho_H', 'P_H  ')
         write (lunp, '(a)') '"Enthalpy, J/kg"'
      case ('Rho_S', 'P_S  ')
         write (lunp, '(a)') '"Entropy, J/kg-mol"'
      end select

      write (lunp, '(4a)') &        ! Species names
      (quotes, trim (name_input(i)), ' mass fraction', quotes, i = 1, n_species)

      if (specified_state(1:1) /= 'P') write (lunp, '(a)') '"Pressure, Pa"'

      select case (specified_state)
      case ('Rho_T', 'P_T  ')
                                    ! Suppress a duplicate temperature
      case default
         write (lunp, '(a)') '"Temperature, K"'
      end select

      write (lunp, '(a)') &
         '"rho, kg/m^3"', '"h, J/kg"', '"s, J/kg-mol"', '"Frozen a, m/s"'

!     Single zone header:

      write (lunp, '(a)') 'ZONE T="Center to wall"'
      write (lunp, '(a,i3,a)') 'I=', neval, ', J=1, ZONETYPE=Ordered'
      write (lunp, '(a)') 'DATAPACKING=BLOCK'

      write (lunp, '(1p, 8e14.6)') reval(:)

      do i = 1, 2 + n_species
         write (lunp, '(1p, 8e14.6)') fplot(:,i)
      end do

      i = 3 + n_species
      if (specified_state(1:1) /= 'P') then
         write (lunp, '(1p, 8e14.6)') fplot(:,i)
         i = i + 1
      end if

      select case (specified_state)
      case ('Rho_T', 'P_T  ')       
                                    ! Suppress a duplicate temperature
      case default
         write (lunp, '(1p, 8e14.6)') fplot(:,i)
         i = i + 1
      end select

      do i = nfplot - 3, nfplot     ! rho, h, entropy, and frozen speed of sound
         write (lunp, '(1p, 8e14.6)') fplot(:,i)
      end do

      close (lunp)

      if (integrate) then

         call integrate_bulk_mass_and_enthalpy &
                        (neval, reval, v1eval, v2eval, bulk_mass, bulk_enthalpy)
      end if

      end subroutine save_results

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine starting_guess (mode)

!     mode = 0:  Initialize the first nozzle throat calculation.  Others start
!                from the previous solution.
!     mode = 1:  Set [xa, xb] containing a zero for a eq. (1) of separable case
!     mode = 2:  Set [xa, xb] containing a zero for a eq. (2) of separable case
!     mode = 3:  ? what does the 2-D case need?

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, intent (in) :: mode ! See above

      real :: gammae, term, To

      select case (specified_state)

      case ('Ht_MF', 'Ht_Ru')  ! Tahir might explain some time

         select case (mode)

         case (0) ! Starting guesses for T and P (first point only)

            if (v1 < 15.e+6) then
               To     = 6500.
               gammae = 1.18
            else if (v1 < 20.e+6) then
               To     = 7000.
               gammae = 1.16
            else
               To     = 7500.
               gammae = 1.15
            end if

            term = one / (one + half * (gammae - one) * fMsq)
            T    = term * To
            P    = 3.41831 * term**(gammae / (gammae - one)) * (v1**0.397) * v2

            write (lunlog, '(/, (a, 1p, 3e14.6))') &
               ' Initial v1, v2, gammae: ', v1, v2, gammae, ' P0, T0: ', P, T

            x(1) = P     ! For Newton option
            x(2) = T

         case (1) ! Set the interval containing the T that solves eq. (1)

            xa1  = 0.8 * T
            xb1  = 1.1 * T
            tol1 = tol * T

         case (2) ! Set the interval containing the P that solves eq. (2)

            xa2  = 0.8 * P
            xb2  = 1.1 * P
            tol2 = tol * P

         end select

      end select

      end subroutine starting_guess

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine store_results (ibc, n, nf, k, f)

!     Store the flow variables indicated by ibc into f(k,1:nf) for the current
!     evaluation point k (k being in 1:n).  Also store plottable results.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, intent (in) :: ibc    ! DPLR's code:
                                     ! 60: primitive vars. rho_s,u,v,Tv,T
                                     ! 61: primitive vars. P,c_s[2-ns],u,v,Tv,T
                                     ! 62: conserved vars. rho_s,rhou,rhov,Ev,E
      integer, intent (in) :: n, nf, k
      real, intent (inout) :: f(n,nf)

      integer :: i
      real    :: E

!     Ensure consistency with DPLR-derived quantities such as molecular weights:

      select case (specified_state)

         case ('Ht_MF', 'Ht_Ru')

            rho =    P / (R_bar * T)
            u   = gamf * (R_bar * T) * fMsq  ! u**2
            h   = v1 - half * u
            u   = sqrt (u)                   ! u

         case default  ! Plain Gary cases; reuse part of residual_eq_1

            call frozen_gas_constants (T)

            rho = P / (R_bar * T)
            u   = gamf * (R_bar * T) * fMsq      ! u**2
            h   = Cp_bar * T + ev_bar + h_ref_sum
            ho  = h + half * u
            u   = sqrt (u)

      end select

      select case (ibc)

         case (60)

            f(k,1:n_species) = rho * mass_fraction(:)
            f(k,n_species+1) = u
            f(k,n_species+2) = zero
            f(k,n_species+3) = T
            f(k,n_species+4) = T

         case (61)

            f(k,1) = P
            f(k,2:n_species) = mass_fraction(2:n_species)
            f(k,n_species+1) = u
            f(k,n_species+2) = zero
            f(k,n_species+3) = T
            f(k,n_species+4) = T

         case (62)

            E = rho * h - P
            f(k,1:n_species) = rho * mass_fraction(:)
            f(k,n_species+1) = rho * u
            f(k,n_species+2) = zero
            f(k,n_species+3) = rho * ev_bar
            f(k,n_species+4) = E

      end select

!     Plottable results:

      fplot(k,1)   = v1
      fplot(k,2)   = v2;              i = 2 + n_species
      fplot(k,3:i) = mass_fraction(:)
      fplot(k,i+1) = P
      fplot(k,i+2) = T
      fplot(k,i+3) = rho
      fplot(k,i+4) = h
      fplot(k,i+5) = s
      fplot(k,i+6) = u

      end subroutine store_results

   end program nozzle_throat_conditions

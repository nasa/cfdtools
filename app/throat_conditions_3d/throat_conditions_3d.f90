!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program throat_conditions_3d
!
!  Description:
!  ============
!
!     Calculate boundary conditions at an arc-jet nozzle throat that correspond
!     to a choice of distributions of two real gas flow variables over the 2-D
!     cross-section, which may be rectangular, circular, or semi-elliptic.
!     Results are output in the form of a pointwise boundary condition for the
!     DPLR 3-D flow solver.  This is a generalization of the earlier utility,
!     NOZZLE_THROAT_CONDITIONS, which applies to axisymmetric flow within circ-
!     ular nozzles. The flow is still considered to be 1-dimensional (u,0,0) at
!     any throat point.
!
!  Approach:
!  =========
!
!     The equilibrium gas compositions are calculated with a C implementation
!     of the relevant portions of the CEA program by Gordon and McBride dated
!     May 20, 1998.  Gary Allen wrote the C package at NASA Ames Research Center
!     and continues to maintain it.  Initially, a C driving program was invoked
!     here via a system call whenever a gas composition was required.  This has
!     now been streamlined with a C interface routine (equilibrium_gas.c) even
!     though its inefficiency was not much of an issue.
!
!  Extension to a Bidirectional Boundary:
!
!     1-D profiles of the specified variables V1, V2 are defined in the control
!     file for each of the major (horizontal) and minor (vertical) axes of the
!     nozzle cross-section, on the indicated uniform grids.  The outer boundary
!     inputs are derived by linear interpolations w.r.t. arc length.  (They are
!     probably constant on the outer edge, but need not be.)  Interior uniform
!     grid point values of V1 and V2 are derived from the edge values via inter-
!     polation.  The specified calculation for the complete in-flow may then be
!     applied to every point (i,j) of the uniform grid.  Point (1,1) at the
!     nozzle center uses a starting guess algorithm.  Further points (i,j) are
!     started from the (i-1,j) point solution except that the (1,j) solution is
!     started from the (1,j-1) solution.
!
!     If economy mode is specified (as recommended), only the axis grid points
!     are solved for by the full method.  Interior points are determined by
!     interpolation, with an extra call to the equilibrium compostion routine
!     to ensure valid flow states.  Note that this adjustment affects results
!     slightly when the option to specify bulk flow quantities is used, but the
!     final difference between the requested and achieved bulk quantities is
!     quite insignificant given the uncertainties in arc-jet operations.
!
!     Three choices are available for the interim uniform grid, for historical
!     reasons that are associated with the three choices for interpolation from
!     flow profiles on the major and minor axes to the interior points.  Some
!     combinations are no longer recommended, but have been retained.
!
!     Recommended choices:
!
!        Nozzle Type           Grid Topology           Interpolation Method
!
!        Rectangular           Rectangular             Cartesian product
!        Circular              Polar                   TFI or Cartesian product
!        Semi-ellipse          Specialized algebraic   Scaling of minor profile
!
!     Programmer note:
!
!     Before Cartesian product interpolation was implemented, TFI was found to
!     behave poorly enough for the rectangle case (yet well for the ellipse)
!     that a polar form of grid was implemented for the rectangle nozzle case.
!     This option treats the rectangle just as for the semi-ellipse initially,
!     with a slight adjustment to ensure that one of the polar angles passes a
!     spoke through the outer corner of the associated rectangle.  Then the
!     polar grid (no longer quite uniform in general) is replaced by another
!     grid with the spokes extrapolated to the desired rectangle.  (Slight
!     concavity in flat top surfaces formed by total enthalpy, say, appear
!     unavoidable with TFI on the simple rectangular grid.  However, Cartesian
!     product interpolation has since been found to behave well.)
!
!     If a target grid is specified, bidirectional interpolation within the
!     preliminary uniform grid is performed.  The boundaries of the two grids
!     should match for sensible results.
!
!  Flow Specification Options (one pair for each axis):
!  ===========================
!
!     Gary provides 6 options for the given pair of flow specifications.  These
!     are preserved as options here, but normal usage is likely to specify the
!     total enthalpy/mass flux combination.
!
!     Option   Flow variables specified          Units
!
!     Ht_MF    total stagnation enthalpy         J/kg
!              mass flux                         kg/(m^2.s)
!              (safeguarded Newton iteration)
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
!     Ht_MF Algorithm:
!     ----------------
!
!     This option determines throat conditions at the indicated frozen Mach
!     number (not necessarily 1) for given total stagnation (reservoir) enthalpy
!     and mass flow rate per unit area.  Two nonlinear equations are solved via
!     a safeguarded two-variable Newton iteration with central difference deriv-
!     atives.
!
!     The equations at a point r from the centerline are:
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
!        u       =  Mf af  =  Mf sqrt (gamf Rbar T)
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
!     The two variables to be solved for are P and T, with the mass fractions
!     essentially converging as well during the iterative solution.
!
!     Ht_MF Starting Guesses (for a Profile from the 2-D case):
!     ---------------------------------------------------------
!
!     For a specified mass fraction of argon (say 10%), fractions for N2 and
!     O2 are readily derived from the values 0.767 and 0.233 for pure air, with
!     other mass fractions set to 0.
!
!     For the current point in space, the previous solution serves as a good
!     starting guess.  For the first point r in space, starting guesses for P
!     and T are implemented as follows:
!
!        P  =  C (1 / t1) ^ (gamf / (gamf - 1)) rho u (r) ho(r)^n
!        T  =    (1 / t1) * To
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
!     Control file for THROAT_CONDITIONS_3D (2 flow variables prescribed)
!     ========================================================================
!     Rectangular   ! Nozzle type: Rectangular | Circular | Semi-Elliptic | Pol.
!     Ht_MF         ! V1_V2 spec: Ht_MF, Rho_T, Rho_H, Rho_S, P_T, P_H, or P_S
!     Economy       ! Fill mode: Full | Economy (Axes + Interpoln.) | Quick look
!     Cartesian     ! Interpolation method: TFI | Cartesian[ Product] | Scaling
!     1.0           ! Frozen Mach number
!     ------------------------------------------------------------------------
!     Iteration controls
!     ------------------------------------------------------------------------
!     1             ! Iterate for target bulk enthalpy?   0 = no | 1 = yes
!     2.9E+07       ! Target bulk enthalpy, J/kg
!     1             ! Iterate for target bulk mass flow rate?  0 | 1
!     90.           ! Target bulk mass flow rate, kg/s
!     ------------------------------------------------------------------------
!     V_1 specifications, Horizontal Axis
!     ------------------------------------------------------------------------
!     Linear        ! Uniform, Linear, Parabolic, Nthdegree, Sinusoid, ...
!     0.   3.E+07   ! (r, V1) at nozzle center-line
!     1.   3.E+07   ! (r, V1) at nozzle wall
!     999.   999.   ! "Width" (several profiles) or Width + Steepness (Sigmoid)
!     none          ! Dataset file name (if relevant)
!     ------------------------------------------------------------------------
!     V_2 specifications, Horizontal Axis
!     ------------------------------------------------------------------------
!     Sinusoid      ! Shape (..., Mthdegree, Sigmoid, Gaussian, Lorentz, ...)
!     0.   250.     ! (r, V2) at nozzle center-line
!     1.   300.     ! (r, V2) at nozzle wall
!     2.     999.   ! "Width" or "Steepness"
!     none          ! Dataset file name
!     ------------------------------------------------------------------------
!     V_1 specifications, Vertical Axis
!     ------------------------------------------------------------------------
!     Linear        ! Shape:  Uniform, Linear, Parabolic, Sinusoid, Dataset ..
!     0.   3.E+07   ! (r, V1) at nozzle center-line
!     0.5  3.E+07   ! (r, V1) at nozzle wall
!     999.   999.   ! "Width" and "Steepness"
!     none          ! Dataset file name (if relevant)
!     ------------------------------------------------------------------------
!     V_2 specifications, Vertical Axis
!     ------------------------------------------------------------------------
!     Sinusoid      ! Shape
!     0.   250.     ! (r, V2) at nozzle center-line
!     0.5  300.     ! (r, V2) at nozzle wall
!     2.     999.   ! "Width" and "Steepness"
!     none          ! Dataset file name
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
!     33 17          ! # initial uniform (i,j) pts.; ni < 0 for diagnostics
!     grid.g         ! File containing target grid for interpolating to | none
!     61             ! DPLR's input profile BC code:  60 | 61 | 62
!     conditions.f   ! File name for output throat conditions (*.f | *.pbca)
!     conditions.dat ! Plottable data (Tecplot ASCII)
!     2.             ! Controls blending of polar grids for the rectangle case
!     0 0 1          ! Converts flow solver xyz convention to internal order
!     0 1 0          ! In this example, the flow solver x is downstream, y is
!     1 0 0          ! up, and z is spanwise (right-handed)
!
!  Output format (portion of DPLR *.pbca file):
!  ============================================
!
!     Unless the file extension for the output throat conditions is 'pbca', we
!     simply write a PLOT3D function file that can be checked by plotting.
!     Otherwise, some editing is avoided because this program inserts the
!     lines expected between blocks, but the partial pbca file is not plottable.
!
!  Further Usage Notes:
!  ====================
!
!     (0) Run the code like this:   throat_conditions_3d < xxx.inp > xxx.log
!
!     (1) Internal units have x spanwise/horizontal and y vertically up. Since
!         z is assumed to be constant, it is not relevant except when the
!         interim uniform-grid results are to be interpolated to a CFD grid.
!         The permutation matrix input should convert from the flow solver xyz
!         convention to make +spanwise <-> +x and +vertical <-> +y.  The 3rd CFD
!         coordinate is then used as the constant z during the 3-space flow
!         interpolations to the CFD grid block face at the throat.
!         Since the working uniform grid has z = 0, it makes sense for the
!         target grid to match that (after any permutation), although the
!         projection process will tolerate permuted z = any reasonable constant.
!         The plottable output interpolations contain just [permuted] x and y.
!
!     (2) The plain uniform rectangle grid was initially abandoned when it
!         produced disappointing (non-convex) TFI results.  The rectangle was
!         temporarily treated as an ellipse via a polar grid.  However, when
!         Cartesian products were provided as an alternative to TFI, the
!         plain rectangle grid option was restored.  Now, the nozzle type and
!         the interpolation type are independent controls, with 'R' meaning
!         the plain rectangular uniform grid and 'P' (polar) meaning the
!         polar form of grid for rectangular nozzles (no longer recommended).
!
!     (3) When rectangular nozzles are treated on an interim polar grid with
!         elliptic outer boundary (to exploit the better behavior of TFI on
!         such boundaries), that interim grid is converted to the desired
!         rectangle shape by stretching the spokes.  It remains a polar grid.
!         The plain rectangle form (uniform radial spacing) is blended with the
!         interim elliptic form in order to moderate cell skewness, for all
!         values of control input rectangle_power = p, which is applied to the
!         ratio (i - 1)/(ni - 1) at radial point i.
!
!             0 < p <= 1. = more of plain uniform rectangle grid
!                 p >  1. = more of interim elliptic grid
!
!     (4) The semi-ellipse case with a wall at the major axis needs a vertical
!         profile between center and upper wall that decays at both ends.  Use
!         the MTHDEGREE option in this case (where the "center" control value
!         now applies to midway between r_center and r_wall).  The NTHDEGREE
!         option has the peak at r_center.
!
!     (5) The polar grid that appears natural for the semi-ellipse proves to
!         be inferior to the specialized grid introduced later, consisting of
!         a family of i-lines arranged to  be normal to the flat and outer
!         boundaries at both ends.  The vertical center-line flow profile is
!         simply scaled to fill the interior flow.
!
!     (6) Use of fill_mode = 'economy' means full calculations are performed
!         along the major and minor axes only, followed by interpolations
!         for the interior points.  Results appear almost indistinguishable
!         from the full calculations in the examples compared, although the
!         interpolated flow quantities are adjusted with equilibrium_composition
!         calls to be sure they satisfy the equations of state.
!
!     (7) For a quick look at the input distributions of V1 and V2, use option
!         fill_mode = 'quick_look'.  The program saves these distributions in
!         plottable form then stops.  See (8).
!
!     (8) Before invoking the bulk flow iteration, perform a "quick-look" run
!         to see what the integrated quantities are, so that reasonable target
!         bulk values can be attained knowing that the problem achieves them by
!         scaling the input wall values of the axis profiles.  The range that
!         is searched by the zero-finder is [0.5, 2.0] for each multiplier.
!
!     (9) The integrated bulk quantities displayed are adjusted to refer to the
!         entire nozzle cross-section.  Enter target values accordingly.
!
!    (10) Excessive "step-halving iteration failure" messages usually mean that
!         unreasonable conditions are being sought.  Use of scaled variables and
!         scaled residuals makes it very unlikely that the Newton iterations
!         will fail.  They can be displayed by entering -ni_uniform for the
!         first dimension of the desired uniform internal grid.
!
!    (11) The Lorentz profile tends toward zero less rapidly than the Gaussian,
!         and is more peaky at r = 0 for the same "width" parameter.
!
!  History:
!  ========
!
!     06/09/04  DKP/DAS  Earlier THROAT_CONDITIONS:   Dinesh Prabhu/D. Saunders.
!     12/19/05  TG /DAS  Initial NOZZE_THROAT_CONDITIONS (Ht_Ru option).
!     12/20/05  TG /DAS  Added the Ht_MF alternative (Newton iteration) to the
!                        Ht_Ru option (non-derivative zero finder). Newton wins.
!     11/28/06  DKP/DAS  Lowered the convergence tolerance from 1.e-5 to 1.e-6.
!     11/29/06   "   "   Put the tolerance back to 1.e-5.  It seems the "Gary"
!                        part of the scheme can give different results for
!                        different initial estimates of the mass fractions.
!                        This needs to be investigated.
!     02/21/07-   DAS    THROAT_CONDITIONS_3D derived from NOZZLE_THROAT_CONDS.
!     02/28/07           for rectangular or semi-elliptic nozzles.
!     03/05/07  TG /DAS  Economy mode for the 'Ht_MF' option (semi-axes
!                        followed by TFI of interior P & T, etc.).
!     03/06/07   "   "   Added nth deg. polynomial choice to get flatter tops.
!     03/07/07    DAS    Tried to overcome concavity in rectangle flow by
!                        morphing to a Cartesian grid with elliptic outer
!                        boundary, but if anything it's slightly poorer.
!     03/08/07     "     Try a polar grid for the rectangle case: use an
!                        interim elliptic grid first, then stretch the spokes
!                        to form the desired rectangle.
!     03/10/07     "     Option to blend the stretched rectangle with the
!                        interim ellipse to soften the corner along the
!                        diagonal.
!     03/12/07     "     Interpolation to CFD block face(s) completed.
!       "  "    TG/DAS   Provision for frozen Mach number other than 1.
!       "  "     "  "    Belated realization that the semi-ellipse case with a
!                        wall at the major axis needs a new profile shape -
!                        introduced MTHDEGREE option.
!     03/13/07  D.Prabhu Dinesh suggested the sigmoid shape function, which
!                 DAS    required another control dubbed "steepness".
!     03/16/07     "     Added the 'quick-look' option.  Implemented the
!                        Cartesian product option as an alternative to TFI
!                        to keep Tahir happy (with concavity predicted, at
!                        least on the plain rectangle grid).
!     03/21/07     "     The 'S' option, now distinct from 'E' for the semi-
!                        ellipse, introduces a better grid consisting of a
!                        family of i-lines normal to the major axis and the
!                        outer wall at both ends.  This appears the best
!                        choice, whether the minor profile peak is near the
!                        middle of the vertical axis (wall at the major axis)
!                        or not.
!     03/29/07     "     Incorporated the scaling found to overcome odd
!                        behavior with a low-pressure NOZZLE_THROAT_CONDITIONS
!                        case.
!     04/07/07     "     Bulk mass flow rate and bulk enthalpy options are
!                        functioning.
!     04/09/07     "     Introduced "Circular" as a nozzle-type, where semi-
!                        ellipse had seemed adequate earlier.  A plain polar
!                        grid is implied, while the semi-ellipse case now
!                        implies a solid wall along the major axis, and the
!                        minor axis profiles should be appropriate (MTHDEGREE
!                        or UNIFORM being about the only sensible choices).
!     04/11/07     "     The finite differencing needed more careful starting
!                        guesses for the mass fractions to reduce noise.
!     05/02/07     "     A uniform profile should work with the wall value, not
!                        the center value, because the bulk iterations adjust
!                        wall values.  Added the option to tailor the flow
!                        interpolated to the CFD grid to look more like a
!                        DPLR pbca file if the file name ends in pbca.
!                        Otherwise, it remains a plottable function file.
!     06/14/07     "     Added Gaussian profile option.
!     07/05/07     "     Added Lorentz  profile option.
!     10/31/08     "     Replaced the system call in equilibrium_composition.f90
!                        to Gary's CEA driving program with a call to interface
!                        routine equilibrium_gas.c, without change here except
!                        to update the documentation and to control the extra
!                        printing to standard output.
!     03/25/10     "     Maria Pulsonetti needed CO2, CO, C species in the top
!                        level tables (already handled at the lower levels).
!                        WARNING:  Polyatomic species have more than one
!                        characteristic vibrational temperature, so if there is
!                        any significant CO2 fraction, enthalpy and temperature
!                        calculations are only approximate.
!     09/14/11      "    Transcribed Tahir's more thorough handling of the
!                        CO2 and NCO triatomics with three vibrational modes as
!                        first done in NOZZLE_THROAT_CONDITIONS.
!                        Note that lewis.c does not handle NCO yet, though.
!
!                        Redirecting standard output (unit 6) has been seen to
!                        interact badly with singular matrix diagnostics from
!                        the lewis function.  Therefore, open a log file on
!                        unit 7 and write to that.  Any C code diagnostics
!                        still come to the screen.
!
!                        [Later:  This experience from the axisymmetric solver
!                        remains a mystery.  For now, we stay with unit 6, and
!                        diagnostics from both languages go to the same file
!                        or to the screen if standard output is not redirected.]
!
!                        Also, Gary points out that Ar+ should have been among
!                        the subset of species here all along (at least for
!                        high temperature cases).
!
!  Authors:
!  ========
!
!     Analysis:          Tahir Gokcen,   ELORET Corp./NASA Ames Research Center
!     Implementation:    David Saunders, ELORET Corp./NASA Ames Research Center
!                                   now: ERC, Inc./NASA ARC
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure    ! Needed for surface interpolation package & I/O
   use xyzq_io_module          ! Simplifies reading CFD grid & writing results
   use trigd

   implicit none

!  Local constants:

   integer, parameter ::    &
      itmax       = 30,     &  ! Newton iteration limit
      lunctl      = 5,      &  ! control file (standard input/command line)
      lunlog      = 6,      &  ! 6 here can cause standard output f90/C mix-up
      lunin       = 1,      &  ! for possible datasets specifying V1 and/or V2
      lunoutf     = 2,      &  ! for throat-condition results in PLOT3D form
      lunoutg     = 3,      &  ! for throat-condition uniform grid if not input
      lunplot     = 4,      &  ! for throat-condition uniform grid results
      lunprint    = lunlog, &  ! lunprint = -lunlog suppresses Newton itn. print
      max_options = 7,      &  ! number of specified state cases handled
      max_species = 21,     &  ! length of species names dictionary
      nv          = 2          ! # nonlinear equations (that might be) solved

   real, parameter ::       &
      half   = 0.5,         &
      one    = 1.0,         &
      third  = 1.0 / 3.0,   &
      toler  = 1.E-6,       &  ! Residual tolerance (Newton, scaled x and f)
      two    = 2.0,         &
      zero   = 0.0,         &
      R_hat  = 8.31441e+03     ! Universal gas constant, J/(kmol.K)

   logical, parameter ::    &
      false = .false.,      &
      true  = .true.

   character, parameter ::  &
      dataset * 7 = 'DATASET'  ! Prescribed option for V1 or V2

!  Local variables:  'h' => horizontal axis; 'v' => vertical axis

   integer :: &
      i, ibc, ier, ios, iter, ip(nv), j, j1, jmajor, k, l, n, nb, nblocks,     &
      ndim, nf, nfplot, ni, ni_uniform, nj, nj_uniform, nj_ellipse, n_species, &
      n_v1_shape_h, n_v1_shape_v, n_v2_shape_h, n_v2_shape_v, status_CEA

   real :: &
      a, ae, b, c_over_M, Cp_bar, Cv_bar, e, ev_bar, gamf, gamf0, h, ho,       &
      h_ref_sum, P, P0, pi, r, rectangle_power, rho, R_bar, s, T, T0, t1,      &
      alpha, cube_root_eps, dr, fMsq, fnorm, fnorm0, frozen_Mach,              &
      bulk_enthalpy, target_bulk_enthalpy, bulk_mass, target_bulk_mass,        &
      quadrature_multiplier,                                                   &
      r_center_h, r_center_v, r_wall_h, r_wall_v,  u, v1, v2,                  &
      v1_center_h, v1_center_v, v1_wall_h, v1_wall_v, v1_width_h, v1_width_v,  &
      v2_center_h, v2_center_v, v2_wall_h, v2_wall_v, v2_width_h, v2_width_v,  &
      v1_steepness_h, v1_steepness_v, v1_wall_h0, v1_wall_v0, xyz_switch(3,3), &
      v2_steepness_h, v2_steepness_v, v2_wall_h0, v2_wall_v0

   real :: &
      mol_weights(max_species), ref_enthalpies(max_species),                   &
      species_multipliers(max_species), theta_vibs(max_species),               &
      theta_CO2(3), theta_NCO(3)

   real :: &
      AJ(nv,nv), aitch(nv), dx(nv), f(nv), ffd(-1:1,nv,nv), fscale(nv),        &
      x(nv), xfd(nv), xlast(nv), xscale(nv)

   real, allocatable, dimension (:) :: &
      Cv_hat, ev_hat, hf, Mol, mass_fraction, mass_fraction_0, mass_fraction_1,&
      mole_fractions, theta_v

   real, pointer, dimension (:) :: &  ! Non-analytic profile definitions
      v1_shape_r_h, v1_shape_f_h, v1_shape_r_v, v1_shape_f_v,                  &
      v2_shape_r_h, v2_shape_f_h, v2_shape_r_v, v2_shape_f_v

   real, allocatable, dimension (:,:,:) :: &
      fplot,   &     ! Specified and calculated flow, internal (not DPLR) units
      P_T,     &     ! For economy mode interpolations
      uniform_v1_v2  ! Prescribed variables evaluated on uniform grid

   logical :: &
      bulk_enthalpy_specified, bulk_mass_specified, cell_centered,             &
      diagnostics, interpolate_to_grid, pbca, skip, show_details

   character :: &
      filename * 80, fill_mode * 1, nozzle_type * 1, interp_type * 1,          &
      specified_state * 5, v1_shape_h * 9, v1_shape_v * 9, v2_shape_h * 9,     &
      v2_shape_v * 9

   character :: &
      main_options(max_options) * 5, species_names(max_species) * 4

   character, allocatable, dimension (:) :: &
      name_input * 4, name_gary * 4

   type (grid_type), pointer, dimension (:) :: &
      uniform_grid, &  ! Single-block preliminary grid and calculated flow
                       ! in DPLR units; must be a pointer to match I/O package
      cfd_grid,     &  ! One array element per target grid block: x, y, z and
                       ! interpolated flow in DPLR units
      cfd_plot         ! Plottable form of interpolations to CFD grid

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
      'Ht_MF', 'Rho_T', 'Rho_H', 'Rho_S', 'P_T  ', 'P_H  ', 'P_S  ' /

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


   pi = two * acos (zero)

!! No - stay with redirecting standard output to any file name on command line.
!! open (lunlog, file='throat_conditions_3d.log', status='unknown')

   call control_inputs ()  ! Internal procedure below

   if (ios /= 0) go to 99


!  Set up work-space associated with the indicated species:
!  --------------------------------------------------------

   n = n_species

   allocate (Cv_hat(n), ev_hat(n), hf(n), Mol(n), mole_fractions(n), theta_v(n))

   do i = 1, n_species

      call lookup2 (max_species, species_names, name_input(i), n)

      if (n <= 0) then
         write (lunlog, '(/, 2a)') ' Unhandled species: ', name_input(i)
         go to 99
      end if

      Cv_hat(i)  = R_hat * species_multipliers(n)
      Mol(i)     = mol_weights(n)
      hf(i)      = ref_enthalpies(n)
      theta_v(i) = theta_vibs(n)

   end do


!  Set up a uniform single-block grid for the initial flow evaluations.
!  The scheme used for interpolation to a CFD grid requires 3-space surface
!  data and use of a derived data type.
!  ------------------------------------------------------------------------

   allocate (uniform_grid(1))

   call uniform_grid_gen ()


!  Discretize the indicated profiles for V1 and V2 on this uniform grid:
!  ---------------------------------------------------------------------

   allocate (uniform_v1_v2(2,ni_uniform,nj_ellipse))

   select case (nozzle_type)

      case ('S')  ! Semi-elliptic nozzle with non-polar grid

         call expand_ellipse_profiles (ni_uniform, nj_uniform,                 &
                                       uniform_grid(1)%x, uniform_grid(1)%y,   &
                                       uniform_v1_v2)

      case default  ! Rectangle and polar semi-ellipse

         call expand_axis_profiles (ni_uniform, nj_ellipse,                    &
                                    uniform_grid(1)%x, uniform_grid(1)%y,      &
                                    uniform_v1_v2)
   end select

   if (ios /= 0) go to 99  ! Bad shape spec.?

   ier = 0

!  Option to seek specified bulk flow rate:
!  ----------------------------------------

   if (bulk_mass_specified) call bulk_iteration (1)

   if (ier /= 0) go to 99


!  Option to seek specified bulk enthalpy:
!  ---------------------------------------

   if (bulk_enthalpy_specified) call bulk_iteration (2)

   if (ier /= 0) go to 99


!  Option to check the specified distributions only:
!  -------------------------------------------------

   if (fill_mode == 'Q') then  ! "Quick look"

      if (nozzle_type == 'P') call extend_spokes ()  ! From interim polar form

      call save_v1_v2 ()

      go to 99

   end if


!  Calculate remaining flow conditions on the uniform grid:
!  --------------------------------------------------------

   call q_allocate (uniform_grid(1), nf, ios)  ! For computed flow in DPLR units

!  Store all possible plottable data, but some of it may be suppressed.

   nfplot = n_species + 8  ! v1, v2, P, T, rho, h, s, a

   allocate (fplot(nfplot,ni_uniform,nj_ellipse))  ! Internal units
   allocate (P_T(2,ni_uniform,nj_ellipse))         ! If economy mode is in use

   if (ndim > 0) then

      v1 = uniform_v1_v2(1,1,1)
      v2 = uniform_v1_v2(2,1,1)

      call starting_guess ()  ! For the first point only

   end if


   call flow_conditions ()

   if (ier /= 0) go to 99


!  Convert interim elliptic grid to rectangular?  (Still polar.)
!  -------------------------------------------------------------

   if (nozzle_type == 'P') call extend_spokes ()


!  Interpolate uniform grid results to a CFD boundary grid?
!  --------------------------------------------------------

   if (interpolate_to_grid) call redistrib_2d ()


!  Save results:
!  -------------

   call save_results ()

!  Done.


99 continue  ! Avoid possible system-dependent STOP behavior.

!  Internal procedures for program nozzle_throat_conditions:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine bulk_iteration (item)

!     Solve the 1-variable nonlinear equation that specifying the indicated bulk
!     flow quantity leads to.  The bulk enthalpy is dependent upon the bulk mass
!     flow (but not vice versa), so it should be determined second.
!     The one variable is a multiplier for the input wall level(s) of the two
!     axis profiles for mass flux (item = 1) or total enthalpy (item = 2).
!     The center levels of these quantities are fixed at the input values.

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

      ni     = uniform_grid(1)%ni
      nj     = uniform_grid(1)%nj
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

                  v2_wall_h = v2_wall_h0 * x
                  v2_wall_v = v2_wall_v0 * x

               case (2) ! Total enthalpy profile is being adjusted

                  v1_wall_h = v1_wall_h0 * x
                  v1_wall_v = v1_wall_v0 * x

            end select

!           Discretize the axis profiles for V1 = total enthalpy and V2 = mass
!           flux, and interpolate to the rest of the uniform grid:

            select case (nozzle_type)

               case ('S')  ! Semi-elliptic with non-polar grid

                  call expand_ellipse_profiles (ni_uniform, nj_uniform,        &
                                                uniform_grid(1)%x,             &
                                                uniform_grid(1)%y,             &
                                                uniform_v1_v2)

               case default  ! Rectangle cases R|P or polar Circle/semi-Ellipse

                  call expand_axis_profiles (ni_uniform, nj_ellipse,           &
                                             uniform_grid(1)%x,                &
                                             uniform_grid(1)%y,                &
                                             uniform_v1_v2)
            end select

!           Integrate the bulk quantities:

            call integrate_bulk_mass_and_enthalpy (ni, nj, 2,                  &
                                                   uniform_grid(1)%x,          &
                                                   uniform_grid(1)%y,          &
                                                   uniform_v1_v2,              &
                                                   bulk_mass, bulk_enthalpy)
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

               v2_wall_h = v2_wall_h0 * x
               v2_wall_v = v2_wall_v0 * x

               write (lunlog, '(/, a, 1p, 2e15.7, a)')                         &
          ' Adjusted mass flux profile wall values (major, minor axis):     ', &
                  v2_wall_h, v2_wall_v, ' kg/(m^2.s)'

            case (2) ! Total enthalpy profile is being adjusted

               v1_wall_h = v1_wall_h0 * x
               v1_wall_v = v1_wall_v0 * x

               write (lunlog, '(/, a, 1p, 2e15.7, a)')                         &
          ' Adjusted total enthalpy profile wall values (major, minor axis):', &
                  v1_wall_h, v1_wall_v, ' J/kg'

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

      subroutine cartesian_product (nozzle_type, ni, nj, nf, f)

!     Fill the interior function values, including outer boundary(s), via the
!     Cartesian product method, for a quadrant of a rectangular or elliptic
!     nozzle.  The semi-axis function values are transformed to [0, 1].  Each
!     appropriate product is calculated then transformed back to the axis
!     function range.
!
!     Polar grids are handled differently from plain rectangular grids because
!     the vertical axis is at j = nj, not i = 1.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Argument:

      character, intent (in)                            :: nozzle_type * 1
      integer,   intent (in)                            :: ni, nj, nf
      real,      intent (inout), dimension (nf, ni, nj) :: f

!     Local variables:

      integer :: i, j, n
      real    :: f1, f2, fscale, fshift, r, rm1

!     Execution:

      do n = 1, nf

         f1 = f(n,1,1);  f2 = f1
         do j = 1, nj, nj - 1
            do i = 1, ni
               f1 = min (f1, f(n,i,j))
               f2 = max (f2, f(n,i,j))
            end do
         end do

         if (f1 == f2) then
            fscale = one / f1
            fshift = zero
         else
            fscale = one / (f2 - f1)
            fshift = -f1 * fscale
         end if

         select case (nozzle_type)

            case ('R')  ! Rectangular grid (plain)

               do j = 2, nj - 1
                  f2 = fscale * f(n,1,j) + fshift
                  do i = 2, ni
                     f1 = fscale * f(n,i,1) + fshift
                     f(n,i,j) = (f1 * f2 - fshift) / fscale
                  end do
               end do

            case ('C', 'E', 'P')  ! Polar grid of one type or another
                                  ! Not really a Cartesian product now
               do j = 2, nj - 1
                  r = real (j - 1) / real (nj - 1) ! Since it's uniform in theta
                  rm1 = one - r                    ! Should do all n at once
                  do i = 2, ni                     ! How is comparison with TFI?
                     f(n,i,j) = rm1 * f(n,i,1) + r * f(n,i,nj)
                  end do
               end do

         end select

      end do  ! Next function

      end subroutine cartesian_product

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine control_inputs ()

!     Read the control file.  All variables are inherited from the host.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Echo the control file to the output log (standard output):

      do ! Until EOF
         read (lunctl, '(a)', iostat=ios) filename ! Handy buffer
         if (ios < 0) exit
         write (lunlog, '(1x, a)') trim (filename)
      end do

      ios = 0;  rewind (lunctl)

!!    cube_root_eps = epsilon (cube_root_eps) ** third ! For central differences
      cube_root_eps =                   1.e-7 ** third ! In view of the Gary pkg

!     Read the control file:

      do i = 1, 3
         read (lunctl, *)
      end do

      read (lunctl, *) nozzle_type;  call upcase (nozzle_type)

      select case (nozzle_type)

         case ('E', 'S')  ! Semi-elliptic

            quadrature_multiplier = 2.

         case default

            quadrature_multiplier = 4.

      end select

      read (lunctl, *) specified_state  ! Main option descriptor

      call lookup2 (max_options, main_options, specified_state, n)

      if (n > 2) then  ! Just call Gary's package at each point r
         ndim = 0
!!    else if (n == 2) then
!!       ndim = 1      ! 2 separable equations (nested 1-variable iterations)
      else if (n == 1) then
         ndim = 2      ! 2 nonlinear equations (2-variable Newton iteration)
      else
         write (lunlog, '(/, 2a)') ' Unhandled main option: ', specified_state
         ios = 1;  go to 99
      end if

      read (lunctl, *) fill_mode    ! Full interior | axes+interp (economy mode)
      call upcase (fill_mode)

      read (lunctl, *) interp_type  ! TFI or Cartesian product
      call upcase (interp_type)

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

!     Variable 1 specs, horizontal axis:

      do i = 1, 3
         read (lunctl, *)
      end do

      read (lunctl, *) v1_shape_h
      read (lunctl, *) r_center_h, v1_center_h
      read (lunctl, *) r_wall_h,   v1_wall_h;   v1_wall_h0 = v1_wall_h
      read (lunctl, *) v1_width_h, v1_steepness_h
      read (lunctl, *) filename

      call upcase (v1_shape_h)

      if (trim (v1_shape_h) == dataset) then ! See another internal procedure

         call read_discretized_shape (n_v1_shape_h, v1_shape_r_h, v1_shape_f_h)
         if (ios /= 0) go to 99

      end if

!     Variable 2 specs, horizontal axis:

      read (lunctl, *);  read (lunctl, *);  read (lunctl, *)
      read (lunctl, *) v2_shape_h
      read (lunctl, *) r_center_h, v2_center_h
      read (lunctl, *) r_wall_h,   v2_wall_h;   v2_wall_h0 = v2_wall_h
      read (lunctl, *) v2_width_h, v2_steepness_h
      read (lunctl, *) filename

      call upcase (v2_shape_h)

      if (trim (v2_shape_h) == dataset) then

         call read_discretized_shape (n_v2_shape_h, v2_shape_r_h, v2_shape_f_h)
         if (ios /= 0) go to 99

      end if

!     Variable 1 specs, vertical axis:

      read (lunctl, *);  read (lunctl, *);  read (lunctl, *)
      read (lunctl, *) v1_shape_v
      read (lunctl, *) r_center_v, v1_center_v
      read (lunctl, *) r_wall_v,   v1_wall_v;   v1_wall_v0 = v1_wall_v
      read (lunctl, *) v1_width_v, v1_steepness_v
      read (lunctl, *) filename

      call upcase (v1_shape_v)

      if (trim (v1_shape_v) == dataset) then ! See another internal procedure

         call read_discretized_shape (n_v1_shape_v, v1_shape_r_v, v1_shape_f_v)
         if (ios /= 0) go to 99

      end if

!     Variable 2 specs, vertical axis:

      read (lunctl, *);  read (lunctl, *);  read (lunctl, *)
      read (lunctl, *) v2_shape_v
      read (lunctl, *) r_center_v, v2_center_v
      read (lunctl, *) r_wall_v,   v2_wall_v;   v2_wall_v0 = v2_wall_v
      read (lunctl, *) v2_width_v, v2_steepness_v
      read (lunctl, *) filename

      call upcase (v2_shape_v)

      if (trim (v2_shape_v) == dataset) then

         call read_discretized_shape (n_v2_shape_v, v2_shape_r_v, v2_shape_f_v)
         if (ios /= 0) go to 99

      end if

      a = r_wall_h - r_center_h  ! Major semi-axis
      b = r_wall_v - r_center_v  ! Minor   "   "

!     Starting guesses for the species mass fractions:

      read (lunctl, *);  read (lunctl, *);  read (lunctl, *)
      read (lunctl, *) n_species

      nf = 5 + n_species  ! For 3-D and two temperatures (though Tv = T)

      n = n_species
      allocate (name_input(n))
      allocate (mass_fraction(n), mass_fraction_0(n), mass_fraction_1(n))

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
      read (lunctl, *) ni_uniform, nj_uniform  ! # points for initial calcs.

      show_details = ni_uniform < 0    ! True means all possible printout;
      diagnostics  = show_details      ! false means suppress most but not all;
      ni_uniform   = abs (ni_uniform)  ! the original "diagnostics" control
                                       ! variable is toggled accordingly.

!     Target CFD grid for output results, if any:

      read (lunctl, *) filename

      interpolate_to_grid = filename(1:4) /= 'none'

      if (interpolate_to_grid) then

         l = len_trim (filename)
         open (lunin, file=filename(1:l), status='old', iostat=ios)
         if (ios /= 0) then
            write (lunlog, '(/, 2a)') &
               ' Unable to open CFD grid file ', filename(1:l)
            go to 99
         end if

         call xyzq_read (lunin, -lunin, true, nblocks, nf, cell_centered,      &
                         cfd_grid, ios)

         if (ios /= 0) then
            write (lunlog, '(/, 2a)') &
               ' Trouble reading CFD grid file ', filename(1:l)
            go to 99
         end if

         close (lunin)

      else
         open (lunoutg, file='uniform.g', status='unknown')
      end if

      read (lunctl, *) ibc       ! DPLR input profile BC code
      read (lunctl, *) filename  ! Output flow results (PLOT3D function file)
      l = len_trim (filename)

      if (fill_mode /= 'Q') open (lunoutf, file=filename(1:l), status='unknown')

!     Option to save some editing of the eventual pointwise BC file:

      pbca = filename(l-3:l) == 'pbca'

      read (lunctl, *) filename  ! Plottable results (Tecplot, uniform grid)

      open (lunplot, file=filename, status='unknown')

      read (lunctl, *) rectangle_power  ! For smoothing the stretched spokes

      do i = 1, 3
         read (lunctl, *) xyz_switch(i,:)
      end do

      close (lunctl)

   99 return

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

      subroutine expand_axis_profiles (ni, nj, xu, yu, v1v2)

!     Set up the prescribed two variables on the uniform grid as follows:
!
!     (a) Rectangular nozzle procedures (and early semi-elliptic procedures):
!
!       > Set up the profiles defined for the major (horizontal) and minor
!         (vertical) axes of the nozzle.  The definitions are probably analytic,
!         but a discrete dataset (possibly from arc-jet measurements) is OK too.
!       > Derive outer boundary values (probably constant) via linear interpola-
!         tion vs. arc length.
!       > Fill the interior function values via the indicated interpolation
!         method.  (The interior grid coordinates are already given.)
!
!     (b) Simpler semi-elliptic nozzle scheme (non-polar grid):
!         See subroutine expand_ellipse_profiles.
!
!     The uniform grid is passed as arguments to simplify the nomenclature.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)  :: ni, nj               ! Uniform grid dimensions
      real,    intent (in)  :: xu(ni,nj), yu(ni,nj) ! Uniform grid coordinates
      real,    intent (out) :: v1v2(2,ni,nj)        ! Desired 2-D distribution

!     Local variables:

      integer :: i, n
      real    :: r, rm1
      real, allocatable, dimension (:)   :: r_data, v_data, r_eval, v_eval
      real, allocatable, dimension (:,:) :: xtfi, ytfi

      logical :: analytic

!     Execution:

!     Major axis, v1:
!     ---------------

      n = ni

      allocate (r_eval(n), v_eval(n))

      r_eval(:) = xu(:,1)

      call discretize_shape (n, r_eval, v1_shape_h, v1_center_h, v1_wall_h,    &
                             v1_width_h, v1_steepness_h, v_eval, analytic)
      if (ios /= 0) go to 99  ! Bad shape spec.

      if (.not. analytic) call redistrib_1d &
         (n_v1_shape_h, v1_shape_r_h, v1_shape_f_h, n, r_eval, v_eval)

      v1v2(1,:,1) = v_eval(:)

!     Major axis, v2:
!     ---------------

      call discretize_shape (n, r_eval, v2_shape_h, v2_center_h, v2_wall_h,    &
                             v2_width_h, v2_steepness_h, v_eval, analytic)
      if (ios /= 0) go to 99  ! Bad shape spec.

      if (.not. analytic) call redistrib_1d &
         (n_v2_shape_h, v2_shape_r_h, v2_shape_f_h, n, r_eval, v_eval)

      v1v2(2,:,1) = v_eval(:)

      deallocate (r_eval, v_eval)

!     Minor axis and outer boundary are nozzle-specific:
!     --------------------------------------------------

      ios = 0

      select case (nozzle_type)

         case ('R')  ! Rectangle (plain)

!           Minor axis, v1:
!           ---------------

            n = nj

            allocate (r_eval(n), v_eval(n))

            r_eval(:) = yu(1,:)

            call discretize_shape (n, r_eval, v1_shape_v, v1_center_v,         &
                                   v1_wall_v, v1_width_v, v1_steepness_v,      &
                                   v_eval, analytic)
            if (ios /= 0) go to 99  ! Bad shape spec.

            if (.not. analytic) call redistrib_1d &
               (n_v1_shape_v, v1_shape_r_v, v1_shape_f_v, n, r_eval, v_eval)

            v1v2(1,1,:) = v_eval(:)

!           Minor axis, v2:
!           ---------------

            call discretize_shape (n, r_eval, v2_shape_v, v2_center_v,         &
                                   v2_wall_v, v2_width_v, v2_steepness_v,      &
                                   v_eval, analytic)
            if (ios /= 0) go to 99  ! Bad shape spec.

            if (.not. analytic) call redistrib_1d &
               (n_v2_shape_v, v2_shape_r_v, v2_shape_f_v, n, r_eval, v_eval)

            v1v2(2,1,:) = v_eval(:)

            deallocate (r_eval, v_eval)

!           Outer boundary arc lengths (j = nj and i = ni edges):
!           -----------------------------------------------------

            n = ni + nj - 1;  allocate (r_eval(n))

            r_eval(1) = zero
            do i = 2, ni
               r_eval(i) = r_eval(i-1) + &
                  sqrt ((xu(i,nj) - xu(i-1,nj))**2 + (yu(i,nj) - yu(i-1,nj))**2)
            end do

            j = nj
            do i = ni + 1, n
               j = j - 1
               r_eval(i) = r_eval(i-1) + &
                  sqrt ((xu(ni,j+1) - xu(ni,j))**2 + (yu(ni,j+1) - yu(ni,j))**2)
            end do

!           Arc-length-based linear interpolation between axis end points:

            r = r_eval(n);  r_eval(:) = r_eval(:) / r

            do i = 2, ni
               r = r_eval(i);  rm1 = one - r
               v1v2(:,i,nj) = (one - r) * v1v2(:,1,nj) + r * v1v2(:,ni,1)
            end do

            j = nj
            do i = ni + 1, n - 1
               j = j - 1
               r = r_eval(i);  rm1 = one - r
               v1v2(:,ni,j) = (one - r) * v1v2(:,1,nj) + r * v1v2(:,ni,1)
            end do

         case ('C', 'E', 'P')  ! Polar grids; treat edge i = 1 as collapsed

!           Minor axis, v1:
!           ---------------

            n = ni

            allocate (r_eval(n), v_eval(n))

            r_eval(:) = yu(:,nj)

            call discretize_shape (n, r_eval, v1_shape_v, v1_center_v,         &
                                   v1_wall_v, v1_width_v, v1_steepness_v,      &
                                   v_eval, analytic)
            if (ios /= 0) go to 99  ! Bad shape spec.

            if (.not. analytic) call redistrib_1d &
               (n_v1_shape_v, v1_shape_r_v, v1_shape_f_v, n, r_eval, v_eval)

            v1v2(1,:,nj) = v_eval(:)

!           Minor axis, v2:
!           ---------------

            call discretize_shape (n, r_eval, v2_shape_v, v2_center_v,         &
                                   v2_wall_v, v2_width_v, v2_steepness_v,      &
                                   v_eval, analytic)
            if (ios /= 0) go to 99  ! Bad shape spec.

            if (.not. analytic) call redistrib_1d &
               (n_v2_shape_v, v2_shape_r_v, v2_shape_f_v, n, r_eval, v_eval)

            v1v2(2,:,nj) = v_eval(:)

            deallocate (r_eval, v_eval)

!           Outer boundary arc lengths (i = ni edge only):
!           ----------------------------------------------

            v1v2(1,1,:) = v1v2(1,1,1)  ! Collapsed edge
            v1v2(2,1,:) = v1v2(2,1,1)

            n = nj;  allocate (r_eval(n))

            r_eval(1) = zero

            do j = 2, nj
               r_eval(j) = r_eval(j-1) + &
                  sqrt ((xu(ni,j) - xu(ni,j-1))**2 + (yu(ni,j) - yu(ni,j-1))**2)
            end do

!           Arc-length-based linear interpolation between axis end points:

            r = r_eval(n);  r_eval(:) = r_eval(:) / r

            do j = 2, nj - 1
               r = r_eval(j);  rm1 = one - r
               v1v2(:,ni,j) = (one - r) * v1v2(:,ni,1) + r * v1v2(:,ni,nj)
            end do

         case default

            write (lunlog, '(/, 2a)') 'Unknown nozzle type: ', nozzle_type,    &
               'Options: ', &
               'Rectangle, Polar [form of rectangle], Circular & Semi-Ellipse.'

            ios = 1;  go to 99

      end select

      deallocate (r_eval)

!     Fill the interior distribution of v1, v2:
!     -----------------------------------------

      select case (interp_type)

         case ('T')  ! TFI

            allocate (xtfi(ni,nj), ytfi(ni,nj))

            xtfi(:,:) = xu(:,:)  ! Since TFIF2D fills in (x,y) as well as f.
            ytfi(:,:) = yu(:,:)

            call tfif2d (2, ni, 1, ni, 1, nj, xtfi, ytfi, v1v2)

            deallocate (xtfi, ytfi)

         case ('C')  ! Cartesian product

            call cartesian_product (nozzle_type, ni, nj, 2, v1v2)

      end select

!     Use fill_mode = 'q[uick-look]' to look at the results promptly.

 99   return

      end subroutine expand_axis_profiles

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine expand_ellipse_profiles (ni, nj, xu, yu, v1v2)

!     This variant did not fit very well into the earlier expand_axis_profiles.
!     It assumes the specialized grid of nozzle_type 'S' or 'C', which has the
!     i = ni edge collapsed.
!     At each i-station, just scale the minor axis profile to the range implied
!     by the major profile.
!     The arguments simplify the nomenclature.  Other variables (in particular
!     jmajor) are shared via the host program.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)  :: ni, nj               ! Uniform grid dimensions
      real,    intent (in)  :: xu(ni,nj), yu(ni,nj) ! Uniform grid coordinates
      real,    intent (out) :: v1v2(2,ni,nj)        ! Desired 2-D distribution

!     Local variables:

      integer :: i, n
      logical :: analytic

      real, allocatable, dimension (:) :: r_data, v_data, r_eval, v_eval

!     Execution:

!     Minor axis, v1:
!     ---------------

      n = nj

      allocate (r_eval(n), v_eval(n))

      r_eval(:) = yu(1,:)

      call discretize_shape (n, r_eval, v1_shape_v, v1_center_v,               &
                             v1_wall_v, v1_width_v, v1_steepness_v,            &
                             v_eval, analytic)
      if (ios /= 0) go to 99  ! Bad shape spec.

      if (.not. analytic) call redistrib_1d (n_v1_shape_v, v1_shape_r_v,       &
                                             v1_shape_f_v, n, r_eval, v_eval)
      v1v2(1,1,:) = v_eval(:)

!     Minor axis, v2:
!     ---------------

      call discretize_shape (n, r_eval, v2_shape_v, v2_center_v,               &
                             v2_wall_v, v2_width_v, v2_steepness_v,            &
                             v_eval, analytic)
      if (ios /= 0) go to 99  ! Bad shape spec.

      if (.not. analytic) call redistrib_1d (n_v2_shape_v, v2_shape_r_v,       &
                                             v2_shape_f_v, n, r_eval, v_eval)
      v1v2(2,1,:) = v_eval(:)

      deallocate (r_eval, v_eval)

!     Major axis, v1:
!     ---------------

!     Semi-ellipse has a wall along the major axis:

      jmajor = (nj + 1) / 2

      n = ni

      allocate (r_eval(n), v_eval(n))

      r_eval(:) = xu(:,1)

      call discretize_shape (n, r_eval, v1_shape_h, v1_center_h, v1_wall_h,    &
                             v1_width_h, v1_steepness_h, v_eval, analytic)
      if (ios /= 0) go to 99  ! Bad shape spec.

      if (.not. analytic) call redistrib_1d (n_v1_shape_h, v1_shape_r_h,       &
                                             v1_shape_f_h, n, r_eval, v_eval)

      v1v2(1,:,jmajor) = v_eval(:)

!     Major axis, v2:
!     ---------------

      call discretize_shape (n, r_eval, v2_shape_h, v2_center_h, v2_wall_h,    &
                             v2_width_h, v2_steepness_h, v_eval, analytic)
      if (ios /= 0) go to 99  ! Bad shape spec.

      if (.not. analytic) call redistrib_1d &
         (n_v2_shape_h, v2_shape_r_h, v2_shape_f_h, n, r_eval, v_eval)

      v1v2(2,:,jmajor) = v_eval(:)

      deallocate (r_eval, v_eval)

!     Fill the bulk of the nozzle flow by scaling the minor profile:
!     --------------------------------------------------------------

      call scale_minor_profile (ni, nj, 2, jmajor, v1v2)

   99 return

      end subroutine expand_ellipse_profiles

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine extend_spokes ()

!     Turn the interim polar elliptic grid into the final rectangular grid by
!     extending the spokes and redistributing.  Retain a portion of the ellipse
!     grid to round off the diagonal region a little.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, j
      real    :: fe, fr, slope, xr, yr

      do j = 2, nj_uniform
         slope = (uniform_grid(1)%y(ni_uniform,j,1) - r_center_v) /            &
                 (uniform_grid(1)%x(ni_uniform,j,1) - r_center_h)
         do i = 2, ni_uniform
            xr = uniform_grid(1)%x(i,1,1)
            yr = (xr - r_center_h) * slope + r_center_v
            fr = (real (i - 1) / real (ni_uniform - 1)) ** rectangle_power
            fe = one - fr
            uniform_grid(1)%x(i,j,1) = fr * xr + fe * uniform_grid(1)%x(i,j,1)
            uniform_grid(1)%y(i,j,1) = fr * yr + fe * uniform_grid(1)%y(i,j,1)
         end do
      end do

      uniform_grid(1)%y(ni_uniform,nj_uniform,1) = r_wall_v

      do j = nj_uniform + 1, nj_ellipse - 1
         slope = (uniform_grid(1)%y(ni_uniform,j,1) - r_center_v) /            &
                 (uniform_grid(1)%x(ni_uniform,j,1) - r_center_h)
         do i = 1, ni_uniform
            yr = uniform_grid(1)%y(i,nj_ellipse,1)
            xr = (yr - r_center_v) / slope + r_center_h
            fr = (real (i - 1) / real (ni_uniform - 1)) ** rectangle_power
            fe = one - fr
            uniform_grid(1)%x(i,j,1) = fr * xr + fe * uniform_grid(1)%x(i,j,1)
            uniform_grid(1)%y(i,j,1) = fr * yr + fe * uniform_grid(1)%y(i,j,1)
         end do
      end do

      end subroutine extend_spokes

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine flow_conditions ()

!     Calculate the flow conditions at each point of the grid by the indicated
!     method.  This routine initially formed the bulk of the main program, but
!     it was modularized as part of the option to iterate to achieve specified
!     bulk flow rate and/or bulk enthalpy (Ht_MF case only).
!
!     Note that the starting guess for point (1,1) is set at the higher level.
!     On the first call, point (i-1,j) results are the starting guesses for
!     point (i,j), while point (1,j-1) starts point (1,j).  On further calls
!     to equilibrium_composition, the previous iteration starts point (1,1) but
!     the rest of the points are done the same way as for the first call.
!
!     Economy mode means only the axis points are fully calculated: interior
!     points are derived by interpolation followed by ensuring that the equi-
!     librium composition is consistent with the results for P and T.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: iend, jend
      logical :: scaling

!     Other variables are inherited from the calling program.

!     Execution:

      iend = ni_uniform;  jend = nj_ellipse  ! Unless it's economy mode

!     For each uniform grid point ...

      do j = 1, nj_ellipse

         do i = 1, ni_uniform

!           Skip interior points if economy mode is specified:

            select case (fill_mode)

               case ('E') ! Economy mode: just do the axes

                  select case (nozzle_type)

                     case ('R')  ! Rectangular grid (plain)

                        skip = i > 1 .and. j > 1
                        iend = 1   ! jend = default

                     case ('C', 'E', 'P')  ! Polar grid

                        skip = j > 1 .and. j < nj_ellipse
                        jend = nj_ellipse  ! iend = default

                     case ('S')  ! Semi-ellipse via scaling

                        skip = i > 1 .and. j /= jmajor
                        iend = 1;  jend = jmajor

                  end select

               case default ! Do all interior points

                  skip = false

            end select

            if (skip) cycle

            if (diagnostics) then
               write (lunlog, '(/, a, 2i4)') ' i, j:', i, j
            else
               write (lunlog, '(   a, 2i4)') ' i, j:', i, j
            end if

            v1 = uniform_v1_v2(1,i,j)
            v2 = uniform_v1_v2(2,i,j)

            call mass_mole_conversion (1)  ! Mass fractions to mole fractions

            select case (ndim)

            case (0) ! Straight call to Gary's package

               if (i == 1    .and. j == 1     .or.  &
                   i == iend .and. j == jend) then  ! Force printing the details
                  status_CEA = 1
               else
                  status_CEA = 0
               end if

               call control_outputs ()

               call equilibrium_composition (specified_state, v1, v2,          &
                                             n_species, name_gary,             &
                                             mole_fractions, P, rho, T, h, s,  &
                                             ae, status_CEA)
               if (status_CEA /= 0) then
                  write (lunlog, '(a, i4, a, 1p, 2e13.6)') &
                     ' Bad equilibrium gas composition status:', status_CEA,   &
                     '  v1, v2:', v1, v2
                  write (lunlog, '(3a, 2i4)') &
                  ' Specified state case: ', specified_state, '    Point:', i, j
                  ier = -1
                  go to 99
               end if

               call mass_mole_conversion (2)  ! Mole fractions to mass fractions

            case (1) ! 2 separable eqns. via nested non-deriv. zero finders

!!!!           This option from NOZZLE_THROAT_CONDITIONS has been eliminated.

            case (2) ! 2 nonlinear equations; Newton iteration w/ central diffs.

!              The iteration involves solving  J dx = f, where J is the matrix
!              of partial derivatives df(i)/dx(j), updating x as x - dx, and
!              repeating until || f || becomes arbitrarily small.

!              Safeguarding involves halving each pure Newton step until || f ||
!              is reduced from its value for the previous iteration, meaning the
!              update step is actually x <-- x - alpha dx, where alpha <= 1.

!              Tests for convergence involve a choice of tolerance, which is
!              dependent upon the data scaling and precision.  The tolerance
!              is used relative to the norm of the starting guess in SI units.

!              N.B. Don't confuse the "v"s here (P, T) w/ the specified v1, v2.

               ier    = 0
               iter   = 0
               x(1)   = P         ! From starting_guess for point (1,1);
               x(2)   = T         ! from nearby pt. for all other (i,j)
               dx(:)  = zero      ! For iteration 0 printout
               alpha  = zero      !  "      "     "     "
               fnorm0 = 1.E+30    ! Big to avoid iteration 0 test
               fscale(1) = v1     ! Since the equations have these units
               fscale(2) = v2
               xscale(1) = x(1);  x(1) = one
               xscale(2) = x(2);  x(2) = one
               aitch(:)  = (one + x(:)) * cube_root_eps

               status_CEA = 0     ! Suppress printout unless requested
               call control_outputs ()

               if (diagnostics) then
                  write (lunlog, '(a, 1p, 2e15.6)') &
                     ' Current v1, v2:                         ', v1, v2
                  write (lunlog, '(a, 1p, 2e15.6)') &
                     ' Initial variables P, T, unscaled:       ', P, T
                  write (lunlog, '(a, 1p, 2e15.6)') &
                     ' Initial variables x1, x2, scaled:       ', x(:)
                  write (lunlog, '(a, 1p, 2e15.6)') &
                     ' Central difference steps:               ', aitch
                  write (lunlog, '(a, 1p, e15.6)') &
                     ' Convergence tolerance on residual norm: ', toler
               end if

  30           continue

!                 Evaluate the residuals for the current variables.
!                 Actually, to calculate central derivatives while we're at it,
!                 we evaluate the functions for each x at x - h, x, x + h.

                  mass_fraction_0(:) = mass_fraction(:)  ! Start consistently

                  do l = 1, nv

                     xfd(:) = x(:)

                     do k = -1, 1

                        xfd(l) = x(l) + real (k) * aitch(l)
                        mass_fraction(:) = mass_fraction_0(:)

                        status_CEA = 0      ! Suppress printout unless requested
                        call control_outputs ()

                        call residual_2d (xfd(1)*xscale(1), xfd(2)*xscale(2),  &
                                          ffd(k,1,l), ffd(k,2,l))

                        ffd(k,1,l) = ffd(k,1,l) / fscale(1)
                        ffd(k,2,l) = ffd(k,2,l) / fscale(2)

                     end do

!                    Partial derivatives dfk / dxl:

                     do k = 1, nv
                       AJ(k,l) = (ffd(1,k,l) - ffd(-1,k,l)) / (two * aitch(l))
                     end do

                  end do

!                 Set up the residual vector f of function values as the RHS of
!                 the system J dx = f.  The norm of f should converge to ~0.

                  f(1) = ffd(0,1,1)
                  f(2) = ffd(0,2,2)

                  fnorm = max (abs (f(1)), abs (f(2)))
                  mass_fraction(:) = mass_fraction_0(:)

   40             continue

                     status_CEA = 0        ! Suppress printout unless requested
                     call control_outputs ()

!                    Halve the step till ||f|| is reduced (except for 1st time).

                     if (fnorm > fnorm0) then
                        if (alpha > 0.001) then
                           alpha = half * alpha
                           do l = 1, nv
                              x(l) = xlast(l) - alpha * dx(l)
                           end do
                           mass_fraction(:) = mass_fraction_0(:)

                           call residual_2d (x(1)*xscale(1), x(2)*xscale(2),   &
                                             f(1), f(2))

                           f(:) = f(:) / fscale(:)
                           fnorm = max (abs (f(1)), abs (f(2)))
                           go to 40
                        end if
                        go to 83
                     end if

                  if (diagnostics) then
                     dx(:) = alpha * dx(:)
                     write (lunprint, &
                        '(a, i3, a, 1p, e9.2, a, e9.2, 2(a, 2e16.8))')         &
                        ' Itn.', iter, '  ||f||:', fnorm, '  step:', alpha,    &
                        '  dx:', dx, '  x:', x(:)
                  end if

                  if (fnorm < toler .and. iter > 0) go to 90  ! Success


                  iter = iter + 1
                  if (iter == itmax) go to 81

                  fnorm0 = fnorm

!                 Set up LHS matrix for next iteration.  AJ(i,j) = dfi/dxj.

!!!               AJ(1, 1) = ...  ! Already done during finite differencing
!!!               AJ(2, 1) = ...
!!!               AJ(1, 2) = ...
!!!               AJ(2, 2) = ...

!                 Solve  J dx = f.

                  call lusolve (nv, nv, AJ, f, ier)  ! Solution overwrites RHS

                  if (ier /= 0) go to 84

                  do l = 1, nv
                     dx(l)    = f(l)
                     xlast(l) = x(l)
                     x(l)     = x(l) - dx(l)
                  end do

                  alpha = one

                  go to 30                           ! Do another iteration


!              Error handling.

  81           ier = 1
               write (lunlog, '(/, a)') ' Iteration limit reached.'
               if (fnorm < 10. * toler) go to 90

               ier = 2
!!!            go to 99   ! Bad data scaling can stop the convergence test from
               go to 90   ! being satisfied even when the solution is adequate.

  83           ier = 3
               write (lunlog, '(/, a)') ' Step-halving iteration failed.'
               go to 90

  84           ier = 4
               write (lunlog, '(/, a)') ' Singular matrix.'  ! Show soln. anyway
!!!            go to 99

  90           continue    ! Wrap up, having converged or almost converged.

!              Ensure final values:

               mass_fraction(:) = mass_fraction_0(:)

               if (i == 1    .and. j == 1     .or.  &
                   i == iend .and. j == jend) then  ! Force printing the details
                  status_CEA = 1
               else
                  status_CEA = 0
               end if

               call control_outputs ()

               call residual_2d (x(1)*xscale(1), x(2)*xscale(2), f(1), f(2))

               f(:) = f(:) / fscale(:)

               write (lunlog, &
                 '(a, 1p, 2e14.6, a, 2e14.6, a, 2e14.6)') &
                 ' v1, v2:', v1, v2, '   x1, x2:', x(:), '   Residuals:', f(:)

               x(:) = x(:) * xscale(:)  ! For starting the next case

            end select  ! Over main options for specified 2 vars. & soln. method

            mass_fraction_0(:) = mass_fraction(:)

!           Store the flow variables indicated by ibc for this point.
!           Internal units are also stored in fplot(:,i,j) for plotting.

            call store_results (specified_state, ibc, i, j,                    &
                                ni_uniform, nj_ellipse, nf, uniform_grid(1)%q, &
                                nfplot, fplot)
            P_T(1,i,j) = P
            P_T(2,i,j) = T

            if (i == 1) then ! Save values to start the next row
               P0    = P
               T0    = T
               gamf0 = gamf
               mass_fraction_1(:) = mass_fraction(:)
            end if

         end do ! Next i

         P    = P0     ! Starting guesses for the first point of the next row
         T    = T0
         gamf = gamf0
         mass_fraction(:) = mass_fraction_1(:)

      end do ! Next j


!     Economy mode requires interpolating outer boundary and interior values:

      j1 = 2  ! Except ...

      scaling = nozzle_type == 'S';  if (scaling) j1 = 1  ! Semi-ellipse case

      select case (fill_mode)

         case ('E')  ! Economy mode

!           Interpolate P and T (only):

            if (scaling) then  ! Scale minor axis profile at each i

               call scale_minor_profile (ni_uniform, nj_uniform, 2, jmajor, P_T)

            else  ! TFI, Cartesian prod., or linear combination of axis profiles

               call interp_interior (ni_uniform, nj_ellipse, 2,                &
                                     uniform_grid(1)%x, uniform_grid(1)%y, P_T)
            end if

!           Ensure true equilibrium compositions by invoking the 'P_T' option
!           for all the interpolated points:

            write (lunlog, '(/, a)') ' Adjusting interpolated points.'

            k = 2 + n_species

            mass_fraction(:) = fplot(3:k,1,1)

            do j = j1, nj_ellipse  ! Both nozzle types skip j = 1 except 'S'

               do i = 1, ni_uniform

                  select case (nozzle_type)

                     case ('R')  ! Rectangular grid (plain)

                        skip = i == 1

                     case ('C', 'E', 'P')  ! Polar grid

                        skip = j == nj_ellipse

                     case ('S')  ! Semi-ellipse via scaling

                        skip = i == 1 .or. j == jmajor

                  end select

                  if (skip) cycle

                  v1 = P_T(1,i,j)
                  v2 = P_T(2,i,j)

                  call equilibrium_composition (main_options(5), v1, v2,       &
                                                n_species, name_gary,          &
                                                mole_fractions, P, rho, T,     &
                                                h, s, ae, status_CEA)
                  if (status_CEA /= 0) then
                     write (lunlog, '(a, i4, a, 1p, 2e13.6)') &
                        ' Bad equilibrium composition status::', status_CEA,   &
                        '  v1, v2:', v1, v2
                     write (lunlog, '(a, 2i4)') &
                        ' Specified state case: P_T, economy mode point:', i, j
                     ier = -2
                     go to 99
                  end if

                  call mass_mole_conversion (2)  ! Back to mass fractions

!                 Note the dual use of the specified_state argument here:

                  call store_results ('Econo', ibc, i, j,                      &
                                      ni_uniform, nj_ellipse, &
                                      nf, uniform_grid(1)%q, nfplot, fplot)

                  if (i == 1) mass_fraction_1(:) = mass_fraction(:)

               end do ! Next i

               mass_fraction(:) = mass_fraction_1(:)

            end do ! Next j

         case default

      end select

 99   return

      end subroutine flow_conditions

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
         c_over_M  = mass_fraction(i) / Mol(i)

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

      if (show_details) then
         write (lunlog, '(a, 1p, e14.6)') '   gamf: ',  gamf
         write (lunlog, '(a, 1p, e14.6)') '  R_bar: ',  R_bar
         write (lunlog, '(a, 1p, e14.6)') ' Cp_bar: ', Cp_bar
         write (lunlog, '(a, 1p, e14.6)') ' Cv_bar: ', Cv_bar
         write (lunlog, '(a, 1p, e14.6)') ' ev_bar: ', ev_bar
      end if

      end subroutine frozen_gas_constants

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine integrate_bulk_mass_and_enthalpy (ni, nj, nf, x, y, f,        &
                                                   bulk_mass, bulk_enthalpy)

!     Integrate bulk mass flow and bulk enthalpy for a nozzle cross section.
!     The first two functions are assumed to be total stagnation enthalpy and
!     mass flow rate per unit area, respectively.  The bulk enthalpy is the
!     integral of rho u ho over the integral of rho u (w.r.t. area).

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)  :: ni, nj, nf
      real,    intent (in),    dimension     (ni, nj) :: x, y
      real,    intent (inout), dimension (nf, ni, nj) :: f  ! Temporarily chngd.
      real,    intent (out) :: bulk_mass, bulk_enthalpy

!     Local variables:

      real :: area
      real, allocatable :: ftemp(:,:), quadratures(:)

!     Execution:

      allocate (quadratures(nf), ftemp(ni,nj))  ! Need to integrate rho u ho

      ftemp(:,:) = f(1,:,:)
      f(1,:,:)   = f(1,:,:) * f(2,:,:)

      call trapezoidal_quadrature_2d (ni, nj, nf, x, y, f, area, quadratures)

      bulk_mass     = quadratures(2)
      bulk_enthalpy = quadratures(1) / bulk_mass  ! I rho u ho dA / I rho u dA

!!!   bulk_enthalpy = quadrature_multiplier * bulk_enthalpy  ! No
      bulk_mass     = quadrature_multiplier * bulk_mass  ! x 4 or 2
      area          = quadrature_multiplier * area

      write (lunlog, '(/, (a, 1p, e15.7, a))') &
         ' Approximate nozzle throat area: ', area,          ' m^2',  &
         ' Approximate bulk enthalpy:      ', bulk_enthalpy, ' J/kg', &
         ' Approximate bulk mass flow rate:', bulk_mass,     ' kg/s'

      f(1,:,:) = ftemp(:,:)

      deallocate (quadratures, ftemp)

      end subroutine integrate_bulk_mass_and_enthalpy

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine interp_interior (ni, nj, nf, xu, yu, f)

!     Fill interior uniform grid flow results by the indicated method.
!     TFI requires initializing the outer edge(s) first.
!     Nozzle-type 'S' (distinct from the earlier 'E') is done more simply in
!     subroutine scale_minor_profile along with case 'C'.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)    :: ni, nj, nf
      real,    intent (in)    :: xu(ni,nj), yu(ni,nj)
      real,    intent (inout) :: f(nf,ni,nj)

!     Local variables:

      integer :: i, j, n
      real    :: r, rm1
      real, allocatable :: s(:), xtfi(:,:), ytfi(:,:)

!     Execution:

!     Avoid a lot of indenting by doing the outer edge even for the Cartesian
!     product case.

!     Also, the semi-ellipse case is handled specially by a scheme that doesn't
!     fit very well with the earlier options.

!     Follow the outer boundary procedure in subroutine expand_axis_profiles:

      select case (nozzle_type)

         case ('R')  ! Rectangle (plain)

!           Outer boundary arc lengths (j = nj and i = ni edges):

            n = ni + nj - 1;  allocate (s(n))

            s(1) = zero
            do i = 2, ni
               s(i) = s(i-1) + &
                  sqrt ((xu(i,nj) - xu(i-1,nj))**2 + (yu(i,nj) - yu(i-1,nj))**2)
            end do

            j = nj
            do i = ni + 1, n
               j = j - 1
               s(i) = s(i-1) + &
                  sqrt ((xu(ni,j+1) - xu(ni,j))**2 + (yu(ni,j+1) - yu(ni,j))**2)
            end do

!           Arc-length-based linear interpolation between axis end points:

            r = s(n);  s(:) = s(:) / r

            do i = 2, ni
               r = s(i);  rm1 = one - r
               f(:,i,nj) = (one - r) * f(:,1,nj) + r * f(:,ni,1)
            end do

            j = nj
            do i = ni + 1, n - 1
               j = j - 1
               r = s(i);  rm1 = one - r
               f(:,ni,j) = (one - r) * f(:,1,nj) + r * f(:,ni,1)
            end do

         case ('C', 'E', 'P')  ! Polar grid

            do j = 1, nj
               f(:,1,j) = f(:,1,1)  ! Collapsed edge
            end do

!           Outer boundary arc lengths (i = ni edge only):

            n = nj;  allocate (s(n))

            s(1) = zero

            do j = 2, nj
               s(j) = s(j-1) + &
                  sqrt ((xu(ni,j) - xu(ni,j-1))**2 + (yu(ni,j) - yu(ni,j-1))**2)
            end do

!           Arc-length-based linear interpolation between axis end points:

            r = s(n);  s(:) = s(:) / r

            do j = 2, nj - 1
               r = s(j);  rm1 = one - r
               f(:,ni,j) = (one - r) * f(:,ni,1) + r * f(:,ni,nj)
            end do

      end select

      deallocate (s)

!     Fill the interior calculated flow:

      select case (interp_type)

         case ('T')  ! TFI

            allocate (xtfi(ni,nj), ytfi(ni,nj))

            xtfi(:,:) = xu(:,:)  ! Since TFIF2D fills in (x,y) as well as f.
            ytfi(:,:) = yu(:,:)

            call tfif2d (nf, ni, 1, ni, 1, nj, xtfi, ytfi, f)

            deallocate (xtfi, ytfi)

         case ('C')  ! Cartesian product

            call cartesian_product (nozzle_type, ni, nj, nf, f)

      end select

   99 return

      end subroutine interp_interior

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
               mole_fractions(i) = mass_fraction(i) / Mol(i)
               sum = sum + mole_fractions(i)
            end do

            Mbar = one / sum  ! Mixture molecular weight

            mole_fractions(:) = Mbar * mole_fractions(:)

         case (2) ! Mole fractions to mass fractions

            do i = 1, n_species
               mass_fraction(i) = mole_fractions(i) * Mol(i)
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

      subroutine redistrib_1d (n, r, v, neval, reval, veval)

!     Redistribute a 1-D (r,v) dataset using a monotonic local cubic spline.
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

      end subroutine redistrib_1d

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine redistrib_2d ()

!     Redistribute the uniform surface dataset to the indicated CFD grid.
!     The new coordinates are assumed to be contained within the data range.
!     The interpolation works with (x,y,z) data, but the uniform working grid
!     z coordinates are zero, and the [permuted] CFD grid z coordinates should
!     be at least constant if not zero.  (The searching really projects each
!     target point onto the uniform grid at z = 0.)

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use adt_utilities

!     Local constants:

      integer :: ndebug = 1, idebug = 2, jdebug = 2

!     Local variables:

      integer :: i, ib, ic, ier, iquad, j, jc, m, n, nb_search, ni, ninside,   &
                 nj, noutside, npts, nquad
      real    :: dmax, dmean, dmin, dsqmin, dtolsq, p, pm1, q, qm1,            &
                 interp_xyz(3), target_xyz(3)

      integer, allocatable :: conn(:,:)
      real,    allocatable :: interp_f(:)

!     Execution:

!     Tolerance for search diagnostics (refine later):

      dtolsq = (0.001) ** 2

      nquad = uniform_grid(1)%ni * uniform_grid(1)%nj

      allocate (conn(3,nquad))  ! For patch # and cell (i,j)

      nb_search = 1

      call build_adt (nb_search, uniform_grid, nquad, conn)

      allocate (interp_f(nf + nfplot), cfd_plot(nblocks))

      do ib = 1, nblocks

         ninside  = 0;  dmax  = zero
         noutside = 0;  dmean = zero
         npts     = 0

         ni = cfd_grid(ib)%ni;  cfd_grid(ib)%mi = ni;  cfd_plot(ib)%mi = ni
         nj = cfd_grid(ib)%nj;  cfd_grid(ib)%mj = nj;  cfd_plot(ib)%mj = nj
                                cfd_grid(ib)%mk = 1;   cfd_plot(ib)%mk = 1

         call q_allocate (cfd_grid(ib), nf, ios)

         call q_allocate (cfd_plot(ib), nfplot, ios)

         do j = 1, nj

            do i = 1, ni

!              Convert from the flow solver coordinates to internal coordinates:

               interp_xyz(1) = cfd_grid(ib)%x(i,j,1) ! Need them as a vector
               interp_xyz(2) = cfd_grid(ib)%y(i,j,1)
               interp_xyz(3) = cfd_grid(ib)%z(i,j,1)

               target_xyz(:) = matmul (xyz_switch, interp_xyz)

               cfd_grid(ib)%x(i,j,1) = target_xyz(1) ! For plottable results
               cfd_grid(ib)%y(i,j,1) = target_xyz(2) ! Plots are 2-D (x,y only)
!!!            cfd_grid(ib)%z(i,j,1) = target_xyz(3) ! Assumed constant

               call search_adt (target_xyz, iquad, p, q, dsqmin, true,         &
                               nb_search, uniform_grid, nquad, conn, interp_xyz)

               if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
                  ninside  = ninside + 1
               else
                  noutside = noutside + 1
               end if
               npts = npts + 1

               dmax  = max (dmax, dsqmin)
               dmean = dmean + dsqmin

!              Interpolate the surface flow at this target point:

               n   = conn(1,iquad) ! Block #
               ic  = conn(2,iquad) ! Lower left quad. indices
               jc  = conn(3,iquad)
               pm1 = one - p
               qm1 = one - q

               interp_f(1:nf) = &
                  qm1 * (pm1 * uniform_grid(n)%q(1:nf,ic,  jc,  1)  + &
                           p * uniform_grid(n)%q(1:nf,ic+1,jc,  1)) + &
                    q * (pm1 * uniform_grid(n)%q(1:nf,ic,  jc+1,1)  + &
                           p * uniform_grid(n)%q(1:nf,ic+1,jc+1,1))

               interp_f(nf+1:nf+nfplot) = &
                  qm1 * (pm1 * fplot(:,ic,  jc  )  + &
                           p * fplot(:,ic+1,jc  )) + &
                    q * (pm1 * fplot(:,ic,  jc+1)  + &
                           p * fplot(:,ic+1,jc+1))

!              Debug output:

               if (ib == ndebug) then
                  if (i == idebug .and. j == jdebug) then
                     write (6, '(a, 3f20.8)') ' target xyz: ', target_xyz
                     write (6, '(a, 3f20.8)') ' interp xyz: ', interp_xyz
                     write (6, '(a, 3i5)') &
                        ' Search result: block, ic, jc: ', n, ic, jc
                     write (6, '(a, 2f20.10)') ' p, q: ', p, q
                     dmin = sqrt (dsqmin)
                     write (6, '(a, f20.10)') ' shortest distance: ', dmin
                     write (6, '(a)') ' Cell vertex values, f(:,i:i+1,j:j+1):'
                     write (6, '(i3, 4f20.8)') &
                        (m, uniform_grid(n)%q(m,ic:ic+1,jc:jc+1,1), m = 1, nf)
                     write (6, '(a)') ' Interpolated values:'
                     write (6, '(i3, f20.8)') (m, interp_f(m), m = 1, nf)
                  end if
               end if

               cfd_grid(ib)%q(:,i,j,1) = interp_f(1:nf)
               cfd_plot(ib)%q(:,i,j,1) = interp_f(nf+1:nf+nfplot)

            end do ! Next i for this target block

         end do ! Next j for this target block

         npts = max (npts, 1) ! In case entire patch was outside the tolerance

         write (lunlog, '(a, i5, a, 2i7, a, 1p, 2e12.5)')                      &
            '  Target patch', ib,                                              &
            ':  # points in/outside tolerance:', ninside, noutside,            &
            ';  max/mean devn.:', sqrt (dmax), sqrt (dmean / real (npts))

      end do ! Next target block

      deallocate (conn, interp_f)

      end subroutine redistrib_2d

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine residual_eq_1 (x, fx)

!     Evaluate the residual for equation (1) of the main option (Ht_MF).

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, intent (in)  :: x
      real, intent (out) :: fx

      integer :: k
      real    :: P_specified, T_specified

!     "x" = T in eq. (1)

!     f(1) = v1 - t1 Cpbar T - evbar - sum (ci hfi) = 0   where v1 = ho
!     f(2) = v2 - P sqrt (gamf / (Rbar T))          = 0     "   v2 = rho u

      P_specified = P  ! The pressure associated with the current T = x
      T_specified = x

!     Compute species densities, etc., for the current estimates of P, T:

      call mass_mole_conversion (1)  ! Mass fractions to mole fractions

      if (diagnostics) then
         write (lunlog, '(/, a, 1p, 2e14.6)') &
            ' Residual_eq_1:  P, T inputs to CEA: ', P_specified, T_specified
         status_CEA = 1  ! Force CEA printout
      else
         status_CEA = 0
      end if

      call equilibrium_composition ('P_T  ', P_specified, T_specified,         &
                                    n_species, name_gary, mole_fractions,      &
                                    P, rho, T, h, s, ae, status_CEA)
      if (diagnostics) &
         write (lunlog, '(a, 1p, 3e14.6)') ' Output rho, h, ae: ', rho, h, ae

      if (status_CEA /= 0) then
         write (lunlog, '(a, i4, a, 2i4, /, a, 1p, 2e13.6)')                   &
            ' Bad equilibrium composition status:', status_CEA,                &
            '  Point i,j:', i, j, ' P & T specified: ', P_specified, T_specified
!!!      go to 99
         stop   ! Hard to avoid, but better than keeping going
      end if

      call mass_mole_conversion (2)  ! Mole fractions to mass fractions

      if (diagnostics) &
         write (lunlog, '(a, /, 1p, 12e11.4)') ' Output mass fractions:',      &
            mass_fraction(:)

!     For this equilibrium composition, update the gas constants:

      call frozen_gas_constants (T_specified)

      t1 = one + half * (gamf - one) * fMsq
      fx = v1 - t1 * Cp_bar * T_specified - ev_bar - h_ref_sum

      if (diagnostics) write (lunlog, '(a, 1p, 3e14.6)') &
         ' Eq. (1) P, T, residual:   ', P, T_specified, fx

!! 99 return
      return

      end subroutine residual_eq_1

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine residual_eq_2 (x, fx)

!     Evaluate the residual for equation (2) of the main option (Ht_MF).

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, intent (in)  :: x
      real, intent (out) :: fx

!     "x" = P in eq. (2)

!     f(1) = v1 - t1 Cpbar T - evbar - sum (ci hfi) = 0   where v1 = ho
!     f(2) = v2 - P sqrt (gamf / (Rbar T))          = 0     "   v2 = rho u

      P    = x
      fx   = v2 - frozen_Mach * P * sqrt (gamf / (R_bar * T))

      if (diagnostics) write (lunlog, '(a, 1p, 3e14.6)') &
         ' Eq. (2) P, T, residual:   ', P, T, fx

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

      subroutine save_results ()

!     Write a 3-D multiblock PLOT3D function file on either the preliminary
!     uniform grid or the indicated CFD grid for DPLR, and also a Tecplot file
!     with the uniform grid results.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      character, parameter :: quotes * 1 = '"'

      integer :: i, ib, n, nb
      logical :: integrate

!     Execution:

      nb = 1  ! Can't be a constant because of xyz_header utility (read & write)

      if (.not. interpolate_to_grid) then

         call xyz_write (lunoutg, true, nb, uniform_grid, ios)

         if (ios /= 0) then
            write (lunlog, '(/, a)') ' Trouble writing uniform grid.'
         end if

         close (lunoutg)

         call q_write (lunoutf, true, nb, nf, uniform_grid, ios)

         if (ios /= 0) then
            write (lunlog, '(/, a)') ' Trouble writing uniform results.'
         end if

      else

         if (.not. pbca) then  ! PLOT3D format

            call q_write (lunoutf, true, nblocks, nf, cfd_grid, ios)

            if (ios /= 0) then
               write (lunlog, '(/, a)') ' Trouble writing flow on CFD grid.'
            end if

         else  ! Save some editing of the DPLR pointwise BC file:

            do ib = 1, nblocks
               write (lunoutf, '(a, /, a, i3, /, a)') &
                  '#', '# Inflow profile for zone', ib, '#'
               write (lunoutf, '(1p, 6e18.10)', iostat=ios) &
                  (cfd_grid(ib)%q(n,:,:,1), n = 1, nf)
               if (ios /= 0) then
                  write (lunlog, '(/, a, i3)') &
                     ' Trouble writing flow on CFD grid.  Block #:', ib
                  exit
               end if
            end do

         end if

      end if

      close (lunoutf)

      integrate = false

      write (lunplot, '(a)') &
         'TITLE = "Nozzle throat conditions"', &
         'VARIABLES = "x, m", "y, m"'

      select case (specified_state) ! For the name of v1

         case ('Ht_MF')
            write (lunplot, '(a)') '"Total enthalpy, J/kg"'
         case ('Rho_T', 'Rho_H', 'Rho_S')
            write (lunplot, '(a)') '"Density, kg/m^3"'
         case ('P_T  ', 'P_H  ', 'P_S  ')
            write (lunplot, '(a)') '"Pressure, Pa"'

      end select

      select case (specified_state) ! For the name of v2

         case ('Ht_MF')
            write (lunplot, '(a)') '"Mass flux, kg/(m^2.s)"'
            integrate = true
         case ('Rho_T', 'P_T  ')
            write (lunplot, '(a)') '"Temperature, K"'
         case ('Rho_H', 'P_H  ')
            write (lunplot, '(a)') '"Enthalpy, J/kg"'
         case ('Rho_S', 'P_S  ')
            write (lunplot, '(a)') '"Entropy, J/kg-mol"'

      end select

      write (lunplot, '(4a)') &        ! Species names
      (quotes, trim (name_input(i)), ' mass fraction', quotes, i = 1, n_species)

      if (specified_state(1:1) /= 'P') write (lunplot, '(a)') '"Pressure, Pa"'

      select case (specified_state)

         case ('Rho_T', 'P_T  ') ! Suppress a duplicate temperature

         case default
            write (lunplot, '(a)') '"Temperature, K"'

      end select

      write (lunplot, '(a)') &
         '"rho, kg/m^3"', '"h, J/kg"', '"s, J/kg-mol"', '"Frozen a, m/s"'

      if (interpolate_to_grid) then

         do ib = 1, nblocks

            write (lunplot, '(a, i2, a)') 'ZONE T="Zone', ib, '"'
            write (lunplot, '(a,i3,a,i3,a)') &
            'I=', cfd_plot(ib)%mi, ', J=', cfd_plot(ib)%mj, ', ZONETYPE=Ordered'
            write (lunplot, '(a)') 'DATAPACKING=BLOCK'

            write (lunplot, '(1p, 8e14.6)') cfd_grid(ib)%x
            write (lunplot, '(1p, 8e14.6)') cfd_grid(ib)%y

            do i = 1, 2 + n_species
               write (lunplot, '(1p, 8e14.6)') cfd_plot(ib)%q(i,:,:,1)
            end do

            i = 3 + n_species
            if (specified_state(1:1) /= 'P') then
               write (lunplot, '(1p, 8e14.6)') cfd_plot(ib)%q(i,:,:,1)
               i = i + 1
            end if

            select case (specified_state)

               case ('Rho_T', 'P_T  ') ! Suppress a duplicate temperature

               case default
                  write (lunplot, '(1p, 8e14.6)') cfd_plot(ib)%q(i,:,:,1)
                  i = i + 1

            end select

            do i = nfplot - 3, nfplot  ! rho, h, entropy, frozen speed of sound
               write (lunplot, '(1p, 8e14.6)') cfd_plot(ib)%q(i,:,:,1)
            end do

         end do ! Next block

      end if

!     Uniform grid result, internal units (added as an extra block if a CFD
!     grid is present):

      ib = 1;  if (interpolate_to_grid) ib = nblocks + 1
      write (lunplot, '(a, i2, a)') 'ZONE T="Zone', ib, '"'
      write (lunplot, '(a,i3,a,i3,a)') &
         'I=', uniform_grid(1)%mi, ', J=', uniform_grid(1)%mj,                 &
         ', ZONETYPE=Ordered'
      write (lunplot, '(a)') 'DATAPACKING=BLOCK'

      write (lunplot, '(1p, 8e14.6)') uniform_grid(1)%x
      write (lunplot, '(1p, 8e14.6)') uniform_grid(1)%y

      do i = 1, 2 + n_species
         write (lunplot, '(1p, 8e14.6)') fplot(i,:,:)
      end do

      i = 3 + n_species
      if (specified_state(1:1) /= 'P') then
         write (lunplot, '(1p, 8e14.6)') fplot(i,:,:)
         i = i + 1
      end if

      select case (specified_state)

         case ('Rho_T', 'P_T  ') ! Suppress a duplicate temperature

         case default
            write (lunplot, '(1p, 8e14.6)') fplot(i,:,:)
            i = i + 1

      end select

      do i = nfplot - 3, nfplot  ! rho, h, entropy, and frozen speed of sound
         write (lunplot, '(1p, 8e14.6)') fplot(i,:,:)
      end do

      close (lunplot)

      if (integrate) then

         ni = uniform_grid(1)%ni
         nj = uniform_grid(1)%nj

         call integrate_bulk_mass_and_enthalpy (ni, nj, nfplot,                &
                                                uniform_grid(1)%x,             &
                                                uniform_grid(1)%y, fplot,      &
                                                bulk_mass, bulk_enthalpy)
      end if

      end subroutine save_results

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine save_v1_v2 ()

!     Save just the specified throat distributions for checking.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      logical :: integrate
      real    :: area, quadratures(2)

      integrate = false

      write (lunplot, '(a)') &
         'TITLE = "Nozzle throat input conditions"', &
         'VARIABLES = "x, m", "y, m"'

      select case (specified_state) ! For the name of v1

         case ('Ht_MF')
            write (lunplot, '(a)') '"Total enthalpy, J/kg"'
         case ('Rho_T', 'Rho_H', 'Rho_S')
            write (lunplot, '(a)') '"Density, kg/m^3"'
         case ('P_T  ', 'P_H  ', 'P_S  ')
            write (lunplot, '(a)') '"Pressure, Pa"'

      end select

      select case (specified_state) ! For the name of v2

         case ('Ht_MF')
            write (lunplot, '(a)') '"Mass flux, kg/(m^2.s)"'
            integrate = true
         case ('Rho_T', 'P_T  ')
            write (lunplot, '(a)') '"Temperature, K"'
         case ('Rho_H', 'P_H  ')
            write (lunplot, '(a)') '"Enthalpy, J/kg"'
         case ('Rho_S', 'P_S  ')
            write (lunplot, '(a)') '"Entropy, J/kg-mol"'

      end select

      write (lunplot, '(a, i2, a)') 'ZONE T="Uniform grid"'
      write (lunplot, '(a,i3,a,i3,a)') &
         'I=', uniform_grid(1)%mi, ', J=', uniform_grid(1)%mj,                 &
         ', ZONETYPE=Ordered'
      write (lunplot, '(a)') 'DATAPACKING=BLOCK'

      write (lunplot, '(1p, 8e14.6)') uniform_grid(1)%x
      write (lunplot, '(1p, 8e14.6)') uniform_grid(1)%y
      write (lunplot, '(1p, 8e14.6)') uniform_v1_v2(1,:,:)
      write (lunplot, '(1p, 8e14.6)') uniform_v1_v2(2,:,:)

      close (lunplot)

      if (integrate) then

         ni = uniform_grid(1)%ni
         nj = uniform_grid(1)%nj

         call integrate_bulk_mass_and_enthalpy (ni, nj, 2,                     &
                                                uniform_grid(1)%x,             &
                                                uniform_grid(1)%y,             &
                                                uniform_v1_v2,                 &
                                                bulk_mass, bulk_enthalpy)
      end if

      end subroutine save_v1_v2

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine scale_minor_profile (ni, nj, nf, jmajor, f)

!     Impose a scaled form of each function's minor profile (i = 1) on the rest
!     of the i stations.  The major profile defines the scaling, and may be at
!     either the major semi-axis or at the middle j station (for the case of a
!     wall along the major axis).

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)    :: ni, nj
      integer, intent (in)    :: nf
      integer, intent (in)    :: jmajor  ! Index of major profiles in-place in f
      real,    intent (inout) :: f(nf,ni,nj)

!     Local variables:

      integer :: i, j, jlow, n
      real    :: fscale, fshift

!     Execution:

      do n = 1, nf

         jlow = nj
         do j = nj - 1, 1, -1
            if (f(n,1,j) < f(n,1,jlow)) jlow = j
         end do

         do i = 2, ni

            if (f(n,1,jlow) == f(n,1,jmajor)) then
               fscale = one
               fshift = zero
            else
               fscale = &
                  (f(n,i,jmajor) - f(n,1,jlow)) / (f(n,1,jmajor) - f(n,1,jlow))
               fshift = &
                  (f(n,1,jmajor) * f(n,1,jlow) - f(n,1,jlow) * f(n,i,jmajor)) /&
                                                (f(n,1,jmajor) - f(n,1,jlow))
            end if

            f(n,i,:) = fscale * f(n,1,:) + fshift

         end do

      end do

      end subroutine scale_minor_profile

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine starting_guess ()

!     Initialize the first nozzle throat calculation.  Others start from a
!     previous solution.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real :: gammae, term, To

!     Starting guesses for T and P (first point only)

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

      end subroutine starting_guess

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine store_results (specified_state, ibc, i, j, ni, nj, nf, f,     &
                                nfplot, fplot)

!     Store the flow variables indicated by ibc into f(1:nf,i,j) for the current
!     evaluation point (i,j) and store plottable results in fplot(1:nfplot,i,j).
!     The arguments simplify nomenclature (only), except 'specified_state'
!     allows economy mode to switch to 'P_T  ' for interpolated points.
!     Other variables are inherited from the host.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      character, intent (in) :: specified_state * 5  ! Made an argument when
                                                     ! economy mode was added
      integer, intent (in) :: ibc    ! DPLR's code:
                                     ! 60: primitive vars. rho_s,u,v,Tv,T
                                     ! 61: primitive vars. P,c_s[2-ns],u,v,Tv,T
                                     ! 62: conserved vars. rho_s,rhou,rhov,Ev,E

      integer, intent (in) :: i, j, nf, ni, nj, nfplot

      real, intent (inout) :: f(nf,ni,nj), fplot(nfplot,ni,nj)

!     Local variables:

      integer :: k
      real    :: E

!     Execution:

!     Ensure consistency with DPLR-derived quantities such as molecular weights:

      select case (specified_state)

         case ('Ht_MF')

            rho =    P / (R_bar * T)
            u   = gamf * (R_bar * T) * fMsq  ! u**2
            h   = v1 - half * u
            u   = sqrt (u)                   ! u

         case default  ! Plain Gary cases; reuse part of residual_eq_1

            call frozen_gas_constants (T)

            rho =    P / (R_bar * T)
            u   = gamf * (R_bar * T) * fMsq      ! u**2
            h   = Cp_bar * T + ev_bar + h_ref_sum
            ho  = h + half * u
            u   = sqrt (u)

      end select

      select case (ibc)

         case (60)

            f(1:n_species,i,j) = rho * mass_fraction(:)
            f(n_species+1,i,j) = u
            f(n_species+2,i,j) = zero
            f(n_species+3,i,j) = zero
            f(n_species+4,i,j) = T
            f(n_species+5,i,j) = T

         case (61)

            f(1,i,j) = P
            f(2:n_species,i,j) = mass_fraction(2:n_species)
            f(n_species+1,i,j) = u
            f(n_species+2,i,j) = zero
            f(n_species+3,i,j) = zero
            f(n_species+4,i,j) = T
            f(n_species+5,i,j) = T

         case (62)

            E = rho * h - P
            f(1:n_species,i,j) = rho * mass_fraction(:)
            f(n_species+1,i,j) = rho * u
            f(n_species+2,i,j) = zero
            f(n_species+3,i,j) = zero
            f(n_species+4,i,j) = rho * ev_bar
            f(n_species+5,i,j) = E

      end select

!     Plottable results:

      if (specified_state == 'Econo') then
         v1 = ho        ! Total enthalpy
         v2 = rho * u   ! Mass flux
      end if

      fplot(1,i,j)   = v1
      fplot(2,i,j)   = v2;              k = 2 + n_species
      fplot(3:k,i,j) = mass_fraction(:)
      fplot(k+1,i,j) = P
      fplot(k+2,i,j) = T
      fplot(k+3,i,j) = rho
      fplot(k+4,i,j) = h
      fplot(k+5,i,j) = s
      fplot(k+6,i,j) = u

      end subroutine store_results

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine uniform_grid_gen ()

!     Generate a single-block uniform grid for initial in-flow calculations.
!     Actually, not all cases lead to truly uniform spacing, but all that is
!     needed is an interim grid from which results can be interpolated to the
!     desired CFD grid later if requested.
!
!     Three types of elliptic grid appear here for historical reasons.  The most
!     obvious polar grid is suited more to a circular nozzle (for a 3-D test
!     article, presumably) or quadrant of an elliptic nozzle than to the semi-
!     ellipse with flat floor that is now better handled specially with a family
!     of i-lines allowing the major profile to bisect the space between the two
!     walls.
!
!     The third type of ellipse is a kludge to help the rectangle case, since
!     TFI behaves better then.  However, the Cartesian product approach is now
!     preferable for the rectangle case.  (The interim ellipse approach simply
!     extends the spokes to form a rectangle just before saving results, thus
!     avoiding concavity observed with flat-top flow distributions on Cartesian
!     rectangular grids.  But to repeat: plain rectangle + Cartesian product
!     appears preferable now.)

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      real, parameter :: ninety = 90.
      real, parameter :: z = 0.  ! Reusable procedures assume 3-space grids

!     Local variables:

      integer :: i, j, jj
      real    :: ai, bi, cs, dtheta, dx, dy, phi, ratio, sn, theta, total, x, y
      real    :: s1(3), s2(3)
      real, allocatable, dimension (:) :: s, x0, y0, z0, xn, yn, zn

!     Execution:

      nj_ellipse = nj_uniform  ! For the rectangle or semi-ellipse, but ...
      if (nozzle_type == 'P') nj_ellipse = ni_uniform + nj_uniform - 1

      uniform_grid(1)%ni = ni_uniform;  uniform_grid(1)%mi = ni_uniform
      uniform_grid(1)%nj = nj_ellipse;  uniform_grid(1)%mj = nj_ellipse
      uniform_grid(1)%nk = 1;           uniform_grid(1)%mk = 1

      call xyz_allocate (uniform_grid(1), ios)

      ratio  =  b / a   ! Minor to major axis lengths - see control_inputs
      dx     =  a / real (ni_uniform - 1)
      dy     =  b / real (nj_uniform - 1)

      select case (nozzle_type)

         case ('R')  ! Rectangular grid (plain)

            do j = 1, nj_uniform
               y = r_center_v + dy * real (j - 1)
               do i = 1, ni_uniform
                  x = r_center_h + dx * real (i - 1)
                  uniform_grid(1)%x(i,j,1) = x
                  uniform_grid(1)%y(i,j,1) = y
                  uniform_grid(1)%z(i,j,1) = z
               end do
               uniform_grid(1)%x(ni_uniform,j,1) = r_wall_h
            end do
            uniform_grid(1)%y(:,nj_uniform,1) = r_wall_v

         case ('C', 'E')  ! Circle | semi-ellipse: nj_uniform = # uniform angles

            dtheta = ninety / real (nj_uniform - 1)

            do j = 1, nj_uniform
               theta = dtheta * real (j - 1)
               cs = cosd (theta)
               sn = sind (theta)
               do i = 1, ni_uniform
                  ai = dx * real (i - 1)
                  bi = ai * ratio
                  uniform_grid(1)%x(i,j,1) = r_center_h + ai * cs
                  uniform_grid(1)%y(i,j,1) = r_center_v + bi * sn
                  uniform_grid(1)%z(i,j,1) = z
               end do
            end do

         case ('P')  ! Polar grid/rectangular nozzle; nj = ni_uni + nj_uni - 1

!           We use a polar elliptic grid that captures the outer corner of the
!           eventual rectangle, so generate it in two stages.

            theta  = 45.;  phi = theta
            dtheta = theta / real (nj_uniform - 1)

            do j = 1, nj_uniform
               theta = dtheta * real (j - 1)
               cs = cosd (theta)
               sn = sind (theta)
               do i = 1, ni_uniform
                  ai = dx * real (i - 1)
                  bi = ai * ratio
                  uniform_grid(1)%x(i,j,1) = r_center_h + ai * cs
                  uniform_grid(1)%y(i,j,1) = r_center_v + bi * sn
                  uniform_grid(1)%z(i,j,1) = z
               end do
            end do

            dtheta = phi / real (ni_uniform - 1)
            jj = nj_ellipse

            do j = 1, ni_uniform - 1
               theta = ninety - dtheta * real (j - 1)
               cs = cosd (theta)
               sn = sind (theta)
               do i = 1, ni_uniform
                  ai = dx * real (i - 1)
                  bi = ai * ratio
                  uniform_grid(1)%x(i,jj,1) = r_center_h + ai * cs
                  uniform_grid(1)%y(i,jj,1) = r_center_v + bi * sn
                  uniform_grid(1)%z(i,jj,1) = z
               end do
               jj = jj - 1
            end do

         case ('S')  ! Semi-ellipse nozzle via specialized grid

!           Strategy:
!
!           Construct i-lines (in the j direction) the mid-points of which will
!           lie along the bisector of the two wall boundaries.  This means nj
!           should be odd.  The profile specified for the major semi-axis will
!           be imposed along this bisector, not along the floor, where the
!           (constant) conditions are implied by the minor semi-axis profile.
!           [Actually, if the minor profile peaks at the center line, the major
!           profile will be imposed along the major axis.  This suits the case
!           of a circular nozzle with a non-axisymmetric test article.]

!           Since arc lengths along most curves are difficult to calculate, we
!           go with the uniform angle distribution for the outer boundary points
!           and use the same relative distribution for the major semi-axis.

!           The outer boundary is the j = nj line; i = ni is a collapsed edge.

            dtheta = ninety / real (ni_uniform - 1)
            j = nj_uniform

            do i = 1, ni_uniform
               theta = ninety - dtheta * real (i - 1)
               cs = cosd (theta)
               sn = sind (theta)
               uniform_grid(1)%x(i,j,1) = r_center_h + a * cs
               uniform_grid(1)%y(i,j,1) = r_center_v + b * sn
               uniform_grid(1)%z(i,j,1) = z
            end do

            allocate (s(ni_uniform))

            call chords2d (ni_uniform, uniform_grid(1)%x(1,j,1),               &
                           uniform_grid(1)%y(1,j,1), true, total, s)

!           Major semi-axis:

            do i = 1, ni_uniform
               uniform_grid(1)%x(i,1,1) = r_center_h + s(i) * a
               uniform_grid(1)%y(i,1,1) = r_center_v
               uniform_grid(1)%z(i,1,1) = z
            end do

!           Minor semi-axis and collapsed outer edge:

            i = ni_uniform
            do j = 1, nj_uniform
               uniform_grid(1)%x(1,j,1) = r_center_h
               uniform_grid(1)%y(1,j,1) = r_center_v + dy * real (j - 1)
               uniform_grid(1)%z(1,j,1) = z
               uniform_grid(1)%x(i,j,1) = r_center_h + a
               uniform_grid(1)%y(i,j,1) = r_center_v
               uniform_grid(1)%z(i,j,1) = z
            end do

!           Morph the minor axis grid line to each i station, controlling the
!           end point slopes as we go:

            j = nj_uniform

            allocate (x0(j), y0(j), z0(j), xn(j), yn(j), zn(j))

            do j = 1, nj_uniform
               x0(j) = r_center_h
               y0(j) = uniform_grid(1)%y(1,j,1)
               z0(j) = z
            end do

            yn(1) = r_center_v
            zn(1) = z
            zn(nj_uniform) = z

            s1(1) = zero   ! Unit vector components at end points
            s1(2) = one
            s1(3) = zero;  s2(3) = zero

            total = sqrt (a*a + b*b)

            do i = 2, ni_uniform - 1

               xn(1) = uniform_grid(1)%x(i,1,1)

               j     = nj_uniform
               xn(j) = uniform_grid(1)%x(i,j,1)
               yn(j) = uniform_grid(1)%y(i,j,1)

               theta = ninety - dtheta * real (i - 1)
               cs = cosd (theta)
               sn = sind (theta)

               s2(1) = b * cs / total  ! Normal to ellipse derived from tangent
               s2(2) = a * sn / total

               call morph_line_3d (1, j, s1, s2, x0, y0, z0, xn, yn, zn)

               do j = 2, nj_uniform - 1
                  uniform_grid(1)%x(i,j,1) = xn(j)
                  uniform_grid(1)%y(i,j,1) = yn(j)
                  uniform_grid(1)%z(i,j,1) = z
               end do

            end do

            deallocate (s, x0, y0, z0, xn, yn, zn)

      end select

      end subroutine uniform_grid_gen

   end program throat_conditions_3d

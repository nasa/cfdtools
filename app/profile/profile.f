C+------------------------------------------------------------------------------
C
      PROGRAM PROFILE
C
C  PURPOSE:
C           PROFILE is a utility for manipulating and/or displaying the
C           coordinates and other properties of airfoils.  It processes
C           one or more profiles at a time in various ways (one way per
C           run).   The input and output profiles may be in one of four
C           formats described below (basically as separate surfaces, as
C           a wrap-around surface, or in three-column form).
C
C           If plotting of the airfoil geometry is  requested,  all  of
C           the profiles  in  the input  dataset will be plotted on the
C           same pair of axes unless the  "THREED"  option is selected.
C           The plots  can  be  of  the  original input data,  or  data
C           derived by PROFILE,  or both.   Curvature distributions may
C           also be plotted as may the optional pressure distributions.
C
C           Plotting of the airfoil geometries, curvature distributions
C           or pressure distributions is handled by a separate program,
C           QPLOT, which should accompany PROFILE.    Users' guides are
C           available for PROFILE and QPLOT.
C
C           The  four  main choices for manipulating input profiles are
C           to REDISTRIBUTE the abscissas using conventional splines or
C           parametric splines depending on whether leading or trailing
C           edges are rounded;  to MODIFY or perturb the  ordinates  in
C           some way according to shape functions added  interactively;
C           to REFINE the ordinates by manipulating the curvature (act-
C           ually 2nd derivative) distribution while seeking or retain-
C           ing some maximum thickness value; and to OPTIMIZE one surf-
C           ace of one airfoil automatically, using a predetermined set
C           of shape functions with some of their parameters  variable,
C           and a target curvature distribution to  be  matched  in the
C           least squares sense.
C
C           Lesser options permit the user to RECTIFY profiles which do
C           not have the point common to the two surfaces in the  usual
C           place, and to NORMALIZE or DENORMALIZE profiles. TRANSFORM-
C           ing between upper/lower surface representation and  camber/
C           thickness representation (either way) is provided for, with
C           decambering as an option.  Applying twist is available from
C           the ROTATE option.
C
C           Two options involve a "secondary" profile  (prompted for at
C           a lower level): an option to COMBINE profiles (one of which
C           may be boundary layer displacement thickness);  and  an op-
C           tion to LOFT linearly between two profiles.
C
C           A  "nose-job"  option permits ROUNDing or SHARPENing of the
C           leading edge region  -  operations which have been made  as
C           reversible as possible.   More generally, the SMOOTH option
C           smooths the entire airfoil (or just one surface) by fitting
C           combinations of "Wagner" functions, which are also employed
C           by the OPTIMIZE and MODIFY options,  or  by implicit and/or
C           explicit smoothing (possibly weighted nonuniformly).
C
C           Tabulation of coordinates along with derivatives and curva-
C           ture is provided, as well as saving the manipulated profile
C           as a new output file.
C
C           Saving of y" distributions is also an option,  for possible
C           editing and reuse in REFINE mode.
C
C           Spreadsheet-compatible output of all likely tabular quanti-
C           ties is also provided for the simpler operating modes. This
C           requires the two surfaces to have common abscissas.
C
C  NOTES:
C           PROFILE  has evolved considerably since its inception as  a
C           basic redistribute-and/or-plot program. Provision for arbi-
C           trary perturbations to the input geometry,  with tabulation
C           of the resultant coordinates and derivatives, required some
C           reorganization, but the structure should now serve for most
C           likely purposes. Some implementation considerations follow.
C
C        *  The case of a 2-element airfoil forced the decision to plot
C           all input profiles on the same page  (although the "THREED"
C           option has since been introduced).    Normalization of such
C           combinations proved an awkward option to provide.  The user
C           should  use  the normalization option carefully.  Normaliz-
C           ation of 3-D wings is not available.
C
C        *  The multiple-profile case also prevented plotting  of  more
C           than one frame type (such as curvature distributions in ad-
C           dition  to  the  airfoils)  -  hence the saving of separate
C           files for later plotting.
C
C        *  Large-scale plots are feasible  (with optional  windowing),
C           but exact scaling cannot be guaranteed because of the vari-
C           ability  of  the  output devices  available.  (Plots of the
C           same data on the same device can vary slightly from plot to
C           plot.)
C
C        *  Derivatives and curvature values are normally estimated  by
C           finite differences for consistency with REFINE and OPTIMIZE
C           modes.    It is well known that these approximations can be
C           poor in the presence of very small X increments and limited
C           precision in the Ys.
C
C        *  An option to plot the full wrap-around curvature  distribu-
C           tion using parametric spline derivatives has been  provided
C           for a proper look at the leading edge region.  But the .ypp
C           file of 2nd derivatives is suppressed in this case to avoid
C           inappropriate use with the REFINE mode.
C
C        *  For simplicity, each of the MODIFY,  REFINE,  and  OPTIMIZE
C           options assumes that the coordinates have been normalized.
C
C  METHOD:
C
C           The basic steps are as follows:
C
C        *  Prompt for mode of operation and the input profile file name.
C
C        *  Set defaults for user inputs and use an input control file
C           to override some of them if necessary.
C
C        *  Scan all of the input profiles, for scaling and normalizing
C           purposes. Use EOF to handle the unknown number of profiles.
C
C        *  Rewind the file and process the (now known number of) profiles
C           one at a time, according to the selected mode.
C
C        *  Write the following output files as requested:
C
C              Revised profile coordinates in one of 4 formats
C
C              Original and/or revised airfoil geometry for plotting
C              (a QPLOT file)
C
C              Tabulated coordinates with derivatives and curvatures
C
C              A more complete, spreadsheet-compatible file
C
C              Second derivatives for possible reuse by REFINE mode
C
C              Original and revised curvature data, including target
C              curvature data for OPTIMIZE mode (another QPLOT file)
C
C              Cps estimated for original and revised airfoil (QPLOT
C              format)
C
C
C  MODES OF OPERATION:
C
C           MODE 0:  "Display" mode - no modifications involved.   Gen-
C                    erate requested output files,  which could include
C                    saving the coordinates in a different format. MODE
C                    <=3 is required for spreadsheet-compatible output.
C
C           MODE 1:  Rearrange or rectify the geometry data so that the
C                    common leading-edge point is indeed the  one  with
C                    minimum abscissa and shift ordinates by  an  input
C                    input amount if required.   Only  the revised pro-
C                    file may be tabulated/plotted in this case. A ver-
C                    tical shift option is also provided.
C
C           MODE 2:  Normalize profile(s) according to the total  range
C                    of x or by input chord & leading edge coordinates.
C                    A negative chord value will denormalize.  The same
C                    input values are used for each element of a multi-
C                    element airfoil.
C
C           MODE 3:  Redistribute the abscissas  and derive correspond-
C                    ing ordinates.   Conventional or parametric spline
C                    techniques are used depending on whether the lead-
C                    ing edge is sharp or rounded.  Distributions along
C                    the arc (in T rather than X) are an option.  Menu:
C
C                    -1 = Read new Xs (or Ts) from a file in standard
C                         PROFILE format (though y coordinates may be
C                         omitted if desired).
C                     0 = Distribute points uniformly.
C                     1 = Distribute points sinusoidally bunched near
C                         the leading edge.
C                     2 = Distribute  points  sinusoidally, near both
C                         the leading and trailing edges.
C                     3 = Sinusoidal bunching around an internal pt.
C                     4 = Vinokur distribution (first, last increments
C                         increments specirfied.
C
C                    A prompt  will  also  be  issued for the number of
C                    points to be generated on each surface.
C
C           MODE 4:  Perturb geometry data  according  to user-selected
C                    shape functions (interactive).
C
C           MODE 5:  Refine the airfoil,  typically modifying  (or  re-
C                    taining) its thickness while retaining (or modify-
C                    ing) its curvature distribution(s).  Numerous user
C                    inputs are prompted for in this case,  since there
C                    are several likely ways provided for  manipulating
C                    the curvature distributions. Defaults are provided
C                    where possible.  The 4 main choices:
C
C                    (1) Leave a surface unchanged  (while presumably
C                        refining the other);
C                    (2) Change the thickness with minimal changes to
C                        the existing curvature  (none of the follow-
C                        ing constraints on the  original y" values);
C                    (3) Impose a constant-2nd-derivative  constraint
C                        in some single region (to remove a bump or a
C                        spike in the curvature distribution, or per-
C                        haps  to  modify regions of flatness by con-
C                        straining second derivatives - and hence the
C                        curvature - away from zero);
C                    (4) Constrain the curvature via an input file of
C                        2nd derivative values (possibly derived from
C                        an earlier run  of  PROFILE,  or prepared by
C                        hand).  The table does not have to cover the
C                        whole abscissa range; linearly interpolating
C                        table look-ups are used.
C
C                    Brief descriptions of the inputs prompted  for  in
C                    "refine" mode follow:
C
C             *  Desired % thickness:  <CR>  retains present thickness.
C
C             *  Width param. for y:  Affects  the  nonuniform  scaling
C                                     applied  to  the ordinates  (both
C                    surfaces).   The default is 2.0.  Larger (3.0-4.0)
C                    tends to retain leading/trailing  edge  shape more
C                    while 1.0 would constrain fore and aft less.
C
C             *  Input y" table:      <CR>  means there is  none,  else
C                                     the file name is  entered.   This
C                    file should be in the standard  "PROFILE"  format.
C                    It can  cover  any  range  of  abscissas.  (Linear
C                    interpolation  is  used.)   It  may  be an  edited
C                    version  of  the  file  from  a  previous  run  of
C                    PROFILE, or it may be much cruder.  The 2nd deriv-
C                    ative values  entered  act  as constraints on  the
C                    curvature since curvature and y" are related if y'
C                    is not large.
C
C             *  Constant y" value:   <CR>  means no such constraint  -
C                                     retain existing curvature  values
C                    as much as possible.   Otherwise,  a  value  of y"
C                    entered will be sought in the  abscissa range that
C                    is prompted for next.
C
C             *  Corresp. x range:     Enter low and high  x  values on
C                                      the same line.   Allow  for  the
C                    fact that strict inequalities are  used  when  the
C                    program tests for being within this range.   E.g.:
C                    Enter  .39 .61  or  .39,.61  if you intend for the
C                    constraint to apply in [0.4,0.6].
C
C             *  Width param. for y":  Default is  3.0.   Affects  non-
C                                      uniform weighting  of  the equa-
C                    tions representing 2nd derivative  constraints  in
C                    the overdetermined system being solved.  Since the
C                    actual values of the  2nd derivatives being sought
C                    also act in  a weighting sense,  effects  of  this
C                    variable are not easy to predict. Values of 2.0 or
C                    1.0 should tend to let y" change more.
C
C             *  x for peak y" weight: The absolute values  of  Y"  are
C                                      so much bigger than those  of  Y
C                    that they all need to be scaled down in the system
C                    being solved.   If you are trying  to  flatten the
C                    curvature plot in some region,  pick the center of
C                    the region for this input. Otherwise, use the mid-
C                    chord value.
C
C             *  y" weights, x/c=0,1:  Default is 0.004. See next item.
C
C             *  peak y" weight:       Default is 0.04.  These  provide
C                                      for the fact that  the  absolute
C                    values of y" are  typically  smaller in  the  mid-
C                    section than near the  leading/trailing  edges, so
C                    they should be weighted more,  especially  in view
C                    of the fact that any  y"  constraints applied  are
C                    typically in the mid-section.  See above.
C
C           MODE 6:  Optimize one surface of one profile  using  a pre-
C                    determined set of shape functions, some parameters
C                    of which are automatically varied so as to achieve
C                    a curvature distribution matching some target cur-
C                    vatures in the least squares sense.
C
C           MODE 7:  Transform  representation  of  profile(s)  between
C                    upper/lower surface  and  camber/thickness (either
C                    way - the user is prompted for the direction).  An
C                    option to decamber a section is also provided.
C
C           MODE 8:  Rotate a profile about some point to apply twist.
C
C           MODE 9:  Combine primary profile with a  secondary  profile
C                    (read from a separate file).  This was prompted by
C                    a need to add or remove a boundary layer displace-
C                    ment thickness  (positive by definition)  but  has
C                    been arranged to handle true  addition/subtraction
C                    of distinct profiles as well.
C
C           MODE 10: Loft between primary and secondary profiles.
C
C           MODE 11: Nose-job option: Round or sharpen the leading edge
C                    region via splines.
C
C           MODE 12: Smooth either surface or both surfaces using least
C                    squares techniques (linear combination of n Wagner
C                    functions plus a "ramp" for thick trailing edges),
C                    or by implicit and/or explicit methods involving a
C                    (possibly nonuniform) weighting of y".
C
C
C  GEOMETRY INPUT:
C
C           Standard PROFILE format is shown below. The lower surface
C           is optional, but a zero must be read for NL if no lower
C           surface is included unless this is the last airfoil in the
C           file (meaning EOF can be used to indicate NL=0).  In this
C           case, a symmetrical airfoil is assumed.
C
C               TITLE                   <CHARACTER*80>
C               NU   Upper surface      <Integer # pts., first token>
C               X         Y             <Reals, first two tokens>
C               X         Y                 :
C               :         :                 :     (May be X/C, Y/C;
C               :         :                 :      Xs are increasing)
C               :         :                 :
C               NL   Lower surface      <Integer, first token> <may be 0>
C               X         Y             <Reals, first two tokens>
C               X         Y                 :
C               X         Y  ! Trailing comments are permitted
C               :         :                 :
C               :         :                 :     (Xs are increasing)
C               :         :                 :
C               ! X         Y           <Point suppressed; NL must be adjusted>
C               :         :                 :
C               :         :                 :
C
C    NOTE:  For standard format, if both surfaces are present, PROFILE
C           expects them to have the same leading edge point.  The trailing
C           edge points may differ.
C
C           The next two formats are wrap-around clockwise and wrap-around
C           counterclockwise, where the coordinates begin at the trailing
C           edge, wrap around the leading edge, and end at the trailing edge.
C           The clockwise case begins with the lower surface, and the counter-
C           clockwise case begins with the upper surface.  The format shown
C           below is essentially the same for both cases.  NPTS is the total
C           number of points on the airfoil.
C
C               TITLE                   <CHARACTER*80>
C               NPTS                    <Integer, first token>
C               X         Y             <Reals, first two tokens>
C               X         Y                 :
C               :         :                 :     (May be X/C, Y/C;
C               :         :                 :      Xs are decreasing
C               :         :                 :      until the leading
C               :         :                 :      edge, then increasing)
C
C    NOTE:  Wrap-around formats do NOT have duplicate leading edge points.
C
C           The fourth format is three-column format.  The airfoil has
C           the same abscissas for both surfaces in the 1st column and
C           ordinates for the upper and lower surfaces in the 2nd and 3rd
C           columns respectively.  Abscissas are increasing as with standard
C           format.  Here NPTS is the number of points on either surface.
C
C               TITLE                           <CHARACTER*80>
C               NPTS                            <Integer, first token>
C               X         YU        YL          <Reals, first 3 tokens>
C               X         YU        YL              :
C               :         :         :               :   (May be X/C, Y/C;
C               :         :         :               :    Xs are increasing)
C               :         :         :               :
C               :         :         :               :
C
C
C  CONTROL INPUTS:
C
C           A file containing keyword inputs and values may be used to
C           override default options.  In general, the keywords refer to
C           the airfoil plot file and other output options, and apply
C           to all modes.   Prompts are issued for inputs needed by a
C           particular mode.
C
C  KEYWORD GUIDELINES AND DEFINITIONS:
C
C           Keyword/value pairs may appear with more than one pair on a
C           line.  However, the multivalued keywords PLTLINE, CPSLINE,
C           CRVLINE, and NOFILE must not appear with other keywords on
C           the same line.
C
C           The default value in each case appears in square brackets.
C
C  KEYWORD  VALUES and synonyms     DESCRIPTION
C  -------  -------------------     -----------
C
C  FORMAT   [SAME]                  One of four formats for output profile.
C           PROFILE or STANDARD     May be in standard PROFILE format (ab-
C           CLOCKWISE or WRAPAROUND scissas increasing), clockwise wrap-around
C           COUNTERCLOCKWISE        format, counterclockwise wrap-around for-
C           THREE-COLUMN or         mat, or 3-column format.  SAME means the
C           THREE_COLUMN or         same format as the input profile.  NOTE:
C           THREECOLUMN  or         To allow easily for several synonyms for
C           TABLE                   for the THREE-COLUMN value, only the first
C                                   5 characters of the value are checked.
C
C  PLTLINE  [DEFAULT]               Controls line types of curves on profile
C           LINE                    plots.  One value may be included for
C           DASH                    each curve on the plot.  The default is
C           DOT                     symbols connected by a solid line, with
C           CHAINDASH               a different symbol type for each succes-
C           CHAINDOT                sive curve.  The first curve typically
C           THICK                   represents the original profile; the
C           SYMBOLS                 second curve represents the revised one.
C                                   Overriding the default might be desirable
C                                   when plotting multi-element airfoils or
C                                   when lines without symbols are required.
C                                   At most 20 curves are provided for.  Note:
C                                   All the line types in QPLOT are available.
C                                   SYMBOLS refers to symbols with no line
C                                   connecting them.
C
C  CPSLINE  [see PLTLINE above]     Controls line types on Cps plots in the
C                                   same manner as PLTLINE above.  One value
C                                   per curve may be included, chosen from
C                                   the same list of values as those shown
C                                   for PLTLINE.
C
C  CRVLINE  [see PLTLINE above]     Controls line types on curvature plots
C                                   in the same way as PLTLINE and CPSLINE.
C
C  CURVATURE or [NONPARAMETRIC] or  CURVATURE and DERIVATIVES are synonymous
C  DERIVATIVES  [FINITE_DIFFERENCE] controls for the type of calculations
C               SPLINE     or       used for derivatives and hence curvature.
C               PARAMETRIC or       The default is separate-surface treatment
C               WRAPAROUND          using finite differences, as needed for
C                                   consistency with PROFILE's REFINE and
C                                   OPTIMIZE options.  The two surfaces appear
C                                   as separate frames in the curvature plot.
C                                   Otherwise, the full wrap-around curvature
C                                   distribution is calculated using a para-
C                                   metric spline and plotted on a single frame.
C
C                                   The default normally suffices except if the
C                                   region of interest is very near a rounded
C                                   leading edge.  Note that not all of the
C                                   possibilities are provided for, such as
C                                   parametric finite differences.
C
C  MINCURVATURE   [-5.]             Cutoff values for plotted curvatures.
C  MAXCURVATURE   [+5.]             Practice shows that +/-5. give useful
C                                   plot scaling by ignoring the high curv-
C                                   ature values near the leading edge.  On
C                                   the other hand, it may well be desired
C                                   to focus on the leading edge region.  Set
C                                   both to 999. to obtain the full range.
C                                   See the CURVATURE/DERIVATIVES control.
C
C  NOFILE   [NONE]                  Used to suppress any combination of the
C           DAT                     seven output files generated by PROFILE.
C           PLT                     The values correspond to the extensions
C           TAB                     of the file names.  See elsewhere for a
C           CRV                     complete description of file contents.
C           YPP                     NONE serves only to assist leaving the
C           CPS                     NOFILE control word in the input file
C           SPREAD                  even if all outputs are desired.
C
C  PLOT     [BOTH]                  Controls plotting of original OR revised
C           ORIGINAL                profile.  The default is to plot both
C           REVISED                 original and revised (if one exists).
C
C  PRECISION   [FULL]               Controls number of digits in output
C           ENGINEERING             airfoil coordinates.  FULL gives F11.8
C                                   if possible, or E15.8 if any X >=10.
C                                   ENGINEERING gives the traditional F10.6
C                                   common to many flow solvers.
C
C  THREED   [FALSE] or [NO]         For plotting of multiple stations from
C           TRUE or YES             a 3-D wing. The default is the 2-D case.
C
C  XAXIS    [6.4]                   Length of x-axis in inches.  The default
C                                   is normally appropriate for an 8.5 x 11
C                                   page in portrait mode.
C
C           The following four keywords apply to windowing.  Any or none
C           of them may be used.
C
C  XMIN     [minima and             Minimum abscissa for desired window
C  XMAX     maxima of               Maximum abscissa for desired window
C  YMIN     the input               Minimum ordinate for desired window
C  YMAX     coordinates]            Maximum ordinate for desired window
C
C
C  SAMPLE CONTROL FILE:
C
C           A sample input file follows.  Note that keywords and values
C           may be separated with blanks, commas, colons, equals signs,
C           pr tabs. Remember, keywords with more than one value should
C           appear on separate lines.  Any keyword or text value may be
C           truncated to unambiguous leading characters.   Blank  lines
C           and trailing ! comments are ignored.
C
C
C           FORMAT = STANDARD   PRECISION = FULL
C           PLOT BOTH  THREED:NO
C           PLTLINE = SOLID, SOLID
C           CPSLINE = DOT, SYMBOLS
C           CRVLINE = SOLID, DASH, CHAINDOT
C           XAXIS = 20.
C           XMIN = 0.  XMAX 0.1
C           MAXCURVATURE = 999.   ! Both 999. means plot the full
C           MINCURVATURE = 999.   ! curvature range
C           DERIVATIVES = PARAMETRIC
C           NOFILE: YPP SPREAD
C
C
C  OUTPUT FILES:
C
C           The following seven output files are generated by  PROFILE.
C           Any  of  the  files  may  be  suppressed  using the keyword
C           NOFILE.  The user is prompted for an identifier, which will
C           become  the first part of each file name.   File extensions
C           are fixed.
C
C  <identifier>.DAT   Contains airfoil coordinates  that  have been re-
C                     vised in some way. May be in one of four formats;
C                     default is the same format as input coordinates.
C
C  <identifier>.PLT   Contains airfoil geometry coordinates  and  other
C                     information necessary for later QPLOTing.  May be
C                     a plot of original profile,  revised profile,  or
C                     both.  (Default is both, superimposed.)
C
C  <identifier>.TAB   Contains tabulated coordinates,  first and second
C                     derivatives and curvatures for the  original  and
C                     revised profile  (if  a  revised profile exists).
C                     Other diagnostics may be written here,  including
C                     a record of selected  shape  functions  from  the
C                     MODIFY option.
C
C  <identifier>.CRV   Contains curvatures of the  original  and revised
C                     profiles (if a revised profile exists)  for later
C                     QPLOTing.  Also to be used as the basis of target
C                     curvature data when OPTIMIZE option is used.  The
C                     MAXCURVATURE and MINCURVATURE keywords  determine
C                     the plot axis range, but the full surfaces appear
C                     in the file (except for the first and last pts.).
C                     See the description of these keywords for how  to
C                     obtain a full wrap-around curvature distribution.
C
C  <identifier>.YPP   Contains second derivatives  in  standard PROFILE
C                     format,  for  possible  reuse  in  "refine" mode.
C                     When a profile has been  revised,  the  file will
C                     contain  only  the  second  derivatives  of   the
C                     revised profile.   Otherwise,  the file will con-
C                     tain only the second derivatives of the  original
C                     profile.
C
C  <identifier>.CPS   Contains QPLOTable estimates of Cps for user-sup-
C                     plied alpha and free stream Mach number  (revised
C                     and/or original).
C
C  <ident>.SPREAD     Contains  spreadsheet-compatible  (tab-delimited)
C                     tabular data.   Available only if MODE <= 3,  and
C                     only if both surfaces have common abscissas.
C
C
C  LOGICAL UNIT NUMBERS:
C
C    LUNCPS      For output Cp estimates in QPLOT format.
C    LUNCRT      For prompts and diagnostics.
C    LUNCRV      For output curvature data in QPLOT format.
C    LUNDAT      For input file of geometry data.
C    LUNINP      For keyword control inputs (disk file).
C    LUNKBD      For responses to prompts.
C    LUNPLT      For airfoil geometry in QPLOT format.
C    LUNREV      For output file of revised geometry data.
C    LUNSPR      For spreadsheet-compatible output.
C    LUNTAB      For tabulations and diagnostics.
C    LUNTBL      For optional table of y" data ("refine" mode; PROFILE format).
C                Alternatively, target curvature distribution ("optimize" mode,
C                QPLOT format).
C    LUNXGR      For file for reading abscissas ("redistribute" mode),
C                or for file containing bump function info. ("optimize" mode,
C                keyword format), or for secondary profile ("combine" mode).
C    LUNYPP      For output y" data in PROFILE format.
C
C
C  EXTERNAL REFERENCES:
C
C    AFGEOM      Calculates various airfoil geometric properties
C    BOUNDS      Finds the max/min values in array(s)
C    COMBINE     Combines primary and secondary profiles
C    COPY,RVERSE Utilities for transferring data
C    CURV2D      Curvature via parametric derivatives
C    FD12K       Finite difference approx. to derivatives/curvature
C    GETCL       Computes aerodynamic coefficients
C    GETCPS      Cheap estimate of Cp distributions
C    LOFT        Lofts linearly between two profiles
C    LSTFIL      Copies a formatted file to screen or disk
C    MAXMIN      Finds maximum and minimum ordinates and abscissas
C                over upper and lower surfaces of all profiles
C    MODIFY      Interactive subroutine to modify  the  geometry  by
C                applying "bump" functions
C    NOSEJOB     Modifies the leading edge region in various ways
C    NRMLIZ      Normalizes coordinates
C    NRMSET      Prompts for the details needed by NRMLIZ
C    ROTATE      Rotates coordinates with option to renormalize
C    OPENER      File opening utility
C    OPTIMIZE    Optimizes one surface using bumps/target curvatures
C    PROTECT     Checks for duplicate points and monotonicity
C    PRREAD      Reads one airfoil profile
C    PRTAB       Tabulates coordinates, derivatives,  and curvatures
C                (which it calculates from the given derivatives)
C    PRWRIT      Saves the profile geometry data on disk in the format
C                described in PRREAD
C    PSFIT       Parametric spline fit and evaluation, used here for
C    PSTVAL      wrap-around curvature calculations.
C    QPLDAT      Entry points QPLFRM and QPLCRV permit saving of
C                QPLOTable data - used here for 2nd derivatives and Cps
C    RDKEYS      Reads keyword input control file
C    READER      Prompting utility
C    RECTIFY     Reorganizes the input geometry so that  the  common
C                leading edge point has in fact the minimum abscissa.
C                Also shifts ordinates by a user-specified amount.
C    REDISTRIB   Redistributes data points using splines
C    REFINE      Refines thickness and/or curvature distribution(s)
C    SCAN2       Scans a string for the next token
C    SELECT      Menu selection by item number or by name
C    TRANSFORM   Transforms representation of one profile, or decambers it.
C    XDERIVS     Forms X derivatives from parametric derivatives
C    XFORM       Called by TRANSFORM, but also here for spreadsheet info.
C
C  ERROR HANDLING:
C    This is mostly confined to reading the input files.  See subroutines
C    PRREAD and RDKEYS.
C
C  ENVIRONMENT:
C    DEC VMS; SGI IRIX; Fortran 90
C
C  HISTORY:
C
C  Dec. '82  DAS/RAK  Original design.
C  12/09/82    LJC    Original coding (plot/redistribute data points).
C  04/29/83    LJC    Added option to "rectify" input data points.
C  July '83    LJC    Added a 3-D capability and an option for reading
C                     a file of new abscissas in "redistribute"  mode.
C  09/27/83    LJC    Interactive MODIFY routine incorporated.
C  Oct. '83    DAS    Integrated alternative version of MODIFY  as the
C                     REFINE option; provided for saving curvature and
C                     y" values; removed REDISTRIB function from main.
C  11/09/83    LJC    Added de-normalizing to normalizing option;  in-
C                     cluded these as a separate MODE; reordered MODEs
C                     in order of complexity.
C  01/20/84    DAS    Incorporated OPTIMIZE mode.
C  04/12/84    DAS    Took advantage of QPLOT's new legend capability.
C  July '84    LJC    Added reading and writing of wraparound formats.
C  07/24/84    LJC    Added calculation of thickness for all modes.
C  Aug. '84    LJC    Changed from namelist to keyword inputs.
C  09/14/84    LJC    Changed legend entry to be read from dataset and
C                     added prompt for title.  Formerly, the title was
C                     read from dataset and the legend was hard-coded.
C  Oct. '84    LJC    Arranged for all plotting to be done  outside of
C                     PROFILE  using  QPLOT.  Formerly,  much  of  the
C                     program  was  devoted  to  plotting  the airfoil
C                     geometry with DISSPLA.
C  Dec. '84    DAS    Incorporated cheap estimates of Cp distributions
C                     using algorithm supplied by Ilan Kroo of the RAC
C                     Branch at NASA Ames.  Added Cl, Cm calculation.
C  01/24/85    LJC    Allowed for original and revised plots with MODE
C                     = 3.  Also modified to take advantage of QPLOT's
C                     new equal axis scaling option.
C  02/12/84    DAS    Added TRANSFORM option.
C  02/13/85    LJC    Added shifting of ordinates to RECTIFY mode.
C  02/19/85    DAS    Fixed REDISTRIB to handle rounded trailing edges
C                     properly; took out "BLUNT" input parameter.
C  02/28/85    LJC    Added 3-column format to PRREAD and PRWRIT.
C  03/22/85    LJC    Allowed for absent  input control file  (meaning
C                     an empty file is not needed). Also added two new
C                     keywords to RDKEYS for controlling line types on
C                     Cps and curvature plots.
C  06/17/85    DAS    Provided generalized distributions in REDISTRIB,
C                     via GETDIS/DSTRIB in place of XGRID.
C  09/05/85    DAS    Fixed minor error in QPLOTable airfoil data file.
C  09/30/85    DAS    Added COMBINE option.
C  10/09/85    DAS    Suppressed plot windowing control values from
C                     QPLOTable file if they are not in use.
C  10/21/85    DAS    Introduced LSTFIL to echo the input control file
C                     to the .TAB file for future reference.
C  11/04/85    DAS    Mixup with END CURVE/END FRAME for OPTIMIZE mode
C                     target curvatures crept in (how? used to be OK).
C  12/30/85    DAS    Added ROUND option (tentatively).
C  04/24/86    DAS    Cp at leading edge was wrong.
C  08/11/86    DAS    Added SMOOTH option (fitting of Wagner functions
C                     in linear least squares sense).
C  09/24/86    DAS    AFGEOM in place of CALCTHICK in main program;
C                     menu options by name now, as well as by number.
C  10/20/86    RAK    AFGEOM returns CAMBER and THICK as fractions; must
C                     convert to % chord at calling level.
C  02/11/87    DAS    Functionality of BOUNDS changed somewhat.
C  03/11/87    DAS    Made use of PROTECT utility (prompted by duplicate
C                     point case that dies ungracefully otherwise).
C  04/24/87    DAS    Menu items clarified and made more unique (because
C                     of SELECT's option for choosing by name or number).
C                     Greater precision in output coordinates now by default.
C                     Also: traps bad IDENT now (more than one token is
C                     presumed to be a goof).
C  04/27/87    DAS    Introduced PRECISION argument in PRWRIT.
C  08/31/87    DAS    Trailing edge angle added as AFGEOM output.
C  09/23/87    DAS    If MODE=0 but PRECISION is not the default, assume
C                     that a revised airfoil dataset is required.
C  04/29/88    DAS    Added spreadsheet-compatible output file.
C  12/01/88    DAS    Turned on spreadsheet file for MODE <= 3 now.
C  11/29/89  DAS/RGL  Plots for multi-section 3D cases with blank subtitle from
C                     the initial prompt now use the case titles as subtitles.
C                     'Original' suppressed from headers if 'Revised' does not
C                     apply.
C  02/14/90    DAS    Installed OPENER in 9 places.  File names are now in
C                     lower case to indulge the Unix community.
C  03/13/90     "     Raised MXPTS from 201 to 300.
C  06/28/90     "     Added the LOFT option.
C  06/29/90     "     Added the ROTATE option.
C  04/22/91     "     Indicated "current axes" for thickness/camber printout.
C  10/22/91     "     Replaced the ROUND option with the "nose-job" option
C                     (round or sharpen).  Providing for reversing the operation
C                     can mean letting through a non-rectified airfoil.
C  10/24/91     "     Provided for full wrap-around curvature plot (parametric
C                     derivatives).  Adjusted the tabulations accordingly.
C                     Added the control scheme description above (adapted from
C                     RDKEYS) and the geometry data description (adapted from
C                     PRREAD).
C  01/07/95     "     Original chord and leading edge are now tabulated;
C                     "normalize" mode *.crv plot file had mismatched X min/max.
C  06/13/95     "     Output of the *.ypp file is now as advertised above.
C  11/21/95     "     TRANSFORM mode now has an option to decamber a section.
C  12/19/96     "     Added SINF and four "COS" shape functions.
C  12/20/96     "     Extended SMOOTH option to allow y"-based implicit and/or
C                     explicit smoothing as an alternative to Wagner fitting.
C  10/20/99     "     Fortran 90 upgrade, mainly to eliminate use of '0' for
C                     carriage control.
C
C  AUTHORS:
C    Leslie Collins, David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C  ACKNOWLEDGMENTS:
C    Robert Kennelly (ex-Sterling, now NASA Ames) provided the inspiration and
C    many key ideas.
C
C    The Aerodynamics Division at NASA Ames funded this software under contract
C    to Sterling Software (previously known as Informatics, Inc.).
C
C-------------------------------------------------------------------------------


C     Declarations:
C     -------------

      IMPLICIT NONE

C  *  Constants:

      REAL, PARAMETER ::
     >   HALF = 0.5E+0, ONE = 1.E+0, UNDEF = 999.E+0, ZERO = 0.E+0

      CHARACTER, PARAMETER ::
     >   BLANK * 1      = ' ',
     >   CURVATURE * 10 = ' curvature',
     >   ENDCRV * 9     = 'END CURVE',
     >   ENDFRM * 9     = 'END FRAME',
     >   LO * 6         = ' lower',
     >   NEW * 3        = 'NEW',
     >   OLD * 3        = 'OLD',
     >   ORIGINAL * 9   = 'Original ',
     >   REVISED * 9    = 'Revised  ',
     >   UP * 6         = ' upper'

      INTEGER, PARAMETER ::
     >   LUNCPS=9, LUNCRT=6, LUNCRV=7,  LUNDAT=2, LUNINP=1,  LUNKBD=5,
     >   LUNPLT=8, LUNREV=3, LUNSPR=13, LUNTAB=4, LUNTBL=10, LUNXGR=11,
     >   LUNYPP=12, MAXPTS=500, MXCRVS=20, MXMENU=12, MXWRAP=2*MAXPTS+1

C     MAXPTS = Maximum number of points allowed for on one surface.
C     MXCRVS = Maximum number of curves allowed for on one plot.
C     MXMENU = Number of options on the main menu.
C     MXWRAP = 2*MAXPTS+1 (for estimation of force coefficients).

C  *  Variables:

      INTEGER
     >   CAMBCASE, FIRST, FORMAT, I, IER, INFORM, IPR, J, LAST, LUNIN,
     >   MARK, MODE, NINC, NL, NLORIG, NOSEMODE, NPR, NPTS, NU, NUORIG,
     >   NWRAP, NXTARG, PRECISION

      REAL
     >   ALFRAD, ALPHA, AREA, ARROW, CAMBER, CD, CENTRD, CHORD, CNORML,
     >   COEFS (MAXPTS, 3), CL, CM, CPB, CPT, CPBOT, CPTOP, CPSTAR,
     >   CPL (MAXPTS), CPU (MAXPTS), CRVMAX, CRVMIN, CRVNEW (MXWRAP),
     >   CRVORIG (MXWRAP), DUMMY, FSMACH, GAMMA, GOGM1,
     >   GP1GM1, MOMENT, PWRAP (MXWRAP), SCALE, TARGCRV (MAXPTS),
     >   TARGX (MAXPTS), TEANGLE, THICK, TOGM1,
     >   X (MAXPTS*2), XAXIS, XC (MAXPTS), XCAM,
     >   XEVAL (MAXPTS), XL (MAXPTS), XLE, XLEFT, XLORIG (MAXPTS),
     >   XMAXSV, XMAX, XMIN, XMINSV, XP (MXWRAP), XPP (MXWRAP), XRIGHT,
     >   XTH, XU (MAXPTS), XUORIG (MAXPTS), XWRAP (MXWRAP),
     >   Y (MAXPTS*2), YEVAL (MAXPTS), YKL (MAXPTS), YKLORG (MAXPTS),
     >   YKU (MAXPTS), YKUORG (MAXPTS), YL (MAXPTS), YLE,
     >   YLORIG (MAXPTS), YMAX, YMAXSV, YMIN, YMINSV, YP (MXWRAP),
     >   YPL (MAXPTS), YPP (MXWRAP), YPPL (MAXPTS), YPPU (MAXPTS),
     >   YPU (MAXPTS), YU (MAXPTS), YUORIG (MAXPTS), YWRAP (MXWRAP)

      LOGICAL
     >   CPSFIL, CRVFIL, DATFIL, DEFAULT, DISTINCT, NOMODS, NOSHOW,
     >   PLTARG, PLTFIL, PLTORIG, PLTREV, QUIT, SAVREV, SPREADFIL,
     >   TABFIL, THREED, UPPER, WRAPCRV, YPPFIL

      CHARACTER
     >   CHOICE * 16, CPSLINE (MXCRVS) * 9, CRVLINE (MXCRVS) * 9,
     >   DATAFILE * 40, ENDSTR * 9, FSTATUS * 9, HTAB * 1, IDENT * 15,
     >   INPUTFILE * 32, LEG * 8, LEGEND * 80, LEGND (MXCRVS) * 80,
     >   MENU (0:MXMENU) * 60, MNTITL * 80, PLTLINE (MXCRVS) * 9,
     >   SBTITL * 80, SBTIT2 * 80, SUBTIT (2) * 24, SURFACE (2) * 13,
     >   XLABEL * 3, YLABEL * 3

C  *  Procedures:

      EXTERNAL
     >   AFGEOM, BOUNDS, COMBINE, COPY, CURV2D, FD12K, GETCL, GETCPS,
     >   LOFT, LSTFIL, MAXMIN, MODIFY, NOSEJOB, NRMLIZ, NRMSET, OPENER,
     >   OPTIMIZE, PROTECT, PRREAD, PRTAB, PRWRIT, PSFIT, PSTVAL,
     >   QPLDAT, RDKEYS, READR, READS, READY, RECTIFY, REDISTRIB,
     >   REFINE, RVERSE, ROTATE, SCAN2, SELECT, SMOOTH, TRANSFORM,
     >   XDERIVS, XFORM

C  *  Storage:

      DATA
     >   SURFACE /'upper surface', 'lower surface'/

      DATA
     >   MENU /
     >   '  0 = DISPLAY (plot/tabulate only or alter format/precision)',
     >   '  1 = RECTIFY leading edge definition; allows vertical shift',
     >   '  2 = NORMALIZE or denormalize coordinates',
     >   '  3 = REDISTRIBUTE the abscissas',
     >   '  4 = MODIFY either surface or both (apply shape functions)',
     >   '  5 = REFINE thickness/curvature',
     >   '  6 = OPTIMIZE one surface ("bumps" + target curvatures)',
     >   '  7 = TRANSFORM YU/YL to/from camber/thickness, or decamber',
     >   '  8 = ROTATE coordinates, with option to renormalize',
     >   '  9 = COMBINE option (add or subtract profiles)',
     >   ' 10 = LOFT linearly between primary and secondary profiles',
     >   ' 11 = NOSE-JOB option: round or sharpen the leading edge',
     >   ' 12 = SMOOTH YU and/or YL: Wagner fn. fits or [im|ex]plicit' /


C     Execution:
C     ----------

      HTAB = CHAR (9)            ! Used to be a PARAMETER (non-standard).

C  *  Select mode of operation:

      WRITE (LUNCRT, 1000)
     >   BLANK,
     >   'Welcome to PROFILE from the Aerodynamics Division, NASA Ames.'

      MODE = 0
      CHOICE = 'DISPLAY'
      NOSHOW = .FALSE.

      CALL SELECT ('Select operating mode.', MXMENU+1, MENU, NOSHOW,
     >             LUNCRT, LUNKBD, MODE, CHOICE, QUIT)
      IF (QUIT) GO TO 810

C-----------------------------------------------------------------------
C         Handle input files (control file; airfoil data file).
C-----------------------------------------------------------------------

      DATAFILE = 'naca0012.dat'
      CALL OPENER (LUNCRT,
     >   'Input airfoil file? <CR>=naca0012.dat: ',
     >   LUNKBD, DATAFILE, LUNDAT, OLD)

      INPUTFILE = 'profile.inp'
      FSTATUS = 'IfPresent'
      LUNIN = LUNINP
      CALL OPENER (LUNCRT,
     >   'Input control file? <CR>=profile.inp or none: ',
     >   LUNKBD, INPUTFILE, LUNIN, FSTATUS)
      IF (FSTATUS == 'MISSING') THEN
         LUNIN = -LUNIN
         INPUTFILE = BLANK
         WRITE (LUNCRT, 1000)
     >      'No control file found - proceeding ...'
      END IF

C  *  Read input control variables (if any):

      CALL RDKEYS (LUNCRT, LUNIN, MXCRVS, UNDEF, FORMAT, PRECISION,
     >   PLTLINE, CPSLINE, CRVLINE, CRVMAX, CRVMIN, WRAPCRV,
     >   PLTORIG, PLTREV, THREED, XAXIS, XMIN, XMAX, YMIN, YMAX,
     >   DATFIL, PLTFIL, CRVFIL, TABFIL, YPPFIL, CPSFIL, SPREADFIL)

C-----------------------------------------------------------------------
C             Set up or adjust the internal control variables.
C-----------------------------------------------------------------------

      TABFIL  = TABFIL .OR.  MODE >= 4
      PLTORIG = PLTFIL .AND. PLTORIG .AND. (MODE == 0 .OR. MODE >= 3)
      PLTREV  = PLTFIL .AND. PLTREV .AND. MODE > 0
      PLTFIL  = PLTFIL .AND. (PLTORIG .OR. PLTREV)
      PLTARG  = CRVFIL .AND. MODE == 6
      CRVFIL  = CRVFIL .AND. MODE /= 7
      CPSFIL  = CPSFIL .AND. MODE /= 7
      SPREADFIL = SPREADFIL .AND. MODE <= 3
      SAVREV  = DATFIL .AND.
     >          (MODE > 0 .OR. FORMAT > 0 .OR. PRECISION /= 1)
      NOMODS  = MODE == 0 .OR. MODE == 7

      IF (WRAPCRV .AND. YPPFIL) THEN
         WRITE (LUNCRT, 1080)
         YPPFIL = .FALSE.
      END IF

      CAMBCASE = 1  ! Transform from Ys to camber/thickness

C  *  CAMBCASE applies to MODE = 7, but having it undefined for other
C     modes is inconvenient for deciding whether to invoke AFGEOM.


C-----------------------------------------------------------------------
C                        Set up the output files.
C-----------------------------------------------------------------------

C  *  Prompt for the identifier to be used for all files.  Figure that
C     more than one token, or something like 'naca0012.dat' is a goof.

  150 CONTINUE
      IDENT = 'profile'
      CALL READS (LUNCRT,
     >   'Identifier for output files? <CR>="profile": ',
     >   LUNKBD, IDENT, DEFAULT, QUIT)
      IF (QUIT) GO TO 810

      FIRST = 1
      LAST = LEN (IDENT)
      CALL SCAN2 (IDENT, ' .;/', FIRST, LAST, MARK)
      IF (LAST /= MARK) GO TO 150

C     Valid identifier may still cause a conflict with the input file.
C     (Unix systems in mind, here.)

      IF (SAVREV) THEN
         IF (IDENT (FIRST:LAST) == DATAFILE (FIRST:LAST)) THEN
            WRITE (LUNCRT, 1000)
     >         'In/out file name conflict.  Choose another identifier.'
            GO TO 150
         END IF
      END IF

      IF (SAVREV)
     >   CALL OPENER (LUNCRT, BLANK, LUNKBD,
     >                IDENT (FIRST:LAST) // '.dat', LUNREV, NEW)
      IF (PLTFIL)
     >   CALL OPENER (LUNCRT, BLANK, LUNKBD,
     >                IDENT (FIRST:LAST) // '.plt', LUNPLT, NEW)
      IF (TABFIL)
     >   CALL OPENER (LUNCRT, BLANK, LUNKBD,
     >                IDENT (FIRST:LAST) // '.tab', LUNTAB, NEW)
      IF (CRVFIL)
     >   CALL OPENER (LUNCRT, BLANK, LUNKBD,
     >                IDENT (FIRST:LAST) // '.crv', LUNCRV, NEW)
      IF (YPPFIL)
     >   CALL OPENER (LUNCRT, BLANK, LUNKBD,
     >                IDENT (FIRST:LAST) // '.ypp', LUNYPP, NEW)
      IF (CPSFIL)
     >   CALL OPENER (LUNCRT, BLANK, LUNKBD,
     >                IDENT (FIRST:LAST) // '.cps', LUNCPS, NEW)
      IF (SPREADFIL)
     >   CALL OPENER (LUNCRT, BLANK, LUNKBD,
     >                IDENT (FIRST:LAST) // '.spread', LUNSPR,
     >                NEW // ':153') ! Raise max. record length

C  *  Prompt for plot title and subtitle:

      MNTITL = 'PROFILE'
      CALL READS (LUNCRT,
     >   'Plot title line?  <CR> uses "PROFILE" for the title:',
     >   LUNKBD, MNTITL, DEFAULT, QUIT)
      IF (QUIT) GO TO 810

      SBTITL = BLANK
      CALL READS (LUNCRT,
     >   'Plot subtitle line?  <CR> means none:',
     >   LUNKBD, SBTITL, DEFAULT, QUIT)
      IF (QUIT) GO TO 810

C-----------------------------------------------------------------------
C                            Scan all datasets.
C-----------------------------------------------------------------------

C  *  Scan all profiles for overall data range (for axis scaling and
C     possible normalization purposes).  Count number of profiles.

      CALL MAXMIN (LUNDAT, MAXPTS, XMINSV, XMAXSV, YMINSV, YMAXSV,
     >             YLE, NPR, X(1), Y(1), XU, XL, YU, YL, LEGND (1), IER)

      IF (IER == 4) GO TO 900
      IF (IER == 3) GO TO 910
      IF (IER == 2) GO TO 920
      IF (NPR == 0) GO TO 930

C  *  Make sure both surfaces have the same leading edge point, unless
C     this is the rectify option:

      IF ((MODE /= 1) .AND.
     >    (XU (1) /= XL (1) .OR. YU (1) /= YL (1))) GO TO 850

      REWIND LUNDAT

      IF (TABFIL) THEN

         WRITE (LUNTAB, 1010)
     >      ' Program PROFILE (Applied Aerodynamics Branch, NASA Ames',
     >      ' Research Center)'
         WRITE (LUNTAB, 1020) ' Operating mode: ', MODE, ' - ', CHOICE
         WRITE (LUNTAB, 1010) ' Control file  : ', INPUTFILE
         WRITE (LUNTAB, 1010) ' Geometry file : ', DATAFILE
         WRITE (LUNTAB, 1010) ' Output file id: ', IDENT (FIRST:LAST)
         WRITE (LUNTAB, 1025) ' Profiles found: ', NPR

         IF (INPUTFILE /= BLANK) THEN

C  *        Echo control keywords to tabulation file.
C           LEGEND is as good as anything for use as a buffer.

            WRITE (LUNTAB, 1000) BLANK, 'Control file contents:', BLANK

            REWIND LUNINP
            CALL LSTFIL (LUNINP, LUNTAB, LEGEND)

         ELSE
            WRITE (LUNTAB, 1020)
     >         ' Control file was absent - all defaults taken.'
         END IF

      END IF

      IF (CPSFIL) THEN

C  *     Don't attempt to estimate Cps for more than one profile:

         CPSFIL = NPR == 1 .AND. .NOT. THREED

C  *     Save user from editing input control file too often:

         IF (CPSFIL) THEN
            CPSFIL = .FALSE.
            CALL READY (LUNCRT,
     >         'Do you really want Cp estimates? (Y/N; <CR>=No) ',
     >         LUNKBD, CPSFIL, DEFAULT, QUIT)
         END IF

         IF (.NOT. CPSFIL) THEN
            CLOSE (UNIT=LUNCPS, STATUS='DELETE')
         ELSE
            ALPHA = 0.
            CALL READR (LUNCRT,
     >         'Enter Alpha, else <CR> gives 0.: ',
     >         LUNKBD, ALPHA, DEFAULT, QUIT)
            IF (QUIT) GO TO 810

            FSMACH = 0.
            CALL READR (LUNCRT,
     >         'Enter free stream Mach number, else <CR> gives 0.: ',
     >         LUNKBD, FSMACH, DEFAULT, QUIT)
            IF (QUIT) GO TO 810

         END IF
      END IF

C-----------------------------------------------------------------------
C                     Begin main loop over all profiles:
C-----------------------------------------------------------------------

      DO IPR = 1, NPR

C  *     Read one profile.  Ignore error processing this time through.

         CALL PRREAD (LUNDAT, LEGND (IPR), MAXPTS, NU, NL, X, Y,
     >                XU, XL, YU, YL, INFORM, IER)

C  *     Normally no use in proceeding if abscissas are not monotonic:

         CALL PROTECT (NU, XU, YU, ARROW, DISTINCT)
         J = 1
         IF (ARROW /= ONE) THEN
            IF (MODE == 11) THEN     ! May need to reverse a round/sharpen
               WRITE (LUNCRT, 1070)
            ELSE IF (MODE /= 1) THEN
               GO TO 860
            END IF
         END IF

         CALL PROTECT (NL, XL, YL, ARROW, DISTINCT)
         J = 2
         IF (ARROW /= ONE) THEN
            IF (MODE == 11) THEN
               WRITE (LUNCRT, 1070)
            ELSE IF (MODE /= 1) THEN
               GO TO 860
            END IF
         END IF

C  *     Save original data before any modifications:

         CALL COPY (NU, XU, XUORIG)
         CALL COPY (NL, XL, XLORIG)
         CALL COPY (NU, YU, YUORIG)
         CALL COPY (NL, YL, YLORIG)
         NUORIG = NU
         NLORIG = NL

C-----------------------------------------------------------------------
C                           Major options.
C-----------------------------------------------------------------------

         IF (MODE == 1) THEN

C  *        Rearrange upper and lower surface points so that the
C           minimum abscissa is the leading edge point. Shift ordinates
C           by input value if required:

            CALL RECTIFY (NU, XU, YU, NL, XL, YL, X, Y, MAXPTS,
     >                    LUNCRT, LUNKBD)

         ELSE IF (MODE == 2) THEN

            IF (IPR == 1) THEN

C  *           Prompt for chord, etc., to normalize/denormalize by:

               CALL NRMSET (LUNCRT, LUNKBD, XMINSV, XMAXSV,
     >                      CNORML, XLE, YLE)
            END IF

            CALL NRMLIZ (NU, XU, XU, XLE, CNORML)
            CALL NRMLIZ (NL, XL, XL, XLE, CNORML)
            CALL NRMLIZ (NU, YU, YU, YLE, CNORML)
            CALL NRMLIZ (NL, YL, YL, YLE, CNORML)

         ELSE IF (MODE == 3) THEN

C  *        Redistribute abscissas using splines:

            CALL REDISTRIB (NU, XU, YU, NL, XL, YL, X, Y,
     >                      MAXPTS, XEVAL, YEVAL, COEFS,
     >                      LUNCRT, LUNKBD, LUNXGR)

         ELSE IF (MODE == 4) THEN

C  *        Modify profile by applying shape functions interactively:

            CALL MODIFY (NU, XU, YU, NL, XL, YL, LUNCRT, LUNKBD,
     >                   LUNTAB)

         ELSE IF (MODE == 5) THEN

C  *        Refine thickness and/or curvature interactively:

            CALL REFINE (NU, XU, YU, NL, XL, YL, MNTITL,
     >                   LUNCRT, LUNKBD, LUNTAB, LUNTBL)

         ELSE IF (MODE == 6) THEN

C  *        Optimize one surface using bumps + target curvatures:

            CALL OPTIMIZE (NU, XU, YU, NL, XL, YL, MNTITL,
     >                     LUNCRT, LUNKBD, LUNTAB, LUNXGR, LUNTBL,
     >                     UPPER, NXTARG, TARGX, TARGCRV)

         ELSE IF (MODE == 7) THEN

C  *        Transform representation (upper/lower <-> camber/thickness)
C           or zero out the camber.  CAMBCASE is an output.

            CALL TRANSFORM (NU, XU, YU, NL, XL, YL, LUNCRT, LUNKBD,
     >                      CAMBCASE)

         ELSE IF (MODE == 8) THEN

C  *        Rotate coordinates, with option to renormalize:

            CALL ROTATE (NU, XU, YU, NL, XL, YL, LUNCRT, LUNKBD)

         ELSE IF (MODE == 9) THEN

C  *        Combine present profile with another to be prompted for.
C           Use YPU, YPPU as scratch for XU2ND, YU2ND, etc.

            CALL COMBINE (NU, XU, YU, NL, XL, YL, MAXPTS, X, Y,
     >                    YPU, YPPU, YPL, YPPL,
     >                    LUNCRT, LUNKBD, LUNXGR)

         ELSE IF (MODE == 10) THEN

C  *        Loft linearly between present profile and another to be prompted
C           for.  Use YPU, YPPU as scratch for XU2ND, YU2ND, etc.

            CALL LOFT (NU, XU, YU, NL, XL, YL, MAXPTS, X, Y,
     >                 YPU, YPPU, YPL, YPPL,
     >                 LUNCRT, LUNKBD, LUNXGR)

         ELSE IF (MODE == 11) THEN

C  *        Round or sharpen the leading edge region:

            CALL NOSEJOB (NU, XU, YU, NL, XL, YL, MAXPTS, X, Y, COEFS,
     >                    LUNCRT, LUNKBD, LUNTAB, NOSEMODE)

         ELSE IF (MODE == 12) THEN

C  *        Smooth YU and/or YL by fitting Wagner shape functions or
C           via implicit and/or explicit y"-based smoothing:

            CALL SMOOTH (NU, XU, YU, NL, XL, YL, LUNCRT, LUNKBD, LUNTAB)

         END IF


C-----------------------------------------------------------------------
C                         Generate the output files.
C-----------------------------------------------------------------------


C  . . . . . . . . . Revised profile coordinates (*.DAT) . . . . . . . .

         IF (SAVREV) THEN

            IF (MODE == 7) THEN
               IF (CAMBCASE == 1) THEN
                  FORMAT = 6
               END IF
            END IF

            IF (FORMAT == 0) FORMAT = INFORM

            IF (MODE == 0) THEN
               CALL PRWRIT (LUNREV, MAXPTS, LEGND(IPR), NU, NL,
     >                      XU, XL, YU, YL, FORMAT, PRECISION)
            ELSE
               CALL PRWRIT (LUNREV, MAXPTS, MNTITL, NU, NL,
     >                      XU, XL, YU, YL, FORMAT, PRECISION)
            END IF

         END IF


C  . . . . . . . . . .  Plottable coordinates (*.PLT) . . . . . . . . . .

C  *     Multi-dataset plot frames should be self-descriptive:

         IF (NPR > 1 .AND. SBTITL == BLANK) THEN
            SBTIT2 = LEGND (IPR)
         ELSE
            SBTIT2 = SBTITL
         END IF

         IF (ABS (XU (1)  - ZERO) > 1.E-3  .OR.
     >       ABS (XU (NU) - ONE)  > 1.E-3) THEN
            XLABEL = ' x '
            YLABEL = ' y '
         ELSE
            XLABEL = 'x/c'
            YLABEL = 'y/c'
         END IF

         IF (PLTORIG) THEN

C  *        Save original airfoil coordinates in QPLOTable form. First store
C           coordinates in wrap-around format, avoiding duplicate leading edge
C           point:

            CALL RVERSE (NLORIG, XLORIG, X)
            CALL RVERSE (NLORIG, YLORIG, Y)
            CALL COPY   (NUORIG, XUORIG, X (NLORIG))
            CALL COPY   (NUORIG, YUORIG, Y (NLORIG))
            NPTS = NLORIG + NUORIG - 1

            IF (IPR == 1 .OR. THREED)
     >         CALL QPLFRM  (LUNPLT, MNTITL, SBTIT2, XLABEL, YLABEL)

            IF (.NOT. PLTREV  .AND.  IPR == NPR  .OR.  THREED) THEN
               ENDSTR = ENDFRM
            ELSE
               ENDSTR = ENDCRV
            END IF

            CALL QPLCRV (LUNPLT, NPTS, X, Y, ENDSTR, UNDEF,
     >                   'SCALE', XMIN, XMAX, YMIN, YMAX, XAXIS, UNDEF,
     >                   LEGND (IPR), PLTLINE (IPR))

         END IF

         IF (PLTREV) THEN

C  *        Now for the revised airfoil:

            CALL RVERSE (NL, XL, X)
            CALL RVERSE (NL, YL, Y)
            CALL COPY   (NU, XU, X (NL))
            CALL COPY   (NU, YU, Y (NL))
            NPTS = NL + NU - 1

            IF (.NOT. PLTORIG)
     >         CALL QPLFRM (LUNPLT, MNTITL, SBTIT2, XLABEL, YLABEL)

            CALL QPLCRV (LUNPLT, NPTS, X, Y, ENDFRM, UNDEF, 'SCALE',
     >                   XMIN, XMAX, YMIN, YMAX, XAXIS, UNDEF,
     >                   'Revised profile', PLTLINE (2))
         END IF


C  . . . . . . . . . . .  Basic tabulations (*.TAB)  . . . . . . . . . .
C                        and 2nd derivatives (*.YPP)
C                        Curvature is also required.

C  *     Don't try to calculate derivatives for a non-rectified airfoil:

         IF (MODE /= 1) THEN

            IF (.NOT. WRAPCRV) THEN

C  *           Derivatives and curvatures of original profile, using Y vs. X:

               CALL FD12K (NUORIG, XUORIG, YUORIG, YPU, YPPU, YKUORG)
               CALL FD12K (NLORIG, XLORIG, YLORIG, YPL, YPPL, YKLORG)

            ELSE

C  *           Derivatives and curvature via X vs. T, Y vs. T parametric form:

               CALL RVERSE (NLORIG, XLORIG, X)
               CALL RVERSE (NLORIG, YLORIG, Y)
               CALL COPY   (NUORIG, XUORIG, X (NLORIG))
               CALL COPY   (NUORIG, YUORIG, Y (NLORIG))
               NPTS = NLORIG + NUORIG - 1

C              Fit a conventional parametric spline (coefs. stored internally):

               CALL PSFIT (NPTS, X, Y, 'C', .FALSE., IER)
               IF (IER /= 0) GO TO 940

C              -NPTS tells PSTVAL to evaluate at the data points,
C              so there's no need to input the Ts.
C              Reuse X/YWRAP for the (unneeded) X/Y outputs.

               CALL PSTVAL (-NPTS, DUMMY, XWRAP, YWRAP, XP, YP,
     >                      XPP, YPP, NPTS, X, Y)

C              Derive the curvature from the derivatives.

               CALL CURV2D (NPTS, XP, XPP, YP, YPP, CRVORIG)

C  *           Set up the curvature and X derivatives in separate-surface
C              form for tabulation.  First, the lower surface:

               SCALE = XLORIG (NLORIG) - XLORIG (1)  ! For dx/dt ~ 0 test.

               CALL XDERIVS (NLORIG, XP, XPP, YP, YPP, YPL, YPPL, SCALE)
               CALL RVERSE  (NLORIG, YPL, YPL)
               CALL RVERSE  (NLORIG, YPPL, YPPL)
               CALL RVERSE  (NLORIG, CRVORIG, YKLORG)

C              Now the upper surface:

               I = NLORIG
               CALL XDERIVS (NUORIG, XP (I), XPP (I), YP (I), YPP (I),
     >                       YPU, YPPU, SCALE)
               CALL COPY (NUORIG, CRVORIG (I), YKUORG)

            END IF

            IF (TABFIL) THEN

C  *           Tabulate original profile with derivatives and curvatures:

               IF (NOMODS) THEN        ! Suppress 'original'
                  SUBTIT (1) = SURFACE (1)
                  SUBTIT (1) (1:1) = 'U'   ! Not 'u'
                  SUBTIT (2) = SURFACE (2)
                  SUBTIT (2) (1:1) = 'L'
               ELSE
                  SUBTIT (1) = ORIGINAL // SURFACE (1)
                  SUBTIT (2) = ORIGINAL // SURFACE (2)
               END IF

               IF (MODE == 7) THEN
                  IF (CAMBCASE == 2) THEN
                     SUBTIT (1) = 'Mean line'
                     SUBTIT (2) = 'Semi-thickness'
                  END IF
               END IF

               CALL PRTAB (LUNTAB, LEGND (IPR), SUBTIT (1), WRAPCRV,
     >                     NUORIG, XUORIG, YUORIG, YPU, YPPU, YKUORG)
               CALL PRTAB (LUNTAB, LEGND (IPR), SUBTIT (2), WRAPCRV,
     >                     NLORIG, XLORIG, YLORIG, YPL, YPPL, YKLORG)
            END IF

C  *        Calculate geometric properties, unless original profile was
C           camber/thickness:

            IF (.NOT. (MODE == 7 .AND. CAMBCASE == 2)) THEN

C  *           AFGEOM needs upper/lower surface in (*,2) form:

               CALL COPY (NUORIG, XUORIG, X (1))
               CALL COPY (NLORIG, XLORIG, X (1+MAXPTS))
               CALL COPY (NUORIG, YUORIG, Y (1))
               CALL COPY (NLORIG, YLORIG, Y (1+MAXPTS))

C              XEVAL(*) is handy for AFGEOM's YPOWER(*) work-space.

               CALL AFGEOM (NUORIG, NLORIG, MAXPTS, X, Y, AREA, CENTRD,
     >                      MOMENT, THICK, XTH, CAMBER, XCAM, TEANGLE,
     >                      XEVAL, YEVAL, COEFS (1,1), COEFS (1,2),
     >                      COEFS (1,3), IER)

               IF (NOMODS) THEN
                  LEG = BLANK
               ELSE
                  LEG = ORIGINAL
               END IF

               DO J = LUNTAB, LUNCRT, LUNCRT - LUNTAB

C  *              Don't clutter screen with original - just revised, if any:

                  IF ((J == LUNCRT .AND. NOMODS) .OR.
     >                (J == LUNTAB .AND. TABFIL)) THEN

                     IF (IER /= 0) WRITE (J, 1055) ORIGINAL, IER

                     WRITE (J, 1050) LEG, THICK * 100., XTH,
     >                               CAMBER * 100., XCAM, AREA, MOMENT
                  END IF

                  IF (J == LUNTAB .AND. TABFIL) THEN
                     CHORD = MAX (XUORIG (NUORIG), XLORIG (NLORIG)) -
     >                       XUORIG (1)
                     WRITE (J, 1060) XUORIG (1), YUORIG (1), CHORD,
     >                               TEANGLE
                  END IF

               END DO

            END IF

         END IF

         IF (YPPFIL .AND. NOMODS) THEN

C  *        Save 2nd derivatives in standard PROFILE format for possible
C           reuse by REFINE mode.  End values are not needed.

            CALL PRWRIT (LUNYPP, MAXPTS, LEGND (IPR), NUORIG-2,
     >                   NLORIG-2, XUORIG (2), XLORIG (2), YPPU (2),
     >                   YPPL (2), 5, 3)
         END IF

C  *     Repeat for the revised airfoil:
C
         IF (MODE > 0) THEN

            IF (.NOT. WRAPCRV) THEN

C  *           Derivatives and curvatures of revised profile, using Y vs. X:

               CALL FD12K (NU, XU, YU, YPU, YPPU, YKU)
               CALL FD12K (NL, XL, YL, YPL, YPPL, YKL)

            ELSE

C  *           Derivatives and curvature via X vs. T, Y vs. T parametric form:

               CALL RVERSE (NL, XL, X)
               CALL RVERSE (NL, YL, Y)
               CALL COPY   (NU, XU, X (NL))
               CALL COPY   (NU, YU, Y (NL))
               NPTS = NL + NU - 1

               CALL PSFIT (NPTS, X, Y, 'C', .FALSE., IER)
               IF (IER /= 0) GO TO 940

               CALL PSTVAL (-NPTS, DUMMY, XWRAP, YWRAP, XP, YP,
     >                      XPP, YPP, NPTS, X, Y)

               CALL CURV2D (NPTS, XP, XPP, YP, YPP, CRVNEW)

               SCALE = XL (NL) - XL (1)

               CALL XDERIVS (NL, XP, XPP, YP, YPP, YPL, YPPL, SCALE)
               CALL RVERSE  (NL, YPL, YPL)
               CALL RVERSE  (NL, YPPL, YPPL)
               CALL RVERSE  (NL, CRVNEW, YKL)

C              Now the upper surface:

               I = NL
               CALL XDERIVS (NU, XP (I), XPP (I), YP (I), YPP (I),
     >                       YPU, YPPU, SCALE)
               CALL COPY (NU, CRVNEW (I), YKU)

            END IF

            IF (TABFIL) THEN
               SUBTIT (1) = REVISED // SURFACE (1)
               SUBTIT (2) = REVISED // SURFACE (2)

               IF (MODE == 7) THEN
                  IF (CAMBCASE == 1) THEN
                     SUBTIT (1) = 'Mean line'
                     SUBTIT (2) = 'Semi-thickness'
                  END IF
               END IF

               CALL PRTAB (LUNTAB, MNTITL, SUBTIT (1), WRAPCRV,
     >                     NU, XU, YU, YPU, YPPU, YKU)
               CALL PRTAB (LUNTAB, MNTITL, SUBTIT (2), WRAPCRV,
     >                     NL, XL, YL, YPL, YPPL, YKL)
            END IF

C  *        Calculate geometric properties, unless revised profile is
C           camber/thickness:

            IF (.NOT. (MODE == 7 .AND. CAMBCASE == 1)) THEN

               CALL COPY (NU, XU, X (1))
               CALL COPY (NL, XL, X (1+MAXPTS))
               CALL COPY (NU, YU, Y (1))
               CALL COPY (NL, YL, Y (1+MAXPTS))

               CALL AFGEOM (NU, NL, MAXPTS, X, Y, AREA, CENTRD,
     >                      MOMENT, THICK, XTH, CAMBER, XCAM, TEANGLE,
     >                      XEVAL, YEVAL, COEFS (1,1), COEFS (1,2),
     >                      COEFS (1,3), IER)

               DO J = LUNTAB, LUNCRT, LUNCRT - LUNTAB

                  IF (.NOT. (J == LUNTAB .AND. .NOT. TABFIL)) THEN
                     IF (IER /= 0) WRITE (J, 1055) REVISED, IER
                     WRITE (J, 1050) REVISED, THICK * 100., XTH,
     >                               CAMBER * 100., XCAM, AREA, MOMENT
                  END IF

                  IF (J == LUNTAB .AND. TABFIL) THEN
                     CHORD = MAX (XU (NU), XL (NL)) - XU (1)
                     WRITE (J, 1060) XU (1), YU (1), CHORD, TEANGLE
                  END IF

               END DO

            END IF

            IF (YPPFIL) THEN

C  *           Save 2nd derivatives in standard PROFILE format for
C              possible reuse.

               CALL PRWRIT (LUNYPP, MAXPTS, MNTITL, NU-2, NL-2,
     >                      XU (2), XL (2), YPPU (2), YPPL (2), 5, 3)
            END IF
         END IF


C  . . . . . . . . . . . Curvature data file (*.CRV) . . . . . . . . . .

         IF (CRVFIL) THEN

C  *        Save original, revised and target curvature data (if any).
C           If this is the RECTIFY or NORMALIZE option, only revised
C           curvature data are saved.
C
            IF (MODE == 0) THEN
               LEG = BLANK
               ENDSTR = ENDFRM
            ELSE
               LEG = ORIGINAL
               ENDSTR = ENDCRV
            END IF

C  *        Normally, the two surfaces appear on separate frames.

            IF (MODE > 0) THEN   ! Fix the common normalize case
               XLEFT  = XMINSV
               XRIGHT = XMAXSV
               CHORD  = MAX (XU (NU), XL (NL)) - XU (1)

               IF (CHORD == ONE) THEN
                  XLEFT  = ZERO
                  XRIGHT = ONE
               END IF
            END IF

            IF (.NOT. WRAPCRV) THEN

C  *           First, the upper surface:

               CALL QPLFRM (LUNCRV, MNTITL, SBTIT2, XLABEL,
     >                      SURFACE (1) // CURVATURE)

               IF (MODE /= 1 .AND. MODE /= 2)
     >            CALL QPLCRV (LUNCRV, NUORIG-2, XUORIG (2), YKUORG (2),
     >                         ENDSTR, UNDEF, BLANK, XMINSV, XMAXSV,
     >                         CRVMIN, CRVMAX, UNDEF, UNDEF, LEG,
     >                         CRVLINE (1))

               IF (PLTARG .AND. UPPER)
     >            CALL QPLCRV (LUNCRV, NXTARG, TARGX, TARGCRV, ENDCRV,
     >                         UNDEF, BLANK, XMINSV, XMAXSV, CRVMIN,
     >                         CRVMAX, UNDEF, UNDEF, 'Target',
     >                         CRVLINE (3))

               IF (MODE > 0)
     >            CALL QPLCRV (LUNCRV, NU-2, XU (2), YKU (2), ENDFRM,
     >                         UNDEF, BLANK, XLEFT, XRIGHT, CRVMIN,
     >                         CRVMAX, UNDEF, UNDEF, REVISED,
     >                         CRVLINE (2))

C  *           Now for the lower surface:

               CALL QPLFRM (LUNCRV, MNTITL, SBTIT2, XLABEL,
     >                      SURFACE (2) // CURVATURE)

               IF (MODE /= 1 .AND. MODE /= 2)
     >            CALL QPLCRV (LUNCRV, NLORIG-2, XLORIG (2), YKLORG (2),
     >                         ENDSTR, UNDEF, BLANK, XMINSV, XMAXSV,
     >                         CRVMIN, CRVMAX, UNDEF, UNDEF, LEG,
     >                         CRVLINE (1))

               IF (PLTARG .AND. .NOT. UPPER)
     >            CALL QPLCRV (LUNCRV, NXTARG, TARGX, TARGCRV, ENDCRV,
     >                      UNDEF, BLANK, XMINSV, XMAXSV, CRVMIN,
     >                      CRVMAX, UNDEF, UNDEF, 'Target', CRVLINE (3))

               IF (MODE > 0)
     >            CALL QPLCRV (LUNCRV, NL-2, XL (2), YKL (2), ENDFRM,
     >                         UNDEF, BLANK, XLEFT, XRIGHT, CRVMIN,
     >                         CRVMAX, UNDEF, UNDEF, REVISED,
     >                         CRVLINE (2))

            ELSE

C  *           Wrap-around form of curvature requested:

               CALL QPLFRM (LUNCRV, MNTITL, SBTIT2, XLABEL,
     >                      'Curvature')

               IF (MODE /= 1 .AND. MODE /= 2) THEN

                  CALL RVERSE (NLORIG, XLORIG, X)
                  CALL COPY   (NUORIG, XUORIG, X (NLORIG))
                  NPTS = NLORIG + NUORIG - 1

                  WRITE (LUNCRV, 1090)
                  CALL QPLCRV (LUNCRV, NPTS, X, CRVORIG, ENDSTR, UNDEF,
     >                         BLANK, XMINSV, XMAXSV, CRVMIN, CRVMAX,
     >                         UNDEF, UNDEF, LEG, CRVLINE (1))
               END IF

               IF (MODE > 0) THEN     ! As above, for the revised profile.

                  CALL RVERSE (NL, XL, X)
                  CALL COPY   (NU, XU, X (NL))
                  NPTS = NL + NU - 1

                  WRITE (LUNCRV, 1090)
                  CALL QPLCRV (LUNCRV, NPTS, X, CRVNEW, ENDSTR, UNDEF,
     >                         BLANK, XLEFT, XRIGHT, CRVMIN, CRVMAX,
     >                         UNDEF, UNDEF, REVISED, CRVLINE (2))
                END IF
             END IF
          END IF


C  . . . . . . . .  Cp estimates in plottable form (*.CPS) . . . . . . .

         IF (CPSFIL) THEN

C  *        Generate approximate Cp distributions for original and/or
C           revised profile.
C           If this is the RECTIFY or NORMALIZE option, only the revised
C           profile is handled.

C  *        Q:  Why not modularize this to unclutter the main program?
C           A:  Too many upper/lower/original/revised arrays involved.

C  *        First, check that abscissas are same for both surfaces, as
C           required by current GETCPS:

            IF (MODE == 0) THEN

C  *           Set "revised" abscissas so just one form of check will do:

               NU = NUORIG
               NL = NLORIG
               CALL COPY (NL, XLORIG, XL)
               CALL COPY (NU, XUORIG, XU)
            END IF

            DO I = 1, MIN (NU, NL)
               IF (XU (I) /= XL (I)) CPSFIL = .FALSE.
            END DO

C  *        Moment calculations use X/C=0.25, so...

            IF (XL (1) /= ZERO) CPSFIL = .FALSE.
            IF (XL (NL) /= ONE) CPSFIL = .FALSE.

            IF (CPSFIL) THEN

               ALFRAD = ALPHA * ATAN (ONE) / 45.E+0
               IF (FSMACH == ZERO) THEN
                  CPSTAR = -999.E+0
               ELSE
                  GAMMA  = 1.4E+0
                  GOGM1  = GAMMA / (GAMMA - ONE)
                  TOGM1  = 2.E+0 / (GAMMA - ONE)
                  GP1GM1 = (GAMMA + ONE) / (GAMMA - ONE)
                  CPSTAR = (2.E+0 / (GAMMA * FSMACH**2)) *
     >               (((TOGM1 + FSMACH**2) / GP1GM1)**GOGM1 - ONE)
                  CPSTAR = MAX (-999.E+0, CPSTAR)
               END IF

               CALL QPLFRM (LUNCPS, MNTITL, SBTIT2, XLABEL, 'Cp')

               ENDSTR = ENDCRV

               IF (MODE /= 1 .AND. MODE /= 2) THEN

C  *              For the original airfoil...

                  CALL GETCPS (NUORIG, XUORIG, YUORIG, YLORIG,
     >                         ALPHA, FSMACH, XC, CPU, CPL)

C  *              XC(J), CPU(J) here correspond to center of Jth panel.
C                 Integrations require wrap-around form  (lower t.e. to
C                 upper t.e.).  This means inserting points at the true
C                 leading, trailing edges - tedious.   Just average the
C                 upper and lower surface values nearest to the leading
C                 and trailing edges.  (Would extrapolation be better?)

                  XWRAP (1) = ONE
                  YWRAP (1) = YLORIG (NLORIG)
                  PWRAP (1) = (CPL (NLORIG-1) + CPU (NLORIG-1)) * HALF

                  I = NLORIG
                  DO J = 2, NLORIG
                     XWRAP (J) = XC (I-1)
                     YWRAP (J) = (YLORIG (I-1) + YLORIG (I)) * HALF
                     PWRAP (J) = CPL (I-1)
                     I = I - 1
                  END DO

                  J = NLORIG+1
                  XWRAP (J) = ZERO
                  YWRAP (J) = YLORIG (1)
                  PWRAP (J) = (CPL (1) + CPU (1)) * HALF

                  NWRAP = NLORIG + NUORIG + 1
                  I = 2
                  DO J = NLORIG+2, NWRAP-1
                     XWRAP (J) = XC (I-1)
                     YWRAP (J) = (YUORIG (I-1) + YUORIG (I)) * HALF
                     PWRAP (J) = CPU (I-1)
                     I = I + 1
                  END DO

                  XWRAP (NWRAP) = ONE
                  YWRAP (NWRAP) = YUORIG (NUORIG)
                  PWRAP (NWRAP) = PWRAP (1)

C  *              Calculate aerodynamic coefficients.  Ignore CD though.

                  CALL GETCL (NWRAP, XWRAP, YWRAP, PWRAP, ALFRAD,
     >                        CL, CD, CM)

                  IF (NOMODS) THEN
                     LEG = BLANK
                  ELSE
                     LEG = ORIGINAL
                  END IF

                  WRITE (LEGEND, 1120) LEG, CL, CM

C  *              Note handling of Cp data range: most negative Cp
C                 determines "top" of the Cp axis, not "bottom".
C                 Also, rounded values are needed for nice scaling.
C                 Use 0.5 as a sensible increment for most Cp distrbns:

                  CPTOP = PWRAP (1)
                  CPBOT = CPTOP
                  CALL BOUNDS (NWRAP, 1, MXWRAP, PWRAP, CPTOP, CPBOT)

                  NINC  = (-CPTOP + 0.499) / HALF
                  CPTOP = -NINC * HALF
                  NINC  = (CPBOT + 0.499) / HALF
                  CPBOT =  NINC * HALF

                  CALL QPLCRV (LUNCPS, NWRAP, XWRAP, PWRAP,
     >                         ENDSTR, UNDEF, BLANK, UNDEF, UNDEF,
     >                         CPBOT, CPTOP, UNDEF, UNDEF,
     >                         LEGEND, CPSLINE (1))

                  WRITE (LUNCRT, 1150) LEG, CL, CM
                  WRITE (LUNTAB, 1100) LEG, ALPHA, FSMACH, CL, CM
                  WRITE (LUNTAB, 1110) (XC (I), CPU (I), CPL (I),
     >                                  I = 1, NUORIG-1)

               END IF

               IF (MODE > 0) THEN

C  *              For the revised airfoil...

                  CALL GETCPS (NU, XU, YU, YL,
     >                         ALPHA, FSMACH, XC, CPU, CPL)

                  XWRAP (1) = ONE
                  YWRAP (1) = YL (NL)
                  PWRAP (1) = (CPL (NL-1) + CPU (NL-1)) * HALF

                  I = NL
                  DO J = 2, NL
                     XWRAP (J) = XC (I-1)
                     YWRAP (J) = (YL (I-1) + YL (I)) * HALF
                     PWRAP (J) = CPL (I-1)
                     I = I - 1
                  END DO

                  J = NL+1
                  XWRAP (J) = ZERO
                  YWRAP (J) = YL (1)
                  PWRAP (J) = (CPL (1) + CPU (1)) * HALF

                  NWRAP = NL + NU + 1
                  I = 2
                  DO J = NL+2, NWRAP-1
                     XWRAP (J) = XC (I-1)
                     YWRAP (J) = (YU (I-1) + YU (I)) * HALF
                     PWRAP (J) = CPU (I-1)
                     I = I + 1
                  END DO

                  XWRAP (NWRAP) = ONE
                  YWRAP (NWRAP) = YU (NU)
                  PWRAP (NWRAP) = PWRAP (1)

                  CALL GETCL (NWRAP, XWRAP, YWRAP, PWRAP, ALFRAD,
     >                        CL, CD, CM)

                  WRITE (LEGEND, 1120) REVISED, CL, CM

                  CPT = PWRAP (1)
                  CPB = CPT
                  CALL BOUNDS (NWRAP, 1, MXWRAP, PWRAP, CPT, CPB)

                  NINC  = (-CPT + 0.499) / HALF
                  CPT   = MIN (CPTOP, -NINC * HALF)
                  NINC  = (CPB + 0.499) / HALF
                  CPB   = MAX (CPBOT,  NINC * HALF)

                  CALL QPLCRV (LUNCPS, NWRAP, XWRAP, PWRAP,
     >                         ENDSTR, UNDEF, BLANK, UNDEF, UNDEF,
     >                         CPB, CPT, UNDEF, UNDEF, LEGEND,
     >                         CPSLINE (2))

                  LEG = REVISED        ! Because LEG is 1 character shorter
                  WRITE (LUNCRT, 1150) LEG, CL, CM
                  WRITE (LUNTAB, 1100) REVISED, ALPHA, FSMACH, CL, CM
                  WRITE (LUNTAB, 1110) (XC (I), CPU (I), CPL (I),
     >                                  I = 1, NU - 1)
               END IF

C  *           Finally, show Cp* as a curve with a legend.
C              Use legend for Mach and alpha for now - maybe in caption later.

               XC (1)  = -0.1E+0
               XC (2)  =  HALF
               CPU (1) =  CPSTAR
               CPU (2) =  CPSTAR

               IF (CPSTAR /= -999.E+0) THEN
                  WRITE (LEGEND, 1130) CPSTAR, FSMACH, ALPHA
               ELSE
                  WRITE (LEGEND, 1140) FSMACH, ALPHA
               END IF

               CALL QPLCRV (LUNCPS, 2, XC, CPU,
     >                      ENDSTR, UNDEF, BLANK, ZERO, ONE,
     >                      UNDEF, UNDEF, UNDEF, UNDEF,
     >                      LEGEND, 'DASH')

            ELSE  ! An empty Cp file may slip through unless:

               CLOSE (UNIT=LUNCPS, STATUS='DELETE')
               WRITE (LUNCRT, 1000) BLANK,
     >            'Cannot estimate Cps - upper/lower abscissas differ.',
     >            'Use REDISTRIBUTE mode first if you want cheap Cps.'
            END IF

         END IF


C  . . . . . . . Spreadsheet-compatible output (*.SPREAD) . . . . . . .

         IF (SPREADFIL) THEN

C  *        Can't do it if abscissas aren't common to both surfaces:

            DO I = 1, MIN (NU, NL)
               IF (XU (I) /= XL (I)) SPREADFIL = .FALSE.
            END DO

            IF (SPREADFIL) THEN

C  *           We have all outputs except mean-line and thickness.
C              Use X and Y as scratch for XFORM to transform in place:

               CALL COPY (NU, YU, X)
               CALL COPY (NU, YL, Y)
               CALL XFORM (.TRUE., NU, X, Y)

               WRITE (LUNSPR, 1010) MNTITL
               IF (SBTIT2 /= BLANK) WRITE (LUNSPR, 1010) SBTIT2

               WRITE (LUNSPR, 1200) XLABEL,
     >            HTAB // YLABEL // UP, HTAB // YLABEL // LO,
     >            HTAB // 'camber',   HTAB // 'semithickness',
     >            HTAB // 'y''' // UP,  HTAB // 'y''' // LO,
     >            HTAB // 'y"' // UP,   HTAB // 'y"' // LO,
     >            HTAB // 'curvature' // UP, HTAB // 'curvature' // LO

               WRITE (LUNSPR, 1210) (XU (I),
     >            HTAB, YU (I), HTAB, YL (I), HTAB, X (I), HTAB, Y (I),
     >            HTAB, YPU (I), HTAB, YPL (I), HTAB, YPPU (I),
     >            HTAB, YPPL (I), HTAB, YKUORG (I), HTAB, YKLORG (I),
     >            I = 1, NU)

            ELSE

C              An empty spreadsheet-compatible file should be removed:

               CLOSE (UNIT=LUNSPR, STATUS='DELETE')
               WRITE (LUNCRT, 1000) BLANK,
     >           'Cannot give spreadsheet file: upper/lower Xs differ.',
     >           'Use REDISTRIBUTE mode first.'
            END IF
         END IF

      END DO  ! Next profile

C-----------------------------------------------------------------------
C                     End of main loop over all profiles.
C-----------------------------------------------------------------------


C  *  Finally, remind the user of the files generated.

      WRITE (LUNCRT, 1000)
      IF (SAVREV) WRITE (LUNCRT, 1010)
     >   '   Modified airfoil:      ', IDENT (FIRST:LAST), '.dat'
      IF (PLTFIL) WRITE (LUNCRT, 1010)
     >   '   Airfoil plot file:     ', IDENT (FIRST:LAST), '.plt'
      IF (TABFIL) WRITE (LUNCRT, 1010)
     >   '   Tabulated results:     ', IDENT (FIRST:LAST), '.tab'
      IF (CRVFIL) WRITE (LUNCRT, 1010)
     >   '   Curvature QPLOT file:  ', IDENT (FIRST:LAST), '.crv'
      IF (YPPFIL) WRITE (LUNCRT, 1010)
     >   '   2nd derivatives file:  ', IDENT (FIRST:LAST), '.ypp'
      IF (CPSFIL) WRITE (LUNCRT, 1010)
     >   '   Cp distributions file: ', IDENT (FIRST:LAST), '.cps'
      IF (SPREADFIL) WRITE (LUNCRT, 1010)
     >   '   Spreadsheet file:      ', IDENT (FIRST:LAST), '.spread'
      GO TO 999


C-----------------------------------------------------------------------
C                           Error handling.
C-----------------------------------------------------------------------

  810 WRITE (LUNCRT, 1000) 'Stopping as requested.'
      GO TO 999
  850 WRITE (LUNCRT, 1000)
     >   'Error. PROFILE does not handle differing leading edge points.'
      GO TO 999
  860 WRITE (LUNCRT, 1010) ' Error in ', SURFACE (J),
     >   ': abscissas must be monotonic except for RECTIFY mode.'
      GO TO 999
  900 WRITE (LUNCRT, 1045) NPR + 1
      GO TO 999
  910 WRITE (LUNCRT, 1030) NPR
      GO TO 999
  920 WRITE (LUNCRT, 1040) NPR
      GO TO 999
  930 WRITE (LUNCRT, 1000) 'Empty input airfoil data file - abort.'
      GO TO 999
  940 WRITE (LUNCRT, 1020)
     >   ' Bad return from PSFIT during curvature calculation.  IER: ',
     >   IER

  999 WRITE (LUNCRT, 1000)

C *** STOP ' ' ! Avoid system differences


C-----------------------------------------------------------------------
C                              Formats.
C-----------------------------------------------------------------------

 1000 FORMAT (1X, A)
 1010 FORMAT (A, A, A)
 1020 FORMAT (/, A, I2, A3, A)
 1025 FORMAT (A, I2)
 1030 FORMAT (/, ' Too many or too few data pts. found on one surface.',
     >        I3, ' profile(s) processed.')
 1040 FORMAT (/, ' Abnormal EOF or other read error.',
     >        I3, ' profile(s) processed.')
 1045 FORMAT (/, ' Missing coordinate in profile ', I2, '.')
 1050 FORMAT (/, 1X, A,
     >        /'   Maximum thickness/chord (current axes) . .', F9.4,
     >        '% at X/C =', F8.5,
     >        /'   Maximum camber/chord . . . . . . . . . . .', F9.4,
     >        '% at X/C =', F8.5,
     >        /'   Unnormalized area  . . . . . . . . . . . .', E13.6,
C****>        /'   Ordinate of centroid . . . . . . . . . . .', E13.6,
     >        /'   Moment of inertia about X axis . . . . . .', E13.6)
 1055 FORMAT (/, ' *** Geometric properties skipped for ', A,
     >        'airfoil.', /, ' IER from AFGEOM: ', I3)
 1060 FORMAT (/'   L.E.:', 2E14.6, '  Chord:', E14.6,
     >        /'   Mean-line angle at T.E. (deg.) . . . . . .', F13.2)
 1070 FORMAT (/, ' WARNING: Common leading edge point is not foremost.',
     >        /, ' Proceeding, but some results such as derivatives',
     >        ' may be affected.', /)
 1080 FORMAT (/, ' WARNING: The second derivatives file (*.ypp) is',
     >        ' being suppressed because', /,
     >        ' wrap-around (parametric) curvature was specified.', /,
     >        ' REFINE and OPTIMIZE modes require the Y vs. X form of',
     >        ' Y".', /)
 1090 FORMAT ('! The curvature distribution is in wrap-around form, ',
     >        'lower t.e. to upper t.e.', /,
     >        '! Clockwise turns along this arc correspond to ',
     >        'negative curvature.')
 1100 FORMAT (/, '1Cp distribution estimates: ', A9, 'airfoil', //,
     >        ' Alpha:', F6.2, '   Free stream Mach:', F5.2, //,
     >        ' Cl:', F7.4, '   Cm:', F8.4, //,
     >        '      X           Cpu        Cpl')
 1110 FORMAT (1X, F10.5, F11.4, F11.4)
 1120 FORMAT (A, ':  Cl =', F7.4, ',  Cm =', F8.4)
 1130 FORMAT ('Cp* =', F8.2,',  Mach =', F4.2, ',  Alpha =', F5.2)
 1140 FORMAT ('Cp* = Inf., Mach =', F4.2, ',  Alpha =', F5.2)
 1150 FORMAT (1X, A, '  Cl:', F7.4, '   Cm:', F8.4)
 1200 FORMAT (A3, 10A)
 1210 FORMAT (1P, (E13.6, 10(A1, E13.6)))

      END PROGRAM PROFILE
C+----------------------------------------------------------------------
C
      SUBROUTINE ACTIVATE (GATHER, NVARS, ACTIVE, ALLVARS,
     >                     ACTVARS, VSCALES )
C
C PURPOSE:  ACTIVATE either gathers an active set of variables  (as a
C           contiguous subset) from a complete set,  or it "scatters"
C           a given active set back into the complete set.
C
C METHOD:   Both operations are combined here at the expense  of  one
C           more argument as a switch, to keep the module count down.
C           An array of logicals switches each variable on or off; no
C           need is seen to count the active ones here.
C
C ARGUMENTS:
C    ARG       DIM TYPE I/O/S DESCRIPTION
C   GATHER      -    L    I   .TRUE.  means  "gather" or "pack", else
C                             .FALSE. means "scatter" or "unpack".
C   NVARS       -    I    I   No. of variables in the complete set.
C   ACTIVE    NVARS  L    I   ACTIVE(I)=.TRUE. means the Ith variable
C                             is active, else it is inactive.
C   ALLVARS   NVARS  R   I/O  Full set of variables; input if GATHER,
C                             and unchanged,  else updated on return.
C   ACTVARS   NVARS  R   I/O  Active subset, packed; output if GATHER
C                             else input (and unchanged).
C   VSCALES   NVARS  R    I   Multiplicative scale factors applied to
C                             packed output variables if "gather" and
C                             divided out if "scatter".   This  array
C                             should parallel the full set, ALLVARS.
C HISTORY:
C
C   01/17/84   DAS    Initial design and code.
C   01/25/84   DAS    Introduced scaling/unscaling.
C
C AUTHOR: David Saunders, Informatics General, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER   NVARS
      REAL      ACTVARS (NVARS), ALLVARS (NVARS), VSCALES (NVARS)
      LOGICAL   ACTIVE (NVARS), GATHER

C ... Local variables:

      INTEGER   I, J

C ... Execution:

      J = 1

      DO I = 1, NVARS

         IF (ACTIVE (I)) THEN

            IF (GATHER) THEN
               ACTVARS (J) = ALLVARS (I) * VSCALES (I)
            ELSE
               ALLVARS (I) = ACTVARS (J) / VSCALES (I)
            END IF

            J = J + 1
         END IF

      END DO

      END SUBROUTINE ACTIVATE
C+----------------------------------------------------------------------
C
      SUBROUTINE ADDBUMPS (NPTS, X, Y, MXPARM, NBUMPS, BNAMES, PARAMS,
     >                     YNEW )
C
C PURPOSE:  ADDBUMPS adds the effect of a list of "bump" functions to
C           each of the given ordinates.  It may update the ordinates
C           in-place if desired.   It is intended for perturbing air-
C           foils, one surface at a time.
C
C METHOD:   Abscissas are assumed to be in the range [0,1], for effi-
C           ciency reasons.  The original ordinates are copied to the
C           new-ordinate array, then the effect of each bump is added
C           to this output array.  Even simple scaling of Y is imple-
C           mented this way, so it is not essential for a SCALE to be
C           the FIRST bump in the list unless the caller is hoping to
C           pass the same array for YNEW as for Y.
C
C ARGUMENTS:
C    ARG    DIM    TYPE I/O/S DESCRIPTION
C   NPTS     -       I    I   Number of ordinates to be perturbed.
C   X       NPTS     R    I   Abscissas (assumed in [0,1]) and the
C   Y       NPTS     R    I   original ordinates for given surface
C   MXPARM   -       I    I   Max. # parameters  for any one bump.
C   NBUMPS   -       I    I   No. perturbing functions to be used.
C   BNAMES  NBUMPS C*(*)  I   Names of the bumps, as recognized by
C                             subroutine BEVAL (which looks at the
C                             first four characters only).
C   PARAMS  MXPARM,  R    I   PARAMS(*,J)  are parameters defining
C           NBUMPS            the Jth bump.
C   YNEW    NPTS     R    O   The updated ordinates - may  be  the
C                             same array as Y if desired,  but be-
C                             ware if SCALE is among the bump set.
C
C PROCEDURES:
C
C    BEVAL   Evaluates indicated bump function at given normalized
C            abscissas, and adds to current ordinates.
C
C HISTORY:
C
C   01/14/84  DAS  Initial design and code.
C   04/11/84   "   Now uses subroutine BEVAL instead of function BUMP.
C   04/21/97   "   The test for 'SCALE' was failing.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   MXPARM, NBUMPS, NPTS
      REAL
     >   PARAMS (MXPARM, NBUMPS), X (NPTS), Y (NPTS), YNEW (NPTS)
      CHARACTER
     >   BNAMES (NBUMPS) * (*)

C     Local constants:

      LOGICAL, PARAMETER ::
     >   ADD = .TRUE.

      INTEGER
     >   I

C     Procedures:

      EXTERNAL
     >   BEVAL

C     Execution:

      YNEW = Y

      DO I = 1, NBUMPS

         IF (BNAMES (I) (1:4) == 'SCAL') THEN

C           This is a scaling of the ordinates (independent of X(I),
C           but still arranged to be additive, not multiplicative):

            CALL BEVAL (BNAMES (I), MXPARM, PARAMS (1, I), ADD, NPTS, Y,
     >                  YNEW)
         ELSE

            CALL BEVAL (BNAMES (I), MXPARM, PARAMS (1, I), ADD, NPTS, X,
     >                  YNEW)
         END IF

      END DO

      END SUBROUTINE ADDBUMPS
C+----------------------------------------------------------------------
C
      SUBROUTINE AFGEOM (NU, NL, NMAX, X, Y, AREA, CENTRD, MOMENT,
     >                   THICK, XTH, CAMBER, XCAM, TEANGLE, YPOWER,
     >                   YINTRP, B, C, D, IER)
C
C     Description:
C
C           AFGEOM (AirFoil GEOMetry) computes several geometrical quantities
C        based on integrals over an airfoil section, including the moment of
C        inertia about the x-axis, and magnitude and location of maximum
C        thickness and camber.  Each integration is performed by summing the
C        analytically-determined integrals of the cubics defined on each
C        subinterval by the spline representation of the appropriate function.
C
C           This version also computes the angle of the mean line at the
C        trailing edge.
C
C     Arguments:
C
C        Name    Dimension   Type   I/O/S  Description
C        NU                   I     I      Number of upper surface points.
C        NL                   I     I      Number of lower surface points.
C                                          4 <= N <= NMAX for N=NU and NL.
C        NMAX                 I     I      Row dimension of X and Y arrays in
C                                          calling program.
C        X       NMAX,2       R     I      Airfoil abscissae (increasing; not
C                                          necessarily same for both surfaces).
C        Y       NMAX,2       R     I      Airfoil surface ordinates.
C                                          (*,1) = upper surface; (*,2) = lower.
C        AREA                 R       O    Cross-sectional area.
C        CENTRD               R       O    Centroid of y ordinates.
C        MOMENT               R       O    Moment of inertia about x-axis.
C        THICK                R       O    Maximum thickness/chord ratio.
C        XTH                  R       O    X coordinate of maximum thickness.
C        CAMBER               R       O    Maximum camber relative to x-axis,
C                                          as fraction of chord.
C        XCAM                 R       O    X coordinate of maximum camber.
C        TEANGLE              R       O    Trailing edge angle of mean line
C                                          with X-axis (in degrees; negative
C                                          for positive camber).
C        YPOWER  NMAX         R         S  For powers of the surface ordinates.
C        YINTRP  NMAX         R         S  For thickness/camber estimate.
C        B,C,D   NMAX         R         S  For cubic spline coeffs.
C        IER                  I       O    < 0 means NX < 4 and CSFIT will fail;
C                                            1 means non-increasing Xs detected.
C     Notes:
C
C        (1)  There are several two dimensional arrays used.  The second
C             index (*, 1) or (*, 2) means upper or lower surface, resp.
C
C        (2)  Special handling of symmetric airfoils has been omitted for
C             simplicity, while special handling of different sets of
C             abscissae (upper/lower) is included for generality. Parametric
C             splines would provide a better fit at leading edge, but are
C             harder to integrate.
C
C     External references:
C
C        CSEVAL   Needed to evaluate spline for generalized thickness calc.
C        CSFIT    Calculates conventional cubic spline coefficients.
C        CSQUAD   Integrates a conventional cubic spline efficiently.
C        PROTECT  Checks for non-increasing Xs.
C
C
C     Author:  Robert Kennelly, NASA-Ames/Sterling Software
C
C     History:
C
C        15 Dec. 1981   RAK   Initial design and coding.
C        31 Mar. 1982   RAK   Installed in FLO6QNM, with cosmetic changes.
C        18 Apr. 1986   RAK   More cosmetics, including IMPLICIT NONE.
C        22 July 1986   RAK   Only need moment of inertia of area.
C        16 Sep. 1986   DAS   Generalized for different upper/lower Xs;
C                             eliminated symmetric airfoil efficiency;
C                             CSQUAD in place of more obvious quartics.
C        20 Oct. 1986   RAK   Minor mods. for consistency with FLO6QNM.
C                             Repaired flakey INTERP loop.  Return THICK
C                             and CAMBER as fraction of chord, not %.
C        31 Aug. 1987   DAS   Introduced trailing edge angle calculation
C                             - this is the obvious place to do it for
C                             the PROFILE application of AFGEOM.
C        22 Apr. 1991 DAS/RAK Max. camber and location of it and of max.
C                             thickness were not allowing for nonzero
C                             leading-edge coordinates.
C        21 Oct. 1999   DAS   ROTATE option can cause non-monotonic Xs.
C                             Return IER = +1 in this case.
C
C-----------------------------------------------------------------------


C     Declarations:
C     -------------

      IMPLICIT NONE

C     Constants:

      REAL, PARAMETER ::
     >   ZERO   = 0.0E+0,
     >   ONE    = 1.0E+0,
     >   TWO    = 2.0E+0,
     >   HALF   = 1.0E+0 / TWO,
     >   THIRD  = 1.0E+0 / 3.0E+0,
     >   FOURTH = 1.0E+0 / 4.0E+0,
     >   RADDEG = 57.29578E+0

C     Arguments:

      INTEGER
     >   IER, NL, NMAX, NU
      REAL
     >   AREA, B (NMAX), C (NMAX), CAMBER, CENTRD, D (NMAX), MOMENT,
     >   TEANGLE, THICK, X (NMAX, 2), XCAM, XTH, Y (NMAX, 2),
     >   YPOWER (NMAX), YINTRP (NMAX)

C     Local variables:

      INTEGER
     >   I, J, K, L, N (2), NX
      REAL
     >   ARROW, RCHORD, RESULT (3, 2), YLOWER, YUPPER
      LOGICAL
     >   INTERP, DISTINCT

C     Procedures:

      EXTERNAL
     >   CSEVAL, CSFIT, CSQUAD, PROTECT

C     Execution:
C     ----------

C     Initialize outputs in case of an early return:

      AREA    = ZERO
      CENTRD  = ZERO
      MOMENT  = ZERO
      THICK   = ZERO
      XTH     = ZERO
      CAMBER  = ZERO
      XCAM    = ZERO
      TEANGLE = ZERO

      N (1)   = NU
      N (2)   = NL

      DO J = 1, 2

         CALL PROTECT (N (J), X (1, J), Y (1, J), ARROW, DISTINCT)

         IF (ARROW /= ONE) THEN
            IER = 1
            GO TO 90
         END IF

      END DO

      RCHORD  = ONE / (MAX (X (N (1), 1), X (N (2), 2)) - X (1, 1))

C     For successive powers of the ordinates ...

      DO I = 1, 3

C        ... and for upper and lower surfaces, ...

         DO J = 1, 2

C           ... form integrand, ...

            NX = N (J)
            DO K = 1, NX
               YPOWER (K) = Y (K, J) ** I
            END DO

C           ... fit conventional cubic spline, ...

            CALL CSFIT (NX, X (1, J), YPOWER, 0, ZERO, 0, ZERO,
     >                  B, C, D, IER)
            IF (IER /= 0) GO TO 90

C           ... then integrate and save the results.  CSQUAD to overwrites
C           C for work-space reasons, but this means that thickness and
C           and camber must be calculated first to avoid refitting spline.

            IF (I == 1 .AND. J == 2) THEN

C              Need common sets of abscissae to estimate thickness easily.
C              Note that a more precise value could be computed using the
C              spline for both surfaces and a 1-dim. optimization routine.

               INTERP = (N (1) /= NX)
               IF (.NOT. INTERP) THEN
                  DO K = 1, NX
                     IF (X (K, 1) /= X (K, 2)) THEN

C                       The ordinates don't match - we'll have to interpolate.

                        INTERP = .TRUE.
                        EXIT

                     END IF
                  END DO
               END IF

               IF (INTERP) THEN

C                 Evaluate lower surface spline at upper surface Xs.

                  CALL CSEVAL (NX, X (1, 2), YPOWER, N (1), X (1, 1),
     >                         B, C, D, YINTRP)
               ELSE

C                 Copy lower surface to simplify next steps.

                  DO K = 1, NX
                     YINTRP (K) = YPOWER (K)
                  END DO
               END IF

C              The number of abscissae is now NU = N (1) either way.

               THICK  = ZERO
               CAMBER = ZERO
               XTH    = ZERO
               XCAM   = ZERO

               DO K = 1, N (1)
                  YUPPER = Y (K, 1)
                  YLOWER = YINTRP (K)
                  IF ((YUPPER - YLOWER) > THICK) THEN
                     THICK = YUPPER - YLOWER
                     XTH = X (K, 1)
                  END IF
                  IF (ABS (YUPPER + YLOWER) > ABS (CAMBER)) THEN
                     CAMBER = YUPPER + YLOWER
                     XCAM = X (K, 1)
                  END IF
               END DO

               THICK  = THICK  * RCHORD
               CAMBER = (CAMBER * HALF - Y (1, 1)) * RCHORD
               XTH  = ( XTH - X (1, 1)) * RCHORD
               XCAM = (XCAM - X (1, 1)) * RCHORD

C              This is also a convenient place to compute the angle that
C              the mean-line at the trailing edge makes with the X-axis.

               TEANGLE = ATAN2 (HALF * ((Y (NU, 1) + YINTRP (NU)) -
     >            (Y (NU - 1, 1) + YINTRP (NU - 1))),
     >            (X (NU, 1) - X (NU - 1, 1))) * RADDEG
            END IF

C           Integrate spline for Y ** I for surface J.  CSQUAD returns all
C           the subintegrals but they're not used here.

            CALL CSQUAD (NX, X (1, J), YPOWER, ZERO, C, C)
            RESULT (I, J) = C (NX)
         END DO
      END DO

      AREA   = (RESULT (1, 1) - RESULT (1, 2))
      CENTRD = (RESULT (2, 1) - RESULT (2, 2)) * HALF
      MOMENT = (RESULT (3, 1) - RESULT (3, 2)) * THIRD

C     Termination:
C     ------------

   90 RETURN

      END SUBROUTINE AFGEOM
C+----------------------------------------------------------------------
C
      SUBROUTINE BEVAL (BNAME, NP, P, ADD, NX, X, FX)
C
C PURPOSE:  BEVAL evaluates the indicated "bump" (shape function) at the
C           given abscissas, which are assumed to be normalized.  As the
C           bump names suggest,  these  functions  are commonly used for
C           perturbing airfoils.   They originated as the  "Hicks-Henne"
C           shape functions.  Some of them are also handy for generating
C           distributions of nonuniform weighting factors.
C
C           An option is provided to add the evaluations to the existing
C           values in the output array rather than just return values of
C           the current function.  This can save work-space in a calling
C           routine that is accumulating the effect of several bumps.
C
C METHOD:   The bump is selected by name rather than code number.   This
C           poses a problem of variable-length names,  resolved here  by
C           dealing with only the first four characters  -  an arbitrary
C           choice.   UPPER case is assumed  -  it seemed unnecessary to
C           handle lower case as well. Similarly, there is no attempt to
C           handle  unnormalized abscissas.
C
C           It was considered too inefficient, at this level, to attempt
C           handling of a given bump function's parameters by name. Thus
C           the ordering of the elements of P(*) IS important.
C
C           The option to ADD rather than just evaluate  is  implemented
C           by ALWAYS adding, but zeroing out first if necessary.   This
C           relieves the calling program from doing the zeroing.
C
C ARGUMENTS:
C    ARG    DIM  TYPE  I/O/S DESCRIPTION
C   BNAME    -   C*(*)   I   Name of desired  bump function.   Only
C                            the first 4 characters are  looked at,
C                            and they must be UPPER case.   Look at
C                            the code for valid names.
C    NP      -     I     I   Number of parameters (other than X) in
C                            the selected bump expression.
C    P      NP     R     I   The given values of the parameters, in
C                            a definite order.
C    ADD     -     L     I   ADD = .TRUE. means the evaluations for
C                            each X(I) are added into FX(I);
C                            ADD = .FALSE. means  FX(I)  is  zeroed
C                            out before this addition.
C    NX      -     I     I   The number of abscissas where the sel-
C                            ected function is to be evaluated.
C    X      NX     R     I   Abscissas in the range [0,1]. The case
C                            of simple scaling is an exception: the
C                            inputs and outputs involve ordinates.
C    FX     NX     R    I/O  The desired bump function values (pos-
C                            sibly added to input values; see ADD).
C
C NOTES:
C   *   The ordering of the shape function parameters is such that
C       the LAST one (P (NP)) is usually a multiplicative factor.
C
C HISTORY:
C   09/15/83   DAS   Initial design as FUNCTION BUMP.
C   01/14/84   DAS   Made simple scaling additive rather than multipli-
C                    cative: S*Y = Y + (S-1)*Y (simplifies usage).
C   02/22/84   DAS   Eliminated SQRT and SIN bump (redundant).
C   04/09/84   DAS   Adapted as SUBROUTINE BEVAL to get the loop inside
C                    rather than outside. Switched to selection by name
C                    rather than by type code; provided for  adding  as
C                    well as just evaluating; Wagner functions in-line.
C   07/17/86   DAS   Added simple "RAMP" function option.   Made Wagner
C                    functions the first choice.
C   08/11/86   DAS   Added "FLAP" and "SLAT" options (simple shearing).
C   07/26/89   RAK   ALOG changed to generic LOG.
C   11/06/93   RAK   "DROOP" function needed (1 - X) factor.
C   06/19/96   DAS   "EXP" function now has value 1. at specified X, as
C                    recommended by James Reuther;  added the symmetric
C                    forms of the modified sine function (SIN1 & SIN2),
C                    requiring use of leading 4 characters, not 3.
C   12/18/96   DAS   Added SINF, COSL, COSR, LCOS, & RCOS functions.
C   06/03/97   DAS   RADDEG name was misleading - changed it to DEGRAD.
C   10/20/99    "    Added SIN3 (fore/aft symmetry via SIN1 on [0,1])
C                    and   SIN4 (  "   "   "   "   "   "   " each half).
C
C AUTHORS: Leslie Collins, Robert Kennelly, David Saunders (Sterling);
C          Ray Hicks (NASA Ames Research Center)
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   NP, NX
      REAL
     >   FX (NX), P (NP), X (NX)
      LOGICAL
     >   ADD
      CHARACTER
     >   BNAME * (*)

C     Local constants:

      REAL, PARAMETER ::
     >   DEGRAD =  0.017453292519943E+0,
     >   PI     =  3.141592653589793E+0,
     >   PIBY2  =  1.570796326794897E+0,
     >   PT5LOG = -0.693147180559945E+0,
     >   RPI    =  0.318309886183790E+0,
     >   HALF   =  5.E-1,
     >   ONE    =  1.E+0,
     >   TWO    =  2.E+0,
     >   ZERO   =  0.E+0

C     Local variables:

      INTEGER
     >   I
      REAL
     >   AEXP, BEXP, CENTER, CENTER2, N, ONEMC, RN, POWER, POWERL,
     >   POWERR, SINXI, TANGNT, THETA, XI
      CHARACTER
     >   KEY * 4

C     Statement functions:

      REAL
     >   EBUMP, SBUMP, PWR, WIDTH, XNRM

      EBUMP (WIDTH, PWR, XNRM) =
     >   XNRM ** PWR * (ONE - XNRM) * EXP (-WIDTH * XNRM)

      SBUMP (WIDTH, PWR, XNRM) =
     >   (MAX (SIN (PI * XNRM ** PWR), ZERO)) ** WIDTH

C     Execution:

C     Check for just evaluating, not adding:

      IF (.NOT. ADD) THEN
         DO I = 1, NX
            FX (I) = ZERO
         END DO
      END IF

C     Avoid comparison of different-length strings:

      KEY = BNAME (1:4)

      IF (KEY == 'WAGN') THEN   ! Wagner functions.

C        Reference:  Ramamoorthy, P. and Padmavathi, K.  "Airfoil Design
C        by Optimization" in J. Aircraft, Vol. 14 (Feb. 1977), 219-221.

         N = P (1)
         RN = ONE / N

         IF (N == ONE) THEN
            DO I = 1, NX
               THETA = TWO * ASIN (SQRT (X (I)))
               FX (I) = FX (I) + ((THETA + SIN (THETA)) * RPI -
     >            (SIN (HALF * THETA)) ** 2) * P (2)
            END DO
         ELSE
            DO I = 1, NX
               THETA = TWO * ASIN (SQRT (X (I)))
               FX (I) = FX (I) + ((SIN (N * THETA) * RN  +
     >            SIN ((N-ONE) * THETA)) * RPI) * P (2)
            END DO
         END IF

      ELSE IF (KEY == 'SINE') THEN  ! Modified "SINE" (unsymmetric):
                                      ! P1 = "center", P2 = "width"
         POWER = PT5LOG / LOG (P (1))
         DO I = 1, NX
            FX (I) = P (3) * SBUMP (P (2), POWER, X (I)) + FX (I)
         END DO

      ELSE IF (KEY == 'SINF') THEN  ! Flipped "SINE" (unsymmetric):
                                      ! P1 = "center", P2 = "width"
         POWER = PT5LOG / LOG (ONE - P (1))
         DO I = 1, NX
            FX (I) = P (3) * SBUMP (P (2), POWER, ONE - X (I)) + FX (I)
         END DO

      ELSE IF (KEY == 'SIN1') THEN  ! Modified "SINE" (symmetric):
                                      ! P1 = "center", P2 = "width"
         CENTER = P (1)
         IF (CENTER <= HALF) THEN
            POWER = PT5LOG / LOG (CENTER)
            DO I = 1, NX
               FX (I) = P (3) * SBUMP (P (2), POWER, X (I)) + FX (I)
            END DO
         ELSE
            POWER = PT5LOG / LOG (ONE - CENTER)
            DO I = 1, NX
               FX (I) = P (3) * SBUMP (P (2), POWER, ONE - X (I)) +
     >                  FX (I)
            END DO
         END IF

      ELSE IF (KEY == 'SIN2') THEN  ! Modified "SINE" (symmetric):
                                      ! P1 = "center", P2 = "width"
         CENTER = P (1)
         IF (CENTER <= HALF) THEN
            POWER = PT5LOG / LOG (ONE - CENTER)
            DO I = 1, NX
               FX (I) = P (3) * SBUMP (P (2), POWER, ONE - X (I)) +
     >                  FX (I)
            END DO
         ELSE
            POWER = PT5LOG / LOG (CENTER)
            DO I = 1, NX
               FX (I) = P (3) * SBUMP (P (2), POWER, X (I)) + FX (I)
            END DO
         END IF

      ELSE IF (KEY == 'SIN3') THEN ! Fore/aft symmetry via SIN1 on [0, 1]

C        Assume [XA, XB] = [0, 1] and CENTER in (0, 0.5]:

         CENTER = P (1)
         POWER  = PT5LOG / LOG (CENTER)

         DO I = 1, NX
            XI = X (I)
            IF (XI > HALF) XI = ONE - XI
            FX (I) = P (3) * SBUMP (P (2), POWER, XI) + FX (I)
         END DO

      ELSE IF (KEY == 'SIN4') THEN ! Fore/aft symmetry via SIN1 on each half

C        Assume [XA, XB] = [0, 0.5] and CENTER in (0, 1).  Left half:

         CENTER  = P (1)
         POWERL  = PT5LOG / LOG (CENTER)
         CENTER2 = ONE - CENTER
         POWERR  = PT5LOG / LOG (CENTER2)

         IF (CENTER <= HALF) THEN
            DO I = 1, NX
               XI = MAX (ZERO, MIN (X (I) * TWO, ONE))
               FX (I) = P (3) * SBUMP (P (2), POWERL, XI) + FX (I)
            END DO
         ELSE
            DO I = 1, NX
               XI = MAX (ZERO, MIN (X (I) * TWO, ONE))
               FX (I) = P (3) * SBUMP (P (2), POWERR, ONE - XI) + FX (I)
            END DO
         END IF

C        Now the [0.5, 0.1] half with center <-- 1 - center

         IF (CENTER2 <= HALF) THEN
            DO I = 1, NX
               XI = MAX (ZERO, MIN ((X (I) - HALF) * TWO, ONE))
               FX (I) = P (3) * SBUMP (P (2), POWERR, XI) + FX (I)
            END DO
         ELSE
            DO I = 1, NX
               XI = MAX (ZERO, MIN ((X (I) - HALF) * TWO, ONE))
               FX (I) = P (3) * SBUMP (P (2), POWERL, ONE - XI) + FX (I)
            END DO
         END IF

      ELSE IF (KEY == 'COSL') THEN  ! 1/4 cosine (or sine), peak at left

C        SIN is used instead of COS because in the alternative SFEVAL form,
C        where COSL was first installed, PI is adjusted so that SIN (PI) > 0.

         DO I = 1, NX
            SINXI = MAX (SIN (PIBY2 * (X (I) + ONE)), ZERO)
            FX (I) = P (2) * SINXI ** P (1) + FX (I)
         END DO

      ELSE IF (KEY == 'COSR') THEN  ! 1/4 cosine (or sine), peak at right

C        SIN is used instead of COS for consistency with COSL.

         DO I = 1, NX
            SINXI = MAX (SIN (PIBY2 * X (I)), ZERO)
            FX (I) = P (2) * SINXI ** P (1) + FX (I)
         END DO

      ELSE IF (KEY == 'LCOS') THEN  ! Inverted 1/4 (co)sine, peak at left

         DO I = 1, NX
            SINXI = MAX (SIN (PIBY2 * X (I)), ZERO)
            FX (I) = P (2) * (ONE - SINXI ** P (1)) + FX (I)
         END DO

      ELSE IF (KEY == 'RCOS') THEN  ! Inverted 1/4 (co)sine, peak at right

         DO I = 1, NX
            SINXI = MAX (SIN (PIBY2 * (X (I) + ONE)), ZERO)
            FX (I) = P (2) * (ONE - SINXI ** P (1)) + FX (I)
         END DO

      ELSE IF (KEY == 'EXPO') THEN  ! "EXPONENTIAL" (peak height = 1.):
                                      ! P1 = "center", P2 = "width"
         ONEMC = ONE - P (1)
         AEXP  = P (1) * (ONE + ONEMC * P (2)) / ONEMC
         BEXP  = P (3) / EBUMP (P (2), AEXP, P (1))

         DO I = 1, NX
            FX (I) = BEXP * EBUMP (P (2), AEXP, X (I))  +  FX (I)
         END DO

      ELSE IF (KEY == 'DROO') THEN  ! "DROOP":  P1 = "width"

         DO I = 1, NX
            FX (I) = ((ONE - X (I)) * EXP (-P (1) * X (I))) * P (2)  +
     >               FX (I)
         END DO

      ELSE IF (KEY == 'LEAD') THEN  ! "LEADING"-edge:  P1 = "power"

         DO I = 1, NX
            FX (I) = ((ONE - X (I)) ** P (1)) * P (2)  +  FX (I)
         END DO

      ELSE IF (KEY == 'TRAI') THEN  ! "TRAILING"-edge:  P1 = "power"

         DO I = 1, NX
            FX (I) = (X (I) ** P (1)) * P (2)  +  FX (I)
         END DO

      ELSE IF (KEY == 'FLAP') THEN

C        "FLAP"-like function (shearing transformation only):
C        P (1) is the hinge point (fraction of chord from leading edge);
C        P (2) is the flap angle in DEGREEs (positive is "down").

         TANGNT = TAN (DEGRAD * P (2))
         DO I = 1, NX
            IF (X (I) > P (1))
     >         FX (I) = FX (I) - (X (I) - P (1)) * TANGNT
         END DO

      ELSE IF (KEY == 'SLAT') THEN

C        "SLAT"-like function (shearing transformation only):
C        P (1) is the hinge point (fraction of chord from leading edge);
C        P (2) is the flap angle in DEGREEs (positive is "down").

         TANGNT = TAN (DEGRAD * P (2))
         DO I = 1, NX
            IF (X (I) < P (1))
     >         FX (I) = FX (I) - (P (1) - X (I)) * TANGNT
         END DO

      ELSE IF (KEY == 'RAMP') THEN  ! "RAMP":  Y = P(1) * X

C        Commonly used in conjunction with Wagner functions.

         DO I = 1, NX
            FX (I) = P (1) * X (I)  +  FX (I)
         END DO

      ELSE IF (KEY == 'SCAL') THEN

C        Simple scaling - X(I) are probably ordinates.  Note that the
C        effect is arranged to be  additive,  not multiplicative,  so
C        that scaling can be treated just as for the  other functions
C        when it is included in a group for perturbing purposes.

         DO I = 1, NX
            FX (I) = X (I) * (P (1) - ONE)  +  FX (I)
         END DO

      ELSE
         WRITE (*, '(/, A)')  ' BEVAL: Illegal bump name.'
         STOP
      END IF

      END SUBROUTINE BEVAL
C+----------------------------------------------------------------------
C
      SUBROUTINE BEVAL2 (NAME, MODE, NP, P, ADD, NX, X, FX, LUNERR, IER)
C
C ONE-LINER:  Shape function utility, second collection
C
C PURPOSE:
C
C        BEVAL2 supplements the earlier BEVAL's collection of "bump"
C     (shape) functions.  The functions here are sufficiently more
C     elaborate than the standard ones to warrant separating them.
C     In fact they may require some initial calls to establish some
C     of their parameters - hence the extra MODE argument.  Otherwise
C     the usage is as for BEVAL, q.v.
C
C ARGUMENTS:
C     ARG    DIM  TYPE  I/O/S DESCRIPTION
C     NAME    -   C*(*)   I   Name of desired shape function.   Only
C                             the first 6 characters are looked at,
C                             and they must be UPPER case.  Look at
C                             the code for valid names.
C     MODE    -     I     I   MODE = 0 means initialize the function:
C                                    some element(s) of P are returned;
C                             MODE = 1 means evaluate the function at
C                                    the given abscissa(s).
C     NP      -     I     I   Number of parameters (other than X) in
C                             the shape function.
C     P      NP     R    I/O  The parameters, in a definite order - see
C                             MODE, and the NOTES below.
C     ADD     -     L     I   ADD = .TRUE. means the evaluations for
C                                   each X (I) are added into FX (I);
C                             ADD = .FALSE. means  FX (I) is zeroed
C                                   before this addition.
C     NX      -     I     I   The number of target abscissas.  NX >= 1.
C     X      NX     R     I   Abscissas in the range [0, 1].
C     FX     NX     R    I/O  The desired bump function values (possibly
C                             added to input values - see ADD).
C     LUNERR  -     I     I   LUNERR < 0 suppresses iteration printout (MODE=0);
C                             |LUNERR| is used for error messages.
C     IER     -     I     I   IER = 0 means successful initialization;
C                             IER < 0 means a failure in the zero-finding
C                                   or in the quadrature - probably the
C                                   no. of fn. evals. limit was exceeded.
C NOTES:
C     By convention, the ordering of the shape function parameters is
C     such that the LAST one, P (NP), is a multiplicative factor.
C     See the code for the meaning of the others.
C
C PROCEDURES:
C     QUANC8RC  Reverse-communication adaptive quadrature
C     ZERORC       "      "      "    1-D zero-finder
C
C HISTORY:
C   02/04/94  D.A.Saunders  Initial adaptation of BEVAL for some new
C                           functions proposed by James.
C
C AUTHORS: James Reuther, David Saunders, NASA Ames, Mt. View, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   MODE, NP, NX, LUNERR, IER
      REAL
     >   FX (NX), P (NP), X (NX)
      LOGICAL
     >   ADD
      CHARACTER
     >   NAME * (*)

C     Local constants:

      INTEGER, PARAMETER ::
     >   MAXFUNZ = 40

      REAL, PARAMETER ::
     >   HALF    =  5.E-1,
     >   ONE     =  1.E+0,
     >   ZERO    =  0.E+0,
     >   LOGPT5  = -0.693147180559945E+0,
     >   TWOPI   =  6.283185307179586E+0,
     >   ABSERR  = 1.E-6,   ! Quadrature's absolute error tolerance
     >   RELERR  = 1.E-5,   ! Tolerance relative to the integral size
     >   TOL     = ZERO,    ! Zero-finder's tolerance
     >   P3A     = -100.,   ! Search interval for P (3)
     >   P3B     = +100.

      CHARACTER, PARAMETER ::
     >   SUBNAME * 6 = 'BEVAL2'

C     Local variables:

      INTEGER
     >   I, I2, ISTATQ, ISTATZ, NFUNQ, NFUNZ
      REAL
     >   ERRESTQ, FLAGQ, FUNQ, FUNZ, HOLDZ (13), P1, P3, P4,
     >   POWER, SINE1, XQ
      LOGICAL
     >   TYPE1

C     Procedures:

      EXTERNAL
     >   QUANC8RC, ZERORC

C     Execution:

      IER = 0

C     Check for just evaluating, not adding:

      IF (.NOT. ADD) THEN
         IF (MODE /= 0) THEN
            FX = ZERO
         END IF
      END IF

      IF (NAME (1 : 4) == 'SINE') THEN

C        Two variations of area-conserving modified sine function.
C        Each produces a family which is symmetric about P1 = 0.5.
C
C        P1 is the zero-crossing value of X in [0, 1];
C        P2 is the power of the sine function (must be a whole number > 0);
C        P3 is the slope of the linear weight function giving zero area;
C        P4 is the overall multiplier.
C
C        'SINE1' uses P1 <= 0.5 functions (or reflections if P1 > 0.5);
C        'SINE2' uses P1 >= 0.5 functions (or reflections if P1 < 0.5).

         P1 = P (1)
         I2 = P (2)     ! (-x) ** y is undefined unless y is an integer

         IF (I2 < 1) THEN
            WRITE (ABS (LUNERR), '(/, A)') ' BEVAL2: P2 < 1.'
            GO TO 900
         END IF

         TYPE1 = NAME (5 : 5) == '1'         ! Else it is '2'
         TYPE1 = TYPE1 .AND. P1 <= HALF  .OR.
     >     .NOT. TYPE1 .AND. P1 >  HALF      ! Expression type now

         IF (TYPE1) THEN
            POWER = LOGPT5 / LOG (P1)
         ELSE
            POWER = LOGPT5 / LOG (ONE - P1)
         END IF

         IF (MODE == 0) THEN

C           Determine the area-conserving weighting factor (only).

            IF (P1 == HALF) THEN

C              Special case:  P3 = 0 is the exact solution.

               P (3) = ZERO
               GO TO 999      ! Avoids excessive indenting.

            END IF

C           We have a zero-finding iteration wrapped around an
C           adaptive quadrature iteration.  Use reverse-communication
C           to avoid passing data via common blocks.

            ISTATZ = 2        ! Initializes the zero-finder
            NFUNZ  = MAXFUNZ

   20       CONTINUE

               CALL ZERORC (P3A, P3B, P3, FUNZ, TOL, NFUNZ, SUBNAME,
     >                      LUNERR, HOLDZ, ISTATZ)

               IF (ISTATZ < 0) THEN  ! Probable application error
                  WRITE (ABS (LUNERR), 1000) ISTATZ
                  IER = ISTATZ
                  GO TO 999

               ELSE IF (ISTATZ > 0) THEN

C                 Calculate the integral for this P3 (iteratively):

                  ISTATQ = 0  ! Initialize the quadrature
                  XQ = ZERO   ! Left end of interval

   30             CONTINUE

C                    Evaluate the shape function (one of two types).
C                    Arrange for higher powers to take the sign of sin ** 1.

                     IF (TYPE1) THEN
                        SINE1 = SIN (TWOPI * XQ ** POWER)
                        FUNQ = (P3 * (XQ - P1) + ONE) *
     >                         SIGN (SINE1 ** I2, SINE1)
                     ELSE
                        SINE1 = SIN (TWOPI * (ONE-XQ) ** POWER)
                        FUNQ = -(P3 * (XQ - P1) + ONE) *
     >                         SIGN (SINE1 ** I2, SINE1)
                     END IF

                     CALL QUANC8RC (XQ, FUNQ, ZERO, ONE, ABSERR, RELERR,
     >                              FUNZ, ERRESTQ, NFUNQ, FLAGQ, ISTATQ)

                     IF (ISTATQ > 0)
     >            GO TO 30

                  IF (ISTATQ /= 0) THEN
                     WRITE (ABS (LUNERR), 1001) FLAGQ, NFUNQ, FUNZ
                     IER = ISTATQ
                     GO TO 999
                  END IF

C                 Else we have evaluated the integral successfully.

                  GO TO 20  ! Continue the search for a zero integral

C              ELSE         ! ISTATZ = 0 - a zero integral has been found
               END IF

            P (3) = P3


         ELSE    ! MODE = 1 (evaluate the shape function with its multiplier)

            P3 = P (3)
            P4 = P (4)
            IF (TYPE1) THEN
               DO I = 1, NX
                  SINE1 = SIN (TWOPI * X (I) ** POWER)
                  FX (I) = FX (I) + P4 * (P3 * (X (I) - P1) + ONE) *
     >                     SIGN (SINE1 ** I2, SINE1)
               END DO
            ELSE
               DO I = 1, NX
                  SINE1 = SIN (TWOPI * (ONE - X (I)) ** POWER)
                  FX (I) = FX (I) - P4 * (P3 * (X (I) - P1) + ONE) *
     >                     SIGN (SINE1 ** I2, SINE1)
               END DO
            END IF

         END IF

      ELSE
         WRITE (ABS (LUNERR), '(/, A)') ' BEVAL2: Bad function name.'
         GO TO 900
      END IF

      GO TO 999

  900 STOP

  999 RETURN

C     Formats:

 1000 FORMAT (/, ' *** BEVAL2 trouble (ZERORC):  status code = ', I2)
 1001 FORMAT (/, ' *** BEVAL2 trouble (QUANC8RC):', /, 5X,
     >        'FLAG = ', F7.3, '   NFUN = ', I5,  '   Integral = ', 1P,
     >        E12.3)

      END SUBROUTINE BEVAL2
C+----------------------------------------------------------------------
C
      SUBROUTINE BLENDFX (NF, NX, MAXNX, FX, WEIGHT, FBLEND)
C
C  PURPOSE:  BLENDFX combines NF discrete functions of X, presumed to be
C            defined at the same abscissas, using the linear combination
C            defined by WEIGHT(1:NF).  The original intended application
C            is to forming composite airfoils from basis airfoils.
C
C  METHOD:   Simple-minded - assume the weights are suitably normalized.
C            If these are design variables,  any normalization  is  best
C            done by the application program anyway.    If both surfaces
C            of a group of profiles are being blended, then each profile
C            would need to be in wrap-around form (one call to BLENDFX),
C            or the surfaces would have to be stored separately (meaning
C            two calls to BLENDFX).
C
C            No attempt is made to check for mismatched abscissas - they
C            are not involved here - just the ordinates.
C
C  ARGUMENTS:
C   ARG      DIM    TYPE I/O/S DESCRIPTION
C    NF       -       I    I   Number of discrete functions to blend
C    NX       -       I    I   Number of data points (same for all fns.)
C  MAXNX      -       I    I   Row dimension of FX(*,*) in calling prog.
C    FX    MAXNX,NF   R    I   The different functions stored by columns
C  WEIGHT    NF       R    I   WEIGHT(J) is applied to FX(I,J) for all I
C  FBLEND    NX       R    O   Blended function desired
C
C  DEVELOPMENT HISTORY:
C  DATE    INITIALS    DESCRIPTION
C  10/24/85   DAS      Initial implementation - after a comment by
C                      Ray Hicks that linear combinations of airfoils
C                      have been used effectively in airfoil design by
C                      optimization.
C
C  AUTHORS:  David Saunders, Informatics, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER
     >   NF, NX, MAXNX

      REAL
     >   FX (MAXNX, NF), WEIGHT (NF), FBLEND (NX)

C ... Local variables:

      INTEGER
     >   I, J

C ... Execution:

C ... Considerations:

C     (1)  NX is likely to be much greater than NF.
C     (2)  Why bother to zero out the sums?

      DO I = 1, NX
         FBLEND (I) = WEIGHT(1) * FX (I, 1)
      END DO

      DO J = 2, NF
         DO I = 1, NX
            FBLEND (I) = WEIGHT(J) * FX (I, J) + FBLEND (I)
         END DO
      END DO

      END SUBROUTINE BLENDFX
C+----------------------------------------------------------------------
C
      SUBROUTINE CALCTHICK (NL, NU, XL, XU, YL, YU, THICKNESS,
     >                      XATMAX, B, C, D, YEVAL)
C
C  PURPOSE: CALCTHICK  determines the maximum thickness of an airfoil
C           and the corresponding abscissa. The thickness is returned
C           as a percentage of the chord (which may or may not be 1).
C
C  METHOD:  If the  upper and lower surfaces are not  defined at  the
C           same  set  of abscissas,  a spline  is  fit to  the lower
C           surface and evaluated using the upper surface  abscissas.
C           Otherwise the thickness is determined using the  original
C           coordinates.  No attempt is made to estimate the location
C           BETWEEN data points where the thickness appears greatest.
C
C  ARGUMENTS: These are obvious, except for B, C, D, and YEVAL, which
C             are  work-space for spline coefficents for fitting just
C             one surface,  plus  evaluation of the spline at the ab-
C             scissas of the other surface.
C
C  HISTORY:
C     Oct. '83  DAS  Original design and coding
C     01/25/84  LJC  Added spline fitting for unlike upper and lower
C                    surface abscissas
C     03/16/84  DAS  Arranged for it not to matter which surface is
C                    upper and which is lower
C     07/20/84  LJC  Replaced SPLINE and SEVAL with CSFIT and CSEVAL
C     01/29/90  DAS  Removed tabs, underscores, and END DOs
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

C     Arguments:

      INTEGER   NL, NU
      REAL      XATMAX, THICKNESS,
     >          B(NL), C(NL), D(NL), YEVAL(NU),
     >          XL(NL), XU(NU), YL(NL), YU(NU)

C     Local variables:

      INTEGER   I, IER
      REAL      DY
      LOGICAL   FITSPLINE

C     Procedures:

      EXTERNAL  CSEVAL, CSFIT

C     Execution:

      THICKNESS = 0.E+0
      FITSPLINE = NU/=NL

      IF (.NOT. FITSPLINE) THEN
         DO I = 1, NU
            IF (XU(I) /= XL(I)) FITSPLINE = .TRUE.
         END DO
      END IF

      IF (FITSPLINE) THEN

C        Spline lower surface and evaluate at upper surface abscissas:

         CALL CSFIT (NL, XL, YL, 2, 0., 2, 0., B, C, D, IER)
         IF ( IER/=0 ) GO TO 900

         CALL CSEVAL (NL, XL, YL, NU, XU, B, C, D, YEVAL)
      END IF

      DO I = 1, NU
         IF (FITSPLINE) THEN
            DY = YU(I) - YEVAL(I)
         ELSE
            DY = YU(I) - YL(I)
         END IF
         IF (THICKNESS < ABS (DY)) THEN
            THICKNESS = ABS (DY)
            XATMAX = XU(I)
         END IF
      END DO

      THICKNESS = THICKNESS * 100.E+0 / (XU(NU) - XU(1))

      RETURN

  900 WRITE (*, '(/, A)') ' CALCTHICK: Error in fitting spline.'

      END SUBROUTINE CALCTHICK
C+----------------------------------------------------------------------
C
      SUBROUTINE CFDISTRIB (NOPTVARS, OPTVARS, SUMSQS)
C
C ACRONYM: Cf. (compare) distributions
C          --            -------
C PURPOSE: CFDISTRIB  computes  the sum-of-squares type of objective
C          function associated with two 1-dimensional distributions.
C          The argument list is that expected by the QNMDIF optimiz-
C          ing algorithm.   The target and current distributions are
C          NOT assumed to be defined at the same abscissas.  Provis-
C          is made for weighting the elements of the sum of squares.
C
C          This routine is application-dependent,  but the structure
C          may be appropriate in other contexts.
C
C          This version of CFDISTRIB  deals with curvature distribu-
C          tions for one airfoil surface. The curvature distribution
C          is a function of "bump" functions applied to the surface.
C          The active variables/parameters for these  bump functions
C          are optimized by QNMDIF to give a perturbed airfoil surf-
C          ace with curvature as close as possible to a  target dis-
C          tribution - an approach to designing or refining  airfoil
C          profiles.
C
C METHOD:  The target distribution is assumed already available, and
C          may cover any part or all of the range of the  calculated
C          distribution, which is generated here with each call. The
C          bulk of the quantities involved are passed through COMMON
C          because of the restricted calling sequence above. An out-
C          line of the steps required for this application  follows.
C
C          *  "Scatter" (unpack) the given active variables.   (This
C             handles cases where certain shape function  parameters
C             are fixed but must be present during evaluation.)
C
C          *  Apply the shape functions to each of the ordinates de-
C             fining the original surface (NOT in-place).
C
C          *  Calculate the resulting curvature distribution.
C
C          *  Evaluate this calculated distribution at each  of  the
C             target abscissas using table-look-up/linear interpola-
C             tion, and accumulate the desired sum of squares.
C
C          *  If necessary, impose constraints, probably in the form
C             of a penalty function constraining the thickness.
C
C ARGUMENTS:
C    ARG       DIM    TYPE I/O/S DESCRIPTION
C  NOPTVARS     -       I    I   Number of variables being optimized
C                                (passed for compatibility with QNM-
C                                DIF, but not actually used directly
C                                here).
C  OPTVARS   NOPTVARS   R    I   Current values of the variables.
C  SUMSQS       -       R    O   Corresponding value of the function
C                                being minimized.
C
C PARAMETER CONSTANTS:
C
C  LUNERR      Logical unit number for error messages.
C  MXBUMPS     Maximum no. of bump functions provided for.
C  MXOPT       Maximum no. of optimization variables allowed.
C  MXPARM      Maximum no. of parameters associated with any of
C              the bump functions (including a multiplier).
C  MXSURF      Maximum no. of data points handled per surface.
C  MXTARG      Maximum no. of target curvature values handled.
C
C COMMON BLOCKS USED:
C
C   /ACTUAL/  (Current values corresponding to input opt. variables)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C   YCALC     MXSURF     R    S   Perturbed airfoil surface
C  Y1CALC     MXSURF     R    S   1st and 2nd derivatives - by-prod-
C  Y2CALC     MXSURF     R    S   ucts of the curvature calculations
C  YKCALC     MXSURF     R    S   Curvature distribution correspond-
C                                 to the current optimization vars.
C
C   /BCHARS/  (Bump names - wish they could be in /BUMPS/)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C  BNAMES    MXBUMPS   C*11   I   List of bump names
C
C   /BUMPS /  (List of bump function parameters, etc.)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C  NBUMPS       -        I    I   Number of bump functions active
C  PARAMS MXPARM,MXBUMPS R  I/O/S Complete set of parameters for the
C                                 active bump functions,   including
C                                 some possible unused ones  present
C                                 because a 2-dim. array approach is
C                                 used for storing them  (as opposed
C                                 to packing the variable numbers of
C                                 parameters associated with differ-
C                                 ent shape functions).  See ACTIVE.
C  ACTIVE MXPARM,MXBUMPS L    I   Marks the bump function parameters
C                                 as active or inactive.  The unused
C                                 ones must be marked inactive along
C                                 with the fixed (but used) ones.
C  VSCALE MXPARM,NBUMPS  R    I   Scale factors needed by optimizing
C                                 algorithm for active variables.
C
C   /ORIGNL/  (Copies of the original airfoil surface data)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C    NXY        -        I    I   Number of points defining  surface
C   XORIG      NXY       R    I   Abscissas for the surface
C   YORIG      NXY       R    I   Ordinates for the unperturbed srf.
C
C   /TARGET/  (Target distribution quantities)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C   NTARG       -        I    I   Number of target data points
C   XTARG     NTARG      R    I   Abscissas for the target data
C  YKTARG     NTARG      R    I   Target curvature values
C  WEIGHT     NTARG      R    I   Weights to be applied (multiplica-
C                                 tively) to each element of the sum
C                                 of squares.  May be all 1s.
C FILES USED:
C    LUN      I/O/S DESCRIPTION
C   LUNERR      O   Error messages
C
C EXTERNAL REFERENCES:
C  ACTIVATE   For inserting active bump variables into complete set
C  ADDBUMPS   For evaluating and applying the bump functions
C  FD12K      For evaluating the perturbed curvature distribution
C  TABLE1     For matching target and current curvature distributions
C
C HISTORY:
C   01/17/84   DAS   Initial design.
C   03/16/84   DAS   Provided for constraining thickness (penalty fun.)
C   02/07/90   DAS   Removed list-directed I/O.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER  NOPTVARS
      REAL     OPTVARS(NOPTVARS), SUMSQS

C ... Parameter constants and
C ... Common blocks:

C ... Quantities passed to objective function routine via COMMON
C     because of fixed calling sequence expected by QNMDIF.  Any
C     changes affecting these COMMONs should also be made in one
C     other routine - OPTIMIZE.
C
      INTEGER, PARAMETER ::
     >   MXBUMPS = 20,  MXPARM = 3, MXOPT = MXPARM*MXBUMPS,
     >   MXSURF  = 180, MXTARG = MXSURF, LUNERR = 6

      REAL            YCALC, Y1CALC, Y2CALC, YKCALC
      COMMON /ACTUAL/ YCALC(MXSURF), Y1CALC(MXSURF), Y2CALC(MXSURF),
     >                YKCALC(MXSURF)

      CHARACTER*11    BNAMES
      COMMON /BCHARS/ BNAMES(MXBUMPS)

      INTEGER         NBUMPS
      REAL            PARAMS, VSCALE
      LOGICAL         ACTIVE
      COMMON /BUMPS / PARAMS(MXPARM,MXBUMPS), VSCALE(MXOPT),
     >                ACTIVE(MXPARM,MXBUMPS), NBUMPS

      INTEGER         NOTHER, NXY
      REAL            XORIG, XOTHER, YORIG, YOTHER
      COMMON /ORIGNL/ XORIG(MXSURF), YORIG(MXSURF),
     >                XOTHER(MXSURF), YOTHER(MXSURF), NOTHER, NXY

      INTEGER         NTARG
      REAL            OPTWRK, PENLTY, TARGTH, XTARG, YKTARG, WEIGHT
      COMMON /TARGET/ XTARG(MXTARG), YKTARG(MXTARG), WEIGHT(MXTARG),
     >                TARGTH, PENLTY, OPTWRK(4*MXSURF), NTARG

C ... Local variables:

      INTEGER  I, IER, INDEX
      REAL     THICK, XATMAX

C ... Procedures:

      EXTERNAL ACTIVATE, ADDBUMPS, CALCTHICK, FD12K, TABLE1
      REAL     TABLE1


C ... Execution:

C ... Unpack the active (optimization) variables so that the bump
C ... functions can be evaluated. (Some parameters may be fixed.)

      CALL ACTIVATE (.FALSE., MXPARM*NBUMPS, ACTIVE, PARAMS,
     >               OPTVARS, VSCALE)

C ... Evaluate and apply the bump functions to the original surface.

      CALL ADDBUMPS (NXY, XORIG, YORIG, MXPARM,
     >               NBUMPS, BNAMES, PARAMS, YCALC)

C ... Generate the corresponding curvature distribution.

      CALL FD12K (NXY, XORIG, YCALC, Y1CALC, Y2CALC, YKCALC)

C ... Check for constraining the thickness using a penalty function:

      IF (PENLTY > 0.E+0) THEN

         CALL CALCTHICK (NXY, NOTHER, XORIG, XOTHER, YCALC, YOTHER,
     >                   THICK, XATMAX, OPTWRK, OPTWRK(MXSURF+1),
     >                   OPTWRK(2*MXSURF+1), OPTWRK(3*MXSURF+1))

C ...    (Note: CALCTHICK doesn't care which surface is which.)

      ELSE
         THICK = TARGTH
      END IF

C ... Compute the objective function (the penalty part being optional):

      SUMSQS = PENLTY * (THICK - TARGTH) ** 2
      INDEX = 1

      DO I = 1, NTARG

C ...    Interpolate the current distribution to each target abscissa.
C        Assume the target abscissas are ordered (to use updated INDEX).

         SUMSQS = (TABLE1 ( NXY, XORIG, YKCALC, INDEX, XTARG(I), IER)
     >            -  YKTARG(I)) ** 2  * WEIGHT(I)  +  SUMSQS

         IF (IER /= 0) GO TO 900

      END DO

      GO TO 999

C ... Error handling:

  900 WRITE (LUNERR, '(/, A, 2I6)')
     >   ' CFDISTRIB: Table look-up error. IER, I: ', IER, I

  999 RETURN

      END SUBROUTINE CFDISTRIB
C+----------------------------------------------------------------------
C
      SUBROUTINE COMBINE (NU, XU, YU, NL, XL, YL,
     >                    MAXPTS, X, Y, XU2ND, YU2ND, XL2ND, YL2ND,
     >                    LUNCRT, LUNKBD, LUN2ND)
C
C  PURPOSE:  COMBINE combines two "airfoils" by addition or subtraction.
C            Adding or removing boundary layer displacement thickness is
C            the original rationale, but there may be other uses.
C
C  METHOD:   COMBINE is intended to be one of PROFILE's high-level modes
C            of operation.   Therefore a "primary" airfoil is assumed to
C            be read by PROFILE in the usual way,  and a secondary "air-
C            foil" (which may actually represent a boundary layer)  will
C            be prompted for here.   The latter should be passed through
C            PROFILE's REDISTRIBUTE mode first if necessary, so that the
C            abscissas match.   No attempt is made to do the redistribu-
C            tion here - one look at REDISTRIB explains why.
C
C  ARGUMENTS:
C   ARG      DIM    TYPE I/O/S DESCRIPTION
C  NU,NL      -       I   I/O  No. of upper/lower surface pts. (primary)
C  XU,XL    NU,NL     R   I/O  Abscissas, upper/lower (primary)
C  YU,YL    NU,NL     R   I/O  Corresponding ordinates, upper/lower.
C                              Surfaces must have common leading edge.
C  MAXPTS     -       I    I   Max. no. pts. provided for on 1 surface.
C  X,Y     2*MAXPTS   R    S   Buffers for reading secondary airfoil.
C  XU2ND,YU2ND MAXPTS R    S   Secondary abscissas found.
C  XL2ND,YL2ND MAXPTS R    S   Secondary ordinates found.
C  LUNCRT     -       I    I   Logical unit for prompts.
C  LUNKBD     -       I    I   Logical unit for responses.
C  LUN2ND     -       I    I   Logical unit for secondary coordinates
C                              in one of the PROFILE formats.
C
C  PROCEDURES:
C    OPENER   File opening utility
C    PRREAD   For reading secondary profile
C    READER   Prompting utility
C
C  HISTORY:
C  09/30/85   DAS   Adapted from REDISTRIB.
C  01/29/90   DAS   Removed the END DOs; installed OPENER.
C
C  AUTHOR:  David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, LUN2ND, MAXPTS, NL, NU
      REAL
     >   X(MAXPTS*2), XL(NL), XU(NU), XL2ND(MAXPTS), XU2ND(MAXPTS),
     >   Y(MAXPTS*2), YL(NL), YU(NU), YL2ND(MAXPTS), YU2ND(MAXPTS)

C ... Local variables:

      INTEGER
     >   FORMAT, I, IER, NL2ND, NU2ND
      REAL
     >   SIGNL, SIGNU
      LOGICAL
     >   ADDED, BLAYER, DEFAULT, QUIT
      CHARACTER
     >   TITLE*80, XY2ND*48

C ... Procedures:

      EXTERNAL
     >   OPENER, PRREAD, READY

C ... Execution:

C ... Read secondary coordinates (in any of the PROFILE formats):

      CALL OPENER (LUNCRT,
     >             'Enter file name for secondary coordinates: ',
     >             LUNKBD, XY2ND, LUN2ND, 'OLD')

      CALL PRREAD (LUN2ND, TITLE, MAXPTS, NU2ND, NL2ND, X, Y,
     >             XU2ND, XL2ND, YU2ND, YL2ND, FORMAT, IER)
      IF (IER /= 0) GO TO 830

C ... Check for mismatched abscissas:

      IF (NU2ND /= NU) GO TO 840
      IF (NL2ND /= NL) GO TO 840

      DO I = 1, NU2ND
         IF (XU2ND(I) /= XU(I)) GO TO 840
      END DO

      DO I = 1, NL2ND
         IF (XL2ND(I) /= XL(I)) GO TO 840
      END DO

C ... Check for adding or subtracting:
C     (Hard to avoid distinguishing boundary layer operation from
C     two-true-airfoils case here.)

      BLAYER = .TRUE.
      CALL READY (LUNCRT,
     >   'Do secondary coordinates represent displacement thickness? '//
     >   '<CR>=Y: ',
     >   LUNKBD, BLAYER, DEFAULT, QUIT)

      ADDED = .TRUE.
      CALL READY (LUNCRT,
     >   'Are they to be ADDed? <CR>=Y=added; N=subtracted: ',
     >   LUNKBD, ADDED, DEFAULT, QUIT)

      IF (ADDED) THEN
         SIGNU = +1.E+0
      ELSE
         SIGNU = -1.E+0
      END IF

      IF (BLAYER) THEN
         SIGNL = -SIGNU
      ELSE
         SIGNL = +SIGNU
      END IF

      DO I = 1, NU2ND
         YU(I) = YU(I) + SIGNU * YU2ND(I)
      END DO

      DO I = 1, NL2ND
         YL(I) = YL(I) + SIGNL * YL2ND(I)
      END DO

      RETURN

C ... Error handling:

  830 WRITE (LUNCRT, 1001)
     >   ' COMBINE: Error reading secondary coordinates - quitting.'
      GO TO 990

  840 WRITE (LUNCRT, 1001)
     >   ' COMBINE: The two sets of coordinates must have the same Xs',
     >   '          for corresponding surfaces.',
     >   '          Use PROFILE''s REDISTRIBUTE option then try again.'
  990 STOP

C ... Formats:

 1001 FORMAT (/, (A))

      END SUBROUTINE COMBINE
C+----------------------------------------------------------------------
C
      SUBROUTINE DUMMY (IENTRY, N, NLDIM, X, F, G, H, L, D, NFTOTL,
     >                  NITER, NTYPE, CONTIN)
C
C
C     Description and Usage:
C
C           This is a null version of the user-supplied routine that is
C        required by QNMDIF.  A more meaningful routine can be used, if
C        necessary, to monitor the progress of the optimization.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C
C                  (see QNMDIF header...)
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     History:
C
C        22 Mar. 1983    RAK    Initial coding.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

      LOGICAL
     >   CONTIN

      INTEGER
     >   IENTRY, N, NLDIM, NFTOTL, NITER, NTYPE

      REAL
     >   X, F, G, H, L, D

      END SUBROUTINE DUMMY
C+----------------------------------------------------------------------
C
      SUBROUTINE GETBUMPS (LUNCRT, LUNKBD, LUNBMP, MXPARM, MXBUMPS,
     >                     NBUMPS, BNAMES, PARAMS, PNAMES, ACTIVE,
     >                     NACTIVE, SCALES)
C
C PURPOSE: GETBUMPS returns a user-specified "bump" set description  in
C          the form of a list of bump (shape function) names, a corres-
C          ponding list of values for all of the parameters  associated
C          with the set, and an indication of which of these parameters
C          are to be considered variable (active) or fixed  (inactive).
C          This is useful when  numerical  optimization  techniques are
C          employed to perturb airfoil profiles.
C
C METHOD:  EITHER -
C          The first N Wagner functions may be invoked:   N is prompted
C          for; scale factors are assigned empirically, of order 1000.;
C          and the initial multipliers are all set to  zero, meaning no
C          "bumps" on the first function evaluation.
C
C          OR -
C          A previously-prepared file is invoked, for keyword-style in-
C          put of any of the bump functions known to subroutine  BEVAL.
C          The format of this file is something like this:
C
C          BUMP: SINE
C          CENTER: 0.5   STATUS: INACTIVE
C          WIDTH:  3.0   STATUS: INACTIVE
C          MULTIPLIER: 0.  STATUS: ACTIVE  SCALE: 100.
C
C          BUMP = EXP
C          POWER = 15.  STATUS = FIXED
C          WIDTH = 10.  STATUS = FIXED   SCALE = 1.
C          MULT = .001  STATUS = VARIABLE  SCALE: 10.
C
C             .
C             .
C             .
C
C          Points to note:
C
C            *  Free format, with several possible delimiters.
C            *  Blank lines are optional.
C            *  Keywords need to be long enough to be unambiguous.
C            *  Any BUMP keyword must be the first for that bump,  and
C               on its own.
C            *  The ordering of subsequent lines describing  a  bump's
C               variables or parameters is unimportant,  but  if  some
C               are omitted,  they will be detected as undefined,  and
C               execution will halt.  There  is  NO ATTEMPT TO DEFAULT
C               values (yet).
C            *  Variable names for a given bump  (e.g. WIDTH)  must be
C               the FIRST keyword on a line, one per line.
C            *  Either ordering of STATUS  and SCALE  within a line is
C               acceptable.
C            *  SCALE is optional if STATUS=FIXED/INACTIVE/CONSTANT.
C            *  SCALE is defaulted if STATUS  is  ACTIVE/FREE/VARIABLE
C               and no entry is given.  (1.0 is the likely default.)
C            *  EOF is used to signal end of data.   Some special key-
C               word would probably be necessary on a mainframe.
C
C ARGUMENTS:
C    ARG    DIM    TYPE I/O/S DESCRIPTION
C   LUNCRT   -       I    I   Logical unit for screen prompts.
C   LUNKBD   -       I    I   Logical unit for keyboard entries.
C   LUNBMP   -       I    I   Logical unit for reading bump info.
C   MXPARM   -       I    I   Max. # parameters  for any one bump.
C   MXBUMPS  -       I    I   Max. # bumps allowed for.
C   NBUMPS   -       I    O   Actual # bumps returned.
C   BNAMES  NBUMPS  C*(*) O   Bump names found.
C   PARAMS  MXPARM,  R    O   PARAMS(1:?,J) are the parameters defining
C           NBUMPS            the Jth bump.
C   PNAMES  MXPARM, C*(*) O   Names of these parameters.
C           NBUMPS
C   ACTIVE  MXPARM,  L    O   ACTIVE(I,J)=.TRUE. if the Ith parameter of
C           NBUMPS            of the Jth bump is to be treated as active
C                             (variable).  Otherwise, the parameter is
C                             either to remain at the value returned here
C                             (fixed), or bump J has fewer than I parameters.
C   NACTIVE  -       I    O   Number of active parameters found.
C   SCALES  MXPARM,  R    O   Scale factors, to be applied multiplicatively
C           NBUMPS            to the active parameters returned here (and
C                             divided out prior to any evaluation of the bumps).
C
C FILES USED: See argument list.
C
C SIGNIFICANT LOCAL VARIABLES:
C    VAR     DIM    TYPE  DESCRIPTION
C   ATTRS     3     C*6   Dictionary of names of "attributes" of a variable
C   BUMPS   NNAMES  C*11  Dictionary of bump function names
C   DICT      4     C*10  Variable dictionary of valid variable names
C   FACTORS MXWAGNER R    Empirical factors used to derive scale factors
C                         for Wagner functions N=2,3,4,.. from the value
C                         entered by the user for N=1
C   LINE      -     C*80  Buffer for one record of bump-description file
C   STATI   NSTATI  C*17  Dictionary of "STATUS" values (e.g., 'ACTIVE')
C
C EXTERNAL REFERENCES:
C   GETLINE  Reads one line; handles comments.
C   LOOKUP   Looks up given key in given dictionary - allows partial
C            matches, and checks for ambiguities.
C   OPENER   File opening utility.
C   PAIRS    Parses a string of keywords + values text into pair(s).
C   READER   Prompting utility with multiple entry points.
C
C ERROR HANDLING:
C     Execution halts with a diagnostic if more bumps or more parameters
C     are found than have been provided for by the calling program. Com-
C     prehensive handling of keyword-driven bump descriptions includes
C     checking for invalid keywords, ambiguous partial keywords, and
C     undefined parameters.
C
C HISTORY:
C   01/20/84   DAS   Initial implementation (Wagner functions only).
C   01/25/84   DAS   Introduced scaling of the variables.
C   02/23/84   DAS   Introduced keyword entry of bump description.
C   03/16/84   DAS   Scale factors for Wagner functions are now pseudo-
C                    variable.  (One can always resort to the keyword-
C                    driven scheme if greater flexibility is needed.)
C   09/24/86   DAS   Clarified use of '|' dictionary entry; minor clean-up.
C   10/24/88   DAS   Removed a scale-factor prompt that showed up undesirably
C                    in the BPLOT application (Wagner fns.; no great loss).
C   02/14/90   DAS   Installed OPENER.
C   12/03/91   DAS   Installed GETLINE; updated coding style.
C   06/19/96   DAS   Added SIN1 and SIN2 functions; EXPONENTIAL function
C                    now expects CENTER, not POWER.
C   12/19/96   DAS   Added SINF and four "COS" functions.
C   10/20/99   DAS   Added SIN3 and SIN4 functions.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, LUNBMP, MXBUMPS, MXPARM, NACTIVE, NBUMPS

      REAL
     >   PARAMS (MXPARM, MXBUMPS), SCALES (MXPARM, MXBUMPS)

      LOGICAL
     >   ACTIVE (MXPARM, MXBUMPS)

      CHARACTER
     >   BNAMES (MXBUMPS) * (*), PNAMES (MXPARM, MXBUMPS) * (*)

C     Local constants:

      INTEGER, PARAMETER ::
     >   MXWAGNER = 20

      REAL, PARAMETER ::
     >   UNDEF = 999.E+0

      LOGICAL, PARAMETER ::
     >   ORDERED = .TRUE.

      CHARACTER, PARAMETER ::
     >   COMMENT * 1 = '!'

C     Local variables:

      INTEGER
     >   I, IENTRY, IK, IOS, IP, J, LAST, N, NEEDED, NKEYS, NP

      REAL
     >   FACTORS (MXWAGNER)

      LOGICAL
     >   DEFAULT, QUIT

      CHARACTER
     >   DICT (3) * 10, FILENAME * 48, KEYS (3) * 10, LINE * 80,
     >   MULTIP * 10, VALUES (3) * 20

C     Procedures:

      EXTERNAL
     >   GETLINE, LOOKUP, OPENER, PAIRS, READI

C     Local character arrays used as static dictionaries:

      INTEGER, PARAMETER ::
     >   NNAMES = 16, NSTATI = 6

      CHARACTER
     >   ATTRS (2) * 6, BUMPS (NNAMES) * 11, STATI (NSTATI) * 17

C ... Helpful comments for future users of LOOKUP (this being the first
C     application):
C
C     The dictionaries must be CAPITALIZED, and LEFT-JUSTIFIED, with no
C     entry exceeding the length indicated in  the  declaration.  Where
C     possible, entries should be ALPHABETIZED for efficiency.  The use
C     of dictionary lookups permits, among other things, partial match-
C     ing and trapping of ambiguities and invalid keywords in a modular
C     sort of way.  The application code runs the risk of becoming NON-
C     mnemonic - it tends to deal with integer item numbers in the dic-
C     tionary rather than with quoted strings. However, the application
C     programmer has the option of staying with strings for readability
C     if efficiency is not an issue, by comparing  DICT (ENTRY) against
C     the (full) string of interest rather than comparing ENTRY against
C     the corresponding item number of interest.    For dictionaries of
C     known length (as in the applications here), the terminating entry
C     '|' is unnecessary.  See LOOKUP for use of '|' to end searches.
C
C     Synonyms may be handled in two ways.   The way that leads to just
C     one result from  LOOKUP  (regardless of which synonym the key re-
C     presents) is illustrated by the "STATI" dictionary.  Briefly, one
C     keyword is treated as fundamental;  any synonyms are entered with
C     this fundamental key appended using any of the delimiters  recog-
C     nized by SCANNR.   This can mean the fundamental string shows  up
C     multiple times, but handling the return from LOOKUP is facilitat-
C     ed.  (The alternative involves checking more than one possibility
C     on return from LOOKUP.)

      DATA
     >   ATTRS  /'SCALE', 'STATUS'/,
     >   BUMPS  /'COSL', 'COSR', 'DROOP', 'EXPONENTIAL',
     >           'LCOS', 'LEADING', 'RCOS', 'SCALE',
     >           'SIN1', 'SIN2', 'SIN3', 'SIN4', 'SINE', 'SINF',
     >           'TRAILING', 'WAGNER'/,
     >   MULTIP /'MULTIPLIER'/,
     >   STATI  /'ACTIVE',
     >           'CONSTANT:INACTIVE',
     >           'FIXED   :INACTIVE',
     >           'FREE    :ACTIVE',
     >           'INACTIVE',
     >           'VARIABLE:ACTIVE'/

      DATA
     >   FACTORS/5*1.E+0, 5*3.E+0, 5*1.E+1, 5*3.E+1/

C     These empirical factors for the relative scaling of the Wagner
C     functions N = 1:MXWAGNER should be tuned some day.


C     Execution:

C     Initialize all parameters as inactive - some will never be used:

      NACTIVE = 0

      DO J = 1, MXBUMPS
         DO I = 1, MXPARM
            ACTIVE (I, J) = .FALSE.
            PARAMS (I, J) = UNDEF
            SCALES (I, J) = 1.E+0
         END DO
      END DO

C     Provide for using the first N Wagner functions easily:

  200 CALL READI (LUNCRT, 'Enter N to use Wagner functions 1:N, ' //
     >            'or <CR> to use an input file: ',
     >            LUNKBD, N, DEFAULT, QUIT )

      IF (.NOT. DEFAULT) THEN

         IF (N > MIN (MXBUMPS, MXWAGNER) .OR. N < 1) GO TO 200

            DO J = 1, N
               BNAMES (J)   = 'WAGNER'
               PNAMES (1, J) = 'N'
               PNAMES (2, J) = MULTIP
               PARAMS (1, J) = J
               PARAMS (2, J) = 0.E+0
               ACTIVE (2, J) = .TRUE.
               SCALES (2, J) = FACTORS (J)
            END DO

            NACTIVE = N
            NBUMPS  = N

      ELSE

C        Bump set has been previously prepared:

         FILENAME = 'bplot.inp'
         CALL OPENER (LUNCRT,
     >      'Enter file name for bumps. Default is bplot.inp: ',
     >      LUNKBD, FILENAME, LUNBMP, 'OLD')

         NBUMPS = 0
  400    CONTINUE

C           Read a line, expecting either 'BUMP', EOF, or empty line:
C
            CALL GETLINE (LUNBMP, COMMENT, LINE, LAST, IOS)

            IF (IOS  < 0)  GO TO 999     ! EOF
            IF (IOS  > 0)  GO TO 800     ! System-dependent error code
            IF (LAST == 0) GO TO 400     ! Empty line

            WRITE (LUNCRT, 1002) LINE (1 : LAST)

C           Organize the line as a keyword/value pair:

            NKEYS = 1
            CALL PAIRS (LINE (1 : LAST), NKEYS, KEYS, VALUES)

            IF (KEYS (1) /= 'BUMP') GO TO 801

C           Identify the bump function name:

            CALL LOOKUP (NNAMES, BUMPS, ORDERED, VALUES (1), IENTRY)
            IF (IENTRY < 1) GO TO 803

            NBUMPS = NBUMPS + 1
            BNAMES (NBUMPS) = BUMPS (IENTRY)

C           Set up a (short) dictionary of variable names for the
C           appropriate bump.  In order to handle many bumps the same
C           way, the order of this dictionary MUST be the order expected
C           by subroutine BEVAL.  This means the dictionary cannot also
C           be alphabetized - hence the two ways of searching.  Note that
C           the items in the CASE statement must be in the same order
C           as they are in the BUMPS dictionary (which IS alphabetized).

            GO TO (410, 420, 430, 440, 450, 460, 470, 480, 490, 500,
     >             510, 520, 530, 540, 550, 560) IENTRY

  410          CONTINUE  ! 'COSL' bump:  (1/4 cosine, peak at left) ** POWER

               NEEDED   = 2
               DICT (1) = 'POWER'
               DICT (2) = MULTIP
               GO TO 600

  420          CONTINUE  ! 'COSR' bump:  (1/4 cosine, peak at right) ** POWER

               NEEDED   = 2
               DICT (1) = 'POWER'
               DICT (2) = MULTIP
               GO TO 600

  430          CONTINUE  ! 'DROOP' bump:  (1 - x) * exp (-WIDTH*x) * MULT

               NEEDED   = 2
               DICT (1) = 'WIDTH'
               DICT (2) = MULTIP
               GO TO 600

  440          CONTINUE
C              'EXPONENTIAL' bump:  MULT * x**P (1-x) exp (-WIDTH*x) /
C                                   CENTER**P (1-CENTER) exp (-WIDTH*CENTER)
               NEEDED   = 3
               DICT (1) = 'CENTER'
               DICT (2) = 'WIDTH'
               DICT (3) = MULTIP
               GO TO 600

  450          CONTINUE  ! 'LCOS' bump:  (1/4 cosine, peak at left) ** POWER

               NEEDED   = 2
               DICT (1) = 'POWER'
               DICT (2) = MULTIP
               GO TO 600

  460          CONTINUE  ! 'LEADING' bump:  (1-x)**POWER * MULT

               NEEDED   = 2
               DICT (1) = 'POWER'
               DICT (2) = MULTIP
               GO TO 600

  470          CONTINUE  ! 'RCOS' bump:  (1/4 cosine, peak at right) ** POWER

               NEEDED   = 2
               DICT (1) = 'POWER'
               DICT (2) = MULTIP
               GO TO 600

  480          CONTINUE  ! 'SCALE' bump:  y = FACTOR * y

               NEEDED   = 1
               DICT (1) = 'FACTOR'
               GO TO 600

  490          CONTINUE  ! 'SIN1' bump (L/L form of SINE)

               NEEDED   = 3
               DICT (1) = 'CENTER'
               DICT (2) = 'WIDTH'
               DICT (3) = MULTIP
               GO TO 600

  500          CONTINUE  ! 'SIN2' bump (R/R form of SINE)

               NEEDED   = 3
               DICT (1) = 'CENTER'
               DICT (2) = 'WIDTH'
               DICT (3) = MULTIP
               GO TO 600

  510          CONTINUE  ! 'SIN3' bump (fore/aft symmetry)

               NEEDED   = 3
               DICT (1) = 'CENTER'
               DICT (2) = 'WIDTH'
               DICT (3) = MULTIP
               GO TO 600

  520          CONTINUE  ! 'SIN4' bump (fore/aft symmetry)

               NEEDED   = 3
               DICT (1) = 'CENTER'
               DICT (2) = 'WIDTH'
               DICT (3) = MULTIP
               GO TO 600

  530          CONTINUE  ! 'SINE' bump (unsymmetric, L/R):
                         !  SIN (PI*x** (LOG (.5)/LOG (CENTER))**WIDTH * MULT

               NEEDED   = 3
               DICT (1) = 'CENTER'
               DICT (2) = 'WIDTH'
               DICT (3) = MULTIP
               GO TO 600

  540          CONTINUE  ! 'SINF' bump (flipped SINE, unsymmetric, R/L)

               NEEDED   = 3
               DICT (1) = 'CENTER'
               DICT (2) = 'WIDTH'
               DICT (3) = MULTIP
               GO TO 600

  550          CONTINUE  ! 'TRAILING' bump:  x**POWER * MULT

               NEEDED   = 2
               DICT (1) = 'POWER'
               DICT (2) = MULTIP
               GO TO 600

  560          CONTINUE  ! 'WAGNER' function:

               NEEDED   = 2
               DICT (1) = 'N'
               DICT (2) = MULTIP
               GO TO 600

  600       CONTINUE
            NP = 1

  610       CONTINUE

C              Read a line, expecting a variable name as first keyword:

               CALL GETLINE (LUNBMP, COMMENT, LINE, LAST, IOS)

               IF (IOS  < 0)  GO TO 805     ! EOF
               IF (IOS  > 0)  GO TO 800     ! System-dependent error code
               IF (LAST == 0) GO TO 610     ! Empty line

               WRITE (LUNCRT, 1002) LINE (1 : LAST)

               NKEYS = 3
               CALL PAIRS (LINE (1 : LAST), NKEYS, KEYS, VALUES)

               CALL LOOKUP (NEEDED, DICT, .NOT. ORDERED, KEYS (1),
     >                      IENTRY)
               IF (IENTRY < 1) GO TO 807

C              Get the value of the variable into PARAMS (*,*) in the
C              order implied by the dynamic dictionary:

               READ (VALUES (1), 1020, ERR=809) PARAMS (IENTRY, NBUMPS)
               PNAMES (IENTRY, NBUMPS) = KEYS (1)

C              Now process any remaining keywords from this line -
C              normally STATUS and SCALE (either order). Special
C              cases: if STATUS is INACTIVE, SCALE is OPTIONAL, and
C              if STATUS is omitted, it defaults to INACTIVE.

               IP = IENTRY
               IK = 2
  700          IF (IK <= NKEYS) THEN

C                 What "attribute" does the next keyword represent?

                  CALL LOOKUP (2, ATTRS, ORDERED, KEYS (IK), IENTRY)
                  IF (IENTRY < 1) GO TO 811

                  IF (ATTRS (IENTRY) == 'STATUS') THEN

C                    Several synonymous values are permitted:

                     CALL LOOKUP (NSTATI, STATI, ORDERED, VALUES (IK),
     >                            IENTRY)
                     IF (IENTRY < 1) GO TO 813

                     IF (STATI (IENTRY) == 'ACTIVE') THEN
                        ACTIVE (IP, NBUMPS) = .TRUE.
                     ELSE  ! Parameter is fixed - ignore any scale factor:
                        IK = 3
                     END IF

                  ELSE IF (ATTRS (IENTRY) == 'SCALE') THEN

C                    Decode the value of the scale factor:

                     READ (VALUES (IK), 1020, ERR = 815)
     >                  SCALES (IP, NBUMPS)

                  END IF

C                 Check for another keyword on this line (3 at most):

                  IK = IK + 1
                  GO TO 700
               END IF

C              Check for more variables needed:

               NP = NP + 1
               IF (NP <= NEEDED)
     >      GO TO 610

C           Still no guarantee the same variable wasn't entered twice.
C           Check for undefined variables for this bump function, and
C           update the number of active variables identified so far:

            DO NP = 1, NEEDED
               IF (PARAMS (NP, NBUMPS) == UNDEF) GO TO 817
               IF (ACTIVE (NP, NBUMPS)) NACTIVE = NACTIVE + 1
            END DO

C           Look for further bumps:

            IF (NBUMPS < MXBUMPS)
     >   GO TO 400

      END IF

      GO TO 999

C     Error Handling:

  800 WRITE (LUNCRT, 1005) IOS
      GO TO 990
  801 WRITE (LUNCRT, 1003) 'Bad keyword where "BUMP" expected: ',
     >                     KEYS (1)
      GO TO 990
  803 WRITE (LUNCRT, 1003) 'Invalid bump name: ', VALUES (1)
      GO TO 990
  805 WRITE (LUNCRT, 1002) 'Unexpected end of data.'
      GO TO 890
  807 WRITE (LUNCRT, 1003) 'Invalid variable name for this bump: ',
     >                     KEYS (1)
      GO TO 890
  809 WRITE (LUNCRT, 1003) 'Bad value for this variable: ', VALUES (1)
      GO TO 890
  811 WRITE (LUNCRT, 1003) 'Bad keyword: ', KEYS (IK)
      GO TO 890
  813 WRITE (LUNCRT, 1003) 'Bad value for STATUS keyword: ',
     >                     VALUES (IK)
      GO TO 890
  815 WRITE (LUNCRT, 1003) 'Bad value for SCALE: ', VALUES (IK)
      GO TO 890
  817 WRITE (LUNCRT, 1002) 'Undefined variable detected:'
      GO TO 890

  890 WRITE (LUNCRT, 1004) NP, NBUMPS
      WRITE (LUNCRT, 1002) 'Reqd. variables:', (DICT (I), I = 1, NEEDED)

  990 STOP

  999 RETURN

C     Formats:

 1002 FORMAT (1X, A)
 1003 FORMAT (1X, A, A)
 1004 FORMAT (' Parameter line', I2, ' of bump number', I3, '.')
 1005 FORMAT (//' *** GETBUMPS: Read error.  IOS = ', I4)
 1020 FORMAT (BN, F20.0)

      END SUBROUTINE GETBUMPS
C+----------------------------------------------------------------------
C
      SUBROUTINE GETCL (NX, XOC, YOC, CP, ALPHA, CL, CD, CM)
C
C
C     Description and Usage:
C
C           Calculate aerodynamic force coefficients for an airfoil by
C        integrating pressure coefficients using the trapezoid method.
C        Adapted from SUBROUTINE FORCF of FLO6 by A. Jameson.
C
C           In this version, airfoil coordinates (XOC, YOC) are assumed
C        to be normalized to the chord, angle of attack is in radians,
C        and reference for the moment coefficient is the quarter chord.
C
C
C     Arguments:
C
C        Name    Dim.     Type   I/O/S   Description
C        NX       -         I      I     No. of airfoil coordinates
C        XOC      NX        R      I     Normalized coords. in wraparound
C        YOC      "         "      "     form (lower t.e. to upper t.e.)
C        CP       NX        R      I     Corresponding Cp values
C        ALPHA    -         R      I     Associated angle of attack
C        CL       -         R      O     Lift coefficient
C        CD       -         R      O     Drag coefficient
C        CM       -         R      O     Moment coefficient (1/4-chord)
C
C
C     Author:  Robert Kennelly, Informatics Inc., 12 Feb. 1982
C
C-----------------------------------------------------------------------


C     Arguments:
C     ----------

      INTEGER NX
      REAL    XOC(NX), YOC(NX), CP(NX), ALPHA, CL, CD, CM

C     Execution:
C     ----------

      XREF = .25E+0
      CL = 0.E+0
      CD = 0.E+0
      CM = 0.E+0

      DO I = 1, NX - 1
         DX = (XOC(I+1) - XOC(I))
         DY = (YOC(I+1) - YOC(I))
         XA = .5E+0 * (XOC(I+1) + XOC(I))
         YA = .5E+0 * (YOC(I+1) + YOC(I))
         CPA = .5E+0 * (CP(I+1) + CP(I))
         DCL = -CPA * DX
         DCD = CPA * DY
         DCM = DCD * YA - DCL * (XA - XREF)
         CL = CL + DCL
         CD = CD + DCD
         CM = CM + DCM
      END DO

C     Express forces in "lab frame" - rotate by ALPHA.

      CLTEMP = CL * COS (ALPHA) - CD * SIN (ALPHA)
      CD = CL * SIN (ALPHA) + CD * COS (ALPHA)
      CL = CLTEMP

      END SUBROUTINE GETCL
C+----------------------------------------------------------------------
C
      SUBROUTINE GETCPS (NPTS, X, YU, YL, ALPHA, FSMACH, XC, CPU, CPL)
C
C ACRONYM: GET (approximate) Cps for an airfoil
C
C PURPOSE: GETCPS is the subroutine version of program FOIL, implemented
C          by Ilan Kroo (NASA Ames/Stanford U.) as a cheap way of gener-
C          ating Cp estimates for a thick airfoil at an angle of attack,
C          using a vortex and source method.    This version was adapted
C          for use by program PROFILE, and may find other uses.
C
C METHOD:  The same abscissas are expected for upper and lower surfaces,
C          for simplicity here.   The coordinates are not necessarily in
C          normalized form.
C
C          The airfoil is modeled with discrete sources and vortices  on
C          the X-axis at 1/4-panel locations.   The thickness effect  is
C          computed from linearized theory with Riegel's correction.  An
C          N*N dense system (N=NPTS-1) is solved for vortex strengths at
C          each coordinate, from which Cp estimates are derived.   These
C          are then corrected for the given free stream Mach number.
C
C          The number of scratch arrays is large enough that no  attempt
C          has been made to pass work-space through the  argument  list.
C          This version is largely self-contained as a result.
C
C REFERENCES:
C          *  Ilan Kroo, Advanced Aerodynamic Concepts Branch (RAC),
C             NASA Ames Research Center, Moffett Field, CA 94035:
C             Forthcoming note on the basic algorithm (1984)
C          *  Schlichting, H. and Truckenbrodt E.,  "Aerodynamics of the
C             Airplane", 1979, McGraw-Hill: Thickness effects (following
C             Riegels).
C
C ARGUMENTS:
C    ARG    DIM     TYPE I/O/S DESCRIPTION
C
C   NPTS     -        I    I   Number of coordinates for each surface
C     X     NPTS      R    I   Abscissas common to both surfaces
C   YU,YL   NPTS      R    I   Ordinates for upper and lower surfaces
C   ALPHA    -        R    I   Desired angle of attack (degrees)
C   FSMACH   -        R    I   Desired free stream Mach number
C    XC     NPTS-1    R    O   Mid-points of panels
C  CPU,CPL  NPTS-1    R    O   Cp estimates at mid-point locations
C
C PARAMETER CONSTANTS:
C
C   MAXN      I     Maximum size of N = NPTS-1 handled by local arrays
C
C EXTERNAL REFERENCES:
C  DECOMP     LU-decomposition of square matrix
C  SOLVE      Solution of factorized linear system
C
C HISTORY:
C   12/03/84  D.A.Saunders   Adapted as subroutine for use in PROFILE.
C   02/21/85    "     "      Switched to version of DECOMP that estimates
C                            the matrix condition number
C   10/23/91    "     "      Switched back to FORSYTHE version of DECOMP
C                            because of a clash with another PROFILE option.
C
C AUTHOR: Ilan Kroo, NASA Ames/Stanford University, Calif.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C ... Local constants:

      INTEGER, PARAMETER ::
     >   MAXN = 201

      REAL, PARAMETER ::
     >   HALF = 0.5, ONE = 1.0, R14 = 0.25, R34 = 0.75, ZERO = 0.

C ... Arguments:

      INTEGER
     >   NPTS
      REAL
     >   X (NPTS), YU (NPTS), YL (NPTS), ALPHA, FSMACH,
     >   XC (NPTS-1), CPU (NPTS-1), CPL (NPTS-1)

C ... Local variables:

      INTEGER
     >   I, J, N
      INTEGER
     >   IP (MAXN)
      REAL
     >   CONST, COSA, DR, PI, SINA, SUM, VTANL, VTANU,
     >   AIC (MAXN, MAXN), DX (MAXN), DYDXM (MAXN),
     >   DYDXT (MAXN), GAMMA (MAXN), RF (MAXN), VTS (MAXN),
     >   XCTL (MAXN), XVORT (MAXN), YM (MAXN+1), YT (MAXN+1)
      LOGICAL
     >   NONZRO

C ... External references:

      EXTERNAL
     >   DECOMP, SOLVE

C ... Execution:

      N = NPTS - 1
      IF (N > MAXN) THEN
         WRITE (*, '(/, A)')
     >      ' *** GETCPS: Too many airfoil coordinates. ***'
         GO TO 99
      END IF

      PI = 4.E+0 * ATAN (ONE)
      DR = PI / 180.E+0

C ... Generate the vortex locations (1/4-panel) and normal vectors:

      DO I = 1, N
         XVORT (I) = X (I) * R34 + X (I+1) * R14
         XCTL (I)  = X (I) * R14 + X (I+1) * R34
         XC (I)    = X (I) * HALF+ X (I+1) * HALF
         DX (I) =  (X (I+1) - X (I))
         YM (I) =  (YU (I) + YL (I)) * HALF
         YT (I) =  (YU (I) - YL (I)) * HALF
      END DO

      YM (N+1) =  (YU (N+1) + YL (N+1)) * HALF
      YT (N+1) =  (YU (N+1) - YL (N+1)) * HALF

      DO I = 1, N
         DYDXT (I) =  (YT (I+1) - YT (I)) / DX (I)
         DYDXM (I) =  (YM (I+1) - YM (I)) / DX (I)
      END DO

C ... The thickness effect is computed from linearized theory
C     with Riegel's correction:

      DO I = 1, N
         SUM = ZERO
         DO J = 1, N
            IF ((XC (I) - X (J)) * (XC (I) - X (J+1)) <= ZERO)
     >         CYCLE
            SUM = SUM + DYDXT (J) *
     >         LOG ((XC (I) - X (J)) / (XC (I) - X (J+1)))
         END DO
         RF (I)  = ONE / SQRT (DYDXT (I) ** 2 + ONE)
         VTS (I) = (ONE + SUM / PI) * RF (I)
      END DO

C ... Set up the influence coefficients of the vortices.  The
C     flow tangency conditions define the right-hand-side:

      SINA = SIN (ALPHA * DR)
      COSA = COS (ALPHA * DR)
      NONZRO = .FALSE.

      DO I = 1, N
         DO J = 1, N
            AIC (I, J) = HALF / (PI * (XCTL (I) - XVORT (J)))
         END DO
         AIC (I, I) = AIC (I, I) + HALF * DYDXT (I) / DX (I)
         GAMMA (I) = SINA - DYDXM (I) * VTS(I) * COSA
         IF (GAMMA (I) /= ZERO) NONZRO = .TRUE.
      END DO

      IF (NONZRO) THEN

C ...    Solve the system for the vortex strengths (else they are zeros):

         CALL DECOMP (N, MAXN, AIC, IP)

         IF (IP (N) == 0) THEN
            WRITE (*, '(/, A)') ' *** GETCPS: Matrix is singular ***'
            GO TO 99
         END IF

         CALL SOLVE (N, MAXN, AIC, GAMMA, IP)

      END IF

C ... Cp is given by  1 - (Vlocal/U)^2.
C     The Karman-Tsien correction for subsonic Mach numbers is applied.

      CONST = SQRT (ONE - FSMACH ** 2)
      DO I = 1, N
         VTANU = VTS (I) * COSA + HALF * RF (I) * GAMMA (I) / DX (I)
         VTANL = VTS (I) * COSA - HALF * RF (I) * GAMMA (I) / DX (I)
         CPU (I) = ONE - VTANU ** 2
         CPL (I) = ONE - VTANL ** 2
         CPU (I) = CPU (I) / (CONST + HALF * (ONE - CONST) * CPU (I))
         CPL (I) = CPL (I) / (CONST + HALF * (ONE - CONST) * CPL (I))
      END DO

   99 RETURN

      END SUBROUTINE GETCPS
C+----------------------------------------------------------------------
C
      SUBROUTINE GETDIS (LUNCRT, LUNKBD, MODE, NP, P, NPTS)
C
C  PURPOSE:  GETDIS does the prompting needed for most uses of DSTRIB,
C            and possibly other 1-D grid utilities.
C
C  ARGUMENTS:
C    ARG    TYPE  I/O/S   DIM     DESCRIPTION
C    LUNCRT   I     I      -      Logical unit for screen prompts.
C    LUNKBD   I     I      -        "   "   "   "  keyboard responses.
C    MODE     I     O      -      MODE with which DSTRIB is to be used
C                                 (-1, 0, 1, 2, or 3 - see DSTRIB - or
C                                 4 for the Vinokur distribution, or
C                                 5 for the sinusoid + quadratic). If
C                                 MODE = -99, calling program can quit.
C    NP       I    I/O     -      Number of distribution params. needed.
C                                 Input with maximum room provided (at
C                                 least 1);  output with actual length
C                                 of P (*) if MODE > 0, else not set.
C    P        R     O      NP     Parameters prompted for.  See NP.
C    NPTS     I     O      -      No. of points in desired distribution.
C
C  PROCEDURES:
C    READER      Prompting utility
C
C  HISTORY:
C  05/26/85    DAS    Original design and code (DSTRIB options only).
C  06/17/85    DAS    Switched to MODE=-99, not NPTS<0, to mean "quit".
C  05/29/90    DAS    Clarified internal point prompt (X, not I).
C  05/04/92    DAS    Provided for the Vinokur distribution.
C  04/01/95    DAS    Provided for FOILGRID.
C  10/18/96    DAS    Replaced FOILGRID with FOILGRD.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, MODE, NP, NPTS
      REAL
     >   P (NP)

C     Local variables:

      CHARACTER
     >   YESNO * 1
      LOGICAL
     >   DEFAULT, QUIT

C     Procedures:

      EXTERNAL
     >   READI, READR

C     Execution:

      WRITE (LUNCRT, '(A)') ' ',
     >' Options (5 with indicated defaults is suggested for airfoils):',
     >  '   -1:  Xs or Ts are to be read from a file, not generated',
     >  '    0:  Uniform distribution',
     >  '    1:  Sinusoidal bunching towards the lower end',
     >  '    2:  Symmetric sinusoidal bunching towards both ends',
     >  '    3:  Sinusoidal bunching around an internal point',
     >  '    4:  Vinokur distribution (1st & last increment specified)',
     >  '    5:  Linear + Quadratic + Sine + Cosine combination',
     >  ' '

      MODE = 5
      CALL READI (LUNCRT,
     >   'Enter distribution choice. <CR> = 5; ^Z (or ^D) = quit: ',
     >   LUNKBD, MODE, DEFAULT, QUIT)

      IF (QUIT) THEN
         MODE = -99
         GO TO 999
      END IF

      IF (MODE == -1) GO TO 999

      NPTS = 100
      CALL READI (LUNCRT, 'How many points?  <CR> gives 100: ',
     >            LUNKBD, NPTS, DEFAULT, QUIT)

      IF (MODE == 0) THEN

C  *     Uniform distribution - no other parameters needed:

      ELSE IF (MODE >= 1 .AND. MODE <= 3) THEN

C  *     Sinusoidal-type distributions require a "WIDTH"-type exponent:

         WRITE (LUNCRT, '(A)') ' ',
     >      ' Exponents > 1.0 give bunching in further regions;',
     >      ' fractional exponents in the range (0.,1.] do not.',
     >      ' '

         P (1) = 1.E+0
         CALL READR (LUNCRT, 'Enter exponent. (Default is 1.) ',
     >               LUNKBD, P (1), DEFAULT, QUIT)

         IF (MODE == 3) THEN

            NP = 2

C           Default here is misleading - data not necessarily in [0.,1.]

            P (2) = 0.5E+0
            CALL READR (LUNCRT,
     >     'Internal pt. (X or T) about which pts. are to be bunched: ',
     >         LUNKBD, P (2), DEFAULT, QUIT)

         ELSE

            NP = 1

         END IF

      ELSE IF (MODE == 4) THEN   ! Vinokur distribution

         NP = 2
  400    P (1) = -0.2
         WRITE (LUNCRT, '(A)')
     >      ' First increment?  +ve is absolute; -ve is relative;'
         CALL READR (LUNCRT,
     >      '-r means r% of range; default = 0.2% of range: ',
     >      LUNKBD, P (1), DEFAULT, QUIT)
         IF (P (1) == 0.) GO TO 400

         P (2) = P (1) + P (1)
         CALL READR (LUNCRT,
     >      'Last increment?  Default = twice the first: ',
     >      LUNKBD, P (2), DEFAULT, QUIT)
         IF (P (2) == 0.) GO TO 400

      ELSE IF (MODE == 5) THEN   ! Linear + quadratic + sine + cosine

         NP = 4
         P (1) = 0.04
         P (2) = 0.0
         P (3) = 0.3
         P (4) = 0.66

         WRITE (LUNCRT, '(A)')
     >      ' Defaults for L, Q, S, C terms are 0.04, 0.0, 0.3, 0.66.'
         YESNO = 'Y'
         CALL READC (LUNCRT, 'Take the defaults? (Y/N; <CR>=Y): ',
     >      LUNKBD, YESNO, DEFAULT, QUIT)

         IF (YESNO /= 'Y') THEN
            CALL READR (LUNCRT, 'Weight on    LINEAR term?  [0.04]: ',
     >         LUNKBD, P (1), DEFAULT, QUIT)
            CALL READR (LUNCRT, 'Weight on QUADRATIC term?  [0.00]: ',
     >         LUNKBD, P (2), DEFAULT, QUIT)
            CALL READR (LUNCRT, 'Weight on      SINE term?  [0.30]: ',
     >         LUNKBD, P (3), DEFAULT, QUIT)
            CALL READR (LUNCRT, 'Weight on    COSINE term?  [0.66]: ',
     >         LUNKBD, P (4), DEFAULT, QUIT)
         END IF

      END IF

  999 RETURN

      END SUBROUTINE GETDIS
C+----------------------------------------------------------------------
C
      SUBROUTINE GETSHAPE (LUNCRT, LUNKBD, PROMPT, SUPRES, IBUMP, BNAME,
     >                     NPARAM, PARAMS, PNAMES, BMULT)
C
C  PURPOSE:  GETSHAPE interactively defines a Hicks-Henne-type shape
C            function or "bump", and its parameters.  A multiplier is
C            also returned if argument BMULT is input as zero.  (The
C            application may already have an intended multiplier, so
C            the prompt needs to be suppressible.)  A "quit" request
C            is signalled by NPARAM = 0 on output.
C
C  PROCEDURES:
C     READER  Prompting utility
C     SELECT  Menu utility
C
C  HISTORY:
C     12/21/96  DAS  Initial adaptation of PROFILE's MODIFY/MODSRF module,
C                    including introduction of SELECT menu utility, as
C                    needed for the implicit/explicit smoothing option.
C     12/24/96   "   Added PNAMES argument for use by MODIFY/MODSRF.
C     12/26/96   "   MODIFY's choice of "done" as the default forced adding
C                    this option to the menu.
C     12/18/96   "   Added SINF, COSL, COSR, LCOS, & RCOS functions.
C     10/20/99   "   Added SIN3 (fore/aft symmetry via SIN1 on [0,1])
C                    and   SIN4 (  "   "   "   "   "   "   " each half).
C
C  AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C
C ------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER   LUNCRT, LUNKBD  ! I   Screen & keyboard
      CHARACTER PROMPT * (*)    ! I   Prompt to appear first on screen;
                                !     see usage described in SELECT
      LOGICAL   SUPRES          ! I   SUPRES = T means don't display the
                                !     menu (unless the user asks for it)
      INTEGER   IBUMP           ! I/O "Bump" number matching BNAME if
                                !     BNAME is non-blank (reqd. by SELECT
                                !     for defaulting); IBUMP = 0 on output
                                !     means "no more shape functions"
      CHARACTER BNAME * (*)     ! I/O "Bump" name as expected by BEVAL;
                                !     blank on input means no default;
                                !     must match IBUMP if non-blank
      INTEGER   NPARAM          ! O   # shape fn. parameters, including
                                !     a multiplier; NPARAM <= 3 in BEVAL;
                                !     NPARAM = 0 means user aborted
      REAL      PARAMS (*)      ! O   Parameter values; PARAMS (NPARAM)
                                !     is (generally) the multiplier,
                                !     copied from BMULT if BMULT |= 0.
      CHARACTER PNAMES (*) * 10 ! O   Parameter names for the chosen bump
                                !     as needed for writing to a log file
      REAL      BMULT           ! I   Multiplier, prompted for if BMULT is
                                !     0., else copied to PARAMS (NPARAM)
C-------------------------------------------------------------------------

C  *  Local constants:

      INTEGER, PARAMETER ::
     >   MXBUMP = 19,    ! Max. # bump functions available
     >   MXPARM = 3      ! Max. # parameters per function incl. multiplier

      CHARACTER, PARAMETER ::
     >   BLANK * 1 = ' '

C  *  Local variables:

      INTEGER
     >   I, IPARM (MXBUMP, 3), IPNAME, NPARMS (MXBUMP), NPROMPT

      CHARACTER
     >   LPROMPT * 43, MENU (0:MXBUMP) * 63, PARMS (6) * 10

      LOGICAL
     >   CR, QUIT

C  *  Procedures:

      EXTERNAL
     >   READR, SELECT

C  *  Storage:

      DATA
     >   MENU /
     >' 0.  Done.                                                     ',
     >' 1.  SCALE:  Y <- Y * P1                                       ',
     >' 2.  RAMP:   Y <- X * P1                                       ',
     >' 3.  FLAP:   Y <- Y - (X - P1) * TAN (P2 deg.)  for X > P1     ',
     >' 4.  SLAT:   Y <- Y - (P1 - X) * TAN (P2 deg.)  for X < P1     ',
     >' 5.  TRAIL:  X ** P1                                           ',
     >' 6.  DROOP:  (1 - X) * EXP (-P1 * X)                           ',
     >' 7.  LEAD:   (1 - X) ** P1                                     ',
     >' 8.  EXPO:   (X^P (1-X) EXP (-P2 X)) / (P1^P (1-P) EXP (-P2 P))',
     >' 9.  SINE:   SIN ** P2  of  pi * X ** (LOG (0.5) / LOG (P1))   ',
     >'10.  SINF:   Flipped form of 9 [SINE] - left & right swapped   ',
     >'11.  SIN1:   Symmetric form of 9 [SINE] - left half            ',
     >'12.  SIN2:   Symmetric form of 9 [SINE] - right half           ',
     >'13.  SIN3:   Ensures fore/aft symmetry (SIN1 on [0, 1])        ',
     >'14.  SIN4:   Ensures fore/aft symmetry (SIN1 on each half)     ',
     >'15.  COSL:   1/4 cosine, peak at left; power P1                ',
     >'16.  COSR:   1/4 cosine, peak at left; power P1                ',
     >'17.  LCOS:   1 - 1/4 cosine, peak at left; power P1            ',
     >'18.  RCOS:   1 - 1/4 cosine, peak at left; power P1            ',
     >'19.  WAGNER: Wagner function; P1 = order of term N             '/

      DATA
     >   NPARMS
     >      /1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2/

      DATA
     >   IPARM
     >      /6, 6, 1, 1, 2, 3, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 5,
     >       0, 0, 4, 4, 6, 6, 6, 3, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6,
     >       0, 0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0/

      DATA
     >   PARMS
     >      /'CENTER', 'POWER', 'WIDTH', 'ANGLE', 'ORDER', 'MULTIPLIER'/

C     IPARM (I, J) selects PARMS (J) for shape function I.

      SAVE
     >   IPARM, MENU, PARMS

C  *  Execution:

C     Prompt for a shape function by number or name.

      CALL SELECT (PROMPT, MXBUMP + 1, MENU, SUPRES, LUNCRT, LUNKBD,
     >             IBUMP, BNAME, QUIT)

      IF (QUIT) GO TO 900
      IF (IBUMP == 0) GO TO 950

      NPARAM  = NPARMS (IBUMP)
      NPROMPT = NPARAM

      IF (BMULT /= 0.) THEN
         IF (IPARM (IBUMP, NPARAM) == 6) THEN
            PARAMS (NPARAM) = BMULT
            NPROMPT = NPROMPT - 1
         END IF
      END IF

      DO I = 1, NPROMPT

C  *     Prompt for and read shape function parameters required:

         IPNAME = IPARM (IBUMP, I)
  200    WRITE (LPROMPT, '(5A)')
     >      'Shape function ', MENU (IBUMP) (6 : 13),
     >      '  Enter ', PARMS (IPNAME), ': '

         CALL READR (LUNCRT, LPROMPT, LUNKBD, PARAMS (I), CR, QUIT)

         IF (QUIT) GO TO 900
         IF (CR) GO TO 200

         PNAMES (I) = PARMS (IPNAME)
      END DO

      GO TO 999


  900 NPARAM = 0  ! Quit
      GO TO 999

  950 NPARAM = 1  ! Distinguish "quit (start over)" from "done"

  999 RETURN

      END SUBROUTINE GETSHAPE
C+-----------------------------------------------------------------------
C
      SUBROUTINE IMPSMOOTH (NPTS, X, Y, ISRF, LUNCRT, LUNKBD, LUNOUT)
C
C  PURPOSE:  IMPSMOOTH applies implicit and/or explicit smoothing to one
C            surface of an airfoil, which need not be normalized.
C
C  METHOD:   Basically,
C
C               _
C               y = y + dx**2 eps y"  (explicit; dx**2 for stability)  or
C               _       _
C               y - eps y" = y        (implicit/stable)
C
C            where y" is the second derivative w.r.t. arc length, s,
C            and eps may be a function of s, using a choice of standard
C            shape functions.
C
C            Suggested usage:  implicit with eps ~ 0.01  and/or
C                              explicit with eps ~ 0.1 (1+ iterations)
C
C  ARGUMENTS:
C
C  VAR   DIM   I/O/S   DESCRIPTION
C  NPTS   -    I       Number of points on current surface
C  X     NPTS  I       Surface coordinates
C  Y     NPTS  I/O     (not necessarily normalized)
C  ISRF   -    I       1 means uppper surface;
C                      2 means lower surface
C  LUNCRT -    I       Logical units for screen, keyboard, and
C  LUNKBD -    I       log file
C  LUNOUT -    I
C
C  PROCEDURES:
C
C  BEVAL     Evaluates numerous shape functions on [0., 1]
C  CHORDS2D  Normalized chord lengths
C  READER    Prompting utility
C  TRDIAG    Tridiagonal solver
C
C  HISTORY:
C
C  12/20/96  DAS  Initial form of scheme outlined by James Reuther,
C                 adapted from WAGSMOOTH (formerly SMSURF).
C  12/24/96  DAS  GETSHAPE needed IBUMP and PNAMES arguments.
C  01/31/97  DAS  Steve Edwards pointed to the need for dT ~ dS**2 to
C                 stabilize the explicit method.
C
C  AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C
C-------------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   ISRF, LUNCRT, LUNKBD, LUNOUT, NPTS

       REAL
     >   X (NPTS), Y (NPTS)

C  *  Local constants:

      INTEGER, PARAMETER ::
     >   MXPTS = 256   ! Max. # points per surface handled

      REAL, PARAMETER ::
     >   ONE   = 1.E+0,
     >   ZERO  = 0.E+0

      CHARACTER, PARAMETER ::
     >   BLANK * 1 = ' '

C  *  Local variables:

      INTEGER
     >   I, IBUMP, ITER, NITER, NPARAM

      REAL
     >   A (MXPTS), B (MXPTS), C (MXPTS), EPS (MXPTS), S (MXPTS),
     >   PARAMS (3), AI, BI, CI, BMULT, EPSEXP, EPSIMP, HL, HR, TERM,
     >   XLE, XTE, YLE, YTE

      CHARACTER
     >   BNAME * 4, PNAMES (3) * 10, SURFCE (2) * 5

      LOGICAL
     >   CR, EOF, SAME, SUPRES, UNIFORM

C  *  Procedures:

      EXTERNAL
     >   BEVAL, CHORDS2D, READI, READR, READY, TRDIAG

C  *  Storage:

      SAVE
     >   BNAME, EPSEXP, EPSIMP, IBUMP, NITER, NPARAM, PARAMS

      DATA
     >   SURFCE / 'UPPER', 'LOWER' /

C  *  Execution:

      IF (NPTS > MXPTS) THEN
         WRITE (LUNCRT, '(/, (A, I4))') ' IMPSMOOTH: ', NPTS,
     >      ' points exceeds local limit: ', MXPTS
         STOP
      END IF

      IF (ISRF == 1) THEN

         WRITE (LUNCRT, '(/, 2A)')
     >      ' Enter smoothing controls for the UPPER surface.',
     >      ' (0. = leave it alone.)'
         EPSIMP = 0.01
         CALL READR (LUNCRT,
     >      '[Peak] IMPLICIT smoothing coefficient? (>=0.; [0.01]) ',
     >      LUNKBD, EPSIMP, CR, EOF)

         EPSEXP = 0.1
         CALL READR (LUNCRT,
     >   '[Peak] EXPLICIT smoothing coefficient? (>=0.; [0.10*dS**2]) ',
     >      LUNKBD, EPSEXP, CR, EOF)

         NITER = 0
         IF (EPSEXP /= ZERO) THEN
            NITER = 1
            CALL READI (LUNCRT, '# explicit smoothing iterations? [1] ',
     >      LUNKBD, NITER, CR, EOF)
         END IF

         UNIFORM = .FALSE.
         CALL READY (LUNCRT,
     >      'Uniform smoothing along the chord? (Y/N; [N]) ',
     >      LUNKBD, UNIFORM, CR, EOF)

         IF (.NOT. UNIFORM) THEN

            SUPRES = .FALSE.  ! Don't suppress the menu
            BNAME  = BLANK    ! Suppresses defaulting
            IBUMP  = 0        ! Matches BNAME
            BMULT  = ONE      ! Non-zero suppresses prompt for a multiplier

            CALL GETSHAPE (LUNCRT, LUNKBD,
     >         'Choose a weighting function.',
     >         SUPRES, IBUMP, BNAME, NPARAM, PARAMS, PNAMES, BMULT)

         ELSE
            BNAME = 'N/A '
         END IF

      ELSE  ! ISRF = 2

         IF (EPSIMP == ZERO .AND. EPSEXP == ZERO) THEN

            WRITE (LUNCRT, '(/, A)')
     >         ' Enter smoothing controls for the LOWER surface.'
            SAME   = .FALSE.
            SUPRES = .FALSE.
            BNAME  = BLANK
            IBUMP  = 0

         ELSE

            SAME = .TRUE.
            CALL READY (LUNCRT,
     >         'Smooth LOWER surface the same way? (Y/N; [Y]) ',
     >         LUNKBD, SAME, CR, EOF)

            SUPRES = SAME
         END IF

         IF (.NOT. SAME) THEN
            EPSIMP = 0.01
            CALL READR (LUNCRT,
     >         '[Peak] IMPLICIT smoothing coefficient? (>=0.; [0.01]) ',
     >         LUNKBD, EPSIMP, CR, EOF)

            EPSEXP = 0.01
            CALL READR (LUNCRT,
     >         '[Peak] EXPLICIT smoothing coefficient? (>=0.; [0.01]) ',
     >         LUNKBD, EPSEXP, CR, EOF)

            NITER = 0
            IF (EPSEXP /= ZERO) THEN
               NITER = 1
               CALL READI (LUNCRT,
     >            '# explicit smoothing iterations? [1] ',
     >            LUNKBD, NITER, CR, EOF)
            END IF

            UNIFORM = .FALSE.
            CALL READY (LUNCRT,
     >         'Uniform smoothing along the chord? (Y/N; [N]) ',
     >         LUNKBD, UNIFORM, CR, EOF)

            IF (.NOT. UNIFORM) THEN

               BMULT = ONE
               IF (BNAME == 'N/A ') BNAME = BLANK
               IF (BNAME == BLANK)  IBUMP = 0

               CALL GETSHAPE (LUNCRT, LUNKBD,
     >            'Choose a weighting function.',
     >            SUPRES, IBUMP, BNAME, NPARAM, PARAMS, PNAMES, BMULT)

            ELSE
               BNAME = 'N/A '
            END IF

         END IF

      END IF

      IF (EOF .OR. (EPSIMP == ZERO .AND. EPSEXP == ZERO)) GO TO 900


C     Log the details:

      WRITE (LUNOUT, 1100) SURFCE (ISRF), EPSIMP, EPSEXP, NITER

      IF (.NOT. UNIFORM) THEN
         PARAMS (NPARAM) = ONE  ! May be left over as eps from a previous call
         PNAMES (NPARAM) = 'MULTIPLIER'

         WRITE (LUNOUT, 1200)
     >      BNAME, (PNAMES (I), PARAMS (I), I = 1, NPARAM)
      END IF


C     Calculate normalized chord-lengths.  TERM = total arc is not used.

      CALL CHORDS2D (NPTS, X, Y, .TRUE., TERM, S)

      XLE = X (1)
      YLE = Y (1)
      XTE = X (NPTS)
      YTE = Y (NPTS)

C     Implicit smoothing (for the lower frequencies):

      IF (EPSIMP /= ZERO) THEN

         IF (UNIFORM) THEN
            DO I = 1, NPTS
               EPS (I) = EPSIMP
            END DO
         ELSE
            PARAMS (NPARAM) = EPSIMP
            CALL BEVAL (BNAME, NPARAM, PARAMS, .FALSE., NPTS, S, EPS)
         END IF

         A (1) = ZERO
         B (1) = ONE
         C (1) = ZERO

         DO I = 2, NPTS - 1
            TERM  = -2. * EPS (I) / (S (I + 1) - S (I - 1))
            A (I) = TERM / (S (I) - S (I - 1))
            C (I) = TERM / (S (I + 1) - S (I))
            B (I) = ONE - (A (I) + C (I))
         END DO

         A (NPTS) = ZERO
         B (NPTS) = ONE
         C (NPTS) = ZERO

         CALL TRDIAG (A, B, C, Y, Y, NPTS)

      END IF


C     Explicit smoothing (for the higher frequencies):

      IF (EPSEXP /= ZERO) THEN

         IF (UNIFORM) THEN
            DO I = 1, NPTS
               EPS (I) = EPSEXP
            END DO
         ELSE
            PARAMS (NPARAM) = EPSEXP
            CALL BEVAL (BNAME, NPARAM, PARAMS, .FALSE., NPTS, S, EPS)
         END IF

         DO ITER = 1, NITER

            DO I = 1, NPTS
               A (I) = Y (I)
            END DO

            DO I = 2, NPTS - 1
               HL    = S (I) - S (I - 1)
               HR    = S (I + 1) - S (I)
               TERM  = 2. * EPS (I) * ((MIN (HL, HR) ** 2) / (HL + HR))
               AI    = TERM / HL
               CI    = TERM / HR
               BI    = ONE - (AI + CI)
               Y (I) = AI * A (I - 1) + BI * A (I) + CI * A (I + 1)
            END DO

         END DO

      END IF

      X (1)    = XLE
      X (NPTS) = XTE
      Y (1)    = YLE
      Y (NPTS) = YTE

  900 RETURN

C     Formats:

 1100 FORMAT (//, ' Details of smoothing for ', A5, ' surface:',
     >        //, ' [Peak] implicit coefficient:', F8.3,
     >        /,  ' [Peak] explicit coefficient:', F8.3,
     >        /,  ' Number of explicit iterations: ', I3)
 1200 FORMAT (    ' Nonuniform weighting function: ', A,
     >        /,  ' Weighting function parameters:',
     >        /,  (4X, A10, F10.6))

      END SUBROUTINE IMPSMOOTH
C+----------------------------------------------------------------------
C
      SUBROUTINE LOFT (NU, XU, YU, NL, XL, YL,
     >                 MAXPTS, X, Y, XU2ND, YU2ND, XL2ND, YL2ND,
     >                 LUNCRT, LUNKBD, LUN2ND)
C
C  PURPOSE:  LOFT performs the simplest type of lofting between two wing
C            sections by linear interpolation in the third dimension.
C
C  METHOD:   LOFT is intended to be one of PROFILE's high-level modes of
C            operation.   Therefore a "primary" airfoil is assumed to be
C            read by PROFILE in the usual way, and a "secondary" airfoil
C            is prompted for here.   The latter should be passed through
C            PROFILE's REDISTRIBUTE mode first if necessary, so that the
C            numbers of abscissas match the corresponding numbers on the
C            primary section's upper and lower surfaces.   No attempt is
C            made to do the redistribution here - one look at  REDISTRIB
C            explains why.
C
C  ARGUMENTS:
C   ARG      DIM    TYPE I/O/S DESCRIPTION
C  NU,NL      -       I   I/O  No. of upper/lower surface pts. (primary)
C  XU,XL    NU,NL     R   I/O  Abscissas, upper/lower (primary)
C  YU,YL    NU,NL     R   I/O  Corresponding ordinates, upper/lower.
C                              Surfaces must have common leading edge.
C  MAXPTS     -       I    I   Max. no. pts. provided for on 1 surface.
C  X,Y     2*MAXPTS   R    S   Buffers for reading secondary airfoil.
C  XU2ND,YU2ND MAXPTS R    S   Secondary abscissas found.
C  XL2ND,YL2ND MAXPTS R    S   Secondary ordinates found.
C  LUNCRT     -       I    I   Logical unit for prompts.
C  LUNKBD     -       I    I   Logical unit for responses.
C  LUN2ND     -       I    I   Logical unit for secondary coordinates
C                              in one of the PROFILE formats.
C
C  PROCEDURES:
C    LINTRP   Modularizes constant interpolation between multiple pairs
C    NRMLIZ   Normalizes/denormalizes coordinates
C    OPENER   File opening utility
C    PRREAD   For reading secondary profile
C    RDREALS  Prompts for one or more real values
C
C  HISTORY:
C  06/28/90   DAS   Adapted from COMBINE option.
C  10/21/93   DAS   Lofting of twisted airfoils meant to safeguard of
C                   checking for equal relative Xs was not appropriate.
C                   Just require the same numbers of points.
C
C  AUTHOR:  David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, LUN2ND, MAXPTS, NL, NU
      REAL
     >   X (MAXPTS*2), XL (NL), XU (NU), XL2ND (MAXPTS), XU2ND (MAXPTS),
     >   Y (MAXPTS*2), YL (NL), YU (NU), YL2ND (MAXPTS), YU2ND (MAXPTS)

C ... Local constants:

      REAL, PARAMETER ::
     >   ONE  = 1.E+0,
     >   ZERO = 0.E+0,
     >   TOL  = 1.E-6   ! For comparing normalized Xs

C ... Local variables:

      INTEGER
     >   FORMAT, I, IER, NL2ND, NU2ND, NVALS
      REAL
     >   CHORD, XLE, Z (3)
      CHARACTER
     >   TITLE*80, XY2ND*48

C ... Procedures:

      EXTERNAL
     >   LINTRP, NRMLIZ, OPENER, PRREAD, RDREALS

C ... Execution:

C ... Read secondary coordinates (in any of the PROFILE formats):

      CALL OPENER (LUNCRT,
     >   'Enter file name for secondary coordinates: ',
     >   LUNKBD, XY2ND, LUN2ND, 'OLD')

      CALL PRREAD (LUN2ND, TITLE, MAXPTS, NU2ND, NL2ND, X, Y,
     >   XU2ND, XL2ND, YU2ND, YL2ND, FORMAT, IER)
      IF (IER /= 0) GO TO 830

C ... Check for mismatched numbers of points:

      IF (NU2ND /= NU) GO TO 840
      IF (NL2ND /= NL) GO TO 840

C ... Programmer note:  Use local variables NU2ND, NL2ND everywhere now -
C     it may save a little code.

C     The RELATIVE distributions of abscissas should also match.
C*****NO: not for twisted sections.  Abandon this check.
C     Use X (*) for normalizing the primary Xs; Y (*) for the secondary Xs:

C*****CALL NRMLIZ (NU2ND, XU, X, XU (1), XU (NU2ND) - XU (1))
C*****CALL NRMLIZ (NU2ND, XU2ND, Y, XU2ND (1),
C****>             XU2ND (NU2ND) - XU2ND (1))

C*****DO 410, I = 1, NU2ND
C*****   IF ( ABS (X (I) - Y (I)) > TOL) GO TO 840
C*410 CONTINUE

C*****CALL NRMLIZ (NL2ND, XL, X, XL (1), XL (NL2ND) - XL (1))
C*****CALL NRMLIZ (NL2ND, XL2ND, Y, XL2ND (1),
C****>             XL2ND (NL2ND) - XL2ND (1))

C*****DO 420, I = 1, NL2ND
C*****   IF (ABS (X (I) - Y (I)) > TOL) GO TO 840
C*420 CONTINUE

C ... Allow for the possibility of starting with normalized sections,
C     which may or may not need to be denormalized:

      IF (XU (1) == ZERO  .AND.  YU (1) == ZERO  .AND.
     >    XU (NU2ND) == ONE) THEN

  510    NVALS = 3
         WRITE (LUNCRT, 1001)
     >      ' Enter primary section''s unnormalized (Xle, Yle)' //
     >      ' and chord.'
         CALL RDREALS (LUNCRT, '$(3 values; <CR> = leave normalized): ',
     >      LUNKBD, NVALS, Z)
         IF (NVALS == 1 .OR. NVALS == 2) GO TO 510

         IF (NVALS == 3) THEN
            CALL NRMLIZ (NU2ND, XU, XU, Z (1), -Z (3))
            CALL NRMLIZ (NU2ND, YU, YU, Z (2), -Z (3))
            CALL NRMLIZ (NL2ND, XL, XL, Z (1), -Z (3))
            CALL NRMLIZ (NL2ND, YL, YL, Z (2), -Z (3))
         END IF
      END IF

      IF (XU2ND (1) == ZERO  .AND.  YU2ND (1) == ZERO  .AND.
     >    XU2ND (NU2ND) == ONE) THEN

  520    NVALS = 3
         WRITE (LUNCRT, 1001)
     >      ' Enter secondary section''s unnormalized (Xle, Yle)' //
     >      ' and chord.'
         CALL RDREALS (LUNCRT, '$(3 values; <CR> = leave normalized): ',
     >      LUNKBD, NVALS, Z)
         IF (NVALS == 1 .OR. NVALS == 2) GO TO 520

         IF (NVALS == 3) THEN
            CALL NRMLIZ (NU2ND, XU2ND, XU2ND, Z (1), -Z (3))
            CALL NRMLIZ (NU2ND, YU2ND, YU2ND, Z (2), -Z (3))
            CALL NRMLIZ (NL2ND, XL2ND, XL2ND, Z (1), -Z (3))
            CALL NRMLIZ (NL2ND, YL2ND, YL2ND, Z (2), -Z (3))
         END IF
      END IF

C ... Now for the third dimension...

  600 NVALS = 3
      WRITE (LUNCRT, '(A)')
      CALL RDREALS (LUNCRT,
     >   ' Enter span stations of primary, secondary, and desired ' //
     >   'sections: ', LUNKBD, NVALS, Z)
      IF (NVALS /= 3) GO TO 600

C ... Do the spanwise linear interpolation for X, Y on each surface,
C     overwriting the primary section data with the results:

      CALL LINTRP (NU2ND, XU, XU2ND, Z (1), Z (2), Z (3), XU)
      CALL LINTRP (NU2ND, YU, YU2ND, Z (1), Z (2), Z (3), YU)
      CALL LINTRP (NL2ND, XL, XL2ND, Z (1), Z (2), Z (3), XL)
      CALL LINTRP (NL2ND, YL, YL2ND, Z (1), Z (2), Z (3), YL)

      RETURN

C ... Error handling:

  830 WRITE (LUNCRT, 1001)
     >   ' LOFT: Error reading secondary coordinates - quitting.'
      GO TO 990

  840 WRITE (LUNCRT, 1001)
     >   ' LOFT: The two sections must have the same RELATIVE point',
     >   '       distributions on corresponding surfaces.',
     >   '       PROFILE''s REDISTRIBUTE and NORMALIZE/DENORMALIZE',
     >   '       options can help.  Quitting...'

  990 STOP

C ... Formats:

 1001 FORMAT (/, (A))

      END SUBROUTINE LOFT
C+---------------------------------------------------------------------
C
      SUBROUTINE MAXMIN (LUNRD, MAXPTS, XMIN, XMAX, YMIN, YMAX,
     >                   YLE, NPR, X, Y, XU, XL, YU, YL, LEGND, IER)
C
C  PURPOSE:  MAXMIN calls PRREAD to scan all of the profiles in the
C            given file.  It counts the number of profiles read and
C            finds the maximum and minimum abscissas and  ordinates
C            over all profiles.  (There may be only one profile.)
C
C  ARGUMENTS:
C   ARG    TYPE   I/O/S    DIM     DESCRIPTION
C  LUNRD     I      I       -      Logical unit used in PRREAD
C  MAXPTS    I      I       -      Maximum number of points allowed
C                                  for on a surface
C  XMIN      R      O       -      Minimum abscissa of profile(s)
C  XMAX      R      O       -      Maximum abscissa of profile(s)
C  YMIN      R      O       -      Minimum ordinate of profile(s)
C  YMAX      R      O       -      Maximum ordinate of profile(s)
C  NPR       I      O       -      Number of profiles read
C  XU        R      S     MAXPTS   Upper surface abscissas
C  XL        R      S     MAXPTS   Lower surface abscissas (if any)
C  YU        R      S     MAXPTS   Upper surface ordinates
C  YL        R      S     MAXPTS   Lower surface ordinates (if any)
C  YLE       R      O       -      Ordinate of leading edge point
C  LEGND     C     S/O      -      Legend entry associated with profile(s).
C                                  The last one found may be used in
C                                  the initial plot setup unless the
C                                  run is a "THREED" case.
C  IER       I      O       -      Error return code.  See PRREAD.
C
C  EXTERNAL REFERENCES:
C  MODULE    DESCRIPTION
C  BOUNDS    Determines maximum and minimum values of array(s)
C  PRREAD    Reads one profile
C
C  AUTHOR: Leslie Collins, Informatics General Corporation, Palo Alto, CA
C
C  HISTORY:
C
C    02/15/83   LJC   Original design and coding
C    03/03/83   LJC   Included indefinite scanning loop to determine
C                     overall data range
C    08/03/84   LJC   Added YLE as an argument (previously determined
C                     in main program)
C    02/11/87   DAS   Functionality of BOUNDS changed somewhat; took
C                     out nonstandard DO WHILE.
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   I, IER, LUNRD, MAXPTS, NPR

       REAL
     >    X, XMAX, XMIN, XL(MAXPTS), XU(MAXPTS), Y, YLE, YMIN, YMAX,
     >    YL(MAXPTS), YU(MAXPTS)

      CHARACTER
     >   LEGND * 80

C  *  Local variables:

      INTEGER
     >   FORMAT, NL, NU

C  *  Set up indefinite loop to scan all profiles for overall data range
C     (for axis scaling and possible normalization purposes).  Count the
C     number of profiles as we go:

C  *  Initialize saved minimums and maximums for first comparison:

      XMIN =  1.E10
      XMAX = -XMIN
      YMIN =  XMIN
      YMAX =  XMAX
      NPR  =  0

  200 CONTINUE

C  *     Read a profile:

         CALL PRREAD (LUNRD, LEGND, MAXPTS, NU, NL, [X], [Y], XU, XL,
     >                YU, YL, FORMAT, IER)
         IF (IER /= 0) GO TO 999

         NPR = NPR+1

C  *     Find minimum and maximum X and Y so far.
C        Note that BOUNDS expects input MAX/MINs.

         CALL BOUNDS (NU, 1, MAXPTS, XU, XMIN, XMAX)
         CALL BOUNDS (NU, 1, MAXPTS, YU, YMIN, YMAX)

         IF (NL /= 0) THEN
            CALL BOUNDS (NL, 1, MAXPTS, XL, XMIN, XMAX)
            CALL BOUNDS (NL, 1, MAXPTS, YL, YMIN, YMAX)
         END IF

C  *     Retrieve Y coordinate of point with minimum X (clumsy):

         DO I = 1, NU
            IF (XU(I) == XMIN) THEN
               YLE = YU(I)
               GO TO 900
            END IF
         END DO

         DO I = 1, NL
            IF (XL(I) == XMIN) THEN
               YLE = YL(I)
               GO TO 900
            END IF
         END DO

  900    CONTINUE

C  *     Look for another profile:

      GO TO 200

  999 RETURN

      END SUBROUTINE MAXMIN
C+----------------------------------------------------------------------
C
      SUBROUTINE MODIFY (NU, XU, YU, NL, XL, YL, LUNCRT, LUNKBD, LUNOUT)
C
C  PURPOSE:  MODIFY drives the interactive application of a variety of
C            shape functions to airfoil surfaces.    It simply invokes
C            MODSRF once for each surface.   (Other versions of MODIFY
C            may need to treat the profile as a whole; hence the above
C            calling sequence.)
C
C  ARGUMENTS:
C     VAR   DIM   I/O/S   DESCRIPTION
C     NU     -      I     Number of points on upper surface
C     XU     NU     I     Abscissas for upper surface
C     YU     NU   I/O     Ordinates for upper surface
C     NL     -      I     Number of points on lower surface
C     XL     NL     I     Abscissas for lower surface
C     YL     NL   I/O     Ordinates for lower surface
C     LUNCRT -      I     Logical unit for screen
C     LUNKBD -      I     Logical unit for keyboard
C     LUNOUT -      I     Logical unit for printed output
C
C  PROCEDURES:
C
C     MODSRF   Perturbs an airfoil surface using shape functions
C
C  HISTORY:
C     09/27/83   LJC   Initial coding
C     10/28/83   DAS   Added LUNs as arguments, for consistency
C
C  AUTHOR: Leslie Collins, Informatics, Palo Alto, CA
C
C-------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER  LUNCRT, LUNKBD, LUNOUT, NU, NL
      REAL     XU (NU), XL (NL), YU (NU), YL (NL)
      EXTERNAL MODSRF

      CALL MODSRF (NU, XU, YU, 1, LUNCRT, LUNKBD, LUNOUT)
      CALL MODSRF (NL, XL, YL, 2, LUNCRT, LUNKBD, LUNOUT)

      END SUBROUTINE MODIFY
C+-----------------------------------------------------------------------
C
      SUBROUTINE MODSRF (NPTS, X, Y, ISRF, LUNCRT, LUNKBD, LUNOUT)
C
C  PURPOSE:  MODSRF allows the user to perturb one airfoil surface by
C            selecting shape functions interactively.  The airfoil is
C            assumed to be normalized.
C
C  METHOD:   DO WHILE more shape functions desired
C               Prompt for and read no. of shape function from menu
C               DO I = 1, number of parameters for this shape function
C                  Prompt for and read function parameters
C               END DO
C               Prompt for and read function multiplier
C            END DO
C
C            Evaluate each function at all (normalized) abscissas,
C            using stored parameters, and add to ordinates.
C
C  ARGUMENTS:
C     VAR   DIM   I/O/S   DESCRIPTION
C     NPTS   -      I     Number of points on current surface
C     X     NPTS    I     Abscissas of current surface (normalized)
C     Y     NPTS  I/O     Ordinates of current surface
C     ISRF   -      I     1 means uppper surface;
C                         2 means lower surface.
C     LUNCRT -      I     Logical unit for screen
C     LUNKBD -      I     Logical unit for keyboard
C     LUNOUT -      I     Logical unit for printed output
C
C  FILES USED:  See arguments.
C
C  PROCEDURES:
C     ADDBUMPS   Applies selected bumps to airfoil ordinates
C     GETSHAPE   Modularization of bump selection (used for smoothing too)
C     READER     Prompting utility
C
C  HISTORY:
C     08/05/83   LJC   Initial design and coding
C     09/07/83   LJC   Added echoing of inputs
C     09/27/83   LJC   Replaced statement functions with subroutine BUMP
C     10/04/83   LJC   Introduced READER routine for accepting inputs
C     01/14/84   DAS   Introduced ADDBUMPS in place of in-line code
C     02/22/84   DAS   Eliminated SQRT and SIN bumps (redundant)
C     04/13/84   LJC   Modified to handle ADDBUMPS which expects bump
C                      names now, not integer code numbers
C     12/27/85   DAS   Streamlined prompts and other I/O
C     08/12/86   DAS   Added RAMP, FLAP, and SLAT options
C     01/29/90   DAS   Removed END DOs and DO WHILE
C     11/09/93   DAS   RAK noticed the DROOP fn. was missing a (1-X) factor
C     06/19/96   DAS   The EXPONENTIAL function now has peak 1. at
C                      specified X; introduced symmetric forms of
C                      the modified sine function.
C     12/18/96   DAS   Installed SINF, COSL, COSR, LCOS, & RCOS functions;
C                      no more concatenations in I/O lists.
C     12/24/96   DAS   Introduced GETSHAPE after it was needed for doing
C                      implicit/explicit smoothing. This keeps the shape
C                      function menu in one place (apart from GETBUMPS),
C                      and allows selecting bumps by name.
C
C  AUTHOR: Leslie Collins, Informatics/NASA Ames, Mt. View, CA
C
C-------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   ISRF, LUNCRT, LUNKBD, LUNOUT, NPTS
       REAL
     >   X (NPTS), Y (NPTS)

C     Local constants:

      INTEGER, PARAMETER ::
     >   MXNSEL = 20,   ! Max. # function selections per surface
     >   MXPARM = 3     ! Max. # parameters per function incl. multiplier

      CHARACTER, PARAMETER ::
     >   BLANK * 1 = ' '

C     Local variables:

      INTEGER
     >   I, IBUMP, NBUMP, NPARAM
      REAL
     >   BMULT, PARAMS (MXPARM), PARM (MXPARM, MXNSEL)
      CHARACTER
     >   BNAME * 4, BNAMES (MXNSEL) * 4, PNAMES (MXPARM) * 10,
     >   SURFCE (2) * 13
      LOGICAL
     >   CR, EOF, RETRY, SUPRES

C     Procedures:

      EXTERNAL
     >   ADDBUMPS, GETSHAPE, READY

C     Storage:

      DATA
     >   SURFCE /'UPPER surface', 'LOWER surface'/

      SAVE
     >   SURFCE

C     Execution:

      BMULT = 0.     ! Non-zero suppresses prompt for a multiplier
      RETRY = .FALSE.

  100 IF (ISRF == 1 .OR. RETRY) THEN
         SUPRES = .FALSE.  ! Don't suppress the shape function menu
         WRITE (LUNCRT, '(A)')
      END IF

      IBUMP = 0
      BNAME = 'DONE'  ! Default is "no more shape functions"
      NBUMP = 0

      CALL GETSHAPE (LUNCRT, LUNKBD,
     >   'First shape function to add to the ' // SURFCE (ISRF) // '?',
     >   SUPRES, IBUMP, BNAME, NPARAM, PARAMS, PNAMES, BMULT)

      IF (NPARAM == 0) GO TO 900  ! Abort
      IF (IBUMP  == 0) GO TO 999  ! Done

      WRITE (LUNOUT, '(/, 3A)' )
     >   ' Shape functions applied to the ', SURFCE (ISRF), ':'

C     Loop over packing of bumps and prompting for further shape functions:

  200 CONTINUE

         NBUMP = NBUMP + 1
         BNAMES (NBUMP) = BNAME
         DO I = 1, NPARAM
            PARM (I, NBUMP) = PARAMS (I)
         END DO

C        Echo selection to output file:

         WRITE (LUNOUT, '(/, 4X, 2A)') BNAME, ':'
         WRITE (LUNOUT, '(4X, A, F10.6)')
     >      (PNAMES (I), PARAMS (I), I = 1, NPARAM)

         IF (NBUMP < MXNSEL) THEN

C           Prompt for the next shape function:

            SUPRES = .TRUE.  ! Suppress the shape function menu
            IBUMP  = 0
            BNAME  = 'DONE'

            CALL GETSHAPE (LUNCRT, LUNKBD,
     >         'Next shape function for the ' // SURFCE (ISRF) // '?',
     >         SUPRES, IBUMP, BNAME, NPARAM, PARAMS, PNAMES, BMULT)

            IF (NPARAM == 0) GO TO 900  ! Abort
            IF (IBUMP  >  0) GO TO 200

         END IF


C     Apply the selected bumps to this surface, in-place:

      CALL ADDBUMPS (NPTS, X, Y, MXPARM, NBUMP, BNAMES, PARM, Y)

      GO TO 999


  900 CONTINUE  ! Exit or start over:

      RETRY = .TRUE.
      CALL READY (LUNCRT,
     >   'Do you want to start this surface over? (Y/N/EOF; <CR>=Yes) ',
     >   LUNKBD, RETRY, CR, EOF)
      IF (EOF)   GO TO 999
      IF (RETRY) GO TO 100


  999 WRITE (LUNOUT, '(A)')

      END SUBROUTINE MODSRF
C+------------------------------------------------------------------------------
C
      SUBROUTINE NOSEJOB (NU, XU, YU, NL, XL, YL, MAXPTS, X, Y,
     >                    COEFS, LUNCRT, LUNKBD, LUNTAB, MODE)
C
C  PURPOSE:
C
C        NOSEJOB performs various modifications to airfoil leading edges.
C     Initially, it has the following options (others may arise):
C
C        > Round off a sharp leading edge.
C        > Sharpen a rounded leading edge.
C
C        After the leading edge is modified, the original chord and thickness
C     could in principle be retrieved here through in-line reuse of PROFILE's
C     REFINE module - maybe some day.
C
C  METHOD:
C
C        A menu is presented for the various options.  The input coordinates
C     are overwritten by the modified coordinates.  The number of points on
C     each surface is held the same.  In the modified nose region, the relative
C     point distribution in terms of arc length is also held the same.  The
C     rounding and sharpening options are thus as reversible as possible,
C     especially if the option to retrieve original chord and thickness is not
C     used.  The user is asked to specify the points I1, I2 for each surface at
C     which the modified surface should blend with the original.
C
C                                      x
C                      I2     x
C                        *
C                     x+
C                   x +
C                 x  +
C                   x +
C                     x+
C                        *
C                      I1     x
C                                      x
C
C        Rounding is achieved by rotating the X axis to be parallel to the
C     bisector of the angle between the tangents at I1 and I2, then fitting a
C     conventional spline to the region near the nose where the "abscissas"
C     (rotated Y coordinates) are monotonic.  Evaluating this spline at points
C     spaced along the arc in the same relative way as the original points is
C     awkward!  Conservatively careful approximations are achieved by evaluating
C     the rounded surface at ~30 points spaced uniformly in the abscissa,
C     refitting those points parametrically, and working with cumulative chord
C     lengths from there.  Local storage is used for this set of gyrations.
C
C        Sharpening the nose is achieved by determining where the extrapolated
C     surfaces intersect, offering that point to the user as the default for
C     the new leading edge, then allowing an alternative point to be entered
C     interactively.  Much the same tedious steps as above are then taken to
C     preserve the original point distribution along the modified portions of
C     each surface, although local spline LCSFIT is used instead of PSFIT,
C     because it should be quite adequate for sharp noses.
C
C  ARGUMENTS:
C   ARG      DIM    TYPE I/O/S DESCRIPTION
C  NU,NL      -       I    I   Number of upper/lower surface pts.
C  XU,XL    NU,NL     R   I/O  Abscissas, upper/lower.
C  YU,YL    NU,NL     R   I/O  Corresponding ordinates, upper/lower.
C  X,Y     NU+NL-1    R    S   Storage for wraparound form of airfoil.
C  MAXPTS     -       I    I   Max. no. of pts. provided for on 1 surface.
C  COEFS   MAXPTS*3   R    S   Coefficients used by conventional spline.
C  LUNCRT     -       I    I   Logical unit for prompts.
C  LUNKBD     -       I    I   Logical unit for responses.
C  LUNTAB     -       I    I   Logical unit for printable record of the
C                              iterations performed by the REFINE option,
C                              if it is ever installed.
C  MODE       -       I    O   The menu choice, in case it is needed by
C                              the calling program.
C
C  PROCEDURES:
C    CHORD         Chord-length utility
C    COPY, RVERSE  Data transfer utilities
C    CSEVAL        Evaluates the spline previously fit by CSFIT
C    CSFIT         Fits a conventional cubic interpolating spline
C    FDCNTR        1st and 2nd derivative by central differencing
C    INTSEC2       Used to find where extrapolated surfaces meet
C    LCSFIT        Local spline, used for the "sharpen" option
C    PSFIT         Parametric form of CSFIT
C    PSTVAL        Evaluates PSFIT's spline at specified values of T
C    ROTATE2D      Rotates point(s) (X, Y) about (P, Q)
C    READER        Prompting utility
C    TSUBJ         Gives arc-length T associated with Jth data point for PSFIT
C    XGRID         Used to generate a uniform distribution
C
C  HISTORY:
C  10/21/91  DAS  Initial implementation, starting with a copy of REDISTRIB.
C                 No use of REFINE yet.  Difficulties with the leading edge
C                 (no longer necessarily the foremost point) would require
C                 RECTIFY as well, and this conflicts with reversing the
C                 operation.  Pursue it again some day, perhaps.
C  03/29/92  DAS  Use with I1=I2=20 revealed that TORIG(NLOCAL=30) needed to
C                 be TORIG(2*NLOCAL) for the "round" option.
C  08/01/93  DAS  INTSEC2 had CALCT1, CALCT2 arguments added.
C  04/03/95  DAS  INTSEC2 had TOL argument added.
C
C  AUTHORS: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, LUNTAB, MAXPTS, MODE, NL, NU

      REAL
     >   COEFS (MAXPTS, 3), X (MAXPTS * 2), XL (NL), XU (NU),
     >   Y (MAXPTS * 2), YL (NL), YU (NU)

C     Local constants:

      INTEGER
     >   MXMENU, NEND, NLOCAL
      REAL
     >   TOL
      CHARACTER
     >   BLANK * 1
      LOGICAL
     >   CALCT1, CALCT2
      PARAMETER
     >  (BLANK   = ' ',
     >   CALCT1  = .TRUE., ! Switches for INTSEC2
     >   CALCT2  = .TRUE.,
     >   TOL     = 1.E-6,  ! Used by INTSEC2 relative to the curve lengths
     >   MXMENU  = 2,
     >   NEND    = 4,      ! Number of original points included at each
                           ! end of the NLOCAL points used to refit the
                           ! nose region parametrically (round option)
     >   NLOCAL  = 40)     ! Number of evaluations in local storage of
                           ! the X vs. Y spline used for rounding, in
                           ! order to translate to arc-lengths carefully
C     Local variables:

      INTEGER
     >   I, I1, I1R, I2, IER, J, J1, J2, K, K1, K2, M, NMONO, NNOSE,
     >   NPTS, NTEMP

      REAL
     >   DUMMY, S1, S2, SCALE, THETA1, THETA2, THETAM, TOFFSET,
     >   XLE, YLE,
     >   TLOCAL (NLOCAL), TORIG (2*NLOCAL), XLOCAL (2*NLOCAL),
     >   YLOCAL (2*NLOCAL)

      LOGICAL
     >   DEFAULT, QUIT

      CHARACTER
     >   MENU (0 : MXMENU) * 37

C     Procedures:

      REAL
     >   CHORD, TSUBJ

      EXTERNAL
     >   CHORD, COPY, CSEVAL, CSFIT, FDCNTR, INTSEC2, LCSFIT, PSFIT,
     >   PSTVAL, READI, READR, ROTATE2D, RVERSE, TSUBJ, XGRID

C     Storage:

      DATA
     >   MENU /
     >   ' PROFILE''s "nose-job" options are:',
     >   '   (1) Round off a sharp leading edge',
     >   '   (2) Sharpen a rounded leading edge'/


C     Execution:

      WRITE (LUNCRT, '(/, A, /, A, A, /)') MENU
  100 M = 0
      CALL READI (LUNCRT, 'What''ll it be? ', LUNKBD, M, DEFAULT, QUIT)
      IF (QUIT) GO TO 999
      IF (M < 1 .OR. M > MXMENU) GO TO 100
      MODE = M

C     Round and sharpen options both need points specified at which to blend.

      I2 = 0
  110 WRITE (LUNCRT, '(A)')
      CALL READI (LUNCRT,
     >   'Index of first UPPER surface point to retain: ',
     >   LUNKBD, I2, DEFAULT, QUIT)
      IF (QUIT) GO TO 999
      IF (I2 <= 1 .OR. I2 >= NU) GO TO 110

      I1 = 0
  120 CALL READI (LUNCRT,
     >   'Index of first LOWER surface point to retain: ',
     >   LUNKBD, I1, DEFAULT, QUIT)
      IF (QUIT) GO TO 999
      IF (I1 <= 1 .OR. I1 >= NL) GO TO 120

      IF (I1 > NLOCAL .OR. I2 > NLOCAL) GO TO 800


C                              ---------------
C                              Rounding option
C                              ---------------

      IF (M == 1) THEN

C        Determine the original relative locations along the arc:

         TORIG (1) = 0.
         J = I1
         DO I = 2, I1
            TORIG (I) = CHORD (XL, YL, J, J - 1) + TORIG (I - 1)
            J = J - 1
         END DO

         I = I1
         DO J = 2, I2
            I = I + 1
            TORIG (I) = CHORD (XU, YU, J - 1, J) + TORIG (I - 1)
         END DO

         NNOSE = I      ! Number of nose points from I1 to I2 inclusive

C        Determine the gradients at I1 and I2 (finite differencing):

         CALL FDCNTR (I1, XL, YL, S1, S1)
         CALL FDCNTR (I2, XU, YU, S2, S2)

         THETA1 = ATAN (S1)
         THETA2 = ATAN (S2)
         THETAM = (THETA1 + THETA2) * 0.5  ! Mean angle, for rotating X axis to
         THETAM = THETAM * 45.0 / ATAN (1.0)

C        Set up the airfoil as a single wraparound curve, eliminating the
C        indicated leading edge region:

         CALL RVERSE (NL, XL, X)
         CALL RVERSE (NL, YL, Y)

         I1R   = NL - I1 + 1
         NTEMP = NU - I2 + 1
         CALL COPY (NTEMP, XU (I2), X (I1R + 1))
         CALL COPY (NTEMP, YU (I2), Y (I1R + 1))

         NPTS = NU + NL - I1 - I2 + 2   ! No. of points left after removing
                                        ! the points to be changed

C        Rotate the airfoil so the X axis is parallel to the angle bisector.
C        The leading edge will do for the center of rotation (arbitrary).

         CALL ROTATE2D (NPTS, X, Y, -THETAM, XU (1), YU (1))

C        Locate the nose portion where the transformed Ys are monotonic:

         DO I = I1R, 2, -1              ! Lower surface
            IF (Y (I - 1) >= Y (I)) THEN
               J1 = I
               GO TO 220
            END IF
         END DO
         J1 = 1

  220    DO I = I1R + 1, NPTS - 1      ! Upper
            IF (Y (I + 1) <= Y (I)) THEN
               J2 = I
               GO TO 240
            END IF
         END DO
         J2 = NPTS

  240    NMONO = J2 - J1 + 1

C        Spline transformed X vs. transformed Y for the nose region:

         CALL CSFIT (NMONO, Y (J1), X (J1), 0, DUMMY, 0, DUMMY,
     >               COEFS (1, 1), COEFS (1, 2), COEFS (1, 3), IER)
         IF (IER /= 0) GO TO 810

C        In order to work with arc lengths and with decent resolution,
C        evaluate this spline at sufficient points for reliable parametric
C        reinterpolation.  We need to include some of the unchanged points,
C        yet we need to include points I1 and I2 exactly in order to find
C        the new arc length (for redistributing it as for the original points).
C        We can't be certain of the number of points beyond J1, J2 that are
C        monotonic in transformed Y.  This is messy!

         K1 = MAX (I1R - NEND, J1)     ! +/-4 should control ends adequately.
         K2 = MIN (I1R + 1 + NEND, J2) ! There may not always be this many.

         J = 0
         DO I = K1, I1R - 1            ! Low end
            J = J + 1
            YLOCAL (J) = Y (I)
         END DO

         K = NLOCAL + 1
         DO I = K2, I1R + 2, -1        ! High end
            K = K - 1
            YLOCAL (K) = Y (I)
         END DO

C        Now for the middle: uniform in the Y range of transformed points
C        I1R and I1R + 1.

         NTEMP = NLOCAL - J - (NLOCAL + 1 - K)
         J = J + 1    ! These are now the YLOCAL (*) indices corresponding
         K = K - 1    ! to data points I1R, I1R + 1

         CALL XGRID (NTEMP, 0, Y (I1R), Y (I1R + 1), YLOCAL (J))

C        Evaluate the conventional spline at these Ys.

         CALL CSEVAL (NMONO, Y (J1), X (J1), NLOCAL, YLOCAL,
     >                COEFS (1, 1), COEFS (1, 2), COEFS (1, 3), XLOCAL)

C        Now refit the rounded curve parametrically:

         CALL PSFIT (NLOCAL, YLOCAL, XLOCAL, 'C', .FALSE., IER)
         IF (IER /= 0) GO TO 820

C        ... and evaluate it at the same relative arc lengths as the originals:

         TOFFSET = TSUBJ (J)  ! TSUBJ extracts an arc length from PSFIT's COMMON

         SCALE = (TSUBJ (K) - TOFFSET) / TORIG (NNOSE)
         DO I = 2, NNOSE
            TORIG (I) = TORIG (I) * SCALE + TOFFSET
         END DO

         CALL PSTVAL (NNOSE, TORIG, Y (I1R), X (I1R), Y (I1R), X (I1R),
     >                Y (I1R), X (I1R), NLOCAL, YLOCAL, XLOCAL)

C        Rotate the new points back to the original coordinate system:

         CALL ROTATE2D (NNOSE, X (I1R), Y (I1R), THETAM, XU (1), YU (1))

C        Finally, overwrite the original points with the rounded ones:

         J = I1R
         DO I = I1 - 1, 1, -1
            J = J + 1
            XL (I) = X (J)
            YL (I) = Y (J)
         END DO

         DO I = 1, I2 - 1
            XU (I) = X (J)
            YU (I) = Y (J)
            J = J + 1
         END DO

         IF (XL (2) <= XL (1) .OR. XU (2) <= XU (1)) THEN
            WRITE (LUNCRT, 1005)
         END IF

C                             -----------------
C                             Sharpening option
C                             -----------------

      ELSE IF (M == 2) THEN

C        Determine where the extrapolated surfaces meet.  INTSEC2's use
C        of local splines means the first four points per curve suffice.

         J1 = 1  ! Estimate of index near point of intersection
         J2 = 1
         CALL INTSEC2 (4, XL (I1), YL (I1), TORIG (1), J1, CALCT1,
     >                 4, XU (I2), YU (I2), TORIG (5), J2, CALCT2, TOL,
     >                 XLE, YLE, -LUNCRT, IER)
         IF (IER /= 0) GO TO 830

         WRITE (LUNCRT, 1004) XLE, YLE
         CALL READR (LUNCRT, 'Different X?  <CR> = above value: ',
     >               LUNKBD, XLE, DEFAULT, QUIT)
         IF (QUIT) GO TO 999
         CALL READR (LUNCRT, 'Different Y?  <CR> = above value: ',
     >               LUNKBD, YLE, DEFAULT, QUIT)
         IF (QUIT) GO TO 999


C        Lower surface first:

C        Determine arc-length distribution of points to be replaced.

         TORIG (1) = 0.
         DO I = 2, I1
            TORIG (I) = CHORD (XL, YL, I - 1, I) + TORIG (I - 1)
         END DO

C        The new leading edge and the first 3 retained points suffice
C        for interpolation in the extrapolated region by local methods.

         XLOCAL (1) = XLE
         YLOCAL (1) = YLE
         TLOCAL (1) = 0.
         I = I1
         DO J = 2, 4
            XLOCAL (J) = XL (I)
            YLOCAL (J) = YL (I)
            TLOCAL (J) = CHORD (XLOCAL, YLOCAL, J - 1, J) + TLOCAL (J-1)
            I = I + 1
         END DO

C        For a sharp leading edge, one chord is near enough to the arc length:

         SCALE = TLOCAL (2) / TORIG (I1)
         DO I = 2, I1
            TORIG (I) = TORIG (I) * SCALE
         END DO

C        Interpolate at the equivalent arc lengths:

         CALL LCSFIT (4, TLOCAL, XLOCAL, .TRUE., 'B', I1 - 1, TORIG,
     >                XL, XL)
         CALL LCSFIT (4, TLOCAL, YLOCAL, .TRUE., 'B', I1 - 1, TORIG,
     >                YL, YL)


C        Repeat for upper surface:

C        Determine arc-length distribution of points to be replaced.

         DO I = 2, I2
            TORIG (I) = CHORD (XU, YU, I - 1, I) + TORIG (I - 1)
         END DO

         I = I2
         DO J = 2, 4
            XLOCAL (J) = XU (I)
            YLOCAL (J) = YU (I)
            TLOCAL (J) = CHORD (XLOCAL, YLOCAL, J - 1, J) + TLOCAL (J-1)
            I = I + 1
         END DO

         SCALE = TLOCAL (2) / TORIG (I2)
         DO I = 2, I2
            TORIG (I) = TORIG (I) * SCALE
         END DO

         CALL LCSFIT (4, TLOCAL, XLOCAL, .TRUE., 'B', I2 - 1, TORIG,
     >                XU, XU)
         CALL LCSFIT (4, TLOCAL, YLOCAL, .TRUE., 'B', I2 - 1, TORIG,
     >                YU, YU)

      END IF


      RETURN

C     Error handling:

  800 WRITE (LUNCRT, 1002) I1, I2, NLOCAL
      GO TO 110
  810 WRITE (LUNCRT, 1003) 'CSFIT', IER
      GO TO 999
  820 WRITE (LUNCRT, 1003) 'PSFIT', IER
      GO TO 999
  830 WRITE (LUNCRT, 1003) 'INTSEC2', IER
C*****GO TO 999

  999 WRITE (LUNCRT, 1001) ' Stopping in NOSEJOB.'
      STOP

C     Formats:

 1001 FORMAT (/, A)
 1002 FORMAT (/, ' Index ', I2, ' and/or ', I2,
     >        ' exceeds local limit of ', I2, '.  Try again.')
 1003 FORMAT (/, ' NOSEJOB:  Bad return from ', A, '.  IER: ', I3)
 1004 FORMAT (/, ' Estimate of leading edge by extrapolation:', /,
     >        ' X = ', G13.6, '   Y = ', G13.6, /)
 1005 FORMAT (/, ' WARNING:  The common leading edge point is no',
     >        ' longer foremost.', /,
     >        ' This may be necessary if you intend to reverse the',
     >        ' operation.', /,
     >        ' But some results such as derivatives are affected. ',
     >        ' Proceeding ...')

      END SUBROUTINE NOSEJOB
C+----------------------------------------------------------------------
C
      SUBROUTINE NRMLIZ (NPTS, XYIN, XYOUT, XYLE, CHORD)
C
C  PURPOSE:
C     NRMLIZ normalizes X or Y coordinates (presumably from airfoils)
C     according to the given chord and leading edge coordinate.  It will
C     DE-normalize if the given chord is negative.  In-place is safe.
C     The relevant formulas are indicated by
C
C        X(norm)   = (X - X(l.e.)) / Chord                 and
C        X(denorm) = X(l.e.) + X(norm) * Chord
C
C  ARGUMENTS:
C   ARG    TYPE  I/O/S   DIM     DESCRIPTION
C   NPTS    I      I      -      Number of coordinates (X or Y)
C   XYIN    R      I    NPTS     Data to be normalized/de-normalized
C   XYOUT   R      O    NPTS     Results. May be same array as XYIN.
C   XYLE    R      I      -      X or Y coordinate of leading edge
C   CHORD   R      I      -      CHORD>0 means normalize by this chord;
C                                CHORD<0 means de-normalize by -CHORD.
C
C  AUTHORS: Various (Informatics, 1983).
C
C   10/21/88  D. Saunders  Ensured 1.0 exactly, where appropriate.
C
C-----------------------------------------------------------------------

C     Arguments:

      REAL XYIN(NPTS), XYOUT(NPTS)

C     Constants:

      REAL, PARAMETER ::
     >   ONE = 1., ZERO = 0.

C     Execution:

C  *  Use local copy of leading edge coordinate, to protect use of
C     an XYIN(*) element passed as this argument:

      XYLEAD = XYLE

      IF (CHORD > ZERO) THEN  ! Normalize coordinates:

         RCHORD = ONE / CHORD
         DO I = 1, NPTS
            XYOUT(I) = (XYIN(I) - XYLEAD) * RCHORD
         END DO

         IF (ABS (XYOUT(NPTS) - ONE) < 1.E-6) XYOUT(NPTS) = ONE

      ELSE  ! De-normalize coordinates:

         DO I = 1, NPTS
            XYOUT(I) = XYLEAD - XYIN(I) * CHORD
         END DO

      END IF

      END SUBROUTINE NRMLIZ
C+----------------------------------------------------------------------
C
      SUBROUTINE NRMSET (LUNCRT, LUNKBD, XMINSV, XMAXSV,
     >                   CNORML, XLE, YLE)
C
C  PURPOSE: NRMSET determines values to normalize/denormalize the
C           profile(s) by - introduced to keep the prompting out of
C           the main program.
C
C  ARGUMENTS:
C   ARG    TYPE  I/O/S   DIM  DESCRIPTION
C   LUNCRT  I    I        -   Logical unit for prompts (screen).
C   LUNKBD  I    I        -   Logical unit for responses (keyboard).
C   XMINSV, R    I        -   Minimum and maximum abscissas found over
C   XMAXSV  R    I        -   all profiles.
C   CNORML  R      O      -   Chord; either input at the terminal or
C                             determined by the data range.
C   XLE     R      O      -   Abscissa of leading edge point; either
C                             input from the terminal or determined from
C                             the data range checking.
C   YLE     R    I/O      -   Ordinate of leading edge point; input as
C                             Y corresp. to XMINSV; may be updated here.
C
C  PROCEDURES:
C  READER   Prompts for and reads integer, real, etc.
C
C  AUTHOR: Leslie Collins, Informatics General, Palo Alto, CA
C
C  HISTORY:
C  08/07/84   LJC   Initial design and coding.
C  10/09/85   DAS   Eliminated normalization of plot scale info.
C  02/07/90   DAS   Eliminated list-directed I/O.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   LUNCRT, LUNKBD
      REAL
     >   CNORML, XLE, XMAXSV, XMINSV, YLE

C  *  Local variables:

      LOGICAL
     >   DEFAULT, QUIT
      REAL
     >   CHORD

C  *  Procedures:

      EXTERNAL
     >   READR

C  *  Execution:

      CHORD = XMAXSV - XMINSV

      WRITE (LUNCRT, 1001) 'Current chord: ', CHORD,
     >   'Enter <CR> to normalize by this chord, or enter a'
      CALL READR (LUNCRT,
     >   'different chord; (a negative value de-normalizes): ',
     >   LUNKBD, CNORML, DEFAULT, QUIT)
      IF (QUIT) GO TO 999

      IF (DEFAULT) THEN

         CNORML = CHORD
         XLE = XMINSV

      ELSE ! YLE is already defined by MAXMIN called from main program.

         WRITE (LUNCRT, 1001)
         CALL READR (LUNCRT, 'Enter abscissa of leading edge point: ',
     >      LUNKBD, XLE, DEFAULT, QUIT)
         IF (QUIT) GO TO 999

         CALL READR (LUNCRT, 'Enter ordinate of leading edge point: ',
     >      LUNKBD, YLE, DEFAULT, QUIT)
         IF (QUIT) GO TO 999

      END IF

      RETURN

  999 WRITE (LUNCRT, 1001) 'Stopping as requested.'
      STOP

 1001 FORMAT (1X, A, G13.6)

      END SUBROUTINE NRMSET
C+----------------------------------------------------------------------
C
      SUBROUTINE OPTIMIZE (NU, XU, YU, NL, XL, YL, TITLE,
     +                     LUNCRT, LUNKBD, LUNPRT, LUNBUMPS, LUNTARG,
     +                     UPPER, NXTARG, TARGX, TARGCRV)
C
C  PURPOSE:  OPTIMIZE  perturbs  one  surface of an airfoil by applying
C            given shape functions to the surface and  optimizing  some
C            of the bump parameters so as to match a given target curv-
C            ature distribution in the least squares sense.
C
C  METHOD:   No attempt is made to optimize both surfaces  in  the same
C            run - separate runs are likely to be more manageable, from
C            both the user's standpoint and the programmer's.
C
C            The chosen bump function set is read from an editable disk
C            file rather than entered interactively.   Repeated running
C            would be too busy otherwise.   Values for ALL of the shape
C            function parameters must be in this file, some as starting
C            guesses for the optimizing algorithm, the others inactive.
C
C            Provision is made for using the first N  Wagner  functions
C            more easily than by reading a text file, since this common
C            case is easily generated.
C
C            The basic procedure follows:
C
C            *  Prompt for which surface to optimize.
C
C            *  Prompt for and read the previously-prepared bump set.
C
C            *  Prompt for and read the target curvature  distribution,
C               assumed to be in simple QPLOT format.
C
C            *  Prompt for any additional constraints (thickness?...).
C
C            *  Set up optimizer  QNMDIF  (initial function evaluation,
C               tolerances, and estimation of optimal finite difference
C               intervals, etc.).
C
C            *  Minimize the sum-of-squares-type objective function.
C
C            *  Update the airfoil surface using the optimal bumps, in-
C               place as expected by the calling program which performs
C               any plotting/tabulating/saving of results requested.
C
C  ARGUMENTS:
C  ARG    DIM  TYPE I/O/S DESCRIPTION
C  NU      -    I     I   Number of upper surface points
C  XU     NU    R     I   Upper surface abscissas
C  YU     NU    R    I/O  Lower surface ordinates
C  NL      -    I     I   Number of lower surface points
C  XL     NL    R     I   Lower surface abscissas
C  YL     NL    R    I/O  Lower surface ordinates
C  TITLE   *    C     I   Title for printout generated by this run (LUNPRT)
C  LUNCRT       I     I   Logical unit for screen
C  LUNKBD       I     I   Logical unit for keyboard
C  LUNPRT       I     I   Logical unit for printout from OPTIMIZE & QNMDIF
C  LUNBUMPS     I     I   Logical unit for previously-prepared bump data
C  LUNTARG      I     I   Logical unit for target curvature distribution
C  UPPER   -    L     O   TRUE if target curvature applies to upper surface
C  NXTARG  -    I     O   Number of target curvature points (for plotting)
C  TARGX NXTARG R     O   Target curvature distribution (needed in COMMON,
C  TARGCRV " "  R     O   but returned to calling program for plotting - an
C                         approach that does not clutter the calling program
C                         with the localized COMMON block).
C
C PARAMETER CONSTANTS:
C    PARAM   TYPE   DESCRIPTION
C  LUNERR      I    Logical unit number for error messages.
C  MXBUMPS     I    Maximum no. of bump functions provided for.
C  MXOPT       I    Maximum no. of optimization variables allowed.
C  MXPARM      I    Maximum no. of parameters associated with any of
C                   the bump functions (including a multiplier).
C  MXSURF      I    Maximum no. of data points handled per surface.
C  MXTARG      I    Maximum no. of target curvature values handled.
C
C COMMON BLOCKS USED:
C
C   /ACTUAL/  (Current values corresponding to input opt. variables)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C   YCALC     MXSURF     R    S   Perturbed airfoil surface
C  Y1CALC     MXSURF     R    S   1st and 2nd derivatives - by-prod-
C  Y2CALC     MXSURF     R    S   ucts of the curvature calculations
C  YKCALC     MXSURF     R    S   Curvature distribution correspond-
C                                 ing to the current bump variables
C
C   /BCHARS/  (Bump function names - can't go in /BUMPS/ but should)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C  BNAMES     MXBUMPS  C*11   S   Names of bumps found by GETBUMPS
C
C   /BUMPS /  (List of bump function parameters, etc.)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C  NBUMPS       -        I    S   Number of bump functions active
C  PARAMS MXPARM,MXBUMPS R    S   Complete set of parameters for the
C                                 active bump functions,   including
C                                 some possible unused ones  present
C                                 because a 2-dim. array approach is
C                                 used for storing them  (as opposed
C                                 to packing the variable numbers of
C                                 parameters associated with differ-
C                                 ent shape functions).  See ACTIVE.
C  ACTIVE MXPARM,MXBUMPS L    S   Marks the bump function parameters
C                                 as active or inactive.  The unused
C                                 ones must be marked inactive along
C                                 with the fixed (but used) ones.
C  VSCALE MXPARM,MXBUMPS R    S   Scale factors needed by the optim-
C                                 izing algorithm.
C
C   /ORIGNL/  (Copies of the original airfoil surface data)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C    NXY        -        I    S   Number of points defining  surface
C   XORIG      NXY       R    S   Abscissas for the surface
C   YORIG      NXY       R    S   Ordinates for the unperturbed srf.
C
C   /TARGET/  (Target distribution quantities)
C    VAR       DIM     TYPE I/O/S DESCRIPTION
C   NTARG       -        I    S   Number of target data points
C   XTARG     NTARG      R    S   Abscissas for the target data
C  YKTARG     NTARG      R    S   Target curvature values
C  WEIGHT     NTARG      R    S   Weights to be applied (multiplica-
C                                 tively) to each element of the sum
C                                 of squares.  May be all 1s.
C
C  SIGNIFICANT LOCAL VARIABLES:
C  VAR       DIM     TYPE   DESCRIPTION
C  D        MXOPT      R    Diagonal factor of Hessian approximation
C  G        MXOPT      R    Objective function gradient approximation
C  H        MXOPT      R    Finite difference intervals
C  L MXOPT*(MXOPT-1)/2 R    Lower triangle factor of Hessian approx.
C  OPTVARS  MXOPT      R    Active bump function variables
C  PNAMES   MXPARM*  C*6    Names of bump function parameters found
C           MXBUMPS
C
C  FILES USED:
C  LUN      I/O/S   DESCRIPTION
C  LUNBUMPS   I     Previously-prepared bump function data
C  LUNCRT     O     User prompts
C  LUNKBD     I     User responses
C  LUNPRT     O     Optimizer printout
C  LUNTARG    I     Target curvature distribution (QPLOT format)
C
C  EXTERNAL REFERENCES:
C  ACTIVATE   Packs or unpacks active variables; scales/unscales too.
C  CALCTHICK  Calculates airfoil thickness.
C  CENDIF     Estimates optimal finite difference intervals.
C  CFDISTRIB  Function needed by QNMDIF - compares distributions and
C             returns corresponding sum of squares.
C  COPY       Copies one array to another.
C  DUMMY      "USER" routine called by QNMDIF - do nothing in this case.
C  GETBUMPS   Preprocesses bump function data read from disk file.
C  LSTFIL     Echoes target curvature file to printable file.
C  OPENER     File opening utility.
C  PRTBUMPS   Prints bump function data.
C  PRWRIT     Used to echo target curvature distribution found.
C  QNMDIF     General purpose minimizer not requiring derivatives.
C  READER     Prompting utility.
C  RDQPL      Reads target curvature distribution in QPLOT format.
C
C  HISTORY:
C  01/18/84   DAS   Initial design and coding.
C  01/25/84   DAS   Added scaling, printing of bump data, etc.
C  01/27/84   DAS   Target data now returned as arguments for plotting
C  01/30/84   LJC   Provided for new CALCTHICK (handles NU/=NL)
C  03/16/84   DAS   Provided for constraining thickness (penalty fun.)
C  10/24/85   DAS   Used LSTFIL in place of PRWRIT for echoing target data
C  08/19/86   DAS   Calling sequence to CENDIF was revised
C  01/29/90   DAS   Removed END DOs, underscores, list-directed I/O, and
C                   concatenations in I/O lists (for IRIS 4D purposes);
C                   installed OPENER.
C
C  AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

C ... Arguments:

      INTEGER   LUNBUMPS, LUNCRT, LUNKBD, LUNPRT, LUNTARG, NL, NU,
     >          NXTARG
      REAL      TARGX(*), TARGCRV(*), XL(NL), XU(NU), YL(NL), YU(NU)
      LOGICAL   UPPER
      CHARACTER TITLE*(*)

C ... Parameter constants and global variables:

C ... Quantities passed to objective function routine via COMMON
C     because of fixed calling sequence expected by QNMDIF.  Any
C     changes affecting these COMMONs should also be made in one
C     other routine - CFDISTRIB.

      INTEGER, PARAMETER ::
     >   MXBUMPS = 20,  MXPARM = 3, MXOPT = MXPARM*MXBUMPS,
     >   MXSURF  = 180, MXTARG = MXSURF, LUNERR = 6

      REAL            YCALC, Y1CALC, Y2CALC, YKCALC
      COMMON /ACTUAL/ YCALC(MXSURF), Y1CALC(MXSURF), Y2CALC(MXSURF),
     >                YKCALC(MXSURF)

      CHARACTER*11    BNAMES
      COMMON /BCHARS/ BNAMES(MXBUMPS)

      INTEGER         NBUMPS
      REAL            PARAMS, VSCALE
      LOGICAL         ACTIVE
      COMMON /BUMPS / PARAMS(MXPARM,MXBUMPS), VSCALE(MXOPT),
     >                ACTIVE(MXPARM,MXBUMPS), NBUMPS

      INTEGER         NOTHER, NXY
      REAL            XORIG, XOTHER, YORIG, YOTHER
      COMMON /ORIGNL/ XORIG(MXSURF), YORIG(MXSURF),
     >                XOTHER(MXSURF), YOTHER(MXSURF), NOTHER, NXY

      INTEGER         NTARG
      REAL            OPTWRK, PENLTY, TARGTH, XTARG, YKTARG, WEIGHT
      COMMON /TARGET/ XTARG(MXTARG), YKTARG(MXTARG), WEIGHT(MXTARG),
     >                TARGTH, PENLTY, OPTWRK(4*MXSURF), NTARG

C ... Local variables (more below):

      INTEGER   I, IER, NOPTVARS
      REAL      D(MXOPT), G(MXOPT), H(MXOPT), L( MXOPT*(MXOPT-1)/2 ),
     >          OPTVARS(MXOPT), THICKNESS, XATMAX
      LOGICAL   DEFAULT, QUIT
      CHARACTER FILENAME*50, PNAMES(MXPARM, MXBUMPS)*6,
     >          SURFACE*1, TEXT*80, YESNO*1

C ... The following locals are required by QNMDIF:

      INTEGER   NFCEN, NFTOTL, NITER, NLDIM, NTYPE
      REAL      EPSMCH, EPSOBJ, ETA, SSQ, SSQMIN, STEPMX, TOL
      LOGICAL   UNITL, LOCAL, CONV, CONTIN, RESCUE, PRINT

C ... Procedures:

      EXTERNAL  ACTIVATE, CALCTHICK, CENDIF, CFDISTRIB, COPY,
     >          DUMMY, GETBUMPS, LSTFIL, OPENER, PRTBUMPS, QNMDIF,
     >          RDQPL, READC, READR, READY

C ... Execution:


      IF (MAX (NU, NL) > MXSURF) THEN
         WRITE (LUNCRT, '(/, A)')
     >      ' Not enough work-space. Recompile OPTIMIZE, CFDISTRIB.'
      END IF

      WEIGHT = 1.E+0

C ... Prompt for the surface being perturbed, and copy it to where
C     the objective function routine can get at it:

  200 CALL READC (LUNCRT, 'Enter U(pper) or L(ower) to indicate'//
     +            ' the surface being optimized: ',
     +            LUNKBD, SURFACE, DEFAULT, QUIT)

      IF (QUIT) THEN
         WRITE (LUNCRT, '(/, A)') ' Stopping at user request.'
         STOP
      END IF

      IF (SURFACE == 'L') THEN
         UPPER = .FALSE.
         NXY = NL
         CALL COPY (NL, XL, XORIG)
         CALL COPY (NL, YL, YORIG)
         NOTHER = NU
         CALL COPY (NU, XU, XOTHER)
         CALL COPY (NU, YU, YOTHER)
      ELSE IF (SURFACE == 'U') THEN
         UPPER = .TRUE.
         NXY = NU
         CALL COPY (NU, XU, XORIG)
         CALL COPY (NU, YU, YORIG)
         NOTHER = NL
         CALL COPY (NL, XL, XOTHER)
         CALL COPY (NL, YL, YOTHER)
      ELSE
         WRITE (LUNCRT, 1001) 'Invalid response. Try again.'
         GO TO 200
      END IF

C ... Get the original thickness ratio, and a target thickness if any:

      CALL CALCTHICK (NL, NU, XL, XU, YL, YU, THICKNESS,
     >                XATMAX, OPTWRK, OPTWRK(MXSURF+1),
     >                OPTWRK(2*MXSURF+1), OPTWRK(3*MXSURF+1))

      WRITE (LUNCRT, 1001) 'Original thickness:', THICKNESS

      TARGTH = THICKNESS
      CALL READR (LUNCRT, 'Enter desired % thickness ' //
     >            '(<CR> means keep same or don''t care): ',
     >            LUNKBD, TARGTH, DEFAULT, QUIT)

C ... Is thickness going to be constrained?

      PENLTY = 0.E+0
      CALL READR (LUNCRT, 'Enter thickness penalty parameter ' //
     >            '(<CR> means no constraint): ',
     >            LUNKBD, PENLTY, DEFAULT, QUIT)

      WRITE (LUNPRT, 1002) 'Case: ', TITLE
      WRITE (LUNPRT, 1001) 'Original thickness:', THICKNESS,
     >                     'Corresponding  x/c:', XATMAX,
     >                     'Target thickness  :', TARGTH,
     >                     'Penalty parameter :', PENLTY

C ... Get the complete set of bump functions, previously prepared:

      CALL GETBUMPS (LUNCRT, LUNKBD, LUNBUMPS, MXPARM, MXBUMPS,
     >               NBUMPS, BNAMES, PARAMS, PNAMES,
     >               ACTIVE, NOPTVARS, VSCALE)

      CALL PRTBUMPS (LUNPRT, MXPARM, NBUMPS, BNAMES, PARAMS,
     >               PNAMES, ACTIVE, NOPTVARS, VSCALE)

C ... Get the target data (assumed to be in simple QPLOT format).
C     It is needed in COMMON for CFDISTRIB, but also by the calling
C     program for plotting.  Arguments were added for this purpose,
C     so as to keep the main program free of COMMON blocks.

      CALL OPENER (LUNCRT,
     >             'Enter target curvature distribution file name: ',
     >             LUNKBD, FILENAME, LUNTARG, 'OLD')

      CALL RDQPL (MXSURF, TEXT, NTARG, XTARG, YKTARG, LUNTARG)

      WRITE (LUNPRT, 1003)
     >   ' ', 'Target curvature distribution found:', ' '

      REWIND LUNTARG
      CALL LSTFIL (LUNTARG, LUNPRT, TEXT)

      NXTARG = NTARG
      CALL COPY (NXTARG,  XTARG, TARGX)
      CALL COPY (NXTARG, YKTARG, TARGCRV)

C ... Prepare for the optimizing algorithm.
C     First, extract the active variables as contiguous starting guesses:

      CALL ACTIVATE (.TRUE., MXPARM * NBUMPS, ACTIVE, PARAMS, OPTVARS,
     >               VSCALE)

      ETA = 1.0E-1
      TOL = 1.0E-3
      SSQMIN = 0.0E+0
      STEPMX = 1.0E+10
      EPSMCH = 5.0E-8
      EPSOBJ = 1.0E+2 * EPSMCH
      PRINT  = .TRUE.
      LOCAL  = .FALSE.
      CONTIN = .TRUE.
      RESCUE = .FALSE.
      UNITL  = .FALSE.
      NITER  = 100
      NTYPE  = 1
      NLDIM  = MAX (NOPTVARS * (NOPTVARS - 1) / 2, 1)

      DO I = 1, NLDIM
         L(I) = 0.0E+0
      END DO
      DO I = 1, NOPTVARS
         H(I) = -1.0E-3
      END DO

      PRINT = .FALSE.
      CALL READY (LUNCRT, 'Do you want full optimization printout? '//
     >            '(Y/N; <CR>=No): ',
     >            LUNKBD, PRINT, DEFAULT, QUIT)

C ... Compute initial value of objective function:

      CALL CFDISTRIB (NOPTVARS, OPTVARS, SSQ)

      WRITE (LUNCRT, 1001) 'Initial function value:', SSQ
      WRITE (LUNPRT, 1001)
      WRITE (LUNPRT, 1001) 'Initial function value:', SSQ, ' '

C ... Estimate good finite difference step-sizes:

      CALL CENDIF (NOPTVARS, OPTVARS, SSQ, EPSOBJ, H, G, D, NFCEN,
     >             LUNPRT, CFDISTRIB)
      NFTOTL = 1 + NFCEN

C ... Minimize the objective function:

      CALL QNMDIF (NOPTVARS, NLDIM, NFTOTL, NITER, NTYPE, LUNPRT,
     >             OPTVARS, SSQ, SSQMIN, G, H, L, D, ETA, TOL, STEPMX,
     >             EPSMCH, EPSOBJ, UNITL, LOCAL, CONV, CONTIN, RESCUE,
     >             PRINT, CFDISTRIB, DUMMY)

C ... Regenerate best result found:

      CALL CFDISTRIB (NOPTVARS, OPTVARS, SSQ)

      WRITE (LUNCRT, 1001) '  Final function value:', SSQ
      WRITE (LUNPRT, 1001)
      WRITE (LUNPRT, 1001) 'Repeat of best function evaluation:', SSQ

      CALL PRTBUMPS (LUNPRT, MXPARM, NBUMPS, BNAMES, PARAMS,
     >               PNAMES, ACTIVE, NOPTVARS, VSCALE)

C ... Update the airfoil permanently:

      IF (UPPER) THEN
         CALL COPY (NU, YCALC, YU)
      ELSE
         CALL COPY (NL, YCALC, YL)
      END IF

C ... Compute the modified thickness achieved:

      CALL CALCTHICK (NL, NU, XL, XU, YL, YU, THICKNESS,
     >                XATMAX, OPTWRK, OPTWRK(MXSURF+1),
     >                OPTWRK(2*MXSURF+1), OPTWRK(3*MXSURF+1))

      WRITE (LUNPRT, 1001)
      WRITE (LUNPRT, 1001) 'Modified % thickness:  ', THICKNESS,
     >                     'Corresponding abscissa:', XATMAX
      WRITE (LUNCRT, 1001) 'Modified % thickness:  ', THICKNESS,
     >                     ' Corresponding abscissa:', XATMAX

C ... Formats:

 1001 FORMAT (' ', A, G13.6)
 1002 FORMAT (/, '1', A, A, //)
 1003 FORMAT (' ', A, A)

      END SUBROUTINE OPTIMIZE
C+------------------------------------------------------------------------------
C
      SUBROUTINE PRREAD (LUNRD, TITLE, MAXPTS, NU, NL, X, Y, XU, XL,
     >                   YU, YL, FORMAT, IER)
C
C  ACRONYM: PRofile: READ one
C
C  PURPOSE: PRREAD reads one airfoil profile per call from a file that is
C           assumed to be open (LUNRD).  Four file formats are supported.
C
C           PRREAD returns the data as XU, YU, XL, YL, matching the standard
C           "PROFILE" format and returns a flag indicating which format the
C           values read were found in.  Standard PROFILE format is shown below.
C           The lower surface values are optional, but a zero must be read for
C           NL if no lower surface is included unless this is the last airfoil
C           in the file (meaning EOF can be used to indicate NL = 0).
C           In this case, PRREAD returns a symmetrical airfoil.
C
C               TITLE                   <Variable-length title>
C               NU   Upper surface      <Integer, first token>
C               X         Y             <Two reals>
C               X         Y                 .
C               .         .                 .     (May be X/C, Y/C;
C               .         .                 .      Xs are increasing.)
C               .         .                 .
C               NL   Lower surface      <Integer, first token> <may be 0 or EOF>
C               X         Y             <Two reals>
C               X         Y                 .
C               X         Y  ! Trailing comments are permitted
C               .         .                 .
C               .         .                 .     (Xs are increasing.)
C               .         .                 .
C               ! X         Y           <Point suppressed; NL must be adjusted>
C               .         .                 .
C
C    NOTE:  If FORMAT = 1 and both surfaces are present, PROFILE expects
C           them to have the same leading edge point.  The trailing edge
C           points may differ.  However, checking for the same leading
C           edge point is left to the calling program in case PRREAD is
C           being used to read other upper/lower surface-type distributions.
C
C           The next two formats are wraparound clockwise and wraparound
C           counterclockwise, where the coordinates begin at the trailing
C           edge, wrap around the leading edge, and end at the trailing edge.
C           The clockwise case begins with the lower surface, and the counter-
C           clockwise case begins with the upper surface.  The format shown
C           below is essentially the same for both cases. NPTS is the total
C           number of points on the airfoil.
C
C               TITLE                   <CHARACTER*80>
C               NPTS                    <Integer, first token>
C               X         Y             <Two reals>
C               X         Y                 .
C               .         .                 .     (May be X/C, Y/C;
C               .         .                 .      Xs are decreasing
C               .         .                 .      until the leading
C               .         .                 .      edge, then increasing)
C
C    NOTE:  Wraparound formats do NOT have duplicate leading edge points.
C
C           The fourth format is called 3-column format.  The airfoil is
C           represented in three columns, with the same abscissas for both
C           surfaces in the 1st column and ordinates for the upper and lower
C           surfaces in the 2nd and 3rd columns respectively.  Abscissas are
C           increasing as with standard format.  Here NPTS is the number of
C           points on either surface.  Further columns may be present in
C           this case only (as for Cp distributions - only the first three
C           columns of the data proper will be read here).
C
C               TITLE                           <CHARACTER*80>
C               NPTS                            <Integer, first token>
C               X         YU        YL     [Cp] <Reals, first 3 tokens;
C               X         YU        YL     [Cp]  further columns are ignored.>
C               .         .         .      [..]     .
C               .         .         .               .   (May be X/C, Y/C;
C               .         .         .               .    Xs are increasing.)
C               .         .         .               .
C               .         .         .               .
C
C  ARGUMENTS:
C   ARG     DIM   TYPE  I/O/S   DESCRIPTION
C   LUNRD    -      I     I     Logical unit number for file being read.
C   MAXPTS   -      I     I     Maximum number of coordinates on any one
C                               surface expected by calling program.
C   TITLE    *      C     O     Variable-length title for this profile.
C   NU,NL    -      I     O     Number of data points found for upper and
C                               lower surfaces.  NL is set to NU if NL = 0
C                               (or if EOF is encountered) on the read.
C   X,Y   MAXPTS*2  R    S/O    Buffers for reading coordinates; returned
C                               with the airfoil in clockwise wraparound form
C                               (NU + NL - 1 points) in case that is more
C                               convenient than the form in XU, YU, XL, YL.
C   XU     MAXPTS   R     O     Upper surface abscissas found.
C   XL     MAXPTS   R     O     Lower surface abscissas.
C   YU,YL  NU,NL    R     O     Corresponding ordinates.
C   FORMAT   -      I     O     FORMAT=1 means data found in standard format
C                                     =2 means clockwise wraparound format
C                                     =3 means counterclockwise wraparound
C                                     =4 means 3-column format
C   IER      -      I     O     IER=0 means one profile found normally;
C                                  =1 means EOF encountered on the first
C                                     read -- normal unless this was the
C                                     first call to PRREAD.
C                                  =2 means an unexpected EOF or other
C                                     read error encountered -- fatal.
C                                  =3 means NU or NL  was found out of
C                                     range (> MAXPTS, or too small).
C                                  =4 means a coordinate was missing.
C
C  PROCEDURES:
C   COPY           Copies a vector
C   GETLINE        Reads a line as text; handles suppressed points & comments
C   RVERSE         Reverses a vector
C   TOKENS         Tokenizes a string
C
C  HISTORY:
C  03/19/82    DAS    Original implementation.
C  06/29/84    LJC    Added reading of wraparound formats.
C  09/14/84    LJC    Changed "legend" entry to be read from the dataset
C                     and added a prompt for the plot title at a higher
C                     level. (Formerly the title was read from the
C                     dataset and the legend (TITLE here) was hard-coded.)
C  10/18/84    DAS    Took out check for duplicate leading edge point -
C                     not wanted if other distributions are being read.
C                     (Left to the calling program where appropriate.)
C  02/27/85    LJC    Added 3-column format.
C  09/18/87    DAS    Handled as few as 2 points on a surface properly;
C                     allowed for EOF on reading NL -- treat as NL = 0.
C  02/06/90  DAS/RAK  GETLINE introduced to permit trailing comments or
C                     suppression of coordinates.
C  10/23/91    DAS    Should use LINE (1:LAST) everywhere, not LINE.
C  09/18/92    DAS    LENTOK = 20 failed on VAX double-precision numbers.
C  11/06/92    DAS    3-column format can now ignore additional columns.
C                     Too few coordinates on a line is now trapped.
C  11/16/93    DAS    Returned additional clockwise wraparound form of
C                     the airfoil in the available X, Y arrays in case
C                     that is more convenient for some applications.
C
C  AUTHOR:     David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   FORMAT, IER, LUNRD, MAXPTS, NU, NL

      CHARACTER
     >   TITLE * 80

      REAL
     >   X (MAXPTS * 2), XU (MAXPTS), XL (MAXPTS), Y (MAXPTS * 2),
     >   YU (MAXPTS), YL (MAXPTS)

C  *  Local constants:

      INTEGER, PARAMETER ::
     >   LENTOK  = 24

      CHARACTER, PARAMETER ::
     >   BLANK * 1   = ' ',
     >   COMMENT * 1 = '!',
     >   IFMT * 8    = '(BN,I24)',
     >   RFMT * 10   = '(BN,F24.0)'

C  *  Local variables:

      INTEGER
     >   I, IOS, LAST, MIDSRF, LE, NPTS, NUMBER

      CHARACTER
     >   LINE * (3 * LENTOK), LIST (3) * (LENTOK)

C  *  Procedures:

      EXTERNAL
     >   COPY, GETLINE, RVERSE, TOKENS

C  *  Execution:

      IER = 0

C  *  Look for descriptive text for a profile.  There may be no more profiles.

      CALL GETLINE (LUNRD, BLANK, TITLE, LAST, IOS)
      IF (IOS  < 0) GO TO 800                    ! Normal EOF
      IF (IOS /= 0) GO TO 900

C  *  Look for a count of the points to follow:

  100 CALL GETLINE (LUNRD, COMMENT, LINE, LAST, IOS)
      IF (IOS  /= 0) GO TO 900
      IF (LAST == 0) GO TO 100    ! Insignificant line

      NUMBER = 1
      CALL TOKENS (LINE (1 : LAST), NUMBER, LIST)

      READ (LIST (1), IFMT, ERR=900) NPTS

      IF (NPTS <= MAXPTS*2 .AND. NPTS > 1) THEN

C  *     Check number of columns in first line of data:

  110    CALL GETLINE (LUNRD, COMMENT, LINE, LAST, IOS)
         IF (IOS  /= 0) GO TO 900
         IF (LAST == 0) GO TO 110
         NUMBER = 3
         CALL TOKENS (LINE (1 : LAST), NUMBER, LIST)

         IF (NUMBER == 2) THEN

C  *        Only two columns of data found. Process one of the first 3 formats.
C           No need to backspace any more.  (But a DO loop is inconvenient.)

            I = 1
  200       CONTINUE
               READ (LIST (1), RFMT, ERR=900) X (I)
               READ (LIST (2), RFMT, ERR=900) Y (I)
               I = I + 1

               IF (I <= NPTS) THEN

  210             CALL GETLINE (LUNRD, COMMENT, LINE, LAST, IOS)

                  IF (IOS  /= 0) GO TO 900
                  IF (LAST == 0) GO TO 210

                  CALL TOKENS (LINE (1 : LAST), NUMBER, LIST)

                  IF (NUMBER < 2) GO TO 910
                  GO TO 200
               END IF

C  *        Distinguish standard and wrap-around formats automatically:

            MIDSRF = MAX (1, NPTS / 4)
            IF (X (MIDSRF) < X (MIDSRF + 1)) THEN

C  *           Abscissas are increasing; "PROFILE" format assumed.
C              NPTS must be less than MAXPTS (not MAXPTS*2) now:

               IF (NPTS <= MAXPTS) THEN
                  FORMAT = 1
                  NU = NPTS

C  *              Copy X and Y into upper surface arrays:

                  CALL COPY (NU, X, XU)
                  CALL COPY (NU, Y, YU)

C  *              Look for a lower surface point count (may be 0 or EOF):

  300             CALL GETLINE (LUNRD, COMMENT, LINE, LAST, IOS)

                  IF (IOS < 0) THEN          ! Normal EOF
                     NL = 0
                  ELSE IF (IOS  /= 0) THEN
                     GO TO 900
                  ELSE IF (LAST == 0) THEN
                     GO TO 300
                  ELSE
                     NUMBER = 1
                     CALL TOKENS (LINE (1 : LAST), NUMBER, LIST)
                     READ (LIST (1), IFMT, ERR=900) NL
                  END IF

                  IF (NL > 0) THEN

                     IF (NL <= MAXPTS) THEN
                        NUMBER = 2
                        DO I = 1, NL
  310                      CALL GETLINE (LUNRD, COMMENT, LINE, LAST,IOS)
                           IF (IOS  /= 0) GO TO 900
                           IF (LAST == 0) GO TO 310
                           CALL TOKENS (LINE (1 : LAST), NUMBER, LIST)
                           IF (NUMBER < 2) GO TO 910
                           READ (LIST (1), RFMT, ERR=900) XL (I)
                           READ (LIST (2), RFMT, ERR=900) YL (I)
                        END DO
                     ELSE
C  *                    Number of lower surface points found out of range:
                        IER = 3
                     END IF

                  ELSE

C  *                 Generate symmetrical lower surface points:

                     DO I = 1, NU
                        XL (I) = XU (I)
                        YL (I) =-YU (I)
                     END DO
                     NL = NU

                  END IF

               ELSE

C  *              Number of upper surface points found out of range:

                  IER = 3

               END IF

            ELSE

C  *           Process wraparound airfoil:

               DO I = 1, NPTS-1  ! Search for leading edge:

                  IF (X (I) < X (I+1)) THEN
                     LE = I
                     GO TO 510
                  END IF
               END DO

C  *           No leading edge found:

               GO TO 900

  510          CONTINUE

C  *           Arrange coordinates in standard PROFILE format:

               IF (Y (MIDSRF) < Y (NPTS - MIDSRF)) THEN

C  *              The lower surface is first; this is the clockwise case:

                  FORMAT = 2
                  NU = NPTS - LE + 1
                  NL = LE
                  CALL RVERSE (NL, X, XL)
                  CALL RVERSE (NL, Y, YL)
                  CALL COPY (NU, X (NL), XU)
                  CALL COPY (NU, Y (NL), YU)

               ELSE  ! Counterclockwise case:

                  FORMAT = 3
                  NL = NPTS - LE + 1
                  NU = LE
                  CALL COPY (NL, X (NU), XL)
                  CALL COPY (NL, Y (NU), YL)
                  CALL RVERSE (NU, X, XU)
                  CALL RVERSE (NU, Y, YU)

               END IF

            END IF

         ELSE IF (NUMBER == 3) THEN  ! Process 3-column format:

            FORMAT = 4
            I = 1
  600       CONTINUE
               READ (LIST (1), RFMT, ERR=900) XU (I)
               READ (LIST (2), RFMT, ERR=900) YU (I)
               READ (LIST (3), RFMT, ERR=900) YL (I)
               I = I + 1
               IF (I <= NPTS) THEN
  610             CALL GETLINE (LUNRD, COMMENT, LINE, LAST, IOS)
                  IF (IOS  /= 0) GO TO 900
                  IF (LAST == 0) GO TO 610
                  CALL TOKENS (LINE (1 : LAST), NUMBER, LIST)
                  IF (NUMBER < 3) GO TO 910
                  GO TO 600
               END IF

            CALL COPY (NPTS, XU, XL)
            NU = NPTS
            NL = NPTS

         ELSE  ! Inappropriate number of columns found:

            GO TO 900

         END IF

      ELSE  ! Number of points found out of range:

         IER = 3

      END IF

C  *  Return a second copy of the coordinates in clockwise wraparound form
C     in case that is more convenient for some application other than PROFILE:

      IF (IER == 0) THEN
         CALL RVERSE (NL, XL, X)
         CALL COPY   (NU, XU, X (NL))
         CALL RVERSE (NL, YL, Y)
         CALL COPY   (NU, YU, Y (NL))
      END IF

      GO TO 999


  800 CONTINUE

C  *  Normal EOF encountered on first read -- no more profiles:

      IER = 1
      GO TO 999

  900 CONTINUE

C  *  Abnormal EOF or other read error -- fatal:

      IER = 2
      GO TO 999

  910 CONTINUE

C  *  Too few values on a line:

      IER = 4
!     GO TO 999

  999 RETURN

      END SUBROUTINE PRREAD
C+------------------------------------------------------------------------
C
      SUBROUTINE PRTAB (LUNTAB, TITLE, SUBTITLE, WRAPCRV,
     >                  N, X, Y, YP, YPP, YK)
C
C  PURPOSE: PRTAB tabulates geometrical properties of an airfoil profile,
C           one surface per call.
C
C  ARGUMENTS:
C   ARG     DIM   TYPE  I/O/S   DESCRIPTION
C   LUNTAB   -      I     I     Logical unit for file being written to
C   TITLE    *      C     I     Main title for tabulation (probably the
C                               title from the file of airfoil coordinates)
C   SUBTITLE *      C     I     Additional subtitle that can help make the
C                               tabulation more self-descriptive
C   WRAPCRV  -      L     I     .TRUE. means the derivatives and curvatures
C                               were calculated parametrically using splines;
C                               .FALSE. means they were finite difference
C                               approximations with surfaces treated separately.
C   N        -      I     I     Number of data points on the surface
C   X        N      R     I     Abscissas for the surface
C   Y        N      R     I     Corresponding ordinates
C   YP       N      R     I     1st derivatives  at given coordinates
C   YPP      N      R     I     2nd derivatives  at given coordinates
C   YK       N      R     I     Curvature values at given coordinates
C
C  HISTORY:
C  01/20/83   LJC   Initial coding
C  04/05/83   LJC   Added calculation of curvatures.
C  04/08/83   LJC   Handle each surface separately.
C  10/17/83   DAS   Removed curvature calculn. - now done by FD12K.
C  10/29/83   DAS   Introduced subtitle.
C  10/28/91   DAS   Introduced WRAPCRV.
C
C  AUTHOR:   Leslie Collins, Informatics, Palo Alto, CA.
C
C----------------------------------------------------------------------

C     Arguments:

      INTEGER       LUNTAB, N
      REAL          X (N), Y (N), YP (N), YPP (N), YK (N)
      LOGICAL       WRAPCRV
      CHARACTER*(*) TITLE, SUBTITLE

C     Local variables:

      INTEGER       I

C     Execution:

      WRITE (LUNTAB, 1001) '1', TITLE, SUBTITLE
      WRITE (LUNTAB, 1002) N
      IF (WRAPCRV) THEN
         WRITE (LUNTAB, 1003) 'parametric spline'
      ELSE
         WRITE (LUNTAB, 1003) 'nonparametric finite differencing'
      END IF
      WRITE (LUNTAB, 1004)
      WRITE (LUNTAB, 1005) (X(I), Y(I), YP(I), YPP(I), YK(I), I = 1, N)

C     Formats:

 1001 FORMAT (/, A1, A, //, 1X, A)
 1002 FORMAT (/, ' Number of points: ', I3)
 1003 FORMAT (/, ' Derivatives and curvature calculated by ',
     >        A, '.')
 1004 FORMAT (//, T9, 'X',  T25, 'Y', T39, 'Y''', T54, 'Y"',
     >        T66, 'CURVATURE' )
 1005 FORMAT (2E15.7, 3E15.6)

      END SUBROUTINE PRTAB
C+----------------------------------------------------------------------
C
      SUBROUTINE PRTBUMPS (LUNPRT, MXPARM, NBUMPS, BNAMES,
     >                     PARAMS, PNAMES, ACTIVE, NACTIVE, SCALES)
C
C PURPOSE: PRTBUMPS prints a description of a given bump set.  It is
C          intended for echoing an initial set prior to optimization
C          of the active variables, and for tabulating the optimized
C          results in unscaled form.
C
C ARGUMENTS:
C    ARG     DIM   TYPE I/O/S DESCRIPTION
C   LUNPRT    -      I    I   Logical unit for tabulations.
C   MXPARM    -      I    I   Max. # parameters  for any one bump.
C   NBUMPS    -      I    I   Number of bump functions involved.
C   BNAMES  NBUMPS  C*(*) I   Names of the bump functions.
C   PARAMS  MXPARM,  R    I   PARAMS(1:?,J) are the parameters de-
C           NBUMPS            fining the Jth bump.
C   PNAMES  MXPARM, C*(*) I   Names of function parameters
C           NBUMPS
C   ACTIVE  MXPARM,  L    I   ACTIVE(I,J)=.TRUE. if the Ith param-
C           NBUMPS            eter of the Jth bump is to be treat-
C                             ed as active (variable).  Otherwise,
C                             the parameter is either to remain at
C                             the value returned here (fixed),  or
C                             bump J has fewer than I parameters.
C   NACTIVE  -       I    I   Number of active parameters.
C   SCALES  MXPARM,  R    I   Scale factors, included here for in-
C           NBUMPS            formation only.
C
C HISTORY:
C   01/25/84   DAS    Initial design and code.
C   02/23/84   DAS    Now prints bump names instead of code numbers.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER   LUNPRT, MXPARM, NACTIVE, NBUMPS
      REAL      PARAMS(MXPARM,NBUMPS), SCALES(MXPARM,NBUMPS)
      LOGICAL   ACTIVE(MXPARM,NBUMPS)
      CHARACTER BNAMES(NBUMPS)*(*), PNAMES(MXPARM,NBUMPS)*(*)

C     Local variables:

      INTEGER   I, J
      REAL      UNDEF

      DATA      UNDEF / 999.E+0 / ! This preset value should match GETBUMPS.

C     Execution:

      WRITE (LUNPRT, '(//, (A))') ' Shape Function Description',
     >                            ' --------------------------'

      DO J = 1, NBUMPS
         WRITE (LUNPRT, '(//, A, I3, 2A, //, A, /)')
     >      ' Function', J, ':  Name: ', BNAMES(J),
     >      ' Parameter    Value       Scale      Active?'

         DO I = 1, MXPARM
            IF (PARAMS(I,J) == UNDEF) CYCLE
               WRITE (LUNPRT, '( 1X, A, F12.8, F12.5, L10 )')
     >            PNAMES(I,J), PARAMS(I,J), SCALES(I,J), ACTIVE(I,J)
         END DO
      END DO

      WRITE (LUNPRT, '(//, A, I4)')
     >   ' Total number of active parameters:', NACTIVE

      END SUBROUTINE PRTBUMPS
C+---------------------------------------------------------------------
C
      SUBROUTINE PRWRIT (LUN, MAXPTS, TITLE, NU, NL, XU, XL, YU, YL,
     >                   FORMAT, PRECISION)
C
C  ACRONYM: Program PRofile: WRITe one airfoil profile
C                   --       ----
C  PURPOSE: PRWRIT writes one airfoil profile to the indicated file, in
C           one of four formats - "standard",  wraparound (either way),
C           or "three-column" - described in PRREAD.
C
C           PRWRIT  may also be used to save camber/thickness distribu-
C           tions,  or to save second derivative information for reuse.
C           Both of these extended uses employ  "standard"  format (not
C           wrap-around), with suitable labeling.
C
C           This version allows for three types of precision,  prompted
C           by the need to retain more digits for manipulating the all-
C           important leading edge region effectively, and to deal with
C           large magnitudes such as when the coordinates are in milli-
C           meters, or when the "ordinates" are really 2nd derivatives.
C
C  METHOD:  Values are passed to PRWRIT as separate surfaces,  and  are
C           returned untouched.  Note that the wrap-around formats omit
C           one of the two leading edge points, which are assumed to be
C           the same point.
C
C  ARGUMENTS:
C   ARG     DIM   TYPE  I/O/S   DESCRIPTION
C   LUN      -      I     I     Logical unit for file being written to.
C   MAXPTS   -      I     I     Max. no. of points on any one surface.
C   TITLE    -     C*(*)  I     Variable-length title for profile.
C   NU,NL    -      I     I     Number of data points for upper and lower
C                               surfaces (NU > 0; NL=0 is OK if FORMAT=1;
C                               NU=NL=no. of pts. in 3-column format and
C                               the camber/thickness distributions - FORMATs
C                               4 & 5).
C   XU     MAXPTS   R     I     Upper surface abscissas.
C   XL     MAXPTS   R     I     Lower surface abscissas (if any).
C   YU,YL  MAXPTS   R     I     Corresponding ordinates, y" values, or
C                               camber/thickness values, accdg. to FORMAT.
C   FORMAT   -      I     I     Requested format for output where:
C                               = 1 means standard PROFILE format (coords.)
C                               = 2 means clockwise wrap-around format
C                               = 3 means counterclockwise wrap-around
C                               = 4 means 3-column format
C                               = 5 means 2nd derivatives (standard format)
C                               = 6 means camber/thickness (standard format)
C   PRECISION -     I     I     Controls number of digits in saved values:
C                               = 1 means "full" (see PURPOSE above);
C                               = 2 means "engineering" or "flow code" (F10.6);
C                               = 3 means "low" - appropriate for y" values.
C  FILES USED:
C   UNIT     I/O/S    DESCRIPTION
C   LUN        O      File (assumed open) to contain one or more datasets
C
C  HISTORY:
C  12/23/82   LJC   Coding adapted from PRREAD.
C  10/30/83   DAS   Introduced alternative E-format so that PRWRIT
C                   can be used to save 2nd derivatives similarly.
C  07/09/84   LJC   Added writing of wraparound formats.
C  02/11/84   DAS   Handled thickness/camber - hard to avoid dupli-
C                   cation of code now.  But PRWRIT is still handy
C                   for getting this stuff out of the main program.
C  02/28/85   LJC   Added 3-column format. (FORMAT=4)
C  04/09/87   DAS   Values except X are written with F10.7 format
C                   instead of F10.6.   (May help difficulties in
C                   critical leading edge region; will not affect
C                   formatted 2-column reads, but may be a little
C                   misleading.)
C  04/23/87   DAS   Erred further in direction of more precision:
C                   use E format except on basically normalized
C                   data; go to F12.8 for normalized airfoils.
C  04/24/87   DAS   Retained F10.6 option after all for old flow
C                   codes (FORMAT=7; standard PROFILE format only.
C  04/27/87   DAS   Reluctantly introduced PRECISION argument after
C                   the above failed to handle FORMAT="flowcode" and
C                   COUNTERCLOCKWISE both.
C  11/03/87   DAS   Switched to G formats for "full" and "low" precision.
C  12/12/91   DAS   Had to "comment out" the text following NU, NL
C                   because of how RDXYZ works now.  (PROFILE uses
C                   PRWRIT to save 2nd derivatives, which are read back
C                   via RDXYZ.)
C
C  AUTHOR:   Leslie Collins, Informatics, Palo Alto, CA.
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      CHARACTER
     >   TITLE * (*)
      INTEGER
     >   FORMAT, LUN, MAXPTS, NL, NU, PRECISION
      REAL
     >   XU (MAXPTS), XL (MAXPTS), YU (MAXPTS), YL (MAXPTS)

C  *  Local variables:

      INTEGER
     >   I
      CHARACTER
     >   FMT * 13

C  *  Execution:

      IF (PRECISION == 1) THEN        ! "Full" precision.
         FMT = '(2(1X,G15.7))'          ! Gw.d needs w >= d + 8 (basically).

      ELSE IF (PRECISION == 2) THEN   ! "Engineering" precision.
         FMT = '(2F10.6)     '          ! Traditional for many flow codes.

      ELSE                              ! "Low" precision.
         FMT = '(2(1X,G13.5))'          ! Appropriate for y" values (FORMAT=5).
      END IF

      IF (FORMAT == 4) FMT (2:2) = '3'

      WRITE (LUN, '(A)') TITLE

      IF (FORMAT == 1) THEN           ! Standard PROFILE format.

         WRITE (LUN, 1001) NU, 'Upper Coordinates'
         WRITE (LUN, FMT) (XU (I), YU (I), I = 1, NU)

         IF (NL > 0) THEN
            WRITE (LUN, 1001) NL, 'Lower Coordinates'
            WRITE (LUN, FMT) (XL (I), YL (I), I = 1, NL)
         END IF

      ELSE IF (FORMAT == 2) THEN      ! Wrap-around (lower surface first).

         WRITE (LUN, 1001) NU + NL - 1, 'Coordinates Clockwise'
         WRITE (LUN, FMT) (XL (I), YL (I), I = NL, 1, -1)
         WRITE (LUN, FMT) (XU (I), YU (I), I = 2, NU)

      ELSE IF (FORMAT == 3) THEN      ! Wrap-around (upper surface first).

         WRITE (LUN, 1001) NU + NL - 1, 'Coordinates Counter-clockwise'
         WRITE (LUN, FMT) (XU (I), YU (I), I = NU, 1, -1)
         WRITE (LUN, FMT) (XL (I), YL (I), I = 2, NL)

      ELSE IF (FORMAT == 4) THEN      ! Three-column format.

         WRITE (LUN, 1001) NU, 'Coordinates per surface'
         WRITE (LUN, FMT) (XU (I), YU (I), YL (I), I = 1, NU)

      ELSE IF (FORMAT == 5) THEN      ! 2nd derivatives (both surfaces).
                                        ! Use standard format.
         WRITE (LUN, 1001) NU, 'Upper 2nd Derivatives'
         WRITE (LUN, FMT) (XU (I), YU (I), I = 1, NU)

         WRITE (LUN, 1001) NL, 'Lower 2nd Derivatives'
         WRITE (LUN, FMT) (XL (I), YL (I), I = 1, NL)

      ELSE  ! FORMAT = 6: Camber and Thickness distributions in Standard format.

         WRITE (LUN, 1001) NU, 'Camber'
         WRITE (LUN, FMT) (XU (I), YU (I), I = 1, NU)

         WRITE (LUN, 1001) NL, 'Thickness'
         WRITE (LUN, FMT) (XL (I), YL (I), I = 1, NL)

      END IF

C  *  Formats:

 1001 FORMAT (I4, ' ! ', A)

      END SUBROUTINE PRWRIT
C+----------------------------------------------------------------------
C
      SUBROUTINE RDKEYS (LUNCRT, LUNIN, MXCRVS, UNDEF, FORMAT,
     >                   PRECISION, PLTLINE, CPSLINE, CRVLINE, MAXCRV,
     >                   MINCRV, WRAPCRV, PLTORIG, PLTREV, THREED,
     >                   XAXIS, XMIN, XMAX, YMIN, YMAX, DATFIL, PLTFIL,
     >                   CRVFIL, TABFIL, YPPFIL, CPSFIL, SPREADFIL)
C
C  PURPOSE:
C
C     RDKEYS reads PROFILE's keyword input control file.   Most of these
C     inputs are common to all modes of operation.  Other mode-dependent
C     inputs are entered interactively in the appropriate subprogram.
C
C  METHOD:
C
C     Read a line of input and break up the line into pairs of keywords and
C     values.  For each pair, look up the keyword in a dictionary and then
C     look up the value in a separate dictionary associated with that keyword.
C     Assign a value to the corresponding internal control variable via the
C     argument list.
C
C     Some keywords may have more than one corresponding value.  In this case,
C     all of the tokens following the first token are treated as values and
C     not pairs.  A full description of the keywords follows.
C
C  KEYWORD GUIDELINES AND DEFINITIONS:
C
C     Keyword/value pairs may appear with more than one pair on a line.
C     However, the multivalued keywords PLTLINE, CPSLINE, CRVLINE, and
C     NOFILE must not appear with other keywords on the same line.
C
C     The default value in each case appears in square brackets.
C
C  KEYWORD  VALUES and synonyms     DESCRIPTION
C  -------  -------------------     -----------
C
C  FORMAT   [SAME]                  One of four formats for output profile.
C           PROFILE or STANDARD     May be in standard PROFILE format (ab-
C           CLOCKWISE or WRAPAROUND scissas increasing), clockwise wraparound
C           COUNTERCLOCKWISE        format, counterclockwise wraparound for-
C           THREE-COLUMN or         mat, or 3-column format.  SAME means the
C           THREE_COLUMN or         same format as the input profile.  NOTE:
C           THREECOLUMN  or         To allow easily for several synonyms for
C           TABLE                   for the THREE-COLUMN value, only the first
C                                   5 characters of the value are checked.
C
C  PRECISION   [FULL]               Controls number of digits in output
C           ENGINEERING             airfoil coordinates.  FULL gives F11.8
C                                   if possible, or E15.8 if any X >=10.
C                                   ENGINEERING gives the traditional F10.6
C                                   common to many flow solvers.
C
C  PLTLINE  [DEFAULT]               Controls line types of curves on profile
C           LINE                    plots.  One value may be included for
C           DASH                    each curve on the plot.  The default is
C           DOT                     symbols connected by a solid line, with
C           CHAINDASH               a different symbol type for each succes-
C           CHAINDOT                sive curve.  The first curve typically
C           THICK                   represents the original profile; the
C           SYMBOLS                 second curve represents the revised one.
C                                   Overriding the default might be desirable
C                                   when plotting multi-element airfoils or
C                                   when lines without symbols are required.
C                                   At most 20 curves are provided for.  Note:
C                                   All the line types in QPLOT are available.
C                                   SYMBOLS refers to symbols with no line
C                                   connecting them.
C
C  CPSLINE  [see PLTLINE above]     Controls line types on Cps plots in the
C                                   same manner as PLTLINE above.  One value
C                                   per curve may be included, chosen from
C                                   the same list of values as those shown
C                                   for PLTLINE.
C
C  CRVLINE  [see PLTLINE above]     Controls line types on curvature plots
C                                   in the same way as PLTLINE and CPSLINE.
C
C  CURVATURE or [NONPARAMETRIC] or  CURVATURE and DERIVATIVES are synonymous
C  DERIVATIVES  [FINITE_DIFFERENCE] controls for the type of calculations
C               SPLINE     or       used for derivatives and hence curvature.
C               PARAMETRIC or       The default is separate-surface treatment
C               WRAPAROUND          using finite differences, as needed for
C                                   consistency with PROFILE's REFINE and
C                                   OPTIMIZE options.  The two surfaces appear
C                                   as separate frames in the curvature plot.
C                                   Otherwise, the full wrap-around curvature
C                                   distribution is calculated using a para-
C                                   metric spline and plotted on a single frame.
C
C                                   The default normally suffices except if
C                                   the region of interest is very near a
C                                   rounded leading edge.  Note that not all
C                                   of the possibilites are provided for, such
C                                   as parametric finite differences.
C
C  MINCURVATURE   [-5.]             Cutoff values for plotted curvatures.
C  MAXCURVATURE   [+5.]             Practice shows that +/-5. give useful
C                                   plot scaling by ignoring the high curv-
C                                   ature values near the leading edge.  On
C                                   the other hand, it may well be desired
C                                   to focus on the leading edge region.  Set
C                                   both to 999. to obtain the full range.
C                                   See the CURVATURE/DERIVATIVES control.
C
C  NOFILE   [NONE]                  Used to suppress any combination of the
C           DAT                     seven output files generated by PROFILE.
C           PLT                     The values correspond to the extensions
C           TAB                     of the file names.  See elsewhere for a
C           CRV                     complete description of file contents.
C           YPP                     NONE serves only to assist leaving the
C           CPS                     NOFILE control word in the input file
C           SPREAD                  even if all outputs are desired.
C
C  PLOT     [BOTH]                  Controls plotting of original OR revised
C           ORIGINAL                profile.  The default is to plot both
C           REVISED                 original and revised (if one exists).
C
C  THREED   [FALSE] or [NO]         For plotting of multiple stations from
C           TRUE or YES             a 3-D wing. The default is the 2-D case.
C
C  XAXIS    [6.4]                   Length of x-axis in inches.  The default
C                                   is normally appropriate for an 8.5 x 11
C                                   page in portrait mode.
C
C           The following four keywords apply to windowing.  Any or none
C           of them may be used.
C
C  XMIN     [minima and             Minimum abscissa for desired window
C  XMAX     maxima of               Maximum abscissa for desired window
C  YMIN     the input               Minimum ordinate for desired window
C  YMAX     coordinates]            Maximum ordinate for desired window
C
C
C  SAMPLE CONTROL FILE:
C
C     A sample input file follows.  Note that keywords and values
C     may be separated with blanks, commas, colons, equals signs,
C     or tabs. Remember, keywords with more than one value should
C     appear on separate lines.  Any keyword or text value may be
C     truncated to unambiguous leading characters.   Blank  lines
C     and trailing ! comments are ignored.
C
C
C           FORMAT = STANDARD   PRECISION = FULL
C           PLOT BOTH  THREED:NO
C           PLTLINE = SOLID, SOLID
C           CPSLINE = DOT, SYMBOLS
C           CRVLINE = SOLID, DASH, CHAINDOT
C           XAXIS = 20.
C           XMIN = 0.  XMAX 0.1
C           MAXCURVATURE = 999.   ! Both 999. means plot the full
C           MINCURVATURE = 999.   ! curvature range
C           DERIVATIVES = PARAMETRIC
C           NOFILE: YPP SPREAD
C
C
C  ERROR HANDLING:
C
C     The dictionary lookup routine returns a pointer indicating which entry
C     matched the key.  When the pointer is negative, meaning no match was
C     found, RDKEYS prints an error message along with the string in question
C     and stops.
C
C  EXTERNAL REFERENCES:
C
C  GETLINE     Reads one line, and handles trailing comments.
C  LOOKUP      Performs dictionary lookups.
C  PAIRS       Breaks up alternating fields into pairs of keywords and values.
C  SCANNR      Looks for non-blank fields in a string.
C  TOKENS      Separates groups of contiguous characters into an array.
C
C  HISTORY:
C
C  08/01/84    LJC    Initial design and coding.
C  12/05/84    DAS    Added CPSFIL handling.
C  02/19/85    DAS    Took out SHARP/BLUNT control - now in REDISTRIB.
C  03/22/85    LJC    Replaced LINE keyword with PLTLINE, CPSLINE and CRVLINE,
C                     to control line types on all plots.
C  10/09/85    DAS    Introduced UNDEF argument because of revised QPLOT.
C  04/24/87    DAS    Added FORMAT=FLOWCODE after PRWRIT was changed to
C                     provide more precision than F10.6.
C  04/27/87    DAS    The above didn't handle FLOWCODE and COUNTERCLOCKWISE
C                     needed by one user - introduced PRECISION instead.
C                     Allowed FORMAT=WRAPAROUND [=CLOCKWISE - arbitrary,
C                     but it had to be one way or the other].
C  05/08/87    DAS    Bugs: VALUES () * 10 is too short for ENGINEERING;
C                     Dimension of LIST () has to be MAX (MXCRVS,NFILDC-1),
C                     except NFILDC is all that will fit on a line anyway
C                     for CPSLINE, CRVLINE, and PLTLINE.
C  04/29/88    DAS    Added spreadsheet-compatible output file control.
C  06/18/91    DAS    Made MXLIST = MXCRVS + 1 = 21 for 3D case of many
C                     sections.  (Plot line types CAN be abbreviated, so
C                     the comment of 05/08/87 is not really right.)
C  10/23/91    DAS    Introduced GETLINE to allow commented control files.
C  10/28/91    DAS    Introduced a means of specifying full wrap-around
C                     curvature distribution (DERIVATIVES/CURVATURE keywords).
C
C  AUTHOR:     Leslie Collins, Informatics General, Palo Alto, CA
C
C------------------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   FORMAT, LUNCRT, LUNIN, MXCRVS, PRECISION
      REAL
     >   MAXCRV, MINCRV, UNDEF, XAXIS, XMAX, XMIN, YMAX, YMIN
      CHARACTER
     >   CPSLINE (MXCRVS) * 9, CRVLINE (MXCRVS) * 9,
     >   PLTLINE (MXCRVS) * 9
      LOGICAL
     >   CPSFIL, CRVFIL, DATFIL, PLTFIL, PLTORIG, PLTREV, SPREADFIL,
     >   TABFIL, THREED, WRAPCRV, YPPFIL

C  *  Local constants:

      INTEGER, PARAMETER ::
     >   MXCHARS = 81, MXLIST = 21, NDERDC = 5, NFILDC = 8, NFMTDC = 8,
     >   NKEYDC  = 17, NLINDC = 9,  NPLTDC = 3, NPREDC = 2, NTHRDC = 4

      CHARACTER, PARAMETER ::
     >   BLANK * 1 = ' '

C  *  Local variables:

      INTEGER
     >   FIRST, I, IOS, J, LAST, MARK, NUMBER, NPAIRS, POINTR
      LOGICAL
     >   CR
      CHARACTER
     >   VAL*3,
     >   BUFFER * (MXCHARS), DERDIC (NDERDC) * 13, KEYDIC (NKEYDC) * 12,
     >   FILDIC (NFILDC) * 6, FMTDIC (NFMTDC) * 11, KEYS (NKEYDC) * 12,
     >   LINDIC (NLINDC) * 9, LIST (MXLIST) * 10, PLTDIC (NPLTDC) * 8,
     >   PREDIC (NPREDC) * 11, THRDIC (NTHRDC) * 8, VALUES (NKEYDC) * 12

C  *  Procedures:

      EXTERNAL
     >   GETLINE, LOOKUP, PAIRS, SCANNR, TOKENS

C  *  Dictionaries:

C     Guidelines:
C     >  All calls to LOOKUP assume that the dictionaries are alphabetical.
C     >  Use of synonyms either leads to long dictionary entries or forces
C        (safe) shortening of the target string to (say) 5 characters.  Any
C        characters beyond that in the other dictionary entries are strictly
C        for programmer readability.

C     DERIVATIVES/CURVATURE dictionary:

      DATA DERDIC
     >   /'FINITE=NONPAR', 'NONPARAMETRIC', 'PARAM=SPLINE',
     >    'SPLINE', 'WRAPA=SPLINE'/

C     Main control keyword dictionary:

      DATA KEYDIC
     >   /'CPSLINE', 'CRVLINE',      'CURVATURE',    'DERIV=CURVA',
     >    'FORMAT',  'MAXCURVATURE', 'MINCURVATURE', 'NOFILE',
     >    'PLOT',    'PLTLINE',      'PRECISION',    'THREED',
     >    'XAXIS',   'XMAX', 'XMIN', 'YMAX', 'YMIN'/

C     NOFILE dictionary:

      DATA FILDIC
     >   /'CPS', 'CRV', 'DAT', 'NONE', 'PLT', 'SPREAD', 'TAB', 'YPP'/

C     FORMAT dictionary:

      DATA FMTDIC
     >   /'CLOCK', 'COUNT', 'PROFI', 'SAME', 'STAND:PROFI',
     >    'TABLE:THREE', 'THREE', 'WRAPA:CLOCK'/

C     Line-type dictionary:

      DATA LINDIC
     >   /'CHAINDASH', 'CHAINDOT', 'DASH', 'DEFAULT', 'DOT',
     >    'LONGDASH', 'SOLID', 'SYMBOLS', 'THICK'/

C     PLOT dictionary:

      DATA PLTDIC
     >   /'BOTH', 'ORIGINAL', 'REVISED'/

C     PRECISION dictionary:

      DATA PREDIC
     >   /'ENGINEERING', 'FULL'/

C     THREED dictionary:

      DATA THRDIC
     >   /'FALSE', 'NO:FALSE', 'TRUE', 'YES:TRUE'/


C  *  Execution:

C  *  Set defaults:

      FORMAT = 0
      PRECISION = 1

      DO J = 1, MXCRVS
         CPSLINE (J) = BLANK
         CRVLINE (J) = BLANK
         PLTLINE (J) = BLANK
      END DO

      MAXCRV = +5.
      MINCRV = -5.
      WRAPCRV = .FALSE.
      PLTORIG = .TRUE.
      PLTREV = .TRUE.
      THREED = .FALSE.
      XAXIS = 6.4
      XMIN = UNDEF
      XMAX = UNDEF
      YMIN = UNDEF
      YMAX = UNDEF
      DATFIL = .TRUE.
      PLTFIL = .TRUE.
      CRVFIL = .TRUE.
      TABFIL = .TRUE.
      YPPFIL = .TRUE.
      CPSFIL = .TRUE.
      SPREADFIL = .TRUE.

C  *  No control file was found if LUNIN < 0.  Return with defaults:

      IF (LUNIN < 0) GO TO 950


    5 CONTINUE

C  *     Read one line of input and break into pairs of keywords and values:

         CALL GETLINE (LUNIN, '!', BUFFER, LAST, IOS)
         IF (IOS  <  0) GO TO 950   ! Normal EOF
         IF (IOS /=  0) GO TO 980   ! Read error
         IF (LAST == 0) GO TO 5

         NPAIRS = NKEYDC
         FIRST = 1

         CALL PAIRS (BUFFER (1 : LAST), NPAIRS, KEYS, VALUES)

         DO I = 1, NPAIRS

C  *        Look for a keyword in the dictionary.

            KEYS (I) (6:) = BLANK        ! Only first 5 are significant
            CALL LOOKUP (NKEYDC, KEYDIC, .TRUE., KEYS (I), POINTR)
            IF (POINTR <= 0) GO TO 960

C  *        Read the accompanying value(s) into an internal control
C           variable/array, or find a text value in the appropriate
C           secondary dictionary and assign accordingly:

            GO TO (10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120,
     >             130, 140, 150, 160, 170) POINTR

   10          CONTINUE  !  --CPSLINE--

C  *           Here's a special case. We need to find the beginning of the
C              first value and separate the rest of the tokens in the buffer.

               CALL SCANNR (BUFFER, FIRST, LAST, MARK)
               FIRST = MARK + 2
               CALL SCANNR (BUFFER, FIRST, LAST, MARK)
               NUMBER = MXLIST
               CALL TOKENS (BUFFER (FIRST:LAST), NUMBER, LIST)

               DO J = 1, NUMBER

C  *              Find each value in the dictionary and assign accordingly:

                  CALL LOOKUP (NLINDC, LINDIC, .TRUE., LIST (J), POINTR)
                  IF (POINTR <= 0) GO TO 970
                  CPSLINE (J) = LIST (J)
               END DO
               GO TO 5

   20          CONTINUE   !  --CRVLINE--

C  *           Another special case. (See CPSLINE above)

               CALL SCANNR (BUFFER, FIRST, LAST, MARK)
               FIRST = MARK + 2
               CALL SCANNR (BUFFER, FIRST, LAST, MARK)
               NUMBER = MXLIST
               CALL TOKENS (BUFFER (FIRST:LAST), NUMBER, LIST)

               DO J = 1, NUMBER

C  *              Find each value in the dictionary and assign accordingly:

                  CALL LOOKUP (NLINDC, LINDIC, .TRUE., LIST (J), POINTR)
                  IF (POINTR <= 0) GO TO 970
                  CRVLINE (J) = LIST (J)
               END DO
               GO TO 5

   30          CONTINUE   !  --CURVATURE--
   40          CONTINUE   !  --DERIVATIVES--

C  *           These two are synonymous controls for derivative estimates.

               VALUES (I) (6:) = BLANK
               CALL LOOKUP (NDERDC, DERDIC, .TRUE., VALUES (I), POINTR)
               IF (POINTR <= 0) GO TO 970

               WRAPCRV = VALUES (I) == 'SPLINE'
               CYCLE


   50          CONTINUE   !  --FORMAT--

C  *           Pass only the first 5 characters of the value in this case.
C              We want to recognize any of a number of synonyms:

               VALUES (I) (6:) = BLANK
               CALL LOOKUP (NFMTDC, FMTDIC, .TRUE., VALUES (I), POINTR)
               IF (POINTR <= 0) GO TO 970
               IF (VALUES (I) == 'PROFI') THEN
                  FORMAT = 1
               ELSE IF (VALUES (I) == 'CLOCK') THEN
                  FORMAT = 2
               ELSE IF (VALUES (I) == 'COUNT') THEN
                  FORMAT = 3
               ELSE IF (VALUES (I) == 'THREE') THEN
                  FORMAT = 4
               ELSE
                  FORMAT = 0
               END IF
               CYCLE

   60          CONTINUE   !  --MAXCURVATURE--
               READ (VALUES (I), 1000) MAXCRV
               CYCLE

   70          CONTINUE   !  --MINCURVATURE--
               READ (VALUES (I), 1000) MINCRV
               CYCLE

   80          CONTINUE   !  --NOFILE--

C  *           Another special case (see CPSLINE):

               CALL SCANNR (BUFFER, FIRST, LAST, MARK)
               FIRST = MARK + 2
               CALL SCANNR (BUFFER, FIRST, LAST, MARK)
               NUMBER = NFILDC
               CALL TOKENS (BUFFER (FIRST:LAST), NUMBER, LIST)

               DO J = 1, NUMBER
                  CALL LOOKUP (NFILDC, FILDIC, .TRUE., LIST (J), POINTR)
                  IF (POINTR <= 0) GO TO 970

                  VAL = LIST (J) (1:3)
                  IF (VAL == 'CPS') THEN
                     CPSFIL = .FALSE.
                  ELSE IF (VAL == 'CRV') THEN
                     CRVFIL = .FALSE.
                  ELSE IF (VAL == 'DAT') THEN
                     DATFIL = .FALSE.
                  ELSE IF (VAL == 'PLT') THEN
                     PLTFIL = .FALSE.
                  ELSE IF (VAL == 'TAB') THEN
                     TABFIL = .FALSE.
                  ELSE IF (VAL == 'YPP') THEN
                     YPPFIL = .FALSE.
                  ELSE IF (VAL == 'SPR') THEN
                     SPREADFIL = .FALSE.
                  ELSE IF (VAL == 'NON') THEN
C                    Do nothing.  'None' allows NOFILE keyword to be present
C                    even if all output files are desired.
                  END IF
               END DO
               GO TO 5

   90          CONTINUE   !  --PLOT--

               CALL LOOKUP (NPLTDC, PLTDIC, .TRUE., VALUES (I), POINTR)
               IF (POINTR <= 0) GO TO 970
               IF (VALUES (I) == 'ORIGINAL') THEN
                  PLTREV = .FALSE.
               ELSE IF (VALUES (I) == 'REVISED') THEN
                  PLTORIG = .FALSE.
               END IF
               CYCLE

  100          CONTINUE   !  --PLTLINE--

C  *           Another special case. (See CPSLINE above.)

               CALL SCANNR (BUFFER, FIRST, LAST, MARK)
               FIRST = MARK + 2
               CALL SCANNR (BUFFER, FIRST, LAST, MARK)
               NUMBER = MXLIST
               CALL TOKENS (BUFFER (FIRST:LAST), NUMBER, LIST)

               DO J = 1, NUMBER

C  *              Find each value in the dictionary and assign accordingly:

                  CALL LOOKUP (NLINDC, LINDIC, .TRUE., LIST (J), POINTR)
                  IF (POINTR <= 0) GO TO 970
                  PLTLINE (J) = LIST (J)
               END DO
               GO TO 5

  110          CONTINUE   !  --PRECISION--

               CALL LOOKUP (NPREDC, PREDIC, .TRUE., VALUES (I), POINTR)
               IF (POINTR <= 0) GO TO 970
               IF (VALUES (I) == 'ENGINEERING') THEN
                  PRECISION = 2
               END IF
               CYCLE

  120          CONTINUE   !  --THREED--

               CALL LOOKUP (NTHRDC, THRDIC, .TRUE., VALUES (I), POINTR)
               IF (POINTR <= 0) GO TO 970
               IF (VALUES (I) == 'TRUE') THEN
                  THREED = .TRUE.
               END IF
               CYCLE

  130          READ (VALUES (I), 1000) XAXIS
               CYCLE

  140          READ (VALUES (I), 1000) XMAX
               CYCLE

  150          READ (VALUES (I), 1000) XMIN
               CYCLE

  160          READ (VALUES (I), 1000) YMAX
               CYCLE

  170          READ (VALUES (I), 1000) YMIN

         END DO ! Next keyword/value pair

C  *     Look for another line of keywords:

      GO TO 5


  950 RETURN


C  *  Error handling:

  960 WRITE (LUNCRT, 1020) 'keyword', BUFFER
      GO TO 999

  970 WRITE (LUNCRT, 1020) 'value', BUFFER
      GO TO 999

  980 WRITE (LUNCRT, 1030) IOS

  999 STOP

C  *  Formats:

 1000 FORMAT (BN, F10.0)
 1020 FORMAT (' Abnormal termination - invalid ', A,
     >         ' in the following line:', /, 1X, A)
 1030 FORMAT (/, ' System error reading a line of keyword text.  IOS: ',
     >        I6)
      END SUBROUTINE RDKEYS
C+-----------------------------------------------------------------------

      SUBROUTINE RECTIFY (NU, XU, YU, NL, XL, YL, X, Y, MAXPTS, LUNCRT,
     >                    LUNKBD )
C
C  PURPOSE:  RECTIFY reorganizes airfoil geometry data so that the common
C            leading edge point is in fact the minimum abscissa.   It can
C            also shift the ordinates to ensure a user-specified leading-
C            edge ordinate, if desired.
C
C  METHOD:
C   *  Merge the separate surfaces as a single wrap-around surface.
C   *  Identify the subscript corresponding to minimum abscissa.
C   *  Separate into two surfaces again, with truly monotonically increasing
C      abscissas.
C   *  Shift ordinates if required.
C
C  ARGUMENTS:
C   ARG    DIM  I/O/S   DESCRIPTION
C   NU      -     I     Number of points on upper surface before and after
C                       rearranging
C   XU   MAXPTS  I/O    Upper surface abscissas in ascending order
C   YU   MAXPTS  I/O    Upper surface ordinates in ascending order
C   NL      -    I      Number of points on lower surface before and after
C                       rearranging
C   XL   MAXPTS  I/O    Lower surface abscissas in ascending order
C   YL   MAXPTS  I/O    Lower surface ordinates in ascending order
C   X   MAXPTS*2    S   Abscissas of both surfaces (wrap-around order)
C   Y   MAXPTS*2    S   Ordinates of both surfaces (wrap-around order)
C  MAXPTS   -    I      Maximum number of points allowed for on a surface
C  LUNCRT   -    I      Logical unit for prompts (screen)
C  LUNKBD   -    I      Logical unit for responses (keyboard)
C
C  EXTERNAL REFERENCES:
C
C  COPY      Copies one array to another
C  RVERSE    Reverses order of lower surface points for wrap-around ordering
C  READER    Prompts for/reads integer, real, etc.; handles <CR>
C
C  HISTORY:
C
C    04/29/83    LJC     Original design and coding
C    02/13/85    LJC     Added shifting of ordinates
C    02/26/85    LJC     Added protection against shifting thicknesses
C
C  AUTHOR: Leslie Collins, Informatics General Corporation, Palo Alto, CA
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, MAXPTS, NL, NU

      REAL
     >   XU(MAXPTS), YU(MAXPTS), XL(MAXPTS), YL(MAXPTS),  X(MAXPTS*2),
     >   Y(MAXPTS*2)

C  *  Local variables:

      INTEGER
     >   I, NPTS

      REAL
     >   LE, SHIFT

      CHARACTER
     >   YESNO

      LOGICAL
     >   DEFAULT, QUIT

C  *  Execution:

C  *  Store coordinates of both surfaces in one wrap-around array:

      CALL RVERSE (NL, XL, X)
      CALL RVERSE (NL, YL, Y)
      CALL COPY   (NU, XU, X(NL))
      CALL COPY   (NU, YU, Y(NL))
      NPTS = NU + NL - 1

C  *  Search for subscript of minimum abscissa:

      DO I = 2, NPTS
         IF (X(I) < X(I-1)) NL = I
      END DO

      NU = NPTS - NL + 1

C  *  Reverse lower surface points so that abscissas are in ascending
C  *  order:

      CALL RVERSE (NL, X, XL)
      CALL RVERSE (NL, Y, YL)

C  *  Store remaining elements of X and Y in XU and YU arrays:

      CALL COPY (NU, X(NL), XU)
      CALL COPY (NU, Y(NL), YU)

C  *  Permit vertical shift if desired:

      WRITE (LUNCRT, 1000)
     >   ' You have an option to shift the airfoil vertically.'
      WRITE (LUNCRT, 1001)
     >   ' The current leading edge ordinate is ', YU(1)
      CALL READR (LUNCRT,
     >   'Enter a new leading edge ordinate or <CR> to leave as is: ',
     >   LUNKBD, LE, DEFAULT, QUIT)

      IF (.NOT. DEFAULT) THEN

         SHIFT = YU(1) - LE
         DO I = 1, NU
            YU(I) = YU(I) - SHIFT
         END DO

C  *     Protect user from shifting a thickness distribution:

         WRITE (LUNCRT, 1000)
     >      ' Is this a camber/thickness distribution? (Y/N, <CR>=N)'
         CALL READC (LUNCRT,
     >      '(If so, the thickness will not be shifted.) ',
     >      LUNKBD, YESNO, DEFAULT, QUIT)

         IF (YESNO == 'N' .OR. DEFAULT) THEN
            DO I = 1, NL
               YL(I) = YL(I) - SHIFT
            END DO
         END IF

      END IF

 1000 FORMAT (/, A)
 1001 FORMAT (A, F10.6)

      END SUBROUTINE RECTIFY
C+------------------------------------------------------------------------------
C
      SUBROUTINE REDISTRIB (NU, XU, YU, NL, XL, YL, X, Y, MAXPTS,
     >                      XEVAL, YEVAL, COEFS, LUNCRT, LUNKBD, LUNXGR)
C
C  PURPOSE:  REDISTRIB redistributes the data points defining a given
C            airfoil profile, using one of several possible methods
C            (uniform, sinusoidal, or Vinokur-type distributions).
C            If the leading edge is rounded, the distributions may be
C            specified as being along the arc if this is preferred over
C            the standard chord-wise distribution.
C
C  METHOD:   A menu is used to determine how the revised coordinates are
C            produced.  The new Xs (or Ts) are either generated over the
C            original data range in a choice of ways, or read from a disk
C            file in standard PROFILE format (TITLE and NU followed by
C            a leading-to-trailing-edge distribution, then either NL
C            and a second leading-to-trailing-edge distribution, or NL=0
C            or, equivalently, EOF).  Thus the disk file may contain points
C            for one or two surfaces.  (A single distribution read will be
C            reused if necessary.)  Any columns after the first are ignored.
C
C            If the airfoil is treated as a single curve, the parametric
C            conventional spline of PSFIT is applied.  Otherwise, several
C            choices of spline are offered.  The original coordinates are
C            overwritten by the redistributed points.
C
C  ARGUMENTS:
C   ARG      DIM    TYPE I/O/S DESCRIPTION
C  NU,NL      -       I   I/O  Number of upper/lower surface pts.
C  XU,XL    NU,NL     R   I/O  Given abscissas, upper/lower.
C  YU,YL    NU,NL     R   I/O  Corresponding ordinates, upper/lower.
C                              The surfaces must have common leading edge
C                              points if they are to be treated as rounded,
C                              and common trailing edge points if they are
C                              to be treated as closed smoothly.
C  X,Y     NU+NL-1    R    S   Workspace for the wraparound form of the
C                              airfoil, if a parametric spline is used.
C  MAXPTS     -       I    I   Max. no. of pts. provided for on 1 surface.
C  XEVAL,YEVAL MAXPTS R    S   Workspace for new abscissas/ordinates.
C  COEFS   MAXPTS*3   R    S   Coefficients used by nonparametric splines;
C                              handy as workspace for the parametric spline.
C  LUNCRT     -       I    I   Logical unit for prompts.
C  LUNKBD     -       I    I   Logical unit for responses.
C  LUNXGR     -       I    I   Logical unit for disk file of input
C                              abscissas or Ts in standard PROFILE format.
C                              (Ordinates are ignored if present.)
C  PROCEDURES:
C    COPY         Copies an array
C    CSEVAL       Evaluates the spline previously fit by CSFIT
C    CSFIT        Fits a conventional cubic interpolating spline
C    DSTRIB       Generates indicated distribution
C    FOILGRD      Linear + Quadratic + Sine + Cosine distribution
C    GETDIS       Prompts for details of generating a distribution
C    GETLINE      File reading utility
C    GETXFORM     Shift/scale factors for transforming [A, B] to [P, Q]
C    HTDIS2       Vinokur's 1-D distribution method
C    LCSFIT       Local cubic spline utility with monotonic and
C                 piecewise linear options as well
C    OPENER       File opening utility
C    PSFIT        Fits an parametric conventional or monotonic cubic spline
C    PSEVAL       Evaluates the PSFIT spline at specified Xs
C    PSTVAL       Evaluates the PSFIT spline at specified Ts
C    READER       Prompting utility
C    RVERSE       Reverses the order of an array
C    TOKENS       Tokenizes a string
C    TSUBJ        Extracts parameter value T associated with data point J
C    USESCALE     Applies the transformation from GETXFORM
C
C  HISTORY:
C  10/28/83  DAS  Code extracted as a module from early PROFILE.
C  07/20/84  LJC  Replaced calls to SPLINE and SEVAL with CSFIT and CSEVAL.
C  08/02/84  LJC  Moved opening of abscissa file from main program and added
C                 prompts for NUEVAL and NLEVAL (originally in namelist).
C  10/20/84  LJC  Modified calls to new versions of PSFIT and PSEVAL.
C  02/19/85  DAS  Handled rounded trailing edges using PSFIT's CLOSED option;
C                 confined rounded leading-edge handling to this routine,
C                 where BLUNT had originally been an input.
C  06/17/85  DAS  Incorporated GETDIS/DSTRIB in place of XGRID and in-line
C                 prompts. (DSTRIB offers generalized distributions -
C                 fractional powers of cosine allowed.)
C  07/23/85  RGL  Allowed for possibility that lower surface only
C                 would be redistributed from a file of abscissas.
C  09/05/85  DAS  Adjusted [T1,T2] definition to reflect modified
C                 parametric spline routines, which use cumulative
C                 arc length now for the parametric variable.
C  08/25/88  DAS  Provided for choosing a monotonic spline in place of
C                 the conventional spline, for nonparametric cases.
C  02/07/90  DAS  Introduced OPENER and GETLINE.
C  03/07/90  DAS  Bug fix: GETDIS return of "quit" was not acted upon.
C  03/23/91  DAS  Replaced MSFIT with LCSFIT in order to provide a
C                 piecewise linear option for wedge-type sections.
C  10/16/91  DAS  PSFIT's calling sequence changed (METHOD); TOL changed
C                 to 1.E-6 instead of 1.E-5 since it is used relatively now.
C  05/04/92  DAS  Provided Vinokur's distribution as an option; provided
C                 for distributions in terms of arc length as well as X;
C                 reused upper surface points read if none found for the
C                 lower surface.
C  08/24/93  DAS  Arc-based redistribution now uses monotonic spline for
C                 X vs. T to avoid possible excursion at the leading edge.
C                 Y vs. T remains conventional.
C  04/01/95  DAS  Installed FOILGRID option.
C  10/18/96  DAS  Replaced FOILGRID with FOILGRD.
C
C  AUTHORS: David Saunders/Leslie Collins, Sterling/NASA Ames, Mt. View, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, LUNXGR, MAXPTS, NL, NU

      REAL
     >   COEFS (MAXPTS, 3), X (MAXPTS*2), XEVAL (MAXPTS), XL (NL),
     >   XU (NU), Y (MAXPTS*2), YEVAL (MAXPTS), YL (NL), YU (NU)

C     Local constants:

      REAL, PARAMETER ::
     >   PERCENT = 1.E-2,
     >   TOL     = 1.E-6,
     >   ZERO    = 0.E+0

      CHARACTER, PARAMETER ::
     >   BLANK * 1   = ' ',
     >   COMMENT * 1 = '!',
     >   IFMT * 8    = '(BN,I20)',
     >   RFMT * 10   = '(BN,F20.0)'

C     Local variables:

      INTEGER
     >   I, IER, IOS, LAST, MODE, NLEVAL, NUEVAL, NUMBER, NP, NPTS

      REAL
     >   A, B, DUMMY, P (4), SCALE, SHIFT, T1, T2

      LOGICAL
     >   ARC, BLUNT, CLOSED, DEFAULT, OPENED, QUIT, REUSE, YES

      CHARACTER
     >   METHOD * 1, TOKEN * 25, XFILE * 50

C     Procedures:

      REAL
     >   TSUBJ

      EXTERNAL
     >   COPY, CSEVAL, CSFIT, DSTRIB, FOILGRD, GETDIS, GETLINE,
     >   GETXFORM, HTDIS2, LCSFIT, OPENER, PSFIT, PSEVAL, PSTVAL,
     >   READC, READY, RVERSE, TOKENS, TSUBJ, USESCALE

C     Execution:

C     Determine whether airfoil is to be treated as one curve or two:

      WRITE (LUNCRT, 1001)
     >   '    PROFILE treats the airfoil as one curve if the leading ',
     >   '    edge is rounded. In this case, if the trailing edge is ',
     >   '    also rounded, the curve is closed smoothly. Otherwise, ',
     >   '    the trailing edge is sharp or open (no continuity).',
     >   BLANK

      BLUNT = .TRUE.
      CALL READY (LUNCRT,
     >   'Is the  LEADing edge rounded? (Y/N; <CR> means Y(es)): ',
     >   LUNKBD, BLUNT, DEFAULT, QUIT)
      IF (QUIT) GO TO 999

      CLOSED = .FALSE.

      IF (BLUNT) THEN
         CALL READY (LUNCRT,
     >      'Is the TRAILing edge rounded? (Y/N; <CR> means NO): ',
     >      LUNKBD, CLOSED, DEFAULT, QUIT)
         IF (QUIT) GO TO 999

      ELSE  ! Don't try to handle a sharp leading edge/rounded trailing edge

C        Non-rounded airfoils may require a choice of interpolation methods:

         WRITE (LUNCRT, 1001)
     >     '    Conventional splines are normally appropriate, but the',
     >     '    monotonic spline or piecewise linear interpolation may',
     >     '    be preferable for angular sections.',
     >     BLANK

  100    METHOD = 'C'
         CALL READC (LUNCRT,
     > 'Conventional spline, Monotonic, or Linear? (C/M/L; <CR> = C): ',
     >      LUNKBD, METHOD, DEFAULT, QUIT)
         IF (QUIT) GO TO 999
         IF (METHOD /= 'C' .AND.
     >       METHOD /= 'M' .AND. METHOD /= 'L') GO TO 100
      END IF

      ARC = .FALSE.          ! For nonparametric cases

      IF (BLUNT) THEN        ! Parametric spline permits two choices:

         METHOD = 'S'
         CALL READC (LUNCRT,
     >      'Distribute along the chord or the arc? (X/S; <CR> = S): ',
     >      LUNKBD, METHOD, DEFAULT, QUIT)
         IF (QUIT) GO TO 999
         ARC = METHOD /= 'X'

C        Set up airfoil as a single wraparound curve, eliminating the
C        duplicate leading edge point.

         IF (XU (1) /= XL (1)  .OR.  YU (1) /= YL (1)) GO TO 830
         IF (CLOSED  .AND.
     >      (XU (NU) /= XL (NL) .OR. YU (NU) /= YL (NL))) GO TO 840

         CALL RVERSE (NL, XL, X)
         CALL RVERSE (NL, YL, Y)
         CALL COPY (NU - 1, XU (2), X (NL + 1))
         CALL COPY (NU - 1, YU (2), Y (NL + 1))
         NPTS = NL + NU - 1

C        Fit a parametric spline, using the "monotonic" option for X vs. T
C        (and a conventional spline for Y vs. T), to avoid possible trouble
C        at the leading edge of blunt airfoils.

         CALL PSFIT (NPTS, X, Y, 'MC', CLOSED, IER)
         IF (IER /= 0) GO TO 900

      END IF

      OPENED = .FALSE.   ! For possible input file of Xs or Ts.


C     Upper surface:
C     -------------

  200 CONTINUE

      YES = .TRUE.
      CALL READY (LUNCRT,
     >   'Do you want to redistribute the UPPER surface? (<CR>=Y): ',
     >   LUNKBD, YES, DEFAULT, QUIT)
      IF (QUIT) GO TO 999

      IF (YES) THEN

         IF (ARC) THEN
            A = TSUBJ (NL)    ! Cumulative chord length limits from PSFIT
            B = TSUBJ (NPTS)  ! This is the range of T to be distributed
         ELSE
            A = XU (1)        ! This is the X range to be distributed
            B = XU (NU)
         END IF

C        Determine the details for this redistribution:

         NP = 2
         CALL GETDIS (LUNCRT, LUNKBD, MODE, NP, P, NUEVAL)

         IF (MODE < -1) THEN

C           User must have changed his mind - quit:

            GO TO 999

         ELSE IF (MODE == -1) THEN

C           Read from file with new Xs (or Ts; standard PROFILE format):

            XFILE = 'xgrid.dat'
            CALL OPENER (LUNCRT,
     >         'Enter filename for new Xs or Ts. <CR> uses xgrid.dat:',
     >         LUNKBD, XFILE, LUNXGR, 'OLD')

            OPENED = .TRUE.

            READ (LUNXGR, 1000, ERR = 870)     ! Skip title

C           Get the no. of pts.  XFILE is a handy character buffer.

            CALL GETLINE (LUNXGR, BLANK, XFILE, LAST, IOS)
            IF (IOS /= 0) GO TO 870
            NUMBER = 1
            CALL TOKENS (XFILE, NUMBER, TOKEN)
            READ (TOKEN, IFMT, ERR = 870) NUEVAL

            IF (NUEVAL > MAXPTS) GO TO 850

C           Pick off the first item on each line:

            DO 320, I = 1, NUEVAL
  310          CALL GETLINE (LUNXGR, COMMENT, XFILE, LAST, IOS)
               IF (IOS /= 0) GO TO 870
               IF (LAST == 0) GO TO 310

               CALL TOKENS (XFILE, NUMBER, TOKEN)
               READ (TOKEN, RFMT, ERR = 870) XEVAL (I)
  320       CONTINUE

            IF (.NOT. ARC) THEN   ! Require a full range of abscissas

               IF (XEVAL (1) /= XU (1) .OR.
     >             XEVAL (NUEVAL) /= XU (NU)) GO TO 880

            ELSE  ! Distribution must be l.e. to t.e. but may be relative

               CALL GETXFORM (XEVAL (1), XEVAL (NUEVAL), A, B,
     >                        SCALE, SHIFT)

               CALL USESCALE ('D', 1, NUEVAL, XEVAL, DUMMY, DUMMY,
     >                        SCALE, SHIFT, IER)

               XEVAL (1) = A      ! Avoid round-off at the end points
               XEVAL (NUEVAL) = B

            END IF

         ELSE     ! Generate the specified distribution

            IF (NUEVAL > MAXPTS) GO TO 850

            IF (MODE <= 3) THEN   ! Uniform or sinusoidal-type

               CALL DSTRIB (MODE, NP, P, NUEVAL, A, B, XEVAL)

            ELSE IF (MODE == 4) THEN   ! Vinokur distribution

C              Relative or absolute?

               IF (P (1) < ZERO) P (1) = -PERCENT * (B - A) * P (1)
               IF (P (2) < ZERO) P (2) = -PERCENT * (B - A) * P (2)

               CALL HTDIS2 (.TRUE., A, B, P (1), P (2), NUEVAL, XEVAL,
     >                      -LUNCRT, IER)
               IF (IER /= 0) GO TO 855

            ELSE IF (MODE == 5) THEN   ! Linear + Quadratic + Sine + Cosine

               CALL FOILGRD (NUEVAL, A, B, P (1), P (2), P (3), P (4),
     >                       XEVAL)
            END IF

         END IF


C        Now evaluate the new upper surface coordinates:

         IF (BLUNT) THEN

            IF (ARC) THEN     ! Evaluate X and Y at the indicated Ts

C              Make use of COEFS (*, 1) for the Ts, and COEFS (*, 2)
C              and COEFS (*, 3) for the unneeded 1st & 2nd derivatives.

               CALL COPY (NUEVAL, XEVAL, COEFS)

               CALL PSTVAL (NUEVAL, COEFS (1, 1), XEVAL, YEVAL,
     >                      COEFS (1, 2), COEFS (1, 3),
     >                      COEFS (1, 2), COEFS (1, 3), NPTS, X, Y)

            ELSE              ! Evaluate Y at the indicated Xs

               T1 = TSUBJ (NL)   ! Define the relevant curve that is
               T2 = TSUBJ (NPTS) ! monotonic in X

               CALL PSEVAL (NUEVAL, XEVAL, YEVAL, COEFS, COEFS, T1,
     >                      T2, NPTS, X, Y, TOL, IER)
               IF (IER /= 0) GO TO 910

            END IF

         ELSE

C           Separate spline for each surface (conventional, monotonic or linear)

            IF (METHOD == 'C') THEN

               CALL CSFIT (NU, XU, YU, 0, DUMMY, 0, DUMMY, COEFS (1,1),
     >                     COEFS (1,2), COEFS (1,3), IER)
               IF (IER /= 0) GO TO 890

               CALL CSEVAL (NU, XU, YU, NUEVAL, XEVAL, COEFS (1,1),
     >                      COEFS (1,2), COEFS (1,3), YEVAL)
            ELSE

               CALL LCSFIT (NU, XU, YU, .TRUE., METHOD, NUEVAL, XEVAL,
     >                      YEVAL, YEVAL)
            END IF

         END IF

C        Overwrite the original coordinates with the redistributed ones:

         CALL COPY (NUEVAL, XEVAL, XU)
         CALL COPY (NUEVAL, YEVAL, YU)
         NU = NUEVAL

      END IF


C     Repeat for lower surface:
C     ------------------------

  400 CONTINUE

      YES = .TRUE.
      CALL READY (LUNCRT,
     >   'Do you want to redistribute the LOWER surface? (<CR>=Y): ',
     >   LUNKBD, YES, DEFAULT, QUIT)
      IF (QUIT) GO TO 999

      IF (YES) THEN

C        Interval to be distributed:

         IF (ARC) THEN
            A = ZERO
            B = TSUBJ (NL)
         ELSE
            A = XL (1)
            B = XL (NL)
         END IF

C        Determine the details for this redistribution:

         NP = 2
         CALL GETDIS (LUNCRT, LUNKBD, MODE, NP, P, NLEVAL)

         IF (MODE < -1) THEN

C           User must have changed his mind - quit:

            GO TO 999

         ELSE IF (MODE == -1) THEN

            IF (.NOT. OPENED) THEN

C              Open file of new Xs or Ts:

               XFILE = 'xgrid.dat'
               CALL OPENER (LUNCRT,
     >            'Enter filename for new abscissas. <CR>=xgrid.dat:',
     >            LUNKBD, XFILE, LUNXGR, 'OLD')

C              Skip upper surface (unused but may be present).
C              HOWEVER: Store it in case only one surface is present.

               READ (LUNXGR, 1000, ERR = 870)  ! Skip title

               CALL GETLINE (LUNXGR, BLANK, XFILE, LAST, IOS)
               IF (IOS /= 0) GO TO 870
               NUMBER = 1
               CALL TOKENS (XFILE, NUMBER, TOKEN)
               READ (TOKEN, IFMT, ERR = 870) NUEVAL

               DO 490, I = 1, NUEVAL
  480             CALL GETLINE (LUNXGR, COMMENT, XFILE, LAST, IOS)
                  IF (IOS /= 0) GO TO 870
                  IF (LAST == 0) GO TO 480

                  CALL TOKENS (XFILE, NUMBER, TOKEN)
                  READ (TOKEN, RFMT, ERR=870) XEVAL (I)
  490          CONTINUE
            END IF

C           Read number of lower surface abscissas (if any):

            CALL GETLINE (LUNXGR, BLANK, XFILE, LAST, IOS)
            IF (IOS < 0) THEN
               REUSE = .TRUE.     ! EOF means same as NL = 0
            ELSE IF (IOS /= 0) THEN
               GO TO 870
            ELSE
               NUMBER = 1
               CALL TOKENS (XFILE, NUMBER, TOKEN)
               READ (TOKEN, IFMT, ERR = 870) NLEVAL

               REUSE = NLEVAL == 0
               IF (NLEVAL > MAXPTS) GO TO 860
            END IF

            IF (REUSE) THEN

               NLEVAL = NUEVAL

               IF (ARC) THEN
                  IF (OPENED) THEN

C                    Recover the relative Ts from the upper surface.
C                    They have been transformed, but transforming then
C                    again is OK apart from round-off.

                     CALL COPY (NLEVAL, COEFS, XEVAL)
                     WRITE (LUNCRT,1001) ' Reusing upper surface Ts ...'
                  END IF
               ELSE IF (OPENED) THEN
                  WRITE (LUNCRT, 1001) ' Reusing upper surface Xs ...'
               END IF

            ELSE  ! Read lower surface Xs (first item on each valid line)

               DO 520, I = 1, NLEVAL
  510             CALL GETLINE (LUNXGR, COMMENT, XFILE, LAST, IOS)
                  IF (IOS /= 0) GO TO 870
                  IF (LAST == 0) GO TO 510

                  CALL TOKENS (XFILE, NUMBER, TOKEN)
                  READ (TOKEN, RFMT, ERR = 870) XEVAL (I)
  520          CONTINUE

            END IF

            IF (.NOT. ARC) THEN   ! Require a full range of abscissas

               IF (XEVAL (1) /= XL (1) .OR.
     >             XEVAL (NLEVAL) /= XL (NL)) GO TO 880

            ELSE  ! Transform relative arc distribution to desired interval

               CALL GETXFORM (XEVAL (1), XEVAL (NLEVAL), A, B,
     >                        SCALE, SHIFT)

               CALL USESCALE ('D', 1, NLEVAL, XEVAL, DUMMY, DUMMY,
     >                        SCALE, SHIFT, IER)

               XEVAL (1) = A      ! Avoid round-off
               XEVAL (NLEVAL) = B

            END IF

         ELSE     ! Generate the specified distribution

            IF (NLEVAL > MAXPTS) GO TO 860

            IF (MODE <= 3) THEN  ! Uniform or sinusoidal-type

               CALL DSTRIB (MODE, NP, P, NLEVAL, A, B, XEVAL)

            ELSE IF (MODE == 4) THEN  ! Vinokur distribution

C              Relative or absolute?

               IF (P (1) < ZERO) P (1) = -PERCENT * (B - A) * P (1)
               IF (P (2) < ZERO) P (2) = -PERCENT * (B - A) * P (2)

               CALL HTDIS2 (.TRUE., A, B, P (1), P (2), NLEVAL, XEVAL,
     >                      -LUNCRT, IER)
               IF (IER /= 0) GO TO 865

            ELSE IF (MODE == 5) THEN   ! Linear + Quadratic + Sine + Cosine

               CALL FOILGRD (NLEVAL, A, B, P (1), P (2), P (3), P (4),
     >                       XEVAL)
            END IF

         END IF


C        Now evaluate the new lower surface coordinates:

         IF (BLUNT) THEN

            IF (ARC) THEN     ! Evaluate X and Y at the indicated Ts

C              [0, T(l.e.)] has been treated as though it were
C              [T(l.e.), T(t.e.)], so the Ts need to be transformed.
C              (We want to come out with X/YEVAL(*) going from l.e. to t.e.)

               CALL GETXFORM (XEVAL (1), XEVAL (NLEVAL), B, A,
     >                        SCALE, SHIFT)

               CALL COPY (NLEVAL, XEVAL, COEFS) ! Use COEFS (*, 1) for Ts

               CALL USESCALE ('D', 1, NLEVAL, COEFS, DUMMY, DUMMY,
     >                        SCALE, SHIFT, IER)

               COEFS (1, 1) = B      ! Avoid round-off
               COEFS (NLEVAL, 1) = A

               CALL PSTVAL (NLEVAL, COEFS (1, 1), XEVAL, YEVAL,
     >                      COEFS (1, 2), COEFS (1, 3),
     >                      COEFS (1, 2), COEFS (1, 3), NPTS, X, Y)

            ELSE              ! Evaluate Y at the indicated Xs

               T1 = ZERO      ! Define monotonic-in-X sub-curve
               T2 = TSUBJ (NL)

               CALL PSEVAL (NLEVAL, XEVAL, YEVAL, COEFS, COEFS, T1,
     >                      T2, NPTS, X, Y, TOL, IER)
               IF (IER /= 0) GO TO 910

            END IF

         ELSE

            IF (METHOD == 'C') THEN

               CALL CSFIT (NL, XL, YL, 0, DUMMY, 0, DUMMY, COEFS (1,1),
     >                     COEFS (1,2), COEFS (1,3), IER)
               IF (IER /= 0) GO TO 890

               CALL CSEVAL (NL, XL, YL, NLEVAL, XEVAL, COEFS (1,1),
     >                      COEFS (1,2), COEFS (1,3), YEVAL)
            ELSE

               CALL LCSFIT (NL, XL, YL, .TRUE., METHOD, NLEVAL, XEVAL,
     >                      YEVAL, YEVAL)
            END IF

         END IF

C        Overwrite the original coordinates with the redistributed ones:

         NL = NLEVAL
         CALL COPY (NLEVAL, XEVAL, XL)
         CALL COPY (NLEVAL, YEVAL, YL)

      END IF

      RETURN


C     Error handling:

  830 WRITE (LUNCRT, 1001)
     >   ' The two surfaces have different leading points - quitting.'
      GO TO 999
  840 WRITE (LUNCRT, 1001)
     >   ' The two surfaces have different trailing points - quitting.'
      GO TO 999
  850 WRITE (LUNCRT, 1002) MAXPTS
      GO TO 200
  855 WRITE (LUNCRT, 1003) 'HTDIS2', IER
      GO TO 200
  860 WRITE (LUNCRT, 1002) MAXPTS
      GO TO 400
  865 WRITE (LUNCRT, 1003) 'HTDIS2', IER
      GO TO 400
  870 WRITE (LUNCRT, 1001) ' Error reading file of input abscissas.'
      GO TO 999
  880 WRITE (LUNCRT, 1001) ' Must supply full range of abscissas.'
      GO TO 999
  890 WRITE (LUNCRT, 1003) 'CSFIT', IER
      GO TO 999
  900 WRITE (LUNCRT, 1003) 'PSFIT', IER
      GO TO 999
  910 WRITE (LUNCRT, 1003) 'PSEVAL', IER
C*****GO TO 999

  999 WRITE (LUNCRT, 1001) ' Stopping in REDISTRIB.'
      STOP ' '

C     Formats:

 1000 FORMAT (A)
 1001 FORMAT (/, (A))
 1002 FORMAT (/, ' Can''t handle so many abscissas. Limit = ', I3)
 1003 FORMAT (/, ' Error in subroutine ', A, '.  IER: ', I3)

      END SUBROUTINE REDISTRIB
C+----------------------------------------------------------------------
C
      SUBROUTINE REFINE (NU, XU, YU, NL, XL, YL, TITLE,
     +                   LUNCRT, LUNKBD, LUNTAB, LUNYPP)
C
C  PURPOSE:  REFINE  perturbs the surface(s) of the given airfoil in  a
C            way that strives for a specified thickness while retaining
C            the original curvature (actually 2nd derivative) distribu-
C            tion as much as possible, specially near the leading edge.
C            It also provides for refining  the curvature distribution,
C            offering a variety of ways to indicate target y" values.
C
C  METHOD:   Solve the overdetermined system obtained from the  simple-
C            minded scaling,
C
C                         ynew(i) = yscale(i) * yold(i)    (i = 2, n-1)
C
C            along with constraints obtained from the finite difference
C            expression for the 2nd derivative at each  interior point:
C
C                         ynew"(i) = old or desired y"(i)  (i = 2, n-1)
C
C            Provision is made for weighting (down) the y" portion.
C
C            Note:   ynew(i) = yscale(i)*yold(i) exactly (i = 1 and n)
C
C            because the leading edge point should come out EXACTLY the
C            same for both surfaces (not just close), and the same  MAY
C            be true of a sharp trailing edge.
C
C            The very sparse system involved was originally solved by a
C            dense method (HDESOL, still shown but commented out); this
C            has since been replaced with a specialized  sparse  method
C            (DTDLSQ).
C
C            Since the desired thickness is not likely to  be  obtained
C            exactly by this linear least squares approach,  the method
C            is iterated until the thickness is arbitrarily close.
C
C            Any target y" values read from a file are assumed to be in
C            the standard PROFILE format for (x,y) pairs, although NU=0
C            or NL=0 is permitted, meaning leave that surface alone.
C
C  ARGUMENTS:
C  VAR   DIM  TYPE I/O/S DESCRIPTION
C  NU     -    I     I   Number of upper surface points
C  XU     NU   R     I   Upper surface abscissas (assumed increasing and
C  YU     NU   R    I/O  Lower surface ordinates             normalized)
C  NL     -    I     I   Number of lower surface points
C  XL     NL   R     I   Lower surface abscissas
C  YL     NL   R    I/O  Lower surface ordinates
C  TITLE  *    C     I   Tabulations title
C  LUNCRT -    I     I   Logical unit for screen
C  LUNKBD -    I     I   Logical unit for keyboard
C  LUNTAB -    I     I   Logical unit for tabulations
C  LUNYPP -    I     I   Logical unit for optional table of y" values
C
C  PARAMETER CONSTANTS:
C  PARAM   TYPE   DESCRIPTION:
C  MAXITER   I    Max. no. iterations for achieving target thickness
C  MAXN      I    MAXSRF-2, because leading and trailing edge points
C                 are not solved for with the interior points
C  MAXSRF    I    Max. no. pts. on one surface allowed by local arrays
C
C  SIGNIFICANT LOCAL VARIABLES:
C  VAR       DIM     TYPE  DESCRIPTION
C**A       2*N,N+1     R   For setting up 2N*N overdetermined system in
C**                        dense form (N=MAXSRF-2)
C**C         2*N       R   Work-space for HDESOL
C    NOTE: A, B, C, D, are of length MAXSRF, not MAXN, because they are
C          reused in the thickness calculation for spline coefs., etc.
C  A       MAXSRF      R   For sparse solution of overdetermined system
C  B       MAXSRF      R    "    "    "    "    "    "    "    "    "
C  C       MAXSRF      R    "    "    "    "    "    "    "    "    "
C  D       MAXSRF      R    "    "    "    "    "    "    "    "    "
C  R          N        R    "    "    "    "    "    "    "    "    "
C  S          N        R    "    "    "    "    "    "    "    "    "
C  YSCALE   MAXSRF     R   Scale factors applied to each y(I), possibly
C                          functions of I
C  YLORIG,  MAXSRF     R   Needed for refining result by repeating the
C  YUORIG                  analysis with an extrapolated target thickness
C
C  FILES USED:
C  LUN       I/O/S     DESCRIPTION
C  LUNCRT      O       User prompts
C  LUNKBD      I       User responses
C  LUNTAB      O       Tabulations
C  LUNYPP      I       Optional input table of y" targets
C
C  PROCEDURES:
C  CALCTHICK   Calculates airfoil thickness
C  COPY        Copies an array
C  LINTRP      Linear interpolation (for thickness iteration)
C  LSTFIL      Echoes the table of y" values to the tabulation file
C  OPENER      File opening utility
C  READER      Prompting utility
C  REFSRF      Performs the above algorithm on one surface at a time
C  RDXYZ       Used to (x, y") pairs - one surface per call
C
C  HISTORY:
C  06/12/95   DAS   Suppressed printing details of linear system in .TAB file.
C                   Thickness iterations still go to .TAB and screen.
C  07/30/93   DAS   Global SAVE introduced.
C  12/12/91   DAS   RDXYZ's calling sequence changed.
C  07/19/90   DAS   No! PRREAD doesn't permit 0 for either surface's
C                   target y" values.  Introduced RDXYZ instead.
C  02/07/90   DAS   Introduced OPENER; removed concatenation from I/O;
C                   revived use of PRREAD to read the y" values.
C  02/10/87   DAS   XATMAXTHICKNESS=1.0 caused failure (YBUMPPARAMS(1)).
C                   Safeguarded this abnormal case. Also safeguarded use
C                   on unnormalized airfoil - shape functions fail then.
C  09/15/86   DAS   PRREAD no longer used for target y" values - it does
C                   not handle NU=0 or NL=0 properly for this context.
C  04/16/86   DAS   Removed YP, YPP, YK arrays and the call to FD12K -
C                   the finite differencing is done in-line in REFSRF.
C  03/19/86   DAS   Streamlined the diagnostic I/O and control of the
C                   prompting done in REFSRF first time in/each surface.
C  10/24/85   DAS   Used LSTFIL in place of PRWRIT for echoing y" table.
C  11/05/84   DAS   Introduced DTDLSQ in place of HDESOL to save storage
C                   and CPU time.
C  10/25/84   DAS   Optional on-line help added; prompts clarified.
C  04/11/84   DAS   Eliminated handling of non-normalized airfoil.
C  01/30/84   LJC   Incorporated new calling sequence for CALCTHICK.
C  12/21/83   DAS   Made use of PRWRIT to print the table of y" values
C                   in the tabulation output (desirable because it is
C                   easy to be working with a table that is something
C                   other than what you thought it was).
C  11/03/83   DAS   Implemented no-y"-constraint case (as opposed to
C                   leave-surface-alone case - shoot for original y").
C  10/28/83   DAS   Renamed REFINE and integrated into PROFILE as an
C                   option distinct from MODIFY.
C  10/26/83   DAS   Incorporated table look-up for y" distributions.
C  10/20/83   DAS   Provided for retaining original curvature (without
C                   knowing what it is ahead of time).
C  10/19/83   DAS   Revised thickness refinement to handle case where
C                   desired thickness = original thickness.
C  10/12/83   DAS   Put y" weighting and constraining into REFSRF so
C                   that surfaces may be treated differently.
C  10/07/83   DAS   Changed back to true y" for each equation in the
C                   overdetermined system (rather than multiplying
C                   throughout by functions of DXI, DXIM1); bounded
C                   RHS y" values in center region (to flatten the
C                   upper curvature distribution somewhat); rethought
C                   the thickness-refinement iteration.
C  09/30/83   DAS   Added refinement step for more precise thickness;
C                   nonlinear weighting of y"; curvature data plotting.
C  09/27/83   DAS   Separated out REFSRF for treating both surfaces.
C  09/26/83   DAS   Put in fancier scaling using sine "bump" (R.Kennelly).
C  09/23/83   DAS   y" before and after from finite differencing now,
C                   for greater consistency.
C  09/20/83   DAS   Initial design and coding (y" from SPLINE; upper
C                   surface only).
C
C  AUTHOR: David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

      SAVE

C ... Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, LUNTAB, LUNYPP, NL, NU

      REAL
     >   XL (NL), XU (NU), YL (NL), YU (NU)

      CHARACTER
     >   TITLE * (*)

C ... Local constants:

      INTEGER, PARAMETER ::
     >   MAXITER   = 10,
     >   MAXSRF    = 301,
     >   MAXN      = MAXSRF - 2

      REAL, PARAMETER ::
     >   ONE       = 1.E+0,
     >   TOLERANCE = 5.E-5,
     >   ZERO      = 0.E+0

      CHARACTER, PARAMETER ::
     >   BLANK * 1 = ' '

C ... Local variables:

      INTEGER
     >   COLUMNS (2), I, IER, ITER, NPTSYPPTABLE (2)

C*****REAL      A (2*MAXN, MAXN+1), C (2*MAXN),  ! For dense proof-of-concept

      REAL
     >   A (MAXSRF), B (MAXSRF), C (MAXSRF), D (MAXSRF), R (MAXN),
     >   S (MAXN), DELTATHICKNESS, DESIREDTHICKNESS, BETTER1, BETTER2,
     >   ORIGINALTHICKNESS, TARGET1, TARGET2, TARGET3, XATMAXTHICKNESS,
     >   XYPPTABLE (MAXSRF, 2), YPPTABLE (MAXSRF, 2), YLORIG (MAXSRF),
     >   YPPWEIGHT (MAXSRF), YSCALE (MAXSRF), YUORIG (MAXSRF),
     >   YBUMPPARAMS (3), YBUMPWIDTH

      LOGICAL
     >   DEFAULT, DETAILS, FINIS, KEEPSAMETHICKNESS, LEAVEALONE (2),
     >   NOYPPTABLE (2), PROMPT, QUIT, YES

      CHARACTER
     >   BUFFER * 80, YPPTABLENAME * 50

C ... Procedures:

      EXTERNAL
     >   CALCTHICK, COPY, LINTRP, LSTFIL, OPENER, RDXYZ, READR, READY,
     >   REFSRF

C ... Execution:

      IF (XU (1) < ZERO .OR. XU (NU) > ONE .OR.
     >    XL (1) < ZERO .OR. XL (NL) > ONE) THEN
         WRITE (LUNCRT, 1005)
     >      ' "Refine" mode requires Xs in [0,1]. Use "Normalize."'
         GO TO 990
      END IF

      IF (MAX (NU, NL) > MAXSRF) THEN
         WRITE (LUNCRT, 1005)
     >      ' Not enough work-space. Recompile module REFINE.'
         GO TO 990
      END IF

C ... The following were moved from REFSRF because of problems with
C     SAVE and DATA statements:

      PROMPT = .TRUE.
      LEAVEALONE (1) = .FALSE.
      LEAVEALONE (2) = .FALSE.

C ... On-line help may encourage first-time users:

      YES = .FALSE.
      CALL READY (LUNCRT,
     >   'Do you want an explanation of the prompts to follow?   '//
     >   '(Y/N; <CR>=No): ',
     >   LUNKBD, YES, DEFAULT, QUIT)

      IF (YES) THEN
         WRITE (LUNCRT, 1001)
         CALL READY (LUNCRT,
     >      '                   <Hit RETURN for more.>',
     >      LUNKBD, YES, DEFAULT, QUIT)
         WRITE (LUNCRT, 1002)
      END IF

      DETAILS = .FALSE.        ! Suppress output to .TAB file
C****      CALL READY (LUNCRT,
C****     >   'Do you want full details of the refinement iterations? '//
C****     >   '(Y/N; <CR>=No): ',
C****     >   LUNKBD, DETAILS, DEFAULT, QUIT)

C ... Calculate the original thickness.  Use available workspace
C     for the spline coefficients and spline evaluation.

      CALL CALCTHICK (NL, NU, XL, XU, YL, YU, ORIGINALTHICKNESS,
     >                XATMAXTHICKNESS,
C****>                A (1, 1), A (1, 2), A (1, 3), A (1, 4))
     >                A, B, C, D)

      WRITE (LUNCRT, 1020) '0Original % thickness:', ORIGINALTHICKNESS

      DESIREDTHICKNESS = ORIGINALTHICKNESS
      CALL READR (LUNCRT,
     >   'Enter desired % thickness, or <CR> to keep same: ',
     >   LUNKBD, DESIREDTHICKNESS, KEEPSAMETHICKNESS, QUIT)

      DELTATHICKNESS = ABS (DESIREDTHICKNESS - ORIGINALTHICKNESS)
      IF (KEEPSAMETHICKNESS .OR. DELTATHICKNESS < 1.E-02) THEN

C ...    Scales on all the ordinates will be essentially 1.0, so
C        pick a reasonable width for the sine bump and don't prompt:

         YBUMPWIDTH = ONE
      ELSE IF (DELTATHICKNESS < 0.2E+0) THEN
         YBUMPWIDTH = 1.5E+0
      ELSE
         YBUMPWIDTH = 2.0E+0
         CALL READR (LUNCRT,
     >      'Enter power of sine to use for y scaling (<CR> = 2.): ',
     >      LUNKBD, YBUMPWIDTH, DEFAULT, QUIT)
      END IF

      WRITE (LUNTAB, 1005)
      WRITE (LUNTAB, 1010) '1Case: ', TITLE, BLANK
      WRITE (LUNTAB, 1020)
     >   ' Original thickness  :', ORIGINALTHICKNESS,
     >   ' Corresponding  x/c  :', XATMAXTHICKNESS,
     >   ' Desired  thickness  :', DESIREDTHICKNESS,
     >   ' y-scaling sine power:', YBUMPWIDTH

      YPPTABLENAME = BLANK
      CALL OPENER (LUNCRT,
     >   'Enter file name for target 2nd derivatives; <CR> if none: ',
     >   LUNKBD, YPPTABLENAME, LUNYPP, 'OLD')

      IF (YPPTABLENAME == BLANK) THEN
         NOYPPTABLE (1) = .TRUE.
         NOYPPTABLE (2) = .TRUE.
      ELSE
         WRITE (LUNTAB, 1010) BLANK, BLANK,
     >      ' Table of y" values used: ', YPPTABLENAME,
     >      ' Values of y" found in table follow:', BLANK, BLANK

         CALL LSTFIL (LUNYPP, LUNTAB, BUFFER)
         REWIND LUNYPP

C ...    Skip the y" file's title:

         READ (LUNYPP, '(A)')

C ...    Look for target y" data one surface at a time:

         COLUMNS (1) = 1
         COLUMNS (2) = 2
         DO I = 1, 2
            CALL RDXYZ (2, LUNYPP, LUNCRT, COLUMNS, MAXSRF,
     >                  NPTSYPPTABLE (I), XYPPTABLE (1, I),
     >                  YPPTABLE (1, I), YPPTABLE (1, I), FINIS, IER)
            IF (IER /= 0) GO TO 920

            NOYPPTABLE (I) = NPTSYPPTABLE (I) == 0
         END DO

      END IF

      CALL COPY (NL, YL, YLORIG)
      CALL COPY (NU, YU, YUORIG)

      YBUMPPARAMS (2) = YBUMPWIDTH
      YBUMPPARAMS (1) = XATMAXTHICKNESS
      IF (XATMAXTHICKNESS == ZERO .OR.
     >    XATMAXTHICKNESS == ONE) YBUMPPARAMS (1) = 0.5E+0

C ... Begin thickness/second-derivative computation - iterated because
C     the desired thickness is not achieved exactly by solution of the
C     overdetermined systems involved.  The target thickness is varied
C     above or below the desired thickness till the computed thickness
C     is arbitrarily close to the desired one.

      ITER = 0
      TARGET3 = DESIREDTHICKNESS

  200 CONTINUE

         YBUMPPARAMS (3) = TARGET3 / ORIGINALTHICKNESS - ONE

C ...    Modify each surface separately:

         CALL REFSRF (1, NU, XU, YU, YSCALE,
     >                YBUMPPARAMS, YPPWEIGHT,
C****>                A, C,
     >                A, B, C, D, R, S,
     >                NOYPPTABLE (1), NPTSYPPTABLE (1),
     >                XYPPTABLE (1, 1), YPPTABLE (1, 1),
     >                DETAILS, PROMPT, LEAVEALONE (1),
     >                LUNCRT, LUNKBD, LUNTAB)

         CALL REFSRF (2, NL, XL, YL, YSCALE,
     >                YBUMPPARAMS, YPPWEIGHT,
C****>                A, C,
     >                A, B, C, D, R, S,
     >                NOYPPTABLE (2), NPTSYPPTABLE (2),
     >                XYPPTABLE (1, 2), YPPTABLE (1, 2),
     >                DETAILS, PROMPT, LEAVEALONE (2),
     >                LUNCRT, LUNKBD, LUNTAB)

C ...    Reestimate the actual thickness achieved:

         IF (ITER == 0) THEN
            WRITE (LUNTAB, 1010)
            WRITE (LUNCRT, 1010)
         ELSE
            BETTER1 = BETTER2
         END IF

         CALL CALCTHICK (NL, NU, XL, XU, YL, YU, BETTER2,
     >                   XATMAXTHICKNESS,
C****>                   A (1, 1), A (1, 2), A (1, 3), A (1, 4))
     >                   A, B, C, D)

         WRITE (LUNTAB, 1040) ' Itn.:', ITER,
     >      'Current % thickness:', BETTER2,
     >      'Corresponding abscissa:', XATMAXTHICKNESS
         WRITE (LUNCRT, 1040) ' Itn.:', ITER,
     >      'Current % thickness:', BETTER2,
     >      'Corresponding x/c:', XATMAXTHICKNESS

C ...    Check for refining this new thickness:

         ITER = ITER + 1
         DELTATHICKNESS = BETTER2 - DESIREDTHICKNESS
         IF (ABS (DELTATHICKNESS) > TOLERANCE) THEN

            IF (ITER < MAXITER) THEN

               IF (ITER == 1) THEN

C ...             Generate second calibration point prior to
C                 linear interpolation:

                  TARGET1 = ORIGINALTHICKNESS
                  TARGET2 = DESIREDTHICKNESS - DELTATHICKNESS
                  TARGET3 = TARGET2

               ELSE IF (ITER >= 3) THEN

C ...             Do a normal iteration, whereas passes 1, 2 were special:

                  TARGET1 = TARGET2
                  TARGET2 = TARGET3
               END IF

               IF (ITER >= 2) THEN
                  CALL LINTRP (1, TARGET1, TARGET2, BETTER1, BETTER2,
     >                         DESIREDTHICKNESS, TARGET3)
               END IF

               CALL COPY (NL, YLORIG, YL)
               CALL COPY (NU, YUORIG, YU)
               GO TO 200

            ELSE
               WRITE (LUNCRT, 1003)
            END IF

         END IF

      GO TO 999

  920 WRITE (LUNCRT, 1030)
     >   ' REFINE: Error reading y" table - aborting. ',
     >   ' IER from RDXYZ: ', IER

  990 STOP

  999 RETURN

C ... Formats:

 1001 FORMAT
     >  (/' REFINE works with  2N  equations in  N  unknowns (the ys):',
     >  //'        N of the form         y (new) = scale * y (old)',
     >   /'    and N of the form   wt * y" (new) = wt * y" (desired)',
     >  //' where  wt  represents weighting of the  second  derivative',
     >   /' equations to equilibrate the two halves of the system.',
     >  //' The first half enables thickness to be changed,  while the',
     >   /' second half enables the second derivatives  (and hence the',
     >   /' curvature distribution) to be smoothed.    The scaling and',
     >   /' the weighting use  "sine"  shape functions which  must  be',
     >   /' controlled by you the user  -  hence the series of prompts',
     >   /' to be discussed next.',//)
 1002 FORMAT
     > (//' Thickness ratio refinement:',
     >  //'    The nonlinear y scaling is intended to preserve as much',
     >   /'    as possible the curvature near the leading and trailing',
     >   /'    edges.    The sine function is centered at the point of',
     >   /'    maximum thickness.    Higher powers of the sine tend to',
     >   /'    preserve the leading/trailing-edge properties better.',
     >  //' Curvature smoothing:',
     >  //'    Typical weighting of the y" equations varies from 0.004',
     >   /'    at the leading and trailing edges to 0.04 at the center',
     >   /'    of the region of interest  (where most of the smoothing',
     >   /'    is sought).  Use a power of 3. or 4. to smooth out some',
     >   /'    NARROW region of noisy curvature,  else  a lesser power',
     >   /'    (1., 1.5, or 2.) if BROAD smoothing is sought (probably',
     >   /'    in conjunction with a table of target 2nd derivatives).',
     >   /)
 1003 FORMAT
     >  (/' REFINE: The T/C refinement iteration has not converged.',
     >   /' You are probably asking for too much.',
     >   /' Either try again, seeking less, or REFINE this result.')

 1005 FORMAT (/, A)
 1010 FORMAT (A, A)
 1020 FORMAT (A, F10.5)
 1030 FORMAT (/, 2A, I5, /)
 1040 FORMAT (A, I3, 4X, A, F10.5, 4X, A, F10.5)

      END SUBROUTINE REFINE
C+----------------------------------------------------------------------
C
      SUBROUTINE REFSRF (ISRF, NSRF, X, Y, YSCALE,
     >                   YBUMPPARAMS, YPPWEIGHT,
C****>                   A, C,
     >                   A, B, C, D, R, S,
     >                   NOYPPTABLE, NPTSYPPTABLE,
     >                   XYPPTABLE, YPPTABLE,
     >                   DETAILS, PROMPT, LEAVEALONE,
     >                   LUNCRT, LUNKBD, LUNTAB)
C
C  PURPOSE:  REFSRF was introduced so that REFINE can treat both surf-
C            aces without repeated code.  See REFINE for more details.
C
C  ARGUMENTS:
C  VAR   DIM   TYPE I/O/S DESCRIPTION
C  ISRF   -     I     I   ISRF=1 means upper surface; 2 means lower
C  NSRF   -     I     I   Number of points on this surface
C  X     NSRF   R     I   Surface abscissas
C  Y     NSRF   R    I/O  Surface ordinates (probably modified on return)
C YSCALE NSRF   R     S   For nonuniform scaling of the ordinates
C YBUMPPARAMS   3  R  I   Parameters for nonuniform y scaling bump
C YPPWEIGHT   NSRF R  S   For nonuniform weighting of y" equations
C**A   2*N,N+1  R     S   For setting up overdetermined system (N=NSRF-2)
C**C     2*N    R     S   For solution of overdetermined system
C  A      N     R     S   For sparse solution of overdetermined system
C  B     N+1    R     S    "    "    "    "    "    "    "    "    "
C  C      N     R     S    "    "    "    "    "    "    "    "    "
C  D      N     R     S    "    "    "    "    "    "    "    "    "
C  R      N     R     S    "    "    "    "    "    "    "    "    "
C  S      N     R     S    "    "    "    "    "    "    "    "    "
C  NOYPPTABLE   L     I   T means no y" table look-up for this surface
C  NPTSYPPTABLE I     I   Obvious input
C  XYPPTABLE, YPPTABLE    Abscissas and ordinates of y" table if present
C  DETAILS -    L     I   Controls detailed writes to LUNTAB
C  PROMPT  -    L    I/O  Since REFSRF is called for either surface and
C                         doesn't know about the outside iteration,
C                         control of the prompts common to both surfaces
C                         is a little awkward.  PROMPT needs to be TRUE
C                         for each surface for the first iteration, and
C                         turned off thereafter.  It also controls some
C                         of the detailed printing to LUNTAB (needed
C                         once for each surface, but not each iteration).
C  LEAVEALONE  L    I/O   Similarly, this logical must be FALSE for the
C                         first time in for each surface, but may be set
C                         TRUE for one of the surfaces at the prompt here.
C  LUNCRT -     I     I   Logical unit for screen
C  LUNKBD -     I     I   Logical unit for keyboard
C  LUNTAB -     I     I   Logical unit for tabulations
C
C  PROCEDURES:
C  BEVAL       For the bumps needed by the y scaling and y" weighting
C  DTDLSQ      Alternative to HDESOL, for special structure involved
C**HDESOL      Linear least squares of dense system by Householder
C**            transformations
C  READER      Prompting utility
C  TABLE1      1-D table look-up (linear interpolation)
C
C  HISTORY:
C  07/30/93   DAS   Introduced global SAVE after all, to be certain.
C  02/07/90   DAS   Removed concatenation from I/O; replaced ('$...')
C                   with (A,$) for IRIS 4D purposes.
C  03/19/86   DAS   Eliminated the END DO VAX dependency for PC compilers
C                   and arranged for suppression of details in .TAB file.
C                   Also avoided problems with SAVE and DATA statements
C                   by making DETAILS, PROMPT, and LEAVEALONE arguments.
C  11/05/84   DAS   Introduced DTDLSQ to avoid large 2-D array.
C                   Left dense-method solution in, but commented out.
C  10/25/84   DAS   Prompts tidied up in view of new on-line help.
C  04/11/84   DAS   Eliminated handling of non-normalized airfoil; put
C                   subroutine BEVAL in place of function BUMP.
C  11/26/83   DAS   Simplified set-up of 2nd derivative constraints  -
C                   prompted by difficulty when I=1.
C  10/28/83   DAS   Revised for integration with REFINE in PROFILE, as
C                   a permanent new capability.
C  10/26/83   DAS   First shot at doing table look-up for  target y"s.
C  10/20/83   DAS   Suppressed prompts for case where a surface is re-
C                   quired to be left unperturbed.
C  10/12/83   DAS   y" weightings and mid-section constant constraints
C                   may be different on each surface now.
C  10/07/83   DAS   Added ISRF to permit different things done to each
C                   surface; enabled bounding of y" in some x/c region
C                   on the upper surface only.
C  09/30/83   DAS   Added saving of plottable curvature data.
C  09/27/83   DAS   Modularized as REFSRF, for treating both surfaces.
C  09/26/83   DAS   Nonuniform y-scaling (modified sine bump, suggest-
C                   by Rob Kennelly).  Both surfaces treated same way.
C  09/23/83   DAS   y"  before and after from finite differencing now,
C                   for greater consistency.
C  09/20/83   DAS   Initial design/coding of special version of MODIFY
C                   (y" from spline; upper surface only).
C
C  AUTHOR: David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

      SAVE

C ... Arguments:

      INTEGER
     >   ISRF, LUNCRT, LUNKBD, LUNTAB, NPTSYPPTABLE, NSRF

C*****REAL A (2*(NSRF-2), NSRF-1), C (2*(NSRF-2)),

      REAL
     >   A (NSRF-2), B (NSRF-1), C (NSRF-2), D (NSRF-2), R (NSRF-2),
     >   S (NSRF-2), YSCALE (NSRF), X (NSRF), Y (NSRF),
     >   YPPWEIGHT (NSRF), XYPPTABLE (*), YPPTABLE (*), YBUMPPARAMS (3)

      LOGICAL
     >   DETAILS, LEAVEALONE, NOYPPTABLE, PROMPT

C ... Local variables:

      INTEGER
     >   I, IER, INDEX, IOS, J, N

      REAL
     >   INNERWEIGHT (2), OUTERWEIGHT (2), SBUMPPARAMS (3),
     >   XATMAXWEIGHT (2), XHI (2), XLO (2), YPPBUMPWIDTH (2),
     >   YPPCONST (2), CI, CIM1, CIP1, DXI, DXIM1, DXSUM, SSQMIN,
     >   WEIGHT, YPPINTERP

      LOGICAL
     >   DEFAULT, NOYPPCONSTRAINT (2), QUIT

      CHARACTER
     >   SRF (2) * 15

C ... Procedures:

      REAL
     >   TABLE1

      EXTERNAL
     >   BEVAL, DTDLSQ, READR, READY, TABLE1

C ... Storage:

      DATA
     >   SRF / ' upper surface ', ' lower surface ' /


C ... Execution:

      IF (LEAVEALONE) GO TO 999

      IF (PROMPT) THEN

         IF (NOYPPTABLE) THEN
            WRITE (LUNCRT, 1010)
     >         ' No target y" values have been read for the',
     >         SRF (ISRF), '-'
            CALL READY (LUNCRT, 'leave the' // SRF (ISRF) //
     >         'unchanged? (Y/N; <CR>=No): ',
     >         LUNKBD, LEAVEALONE, DEFAULT, QUIT)
         END IF

C ...    <Default is FALSE from calling program.>

         IF (LEAVEALONE) THEN
            WRITE (LUNTAB, 1010) ' The', SRF (ISRF), 'is unchanged.'
            PROMPT = ISRF == 1
            GO TO 999
         ELSE
            WRITE (LUNTAB, 1010) ' Controls for', SRF (ISRF), 'follow:'
            WRITE (LUNCRT, 1010) '    Smoothing of', SRF (ISRF),
     >                           'is controlled by the following:'

            IF (NOYPPTABLE) THEN
               CALL READR (LUNCRT,
     > '   Enter y" value sought over some range. <CR>=no constraint: ',
     >            LUNKBD, YPPCONST (ISRF), NOYPPCONSTRAINT (ISRF), QUIT)

               IF (.NOT. NOYPPCONSTRAINT (ISRF)) THEN
C ...             Permit simple flattening of curvature in some interval:

  200             WRITE (LUNCRT, '(A)', ADVANCE = 'NO')
     >              '    Corresp. x/c-range? (2 values):               '
                  READ  (LUNKBD, *, IOSTAT = IOS) XLO (ISRF), XHI (ISRF)
                  IF (IOS /= 0) THEN
                     WRITE (LUNCRT, '(A)')
                     GO TO 200
                  END IF
                  WRITE (LUNTAB, 1020) 'y" constraint and x range:',
     >               YPPCONST (ISRF), XLO (ISRF), XHI (ISRF)
               END IF
            END IF

            XATMAXWEIGHT (ISRF) = 0.5E+0
            CALL READR (LUNCRT,
     >         '   Center of smoothed region? <CR> gives x/c=0.5: ',
     >         LUNKBD, XATMAXWEIGHT (ISRF), DEFAULT, QUIT)

            YPPBUMPWIDTH (ISRF) = 3.E+0
            CALL READR (LUNCRT,
     >         '   Power of sine for y" weights?  <CR> gives 3.0: ',
     >         LUNKBD, YPPBUMPWIDTH (ISRF), DEFAULT, QUIT)

            WRITE (LUNTAB, 1020) ' x/c at peak y" weight    :',
     >                           XATMAXWEIGHT (ISRF)
            WRITE (LUNTAB, 1020) ' Sine power for y" weights:',
     >                           YPPBUMPWIDTH (ISRF)

            OUTERWEIGHT (ISRF) = 0.004E+0
            CALL READR (LUNCRT,
     >         '   Weight at x/c=0,1 for y"? <CR> gives .004: ',
     >         LUNKBD, OUTERWEIGHT (ISRF), DEFAULT, QUIT)
            WRITE (LUNTAB, 1020) ' y" weight at x/c=0 and 1 :',
     >                           OUTERWEIGHT (ISRF)

            INNERWEIGHT (ISRF) = 0.04E+0
            CALL READR (LUNCRT,
     >         '   Peak weight for y"?       <CR> gives .040: ',
     >         LUNKBD, INNERWEIGHT (ISRF), DEFAULT, QUIT)
            WRITE (LUNTAB, 1020) ' Peak y" weight           :',
     >                           INNERWEIGHT (ISRF)
            WRITE (LUNCRT, 1005)
         END IF
      END IF

C ... Set up scale factors for the ordinates (nonlinear thinning):

      DO I = 1, NSRF
         YSCALE (I) = 1.E+0
      END DO

      CALL BEVAL ('SINE', 3, YBUMPPARAMS, .TRUE., NSRF, X, YSCALE)

      IF (DETAILS) THEN
         WRITE (LUNTAB, 1005)
     >      ' ', ' y scale factors sought at each x:', ' '
         WRITE (LUNTAB, 1030) (X (I), YSCALE (I), I = 1, NSRF)
      END IF

C ... Generate the y" weights for this surface:

      SBUMPPARAMS (1) = XATMAXWEIGHT (ISRF)
      SBUMPPARAMS (2) = YPPBUMPWIDTH (ISRF)
      SBUMPPARAMS (3) = INNERWEIGHT (ISRF) - OUTERWEIGHT (ISRF)

      DO I = 1, NSRF
         YPPWEIGHT (I) = OUTERWEIGHT (ISRF)
      END DO

      CALL BEVAL ('SINE', 3, SBUMPPARAMS, .TRUE., NSRF, X, YPPWEIGHT)

      IF (PROMPT) THEN
         PROMPT = ISRF == 1

         IF (DETAILS) THEN
            WRITE (LUNTAB, 1010) ' Weights used for y" at each',
     >                           SRF (ISRF), 'x/c:', ' '
            WRITE (LUNTAB, 1030) (X (I), YPPWEIGHT (I), I = 1, NSRF)
         END IF
      END IF

C ... Set up overdetermined system.
C     The original solution as a dense system is left in comment form
C     as an aid to understanding the solution by a specialized least
C     squares algorithm.
C     First, zero out the system (treated as dense):
C
C     N = NSRF - 2
C
C     DO J = 1, N+1
C        DO I = 1, 2*N
C           A (I, J) = 0.E+0
C        END DO
C     END DO
C
C ... Next, the simple scaling part:
C
C     DO I = 1, N
C        A (I, I)   = 1.E+0
C        A (I, N+1) = Y (I+1)*YSCALE (I+1)
C     END DO
C
C ... Next, the 2nd derivative "constraints":
C
C     DO I = 1, N
C        DXI   = X (I+2) - X (I+1)
C        DXIM1 = X (I+1) - X (I)
C        DXSUM = X (I+2) - X (I)
C        WEIGHT= YPPWEIGHT (I+1)*2.E+0
C        CIM1  = WEIGHT / (DXIM1*DXSUM)
C        CI    =-WEIGHT / (DXIM1*DXI)
C        CIP1  = WEIGHT / (DXI*DXSUM)
C        IF (I>1) A (I+N, I-1) = CIM1
C        A (I+N, I)   = CI
C        A (I+N, I+1) = CIP1
C        A (I+N, N+1) = CIM1*Y (I) + CI*Y (I+1) + CIP1*Y (I+2)
C
C        IF (NOYPPTABLE) THEN
C           IF (.NOT. NOYPPCONSTRAINT (ISRF)) THEN
C ...          Flatten the curvature plot in given x interval:
C
C              IF (X (I+1)>XLO (ISRF) .AND. X (I+1)<XHI (ISRF))
C    >            A (I+N, N+1) = 0.5E+0*WEIGHT*YPPCONST (ISRF)
C           ELSE
C ...          Don't constrain the y" target - shoot for original.
C              CONTINUE
C           END IF
C        ELSE
C           INDEX = 1
C           YPPINTERP = TABLE1 (NPTSYPPTABLE, XYPPTABLE,
C    >                          YPPTABLE, INDEX, X (I+1), IER)
C           IF (IER==0) THEN
C              A (I+N, N+1) = 0.5E+0*WEIGHT*YPPINTERP
C           ELSE IF (IER<=4) THEN
C              WRITE (LUNCRT, *) 'Table look-up error. IER:', IER
C              WRITE (LUNCRT, *) 'NTABLE:', NPTSYPPTABLE,
C    >                           '  ISRF:', ISRF
C              WRITE (LUNCRT, *) 'Abscissa:', X (I+1)
C              STOP
C           ELSE
C ...          Abscissa was out of table range - use original y".
C              CONTINUE
C           END IF
C        END IF
C     END DO
C
C ... Adjust RHS for known leading and trailing edge values:
C
C     Y (1)    = YSCALE (1)*Y (1)
C     Y (NSRF) = YSCALE (NSRF)*Y (NSRF)
C     A (1+N, N+1) = A (1+N, N+1) - 2.E+0*YPPWEIGHT (2)*Y (1) /
C    >                        ( (X (2)-X (1))* (X (3)-X (1)))
C     A (N+N, N+1) = A (N+N, N+1) - 2.E+0*YPPWEIGHT (N)*Y (NSRF) /
C    >                                               (DXI*DXSUM)
C
C     WRITE (LUNTAB, *)
C     WRITE (LUNTAB, *)
C    >   'RHS vector, of length 2N, where N = NSRF-2 =', N, ':'
C     WRITE (LUNTAB, *)
C     WRITE (LUNTAB, ' (1X, 1P, 10E13.4)') (A (I, N+1), I = 1, N)
C     WRITE (LUNTAB, ' (1X, 1P, 10E13.4)') (A (I, N+1), I = N+1, 2*N)
C     WRITE (LUNTAB, *)
C
C ... Solve the system:
C
C     CALL HDESOL (2*N, 2*N, N+1, A, C, SSQMIN)
C
C ... Extract the solution:
C
C     DO I = 1, N
C        Y (I+1) = C (I)
C     END DO


C ... Set up overdetermined system.
C     DTDLSQ was written specially for this diagonal+tridiagonal structure.

      N = NSRF - 2

C ... First, the simple scaling part:

      DO I = 1, N
         D (I) = 1.E+0
         R (I) = Y (I+1) * YSCALE (I+1)
      END DO

C ... Next, the 2nd derivative "constraints":

      DO I = 1, N

         DXI   = X (I+2) - X (I+1)
         DXIM1 = X (I+1) - X (I)
         DXSUM = X (I+2) - X (I)
         WEIGHT= YPPWEIGHT (I+1) * 2.E+0
         CIM1  = WEIGHT / (DXIM1*DXSUM)
         CI    =-WEIGHT / (DXIM1*DXI)
         CIP1  = WEIGHT / (DXI*DXSUM)
         C (I)  = CIM1
         A (I)  = CI
         B (I)  = CIP1
         S (I)  = CIM1 * Y (I) + CI * Y (I+1) + CIP1 * Y (I+2)

         IF (NOYPPTABLE) THEN
            IF (.NOT. NOYPPCONSTRAINT (ISRF)) THEN
C ...          Flatten the curvature plot in given x interval:

               IF (X (I+1)>XLO (ISRF) .AND. X (I+1)<XHI (ISRF))
     >            S (I) = 0.5E+0 * WEIGHT * YPPCONST (ISRF)
            ELSE
C ...          Don't constrain the y" target - shoot for original.
               CONTINUE
            END IF
         ELSE
            INDEX = 1
            YPPINTERP = TABLE1 (NPTSYPPTABLE, XYPPTABLE,
     >                          YPPTABLE, INDEX, X (I+1), IER)
            IF (IER == 0) THEN
               S (I) = 0.5E+0 * WEIGHT * YPPINTERP
            ELSE IF (IER <= 4) THEN
               WRITE (LUNCRT, 1005)
               WRITE (LUNCRT, 1015) ' Table look-up error. IER:', IER,
     >                              ' NTABLE:', NPTSYPPTABLE,
     >                              ' ISRF:', ISRF
               WRITE (LUNCRT, 1020) ' Abscissa:', X (I+1)
               STOP
            ELSE
C ...          Abscissa was out of table range - use original y".
               CONTINUE
            END IF
         END IF

      END DO

C ... Adjust RHS for known leading and trailing edge values:

      Y (1)    = YSCALE (1) * Y (1)
      Y (NSRF) = YSCALE (NSRF) * Y (NSRF)
      S (1) = S (1) - C (1) * Y (1)
      S (N) = S (N) - B (N) * Y (NSRF)

      IF (DETAILS) THEN
         WRITE (LUNTAB, 1005) ' ', ' RHS vector, length 2*(N-2):', ' '
         WRITE (LUNTAB, 1040) R (1 : N)
         WRITE (LUNTAB, 1040) S (1 : N)
      END IF

C ... Solve the system:

      CALL DTDLSQ (N, A, B, C, D, R, S, SSQMIN)

C ... Extract the solution:

      DO I = 1, N
         Y (I+1) = R (I)
      END DO

      WRITE (LUNTAB, 1025) ' Sum of squares obtained  :', SSQMIN
C
  999 RETURN

C ... Formats:

 1005 FORMAT (A)
 1010 FORMAT (/, A, A, A)
 1015 FORMAT (A, I5)
 1020 FORMAT (A, 3F15.6)
 1025 FORMAT (/, A, 1P, E15.6)
 1030 FORMAT (5 (1X, F12.4, F10.5))
 1040 FORMAT (1X, 1P, 10E13.4)

      END SUBROUTINE REFSRF
C+---------------------------------------------------------------------
C
      SUBROUTINE ROTATE (NU, XU, YU, NL, XL, YL, LUNCRT, LUNKBD)
C
C  PURPOSE: ROTATE applies twist to one profile about a variable center
C           of rotation.  The result may be optionally renormalized.
C
C  ARGUMENTS:
C   ARG     DIM   TYPE  I/O/S   DESCRIPTION
C   NU       -      I     I     Upper surface info.
C   XU       NU     R    I/O      "     "     "
C   YU       NU     R    I/O      "     "     "
C   NL       -      I     I     Lower surface info.
C   XL       NL     R     I       "     "     "
C   YL      NPTS    R    I/O      "     "     "
C  LUNCRT    -      I     I     Logical unit for prompting.
C  LUNKBD    -      I     I     Logical unit for responses.
C
C  PROCEDURES:
C    NRMLIZ    Normalizes X or Y for one surface of an airfoil
C    RDREALS   Prompts for more than one value
C    READER    Prompting utility
C    ROTATE2D  Rotates point(s) (X,Y) about (P,Q)
C
C  AUTHOR:     David Saunders, Sterling Software, Palo Alto, CA.
C
C  HISTORY:
C  06/29/90   DAS   Initial implementation.
C  10/21/99    "    Positive rotation should be clockwise (like AOA).
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, NL, NU

      REAL
     >   XL(NL), XU(NU), YL(NL), YU(NU)

C  *  Local constants:

      REAL, PARAMETER ::
     >   FOURTH = 0.25E+0, HALF = 0.5E+0, ONE = 1.E+0, ZERO = 0.E+0

C  *  Local variables:

      INTEGER
     >   CHOICE, NVALS

      REAL
     >   CHORD, TWIST, VALS (2), XTE, YTE

      LOGICAL
     >   DEFAULT, QUIT, RENORM

C  *  Procedures:

      EXTERNAL
     >   NRMLIZ, RDREALS, READI, READR, READY, ROTATE2D

C  *  Execution:

      WRITE (LUNCRT, 1001)
     >   ' Options for applying twist:',
     >   '    1 = Rotate about leading edge',
     >   '    2 = Rotate about trailing edge',
     >   '    3 = Rotate about quarter-chord (on center-line)',
     >   '    4 = Rotate about specified point'

      CHOICE = 1
      CALL READI (LUNCRT, 'Enter selection (<CR> = 1): ',
     >            LUNKBD, CHOICE, DEFAULT, QUIT)

      CHORD = MAX (XU (NU), XL (NL)) - XU (1)
      XTE   = (XU (NU) + XL (NL)) * HALF
      YTE   = (YU (NU) + YL (NL)) * HALF

      IF (CHOICE == 1) THEN
         VALS (1) = XU (1)
         VALS (2) = YU (1)
      ELSE IF (CHOICE == 2) THEN
         VALS (1) = XTE
         VALS (2) = YTE
      ELSE IF (CHOICE == 3) THEN
         VALS (1) = XU (1) + (XTE - XU (1)) * FOURTH
         VALS (2) = YU (1) + (YTE - YU (1)) * FOURTH
      ELSE IF (CHOICE == 4) THEN
  200    NVALS = 2
         CALL RDREALS (LUNCRT, '$Enter center of rotation (2 values): ',
     >      LUNKBD, NVALS, VALS)
         IF (NVALS /= 2) GO TO 200
      END IF

  300 CALL READR (LUNCRT,
     >   'Enter twist (degrees; positive is clockwise): ',
     >   LUNKBD, TWIST, DEFAULT, QUIT)
      IF (DEFAULT) GO TO 300

C  *  Provide for renormalizing (but in general this does not recover
C     the original abscissas):

      IF (XU (1) == ZERO .AND. YU (1) == ZERO .AND. CHORD == ONE) THEN
         RENORM = .TRUE.
         CALL READY (LUNCRT,
     >      'Do you want the result renormalized?  (Y/N; <CR> = Yes): ',
     >       LUNKBD, RENORM, DEFAULT, QUIT)
      ELSE
         RENORM = .FALSE.
      END IF

C  *  Perform the rotation, in-place:

      TWIST = -TWIST ! ROTATE2D thinks in quadrants

      CALL ROTATE2D (NU, XU, YU, TWIST, VALS (1), VALS (2))
      CALL ROTATE2D (NL, XL, YL, TWIST, VALS (1), VALS (2))

C  *  Warn user if leading edge point has been displaced:

      IF (XU (2) < XU (1) .OR. XL (2) < XL (1)) THEN
         WRITE (LUNCRT, 1001)
     >   ' *** WARNING: LE point is no longer the true leading edge.',
     >   '              PROFILE''s "rectify" option may be appropriate.'
         IF (RENORM) WRITE (LUNCRT, 1000)
     >   '              Normalization has been affected also.'
      END IF

      IF (RENORM) THEN
         CHORD = MAX (XU (NU), XL (NL)) - XU (1) ! Not if warning above applied.
         CALL NRMLIZ (NU, XU, XU, XU (1), CHORD)
         CALL NRMLIZ (NU, YU, YU, YU (1), CHORD)
         CALL NRMLIZ (NL, XL, XL, XL (1), CHORD)
         CALL NRMLIZ (NL, YL, YL, YL (1), CHORD)
         WRITE (LUNCRT, 1001)
     >    ' *** Note that the original Xs are NOT retained in general.',
     >    '     Rerun PROFILE in "redistribute" mode to recover them.'
      END IF

 1000 FORMAT (A)
 1001 FORMAT (/, (A))

      END SUBROUTINE ROTATE
C+----------------------------------------------------------------------
C
      SUBROUTINE SMOOTH (NU, XU, YU, NL, XL, YL, LUNCRT, LUNKBD, LUNOUT)
C
C  PURPOSE:  SMOOTH drives smoothing of an airfoil (either surface or
C            both surfaces) by a choice of methods.  The Y coordinates
C            are smoothed in-place.
C
C  ARGUMENTS:
C     VAR   DIM   I/O/S   DESCRIPTION
C     NU     -    I       # points on upper surface
C     XU     NU   I       Upper surface coordinates, leading to trailing
C     YU     NU   I/O     edge (not necessarily normalized)
C     NL     -    I       # points on lower surface
C     XL     NL   I       Lower surface coordinates as for upper
C     YL     NL   I/O
C     LUNCRT -    I       Logical units for screen, keyboard, and
C     LUNKBD -    I       printed coefficients
C     LUNOUT -    I
C
C  PROCEDURES:
C     WAGSMOOTH  Smooths an airfoil surface using Wagner functions
C     IMPSMOOTH  Smooths an airfoil surface via implicit + explicit methods
C
C  HISTORY:
C     08/11/86  DAS  Initial Wagner fn. form (from program SMOOTH).
C     12/20/96  DAS  Added implicit/explicit smoothing option.
C
C  AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C
C-------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER  LUNCRT, LUNKBD, LUNOUT, NU, NL
      REAL     XU (NU), XL (NL), YU (NU), YL (NL)

C     Local variables:

      INTEGER  METHOD
      LOGICAL  CR, EOF

C     Procedures:

      EXTERNAL IMPSMOOTH, READI, WAGSMOOTH

C     Execution:

      WRITE (LUNCRT, '(/, (1X, A))')
     >   'Smoothing methods:',
     >   '   1:  Least squares fitting of Wagner functions 1:N',
     >   '   2:  Implicit + explicit smoothing, variable along chord'

  10  METHOD = 1
      CALL READI (LUNCRT, 'Pick one. Default = 1: ',
     >            LUNKBD, METHOD, CR, EOF)
      IF (EOF) GO TO 99

      IF (METHOD == 1) THEN  ! Wagner function smoothing

         CALL WAGSMOOTH (NU, XU, YU, 1, LUNCRT, LUNKBD, LUNOUT)
         CALL WAGSMOOTH (NL, XL, YL, 2, LUNCRT, LUNKBD, LUNOUT)

      ELSE IF (METHOD == 2) THEN  ! Implicit (+ explicit) smoothing

         CALL IMPSMOOTH (NU, XU, YU, 1, LUNCRT, LUNKBD, LUNOUT)
         CALL IMPSMOOTH (NL, XL, YL, 2, LUNCRT, LUNKBD, LUNOUT)

      ELSE
         GO TO 10
      END IF

   99 RETURN

      END SUBROUTINE SMOOTH
C+---------------------------------------------------------------------
C
      SUBROUTINE TRANSFORM (NU, XU, YU, NL, XL, YL, LUNCRT, LUNKBD,
     >                      ICASE)
C
C  PURPOSE: TRANSFORM  drives transformation of an airfoil from  upper/
C           lower surface representation to thickness/camber represent-
C           ation, or vice versa.   It removes from the calling program
C           such necessities as prompting for which way to go, and dis-
C           allowing it if upper/lower abscissas are not the same, etc.
C           This version provides for zeroing the camber of a section.
C
C  METHOD:  See XFORM for the essentials - separated out in case  it is
C           reusable elsewhere.  Note that ICASE is an OUTput, and that
C           the transformation is done in-place.
C
C           Checking for dissimilar abscissas might save somebody  some
C           grief - quitting here is better than not bothering to check.
C
C  ARGUMENTS:
C   ARG     DIM   TYPE  I/O/S   DESCRIPTION
C   NU       -      I     I     ICASE = 0 or 1: no. of upper surface pts.
C   XU       NU     R     I     and corresp. abscissas, else same as NL,XL.
C   YU       NU     R    I/O    ICASE = 0: upper Ys in/decambered Ys out;
C                               ICASE = 1: upper Ys in/camber out;
C                               ICASE = 2: camber in/upper Ys out.
C   NL       -      I     I     ICASE = 2: no. of values in camber/thickness
C   XL       NL     R     I                distribns. and corresp. abscissas
C                                          else # lower surface pts. and Xs.
C   YL      NPTS    R    I/O    ICASE = 0: lower Ys in/decambered Ys out;
C                               ICASE = 1: lower Ys in/semi-thickness out;
C                               ICASE = 2: semi-thickness in; lower Ys out.
C  LUNCRT    -      I     I     Logical unit for prompting.
C  LUNKBD    -      I     I     Logical unit for responses.
C  ICASE     -      I     O     Lets user say which way to go, and tells the
C                               calling program for tabulation purposes, etc.
C                               ICASE = 0 means zero out the camber;
C                                     = 1 means Y-coordinates in, and
C                                         camber/semi-thickness out;
C                                     = 2 means transform the other way.
C
C  PROCEDURES:
C    READER   Prompting utility.
C    XFORM    Does the actual transformation (in same source module as
C             TRANSFORM because the two go hand in hand, but XFORM may
C             be useful on its own some day).
C
C  HISTORY:
C  02/12/85   DAS   Initial design and code.
C  08/31/87   DAS   'A' format instead of list-directed.
C  11/21/95   DAS   Installed a decambering option.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C----------------------------------------------------------------------
C
      IMPLICIT NONE

C  *  Arguments:

      INTEGER  LUNCRT, LUNKBD, NL, NU, ICASE
      REAL     XL (NL), XU (NU), YL (NL), YU (NU)

C  *  Local variables:

      INTEGER  I
      REAL     YLE
      LOGICAL  DEFAULT, QUIT, SAMEXS, YIN

C  *  Execution:

      WRITE (LUNCRT, 1001)
     >   ' Options for TRANSFORM mode:',
     >   '  0 = Decamber upper/lower surface coordinates',
     >   '  1 = Upper/lower surfaces to camber/thickness distributions',
     >   '  2 = Camber/thickness distributions to upper/lower surfaces',
     >   ' '

      ICASE = 1
      CALL READI (LUNCRT, 'Enter selection (<CR> = 1): ',
     >            LUNKBD, ICASE, DEFAULT, QUIT)

      IF (NU /= NL) THEN
         WRITE (LUNCRT, 1001)
     >      ' Cannot transform - different-sized distributions.',
     >      ' Stopping here.'
         GO TO 900
      END IF

      IF (ICASE /= 2) THEN  ! Allow round-off X differences for ICASE = 2
         SAMEXS = .TRUE.
         DO I = 1, NU
            IF (XU (I) /= XL (I)) SAMEXS = .FALSE.
         END DO

         IF (.NOT. SAMEXS) THEN
            WRITE (LUNCRT, 1001)
     >         ' Cannot transform - different abscissa distributions.',
     >         ' Use REDISTRIBUTE mode and try again.  Stopping here.'
            GO TO 900
         END IF

C        Ensure zero camber at the nose:

         YLE = YU (1)
         IF (YLE /= 0.) THEN
            DO I = 1, NU
               YU (I) = YU (I) - YLE
               YL (I) = YL (I) - YLE
            END DO
         END IF
      END IF

      YIN = ICASE /= 2

      CALL XFORM (YIN, NU, YU, YL)

      IF (ICASE == 0) THEN  ! Zero the camber ...
         DO I = 1, NU
            YU (I) = 0.
         END DO

         CALL XFORM (.FALSE., NU, YU, YL)  ! ... and reverse the transformation
      END IF

      GO TO 999

  900 STOP

  999 RETURN

 1001 FORMAT (/, (A))

      END SUBROUTINE TRANSFORM
C+---------------------------------------------------------------------
C
      SUBROUTINE XFORM (YIN, NPTS, YU, YL)
C
C  PURPOSE: XFORM  transforms one  airfoil  from  upper/lower  surface
C           representation to thickness/camber representation, or vice
C           versa.    Actually, semi-thickness is used, not thickness,
C           while "camber" is more accurately the mean line here:
C
C                 C = (YU + YL)/2               YU = C + T
C                 T = (YU - YL)/2               YL = C - T
C
C  METHOD:  The required output overwrites the input arrays.  Note the
C           assumption of equal numbers of points on both surfaces and
C           matching abscissas.
C
C  ARGUMENTS:
C   ARG     DIM   TYPE  I/O/S   DESCRIPTION
C   YIN      -      L     I     YIN=T means y-coordinates in; camber/semi-
C                                     thickness out;
C                               YIN=F means transform the other way.
C   NPTS     -      I     I     No. of points in each of the distributions
C                               involved.
C   YU      NPTS    R    I/O    If YIN=T: Upper surface in; camber out.
C                               If YIN=F: Camber in; upper surface out.
C   YL      NPTS    R    I/O    If YIN=T: Lower surface in; semi-thickness out.
C                               If YIN=F: Semi-thickness in; lower surface out.
C
C  HISTORY:
C  02/12/85   DAS   Initial design and code.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      LOGICAL  YIN
      INTEGER  NPTS
      REAL     YU (NPTS), YL (NPTS)

C  *  Local variables:

      INTEGER  I
      REAL     FACTOR, TEMP

C  *  Execution:

      IF (YIN) THEN
         FACTOR = 0.5E+0
      ELSE
         FACTOR = 1.E+0
      END IF

      DO I = 1, NPTS
         TEMP   = (YU (I) + YL (I)) * FACTOR
         YL (I) = (YU (I) - YL (I)) * FACTOR
         YU (I) = TEMP
      END DO

      END SUBROUTINE XFORM
C+-----------------------------------------------------------------------
C
      SUBROUTINE WAGSMOOTH (NPTS, X, Y, ISRF, LUNCRT, LUNKBD, LUNOUT)
C
C  PURPOSE:  WAGSMOOTH smooths the given airfoil surface by fitting a
C            linear combination of the first N Wagner shape functions.
C            The airfoil need not be normalized.
C
C  METHOD:   Prompt for whether to smooth this surface. If so, ...
C
C         >  Prompt for number of Wagner functions to use.
C         >  Normalize coordinates: [0.,1.] for X; corresp. scaling for Y.
C         >  If YTE = Y(NPTS) is not zero, use BEVAL's 'RAMP' option to
C            subtract Y = YTE * X at every abscissa, in place.  The
C            Wagner functions will be fitted to the adjusted ordinates.
C         >  Perform the least squares fit.  (See WAGFIT.)
C         >  Evaluate the result, including adding back the ramp term
C            if necessary, all in-place.
C         >  Denormalize.
C
C            Note that local work-space is allocated for WAGFIT -
C            PROFILE's main program does not have enough available,
C            and WAGSMOOTH is somewhat self-contained as a result.
C
C  ARGUMENTS:
C  VAR   DIM   I/O/S   DESCRIPTION
C  NPTS   -    I       Number of points on current surface
C  X     NPTS  I       Abscissas of current surface   (not necessarily
C  Y     NPTS  I/O     Ordinates of current surface        normalized)
C  ISRF   -    I       1 means uppper surface;
C                      2 means lower surface
C  LUNCRT -    I       Logical units for screen, keyboard, and
C  LUNKBD -    I       printed coefficients
C  LUNOUT -    I
C
C  PROCEDURES:
C  READER   Prompting utility
C  BEVAL    Evaluates numerous shape functions at given Xs
C  NRMLIZ   Normalizes/denormalizes coordinates
C  WAGFIT   Sets up and solves the linear least squares problem
C
C  HISTORY:
C  08/11/86   DAS   Initial implementation (originally installed in
C                   program SMOOTH, but PROFILE is the logical place;
C                   expect X in [0.,1.], because PROFILE can normalize).
C  10/22/88   DAS   Introduced normalizing/denormalizing here - too many
C                   steps otherwise for arbitrary coordinates.
C  12/13/96   DAS   Added option to replace points near the trailing edge.
C
C  AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C
C-------------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   ISRF, LUNCRT, LUNKBD, LUNOUT, NPTS

       REAL
     >   X (NPTS), Y (NPTS)

C  *  Local constants:

      INTEGER, PARAMETER ::
     >   MXPTS  = 200,   ! Arbitrary - more pts./surface OK if NWAG < MXWAG
     >   MXWAG  = 24,
     >   MXWORK = (MXWAG + 2) * MXPTS

      REAL, PARAMETER ::
     >   ONE    = 1.E+0,
     >   ZERO   = 0.E+0

      CHARACTER, PARAMETER ::
     >   RAMP * 4 = 'RAMP'

C  *  Local variables:

      INTEGER
     >   I, IER, J, NFIX, NLEFT, NWAG

      REAL
     >   B (2), C (MXWAG), CHORD, RMSDEV, WORK (MXWORK), XLE, XTE, YLE,
     >   YN, YTE

      LOGICAL
     >   ADD, CR, EOF

      CHARACTER
     >   SURFCE (2) * 5

C  *  Procedures:

      EXTERNAL
     >   BEVAL, NRMLIZ, READI, WAGFIT

C  *  Storage:

      SAVE
     >   NWAG
      DATA
     >   SURFCE / 'UPPER', 'LOWER' /

C  *  Execution:

      IF (ISRF == 1) THEN
         NWAG = 10
         WRITE (LUNCRT, '(/, A)')
     >      ' Enter # Wagner functions to use on the UPPER surface.'
         CALL READI (LUNCRT,
     >      '(<CR> = 10; 0 or EOF (^D) = leave surface alone): ',
     >      LUNKBD, NWAG, CR, EOF)
      ELSE  ! ISRF = 2
         IF (NWAG == 0) THEN
            NWAG = 10
            CALL READI (LUNCRT,
     >         ' # Wagner functions to use on the LOWER surface? [10] ',
     >         LUNKBD, NWAG, CR, EOF)
         ELSE
            CALL READI (LUNCRT,
     >   '# Wagner fns. for the LOWER surface? (<CR> = same as upper) ',
     >         LUNKBD, NWAG, CR, EOF)
         END IF
      END IF

      IF (EOF .OR. NWAG <= 0) GO TO 899

      IF (NWAG > MXWAG) THEN
         NWAG = MXWAG
         WRITE (LUNCRT, '(A, I3)') ' Limiting N to maximum of', MXWAG
      END IF

  200 NFIX = 3
      CALL READI (LUNCRT,
     > '# pts. to fix (after smoothing) forward of the TE? (>=0; [3]) ',
     >   LUNKBD, NFIX, CR, EOF)
      IF (EOF) GO TO 899
      IF (NFIX < 0 .OR. NFIX > NPTS / 4) GO TO 200

C  *  Normalize coordinates, since Wagner functions are defined on [0,1]:

      XLE = X (1)
      YLE = Y (1)
      XTE = X (NPTS)
      YTE = Y (NPTS)
      CHORD = XTE - XLE

      CALL NRMLIZ (NPTS, X, X, XLE, CHORD)
      CALL NRMLIZ (NPTS, Y, Y, YLE, CHORD)

C  *  To preserve a thick trailing edge exactly, a "ramp" function
C     must be subtracted first then added back in:

      YN = -Y (NPTS)

      IF (YN /= ZERO) THEN

C  *     Subtract "ramp" function in-place:

         ADD = .TRUE.
         CALL BEVAL (RAMP, 1, [YN], ADD, NPTS, X, Y)
      END IF

C  *  Set up and solve the linear least squares problem:

      CALL WAGFIT (NPTS, X, Y, NWAG, MXWORK, WORK, C (1),
     >             RMSDEV, IER)
      IF (IER /= 0) GO TO 910

      WRITE (LUNOUT, 1120)
     >   SURFCE (ISRF), NWAG, NFIX, RMSDEV, (C (I), I = 1, NWAG)

C  *  Evaluate the fit (in-place):

      ADD = .FALSE.
      DO J = 1, NWAG
         B (1) = J
         B (2) = C (J)

         CALL BEVAL ('WAGNER', 2, B, ADD, NPTS, X, Y)
         ADD = .TRUE.
      END DO

      IF (YN /= ZERO) THEN

C  *     Add back the "ramp" function:

         YN = -YN
         CALL BEVAL (RAMP, 1, [YN], ADD, NPTS, X, Y)
      END IF

C  *  Patch bad points near the trailing edge?

      IF (NFIX > 0) THEN

         DO I = 1, NPTS
            WORK (I) = X (I)
            WORK (I + NPTS) = Y (I)
         END DO

         NLEFT = NPTS - NFIX
         WORK (NLEFT) = WORK (NPTS)  ! Remove NFIX pts
         WORK (NPTS + NLEFT) = WORK (NPTS + NPTS)

         CALL LCSFIT (NLEFT, WORK (1), WORK (1 + NPTS), .TRUE., 'M',
     >                NFIX, X (NLEFT), Y (NLEFT), Y (NLEFT))
      END IF

C  *  Denormalize:

      CALL NRMLIZ (NPTS, X, X, XLE, -CHORD)
      CALL NRMLIZ (NPTS, Y, Y, YLE, -CHORD)

      X (1)    = XLE
      X (NPTS) = XTE
      Y (1)    = YLE
      Y (NPTS) = YTE

  899 RETURN

C  *  Error handling:

  910 WRITE (LUNCRT, 1110) IER
      STOP

C  *  Formats:

 1110 FORMAT (' WAGSMOOTH: Bad return from WAGFIT - abort. IER: ', I1,
     >        /)
 1120 FORMAT (//, ' Details of smoothing for ', A5, ' surface:',
     >        //, ' # Wagner functions:', I3,
     >        /,  ' # patched points forward of TE:', I3,
     >        /,  ' RMS Deviation:', 1P, E14.6,
     >        /,  ' Coefficients:',
     >        /,  (5E14.6))

      END SUBROUTINE WAGSMOOTH

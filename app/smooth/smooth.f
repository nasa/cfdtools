C+------------------------------------------------------------------------------
C
      PROGRAM SMOOTH
C
C  PURPOSE:
C
C     SMOOTH drives a variety of 1-D smoothing/interpolating algorithms,
C     either for evaluating or comparing them, or for generating quanti-
C     ties that may be used elsewhere.   It applies the chosen method(s)
C     to one or more "curves" read from a single file, and saves results
C     in a form suited to plotting by some other program.   On VAXes and
C     IRISes at NASA Ames,  plottable results are tailored to the  QPLOT
C     package.   A Macintosh derivative known as QuickFit is also avail-
C     able for which QuickPlot is appropriate. Results may also be saved
C     in the format used for the input data.
C
C     This version permits repeating the above with further input files.
C     Plots of curves from any one input file will appear on one frame -
C     a new plot frame for a new input file.
C
C     This version also allows selection of X and Y from multiple-column
C     files (same columns for all curves assumed),  and two data formats
C     are supported (the original "NPTS" form and an "indefinite" form).
C
C     Input data file format:
C
C        [TITLE]         <Optional title, up to 80 characters.>
C        [N1]            <# pts. in 1st curve: not needed if "END"
C        X (1)  Y (1)        is present; both are redundant if
C        X (2)  Y (2)        there is only one curve.>
C        X (3)  Y (3)
C        X (4)  Y (4)    <X and Y may be extracted from specified
C        :      :            columns if there are more than two.>
C        :      :
C        # X     Y       <"#" suppresses points or adds comments.>
C        X (N1) Y (N1)
C                        <Blank lines are ignored.>
C        N2  or  END     <Repeat for more curves, omitting title.>
C        X (1)  Y (1)
C        :      :
C        :      :
C
C     Normally,  uniformly-spaced abscissas are generated for evaluating
C     fitted curves in order to plot them, but these may be suppressed.
C
C     Evaluating fitted curves at further abscissas read from  (column 1
C     of) a file is also provided for, with these abscissas being in the
C     same format as above - any further columns are ignored.  
C
C     The original data,  the fitted curve,  and any further evaluations
C     are written to a file of plottable results, named "smooth.plt." In
C     the case of least squares methods, plots of deviations may also be
C     obtained via a separate plot file, "smooth.dev."   (Its scaling is
C     unlikely to be the same as for the original data.)  This file con-
C     tains  X (data), Y (fit), and Y (fit) - Y (data)  in three columns
C     with columns 1 and 3 set ready to plot (but column 2 there in case
C     it is needed).
C
C     Evaluated results may also be saved in "SMOOTH" format as the file
C     "smooth.out" if requested.  (Later: this file name is now prompted
C     for as the first prompt,  which is unintuitive during most  likely
C     uses,  but an artifact of allowing for reading more than one input
C     dataset.)
C
C     Any coefficients computed by a smoothing/fitting method are  saved
C     separately in "smooth.log."
C
C
C  METHOD:  <Outline>
C
C     Open files for plottable and/or reusable results.
C
C     NEW FILE:
C        Prompt for and open an input data file.
C        Read the input data title.
C        Prompt for whether to normalize the data to [0,1].
C        Prompt for plot title and axis labels, and initialize plot frame.
C        Prompt for the number of uniformly-spaced abscissas to be used in
C           evaluating and plotting each fitted curve (may be suppressed).
C        Prompt for the filename of optional user-supplied abscissas for
C           an additional fit evaluation and plot, and for possible reuse.
C        If <user-supplied abscissas> open and read this file.
C        Initialize the tabulated output.
C
C     NEW CURVE FIT:
C        Rewind data file.
C        Select a smoothing routine (else GO TO NEXT).
C        Skip the data file title.
C        Read NX, number of data points in the first dataset (curve).
C
C        WHILE <not EOF> DO
C           Read one dataset (curve).
C           If <first selection>
C              Write plot info for this dataset.
C              If <normalizing> normalize X (*) and Y (*) to [0, 1].
C           If <fit to be plotted> generate uniformly-spaced abscissas
C              in the data range.
C           Prompt for smoothing parameters (if applicable).
C           Perform requested curve fit.
C           If <fit to be plotted> evaluate curve at generated abscissas.
C           If <user-supplied abscissas> evaluate curve at these.
C           If <normalizing> denormalize results.
C           Read NX for next curve.
C        END WHILE
C
C        GO TO NEW CURVE FIT
C
C     NEXT:
C
C        Close input data file for possible reuse of logical unit.
C        If <new file desired> then
C           Output "end-frame" to plottable data file.
C           GO TO NEW FILE
C        Else
C           Close plottable and printable files.
C           Done.
C
C  SIGNIFICANT CONSTANTS:
C
C    MXFCOF   Limit on size of partial Fourier series defined by N.
C    MXFSFIT  Limit on no. of coefs. for FSFIT (expensive in memory).
C    MXLTYP   Limit on no. of line types known to the plot package.
C    MXMENU   Limit on no. of smoothing options (including quit).
C    MXNDEG   Limit on degree of simple polynomial to fit.
C    MXNXK    Limit on no. of knots which ICSFKU can handle.
C    MXPLOT   Limit on no. of points for plotting smoothed curve.
C    MXPTS    Limit on no. of points per input dataset.
C             Actually, MXPTS-1 is the limit - allows for inserting
C             an extra point for closed curves (PSFIT).
C    MXVEC    Limit on no. of discrete vectors handled (VECFIT).
C    MXWAG    Limit on no. of Wagner functions permitted (WAGFIT).
C
C  SIGNIFICANT LOCAL VARIABLES:
C
C    NAME     DIM       DESCRIPTION
C    A      0:MXFCOF    Fourier series coefficients.
C    B      1:MXFCOF    Fourier series coefficients.
C    C      0:MXWAG     1:N+1 used for first N Wagner function coefs
C                             plus ramp function;
C                       0:N   used for polynomial coefs. by PNFIT;
C                       1:N   used for coefs. of N vectors by VECFIT.
C    COEFS   MXPTS,3    Spline coefficients.
C    LEGEND     -       Legend for plot file.
C    TTLDAT     -       Title from dataset.
C    TTLPLT     -       Title for plot file.
C    XLABEL,YLABEL      Axis labels for plot file.
C    WGHTS    MXPTS     Weights needed by spline routine ICSSCU.
C    X        MXPTS     One set of data abscissas.
C    X0       MXPTS     Copy of data abscissas in case of normalization.
C    XPLOT    MXPLOT    Equal-spaced abscissas for plotting smoothed curves.
C    XEVAL    MXPTS     User-supplied abscissas for evaluating   "     "   .
C    XEVAL0     "       Copy of user-supplied abscissas (in case of nrmlzn.)
C    XNORM              Data abscissas normalized to [0,2*pi].
C    XK       MXNXK     Knot locations for ICSFKU, ICSVKU.
C    Y        MXPTS     One set of data ordinates.
C    YDEV     MXPTS     Values of Y (*) - YEVAL (*).
C    YEVAL    MXPTS     Values of smoothed curves at XEVAL (*).
C    YFIT     MXPTS     Values of smoothed curves at X (*).
C    YPLOT    MXPLOT    Values of smoothed curves at XPLOT (*).
C    YSM      MXPTS     Smoothed result from ICSSCU, ICSFKU, ICSVKU.
C    WORK     MXWORK    Other workspace needed by various methods.
C
C  FILES USED:
C
C    LUNCRT     O     Prompts to user, and error messages
C    LUNDAT     I     Data to be smoothed (one or more curves per file)
C    LUNDEV     I     Deviations in plottable form (smooth.dev)
C    LUNKBD     I     User responses to prompts
C    LUNLOG     O     Computed coefficients or diagnostics (smooth.log)
C    LUNOUT     O     Evaluated results in SMOOTH format (default: smooth.out)
C    LUNPLT     O     Original and fitted data in plottable form (smooth.plt)
C    LUNXGR     I     User-supplied abscissas for curve fit evaluation
C
C  PROCEDURES:
C
C    ALPHA    PRMODULES  Identifies non-numeric data.
C    BEVAL    PROFILLIB  Evaluates "bump" functions (WAGFIT).
C    BOUNDS   NUMODULES  Finds max./min. values in an array.
C    COPY     NUMODULES  Copies arrays.
C    CSEVAL   INTERPLIB  Evaluates cubic spline at arbitrary abscissas.
C    CSFIT    INTERPLIB  Fits cubic interpolating spline (many end conds.).
C    CURV1    INTERPLIB  Fits interpolating spline under tension.
C    CURV2    INTERPLIB  Evaluates tension spline fitted by CURV1.
C    FSCOEF   INTERPLIB  Calculates jth Fourier coefficients.
C    FSERIES  INTERPLIB  Fourier-type analysis of irregular data (2 options).
C    FSEVAL   INTERPLIB  Evaluates partial Fourier series.
C    FSFIT    INTERPLIB  Fits Fourier coefficients in least squares sense.
C    FSSUM    INTERPLIB  Evaluates result of FSERIES at given Xs.
C    GETLINE  PRMODULES  Input utility which suppresses ! comments.
C    GETSCALE NUMODULES  Determines shift/scale factors for normalizing X[,Y].
C    GETXFORM NUMODULES  Finds coefs. for transforming interval [a,b] to [p,q].
C    ICSEVU     IMSL     Evaluates given spline.
C    ICSFKU     IMSL     Least squares spline - fixed knots.
C    ICSSCU     IMSL     Simple smoothing spline.
C    ICSVKU     IMSL     Least squares spline - variable knots.
C    IQHSCU     IMSL     Interpolation by quasi-Hermite piecewise polynomials.
C    LCSFIT   INTERPLIB  Local cubic spline method (nonparam. ver. of PLSFIT).
C    LEGENDRE INTERPLIB  Fit/evaluate a linear combination of Legendre polynoms.
C    LINE1D   INTERPLIB  Subroutine form of TABLE1.
C    LSFIT1   INTERPLIB  Local least squares spline method.
C    MSFIT    INTERPLIB  Monotonic cubic spline fit.  See also PLSFIT.
C    OPENER   PRMODULES  Utility for opening files.
C    PLSFIT   INTERPLIB  Parametric local spline method (space-efficient).
C    PNEVAL   INTERPLIB  Evaluates a polynomial at arbitrary abscissas.
C    PNFIT    INTERPLIB  Fits (one) polynomial in least squares sense.
C    PROTECT  NUMODULES  Checks monotonicity and distinctness.
C    PSFIT    INTERPLIB  Parametric interpolating spline package. Includes
C                        PSTVAL and TSUBJ, but PSEVAL is not applied here.
C    QHSFIT   INTERPLIB  Quasi-Hermite spline (reimplementation of IQHSCU).
C    QINTERP  INTERPLIB  Quadratic interpolation.
C    QPLDAT   PRMODULES  Writes data, including namelist, to a QPLOTable
C                        file through entry points QPLCRV, QPLFRM, QPLNXY.
C    RADIAL_BASIS_1D*    INTERPLIB routines for use of radial basis functions.
C    RDLIST   PRMODULES  Reads an indefinite number of integers.
C    RDXYZ2   PRMODULES  Reads X[,Y[,Z]] dataset in SMOOTH-like format.
C    READER   PRMODULES  Prompting utility (ENTRY points READC, READI, ...).
C    SELECT   PRMODULES  Menu utility.
C    SHIFTX   NUMODULES  Adds offset to array elements, in place.
C    SECOND   PRMODULES  Returns CPU time used by calling program so far.
C    SMOOTHXYZ INTERP3D  Smooths a surface or curve using a Gaussian kernel.
C    TABLE1   INTERPLIB  Function for linear interpolation by 1-D table look-up.
C    TOGGLE   PRMODULES  Menu utility for on/off options.
C    TSUBJ    INTERPLIB  See PSFIT.
C    USESCALE NUMODULES  Normalizes/denormalizes X[,Y] using GETSCALE's results.
C    VECFIT   INTERPLIB  Fits linear combination of discrete vectors read from
C                        files prompted for interactively.
C    WAGFIT   INTERPLIB  Fits linear combination of Wagner functions -
C                        appropriate for dealing with airfoils.
C    XGRID    INTERPLIB  Generates abscissas for plotting fitted curves.
C
C  ENVIRONMENT:
C
C    Fortran 90 (OpenVMS, IRIX, Linux)
C
C  HISTORY:
C
C    DAS   06/10/83   Initial design and code (ICSSCU, PNFIT; new plot frame
C                     for each method - multiple curves per frame).
C    CLH   06/16/83   Added ICSFKU option.
C    DAS   06/20/83   Added FSFIT  option for periodic data.
C    DAS   06/22/83   Added FSCOEF option for periodic data.
C    CLH   07/07/83   Added ICSVKU option.
C    CLH   10/11/83   Added ICSMOU option.
C    CLH   12/12/83   Added IQHSCU option.
C    RGL   04/25/84   Updated line type for compatibility with current
C                     QPLOT conventions; added legend for each curve;
C                     installed prompting utility READER.
C    DAS   06/29/84   All curves on a single frame now; QPLOT namelist
C                     introduced; data may be in columns or rows.
C    RGL   07/17/84   Added CSFIT option.
C    RGL   02/15/85   Added option for a file of abscissas to be used
C                     in smoothed curve evaluation and plotting;
C                     reverted to restricting data format to columns;
C                     included on-line description of available methods
C                     and of spline end conditions possible with CSFIT;
C                     incorporated QPLDAT to write the plot file.
C    DAS   03/05/85   Added WAGFIT option.
C    RGL   05/14/85   Added TABLE1 option.
C    RGL   05/23/85   Added printing of calculated coefficients.
C    RGL   07/13/85   Evaluations at uniformly-spaced abscissas (for
C                     probable plotting purposes) may be suppressed now.
C    DAS   12/20/85   Streamlined the terminal I/O.
C    RGL   03/27/86   Suppressed unnecessary output of IMSL error codes.
C    DAS   04/21/86   Removed DO WHILEs and END DOs and use of MAX in
C                     parameter statement in anticipation of move to PC;
C                     added PSFIT option (PSTVAL only - not PSEVAL).
C    DAS   05/27/86   Added VECFIT option.
C    DAS   06/12/86   Added option to scale/unscale data for the case
C                     of polynomial fits (to handle large numbers).
C    DAS   06/16/86   Provided for simple shift of data by the mean,
C                     for polynomial fits.  (These last two are for
C                     comparison purposes; full scaling is probably
C                     unnecessary.)  [P.S.: Introduction of a normalizing
C                     option used by all methods meant giving up the
C                     shift-by-the-mean option for PNFIT later - see 7/89.]
C    DAS   07/17/86   WAGFIT was simplified - needed to handle "ramp"
C                     function differently here.
C    DAS   08/14/86   Installed tension spline (CURV1, CURV2).
C    DAS   12/10/86   Refinements prompted by QuickFit version on Mac:
C                     > go back for another input data file;
C                     > arbitrary units for Fourier methods;
C                     > smarter handling of closed case for PSFIT
C                       (insert a point at I=NX+1 if necessary).
C    DAS   12/23/86   Added MSFIT/CSEVAL combination; made use of SELECT.
C    DAS   01/15/87   PSFIT can now choose MSFIT, which now has optional
C                     periodic boundary conditions.
C    RAK   02/17/87   Installed PLSFIT for testing (the TIGHT option
C                     should give the same answers as MSFIT under PSFIT).
C    DAS   02/20/87   > Installed LINE1D in parallel with TABLE1. (Both
C                       are retained because SMOOTH serves to test them.)
C                     > Provided for plotting the deviations for least
C                       squares methods by setting XEVAL (*) = X (*).
C                     > SELECT now permits no default - handy at "done"
C                       time, when CR = CTRL Z = quit.
C                     > Echoed RMS and maximum deviations to the screen
C                       as well as to the log file.
C    DAS/RAK 3/3/87   > Straightened out some logic errors in the above.
C                     > Added Y (fit) as third column in smooth.dev.
C                     > Introduced PROTECT for checking monotonicity, etc.
C                     > FSFIT now works for N=0.
C    DAS   10/20/87   > Provided for timing (fit + evaluation) for each curve
C                       of each file more precisely.  (Total CPU time is also
C                       measured, but it includes the I/O.)
C                     > CSFIT and MSFIT can handle decreasing abscissas now
C                       because of revised search used by CSEVAL.
C                     > Upped MXPTS to 1001 (but PSFIT/PSEVAL retain internal
C                       COMMON block limit of 500).
C    DAS   11/20/87   Installed LSFIT1 to test it.
C    DAS   06/10/88   Safeguarded Fourier methods better; upped MXFCOF to 180.
C    DAS   10/22/88   Wagner function option now handles arbitrary units,
C                     including descending abscissas.
C    DAS   02/17/89   Handled LSFIT1's option to suppress spline step.
C    DAS   03/30/89   Added prompt for a subtitle on the plot.
C    DAS   05/08/89   LSFIT1 should have NEARPTS >= DEGREE + 2, not + 1.
C    DAS   07/09/89   Introduced GETSCALE/USESCALE for all methods.  This was
C                     prompted by a wind tunnel data case that worked fine in
C                     QPLOT, which normalizes before applying PLSFIT, but was
C                     poorly handled by SMOOTH's PLSFIT option.  (Only PNFIT
C                     had been given scaling/shifting options; now it is clear
C                     that ALL methods are best applied to normalized data.)
C    DAS   07/15/89   Installed FSERIES/FSSUM/FSEVAL2 (mainly to test Macintosh
C                     versions).
C    DAS   08/31/89   Installed LCSFIT.  Small fix to blank out legends after
C                     the first dataset.
C    DAS   09/09/89   EPS needed for WAGFIT [0,1] test; YPLOT(*) was not being
C                     vertically shifted; legend text for LSFIT1 was being
C                     clobbered by deviations info.; shortened it.
C    DAS   05/29/90   "No. of uniform abscissas" prompt was imprecise for
C                     parametric fits.  Normalizing prompt also needs a warning.
C    DAS   07/30/90   Installed RDXYZ which became available from program
C                     PROFILE, allowing commented-out points in the input files
C                     and eliminating list-directed I/O and some error handling.
C    DAS   01/17/91   Replaced XFORMX with GETXFORM/USESCALE.  MXPTS and MXPLOT
C                     raised to 1500.
C    DAS   05/27/91   Results may now be saved in SMOOTH format if desired (as
C                     well as or instead of QPLOT format).  Also handled the
C                     NX=0 dataset case.  [D. Serafini's request in both cases.]
C    DAS   06/07/91   Replaced IMSL's IQHSCU with QHSFIT (same algorithm).
C    DAS   06/17/91   Added date/time stamp as a caption in the QPLOTable file.
C    DAS   06/20/91   LCSFIT now has a periodic end-condition option.
C    DAS   10/16/91   PSFIT and PSTVAL usage changed.
C    DAS   12/11/91   Introduced multi-column capability plus an alternative
C                     "indefinite-number-of-points-per-curve" data format, as
C                     both supported by a revised form of RDXYZ.  The upward
C                     compatibility with SMOOTH's original "NPTS" data format
C                     was achieved at some cost in complexity, and not all cases
C                     of mismatched NPTS values are trapped gracefully.  The
C                     alternative data format is now preferred.
C    DAS   12/31/91   Made TITLE optional for data files and target X file.
C    DAS   02/12/92   TITLE was being skipped at the top of the loop over
C                     datasets, even when none was present.  And the test for
C                     an empty file shouldn't have decremented DATASET if
C                     the (previous) value of NX was not zero.
C    DAS   08/16/97   CSEVAL needs NDATA = -NX for the periodic case now.
C    DAS   10/14/99   Fortran 90 translation (mainly to eliminate 0 carriage
C                     control).  Always display sum of squares & RMS deviation
C                     whether deviations are being plotted or not.
C    DAS   07/07/00   Added QINTERP.
C    DAS   09/27/00   Plugged SMOOTHXYZ into the defunct ICSMOU slot.
C    DAS   01/28/04   Increased LSFIT1's limit to 100 neighboring pts., not 10.
C    DAS   06/14/04   For FSFIT, the transformation to [0, 2pi] was using
C                     X(1) and X(NX) but the abscissas need not be ordered.
C    DAS   11/19/10   Added exponential curve fits AEBXPLUSC and ABXECX.
C    DAS   11/23/10   Added power law curve fit AXBPLUSC.
C    DAS   12/28/10   Exponential & power law fits need not trap unordered data.
C    DAS   06/17/16   Printed coefficients have more digits showing now.
C                     No more use of 1P: ES formatting is preferable.
C    DAS   11/23/16   Make the comment character # instead of !, for Tecplot
C                     compatibility.  This required RDXYZ2 variant with COMMENT
C                     as an argument.  Suppress output of deviations by default.
C    DAS   07/01/20   Added LEGENDRE option.
C    DAS   12/29/22   Added RADIAL_BASIS_1D_WEIGHTS & RADIAL_BASIS_1D_EVAL.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C                  Now with AMA, Inc. at NASA ARC.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Constants:
C     ----------

      CHARACTER * 1, PARAMETER ::
     >   BLANK = ' ', COMMENT = '#', DENRMLIZ = 'D', INDEPEND = 'I',
     >   NRMLIZ = 'N'

      REAL, PARAMETER ::
     >   EIGHT = 8.E+0, EPS = 1.E-6, OMIT = 999.E+0, ONE = 1.E+0,
     >   SMALL = 0.05E+0, ZERO = 0.E+0

      LOGICAL, PARAMETER ::
     >   FALSE = .FALSE., TRUE = .TRUE.

      INTEGER, PARAMETER ::
     >   IMODE = 0, IOMIT = 0, LUNDAT = 1, LUNXGR = 2, LUNPLT = 3,
     >   LUNDEV = 4, LUNKBD = 5, LUNCRT = 6, LUNLOG = 7, LUNOUT = 8,
     >   MXFSFIT = 20, MXLTYP = 6, MXMENU = 24, MXNXK = 28, MXNDEG = 12,
     >   MXLEGENDRE = 20, ! Presumably too-high a degree isn't wise
     >   MXLSFIT1 = 100,  ! Must match hard-coded limit in LSFIT1
     >   MXPTS = 100000, MXPLOT = MXPTS, MXVEC = 10, MXWAG = 20,
     >   MXWORK = 2000000, MXFCOF = (MXPTS-1)/2, NDIM = 2

C     N.B.:  PSFIT is limited by an internal COMMON block which may allow
C     -----  fewer than MXPTS points.
C
C     Notes on work-space allocation:        MXPTS=500    Method
C     -------------------------------          Value
C
C     MXWORK is   MAX (MXPTS * (MXNDEG+4),      9000      PNFIT
C                      MXPTS * (2*MXFSFIT+3),  16500      FSFIT
C                      MXPTS * (MXNXK+6),      17000      ICSFKU, ICSVKU
C                      MXPTS * (MXWAG+2),      11500      WAGFIT
C                      MXPTS * (2*MXVEC+2),    11000      VECFIT
C                      7*MXPTS+14,              3514      ICSSCU
C
C     There has been no attempt to keep storage to an absolute minimum
C     in SMOOTH - reuse of WORK (MXWORK) covers the bulk of this issue.
C

C     Variables:
C     ----------

      REAL
     >   A (0:MXFCOF), B (0:MXFCOF), C (0:MXWAG), COEFS (MXPTS,3),
     >   LEGENDRE_COEFS(0:MXLEGENDRE), ! For linear combination of P0:PN
     >   SCALE (2), SHIFT (2), WGHTS (MXPTS), WORK (MXWORK), X (MXPTS),
     >   X0 (MXPTS), XEVAL (MXPTS), XEVAL0 (MXPTS), XPLOT (MXPLOT),
     >   XK (MXNXK), XNORM (MXPTS), Y (MXPTS), Y0 (MXPTS), YDEV (MXPTS),
     >   YEVAL (MXPTS), YFIT (MXPTS), YPLOT (MXPLOT), YSM (MXPTS)

      INTEGER
     >   I, ICOLUMNS (2), IDATASET, IDATE_TIME (8), IENDL, IENDR, IER,
     >   INCR, INCR0, INDEX, IOS, ITYPE, J, LAST, LEFT, LUN, M, MODE,
     >   NCOL, NDATA, NDEG, NEARPTS, NEVAL, NFCOF, NFCREQ, NIT, NPLOT,
     >   NTPLOT, NUMPLT, NX, NXK

      REAL
     >   AE, ARROW, BE, CE, DELTAT, DERIVL, DERIVR, DEVMAX, DIS, ERROR,
     >   FWHM, RMSDEV, SCR, SHAPE, SM, SSQ, SSQMIN, TBEGIN, TEND, TENSE,
     >   TIMEA, TIMEB, TIME1, TIME2, TOTALT, TWOPI, XAXIS, XMAX, XMEAN,
     >   XMIN, YMAX, YMIN, YNNORM, Y1

      LOGICAL
     >   ADD, ARRHENIUS, AUGMENTED, CLOSED, CYCLIC, DEFAULT, DISTINCT,
     >   FINIS, FIRSTFIL, FIRSTSET, FIRSTTYP, LEASTSQRS, MONO, NMLIZ,
     >   NOABS, NOFIT, ONOFF (1:4), PERIODIC, PLTABS, PLTCRV, PLTDEV,
     >   QUIT, RAMP, SHOWITERS, THRU00, YES, YESOUT, YESPLT, YESTITLE

      CHARACTER
     >   CDATE_TIME (3) * 10, CHOICE * 6, DATAFILE * 60,
     >   DSTRIB * 3, HELP (14) * 58, LEGEND * 91, LINE * 9,
     >   LTYPES (MXLTYP) * 9, MENU (MXMENU) * 63, METHOD * 3,
     >   OPTIONS (1:4) * 34, OUTFILE * 60, REPLY * 1, SHORTLEG * 35,
     >   SUBTTL * 80, TBUF * 18, TTLDAT * 80, TTLOUT * 80,
     >   TTLPLT * 80, XFILE * 60, XLABEL * 60, YLABEL * 60

C     Procedures:
C     -----------

      REAL
     >   CURV2, TABLE1, TSUBJ

      LOGICAL
     >   ALPHA

      EXTERNAL
     >   ALPHA,    AEBXPLUSC, AXBECX,   BEVAL,     BOUNDS,    COPY,
     >   CSEVAL,   CSFIT,     CURV1,    CURV2,     FSCOEF,    FSEVAL,
     >   FSFIT,    GETLINE,   GETSCALE, GETXFORM,  ICSEVU,    ICSFKU,
     >   ICSSCU,   ICSVKU,    LCSFIT,   LEGENDRE,  LINE1D,    LSFIT1,
     >   MSFIT,    OPENER,    PLSFIT,   PNEVAL,    PNFIT,     PROTECT,
     >   PSFIT,    PSTVAL,    QHSFIT,   QINTERP,   QPLCRV,    QPLFRM,
     >   QPLNXY,   RDLIST,    RDXYZ2,   READI,     READR,     READS,
     >   READY,    SECOND,    SELECT,   SHIFTX,    SMOOTHXYZ, TABLE1,
     >   TOGGLE,   TSUBJ,     UGETIO,   USESCALE,  VECFIT,    WAGFIT,
     >   XGRID

      EXTERNAL
     >   RADIAL_BASIS_1D_WEIGHTS, RADIAL_BASIS_1D_EVAL

C     Storage:
C     --------

      DATA WGHTS
     >   /MXPTS * 1.E+0/

      DATA LTYPES
     >   /'SOLID','DASH', 'DOT', 'CHAINDASH', 'CHAINDOT', 'LONGDASH'/

      DATA OPTIONS
     >   /'1:  Explain input data formats?   ',
     >    '2:  Save results in QPLOT format? ',
     >    '3:  Save results in SMOOTH format?',
     >    '4:  Save curve fit deviations?    '/

      DATA HELP
     >   /' [TITLE]         <Optional title, up to 80 characters.>   ',
     >    ' [N1]            <# pts. in 1st curve: not needed if "END"',
     >    ' X (1)  Y (1)        is present; both are redundant if    ',
     >    ' X (2)  Y (2)        there is only one curve.>            ',
     >    ' X (3)  Y (3)                                             ',
     >    ' X (4)  Y (4)    <X and Y may be extracted from specified ',
     >    ' :      :            columns if there are more than two.> ',
     >    ' :      :                                                 ',
     >    ' # X     Y       <"#" suppresses points or adds comments.>',
     >    ' X (N1) Y (N1)                                            ',
     >    '                 <Blank lines are ignored.>               ',
     >    ' N2  or  END     <Repeat for more curves, omitting title.>',
     >    ' X (1)  Y (1)                                             ',
     >    ' :      :                                                 '/

      DATA MENU /
     >'  1:  SMTHXYZ- weighted averaging via a Gaussian kernel        ',
     >'  2:  ICSSCU - cubic spline with variable smoothing parameter  ',
     >'  3:  ICSFKU - least squares cubic spline with fixed knots     ',
     >'  4:  ICSVKU - least squares cubic spline with variable knots  ',
     >'  5:  PNFIT  - least squares polynomial approximation          ',
     >'  6:  FSFIT  - pseudo-Fourier series; least sqrs/irregular data',
     >'  7:  FSERIES- pseudo-Fourier series; trapezoid/irregular data ',
     >'  8:  FSCOEF - Fourier series by recursion; regular data       ',
     >'  9:  QHSFIT - interpolating quasi-Hermite cubic spline        ',
     >' 10:  CSFIT  - interpolating cubic spline/many end conditions  ',
     >' 11:  TABLE1 - or LINE1D - linear interpolation/extrapolation  ',
     >' 12:  PSFIT  - parametric spline: X = X(T), Y = Y(T) via CSFIT ',
     >' 13:  WAGFIT - combination of Wagner functions (least squares) ',
     >' 14:  VECFIT - combination of discrete vectors (least squares) ',
     >' 15:  CURV1  - interpolating cubic spline under tension        ',
     >' 16:  MSFIT  - monotonic interpolating cubic spline            ',
     >' 17:  PLSFIT - parametric local cubic spline (low storage)     ',
     >' 18:  LCSFIT - nonparametric form of PLSFIT (+ linear option)  ',
     >' 19:  LSFIT1 - local lst. sqrs. polynomial/spline hybrid method',
     >' 20:  QINTERP- quadratic interpolation                         ',
     >' 21:  AXBECX - or AEBXPLUSC exponential curve fits (lst. sqrs.)',
     >' 22:  AX^B+C - power law curve fit (least squares)             ',
     >' 23:  LEGENDRE (linear combination of Legendre polynomials 0:n)',
     >' 24:  Radial basis functions (interpolation/no smoothing)      '/


C     Execution
C     ---------

      TWOPI = EIGHT * ATAN (ONE)

      WRITE (LUNCRT, 1003)
     >   ' Program SMOOTH applies 1-D interpolation and',
     >   ' smoothing methods.'

C     One or two formats may be desired for the results:

      ONOFF (1) = FALSE
      ONOFF (2) = FALSE
      ONOFF (3) = TRUE
      ONOFF (4) = FALSE

      CALL TOGGLE (4, OPTIONS, ONOFF, LUNCRT, LUNKBD)

      IF (ONOFF (1)) WRITE (LUNCRT, 1000) BLANK, HELP, BLANK

      YESPLT = ONOFF (2)
      YESOUT = ONOFF (3)
      PLTDEV = ONOFF (4)

      IF (YESPLT) OPEN (UNIT=LUNPLT, FILE='smooth.plt',STATUS='UNKNOWN')

      IF (YESOUT) THEN
         OUTFILE = 'smooth.out'
         CALL OPENER (LUNCRT,
     >      'Enter filename for evaluated results (<CR>=smooth.out): ',
     >      LUNKBD, OUTFILE, LUNOUT, 'UNKNOWN')
      END IF

      IF (PLTDEV) OPEN (UNIT=LUNDEV, FILE='smooth.dev',STATUS='UNKNOWN')

      OPEN (UNIT=LUNLOG, FILE='smooth.log', STATUS='UNKNOWN')

C ... Reassign the logical unit for IMSL error messages:

      CALL UGETIO (3, 0, LUNCRT)

      FIRSTFIL = TRUE

      CALL SECOND (TIME1)

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C                 Begin loop over distinct data files
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  100 CONTINUE

      DATAFILE = 'smooth.dat'
      CALL OPENER (LUNCRT,
     >   'Enter the input data file name (<CR>="smooth.dat"): ',
     >   LUNKBD, DATAFILE, LUNDAT, 'OLD')

C ... Look for the optional title:

  105 CALL GETLINE (LUNDAT, COMMENT, TTLDAT, LAST, IOS)
      IF (IOS  /= 0) GO TO 800
      IF (LAST == 0) GO TO 105

C ... If we have strictly numeric data, assume the title is missing.

      YESTITLE = ALPHA (TTLDAT)

      IF (YESTITLE) THEN
         WRITE (LUNCRT, '(/, A, //, 1X, A, /)')
     >      ' Data file title found:', TTLDAT
      ELSE   ! The line will be reread as part of the data proper.
         TTLDAT = '" "'
      END IF

C ... Prompt for plot title (s) and axis labels:

      IF (YESPLT .OR. PLTDEV) THEN
         TTLPLT = TTLDAT
         CALL READS (LUNCRT,
     >      'Enter plot title (<CR>=file title, or blank if none): ',
     >      LUNKBD, TTLPLT, DEFAULT, QUIT)
         IF (QUIT) GO TO 999
         IF (YESOUT) THEN
            TTLOUT = TTLPLT
            CALL READS (LUNCRT,
     >         'Enter title for evaluated results (<CR>=plot title): ',
     >         LUNKBD, TTLOUT, DEFAULT, QUIT)
            IF (QUIT) GO TO 999
         END IF

         SUBTTL = '" "'
         CALL READS (LUNCRT, 'Enter plot subtitle (<CR>=none): ',
     >      LUNKBD, SUBTTL, DEFAULT, QUIT)
         IF (QUIT) GO TO 999

         XLABEL = 'X'
         CALL READS (LUNCRT, 'Enter abscissa axis label (<CR>="X"): ',
     >      LUNKBD, XLABEL, DEFAULT, QUIT)
         IF (QUIT) GO TO 999

         YLABEL = 'Y'
         CALL READS (LUNCRT, 'Enter ordinate axis label (<CR>="Y"): ',
     >      LUNKBD, YLABEL, DEFAULT, QUIT)
         IF (QUIT) GO TO 999

      ELSE

         IF (YESOUT) THEN
            TTLOUT = TTLDAT
            CALL READS (LUNCRT,
     >    'Enter title for evaluated results (<CR>=input data title): ',
     >                  LUNKBD, TTLOUT, DEFAULT,QUIT)
            IF (QUIT) GO TO 999
         END IF

      END IF

C ... We can't support multi-column files without some kind of prompt:

  110 ICOLUMNS (1) = 1
      ICOLUMNS (2) = 2
      NCOL = 2
      CALL RDLIST (LUNCRT, 'Enter X and Y column numbers. <CR>=1 & 2: ',
     >             LUNKBD, NCOL, ICOLUMNS)
      IF (NCOL  < 0) GO TO 999          ! EOF = Quit
      IF (NCOL == 1) GO TO 110          ! 0 or 2 entries are OK

C ... Provide for normalizing the data:

      WRITE (LUNCRT, 1003)
     >   ' Transforming X and Y to the interval [0,1] benefits some',
     >   ' methods',
     >   ' (but may affect the spacing of evaluated points for',
     >   ' parametric fits).'
      NMLIZ = FALSE
      CALL READY (LUNCRT,
     >   'Do you want to normalize the data? (Yes|No; <CR>=No): ',
     >   LUNKBD, NMLIZ, DEFAULT, QUIT)
      IF (QUIT) GO TO 999

C ... Provide for evaluating the fitted curves on a uniform grid
C     (either for plotting purposes or perhaps for the reusable results):

      WRITE (LUNCRT, 1000)
  130 CONTINUE
         NPLOT = 201
         WRITE (LUNCRT, 1002)
     >      ' Enter the number of uniform Xs',
     >      ' (or points along the arc for parametric fits)'
         CALL READI (LUNCRT,
     >      'at which to evaluate each curve fit. (0=none; <CR>=201): ',
     >      LUNKBD, NPLOT, DEFAULT, QUIT)
         IF (QUIT) GO TO 999
         IF ((NPLOT > 0 .AND. NPLOT < 2) .OR. NPLOT > MXPLOT)
     >GO TO 130

      PLTCRV = NPLOT > 0

C ... Provide for evaluating at abscissas read from a file:

      PLTABS = FALSE
      WRITE (LUNCRT, 1001)
     >   ' Each fit may be evaluated at abscissas read from a file.'
      NEVAL = 0
      XFILE = BLANK
      CALL OPENER (LUNCRT, 'Enter the file name, or <CR> if none: ',
     >   LUNKBD, XFILE, LUNXGR, 'OLD')

      PLTABS = XFILE /= BLANK

      IF (.NOT. (PLTCRV .OR. PLTDEV .OR. PLTABS)) THEN
         WRITE (LUNCRT, 1001)
     >      ' You have suppressed all useful output.  Try again.'
         GO TO 130
      END IF

      IF (PLTABS) THEN

C ...    Read supplied abscissas, ignoring any title and any extra columns:

  140    CALL GETLINE (LUNXGR, COMMENT, TTLDAT, LAST, IOS)
         IF (IOS  /= 0) GO TO 810
         IF (LAST == 0) GO TO 140

C ...    If we have strictly numeric data, assume the title is missing.

         IF (.NOT. ALPHA (TTLDAT)) REWIND LUNXGR

         CALL RDXYZ2 (1, LUNXGR, LUNCRT, COMMENT, 1, MXPTS, NEVAL,
     >                XEVAL, XEVAL, XEVAL, FINIS, IER)
         IF (IER /= 0) GO TO 810   ! RDXYZ explains errors except EOF

         IF (NMLIZ) XEVAL0 (1 : NEVAL) = XEVAL (1 : NEVAL)

         CLOSE (LUNXGR)

      END IF

C ... Initialize the plot frame(s):

C     DATE_AND_TIME returns CCYYMMDD, hhmmss.sss, +/-hhmm, & 8 integers:

      CALL DATE_AND_TIME (CDATE_TIME (1), CDATE_TIME (2),
     >                    CDATE_TIME (3), IDATE_TIME)

      TBUF = 'CCYY-MM-DD @ hh:mm'
      TBUF (1 : 4)  = CDATE_TIME (1) (1 : 4)
      TBUF (6 : 7)  = CDATE_TIME (1) (5 : 6)
      TBUF (9 :10)  = CDATE_TIME (1) (7 : 8)
      TBUF (14:15)  = CDATE_TIME (2) (1 : 2)
      TBUF (17:18)  = CDATE_TIME (2) (3 : 4)

      IF (YESPLT) THEN
         CALL QPLFRM (LUNPLT, TTLPLT, SUBTTL, XLABEL, YLABEL)
         WRITE (LUNPLT, 1000) TBUF
      END IF

      IF (PLTDEV) THEN
         CALL QPLFRM (LUNDEV, TTLPLT, SUBTTL, XLABEL, 'Fit - Data')
         WRITE (LUNDEV, 1000) TBUF
      END IF

C ... Initialize the reusable results file and printable log file:

      IF (FIRSTFIL) THEN
         IF (YESOUT) THEN
            WRITE (LUNOUT, 1000) TRIM (TTLOUT)
         END IF
         WRITE (LUNLOG, 1001)
     >   ' Program SMOOTH:  Application of 1-D Data Fitting Algorithms',
     >   ' :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
      ELSE
         WRITE (LUNLOG, 1000) BLANK, '1', BLANK
      END IF

      WRITE (LUNLOG, 1003) ' Input data file:  ', DATAFILE

      IF (YESPLT) THEN
         WRITE (LUNLOG, '(/, A, /, 1X, A)')
     >      ' Plot title prompted for follows:', TRIM (TTLPLT)
         IF (SUBTTL /= BLANK) WRITE (LUNLOG, '(/, A, /, 1X, A)')
     >      ' Subtitle:', TRIM (SUBTTL)
      END IF

      IF (NMLIZ) THEN
         WRITE (LUNLOG, 1003)
     >      ' Data normalized to [0,1] prior to application of',
     >      ' selected method(s).'
      END IF

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C                Begin loop over curve fits requested
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      FIRSTTYP = TRUE
      NUMPLT = 0

  200 CONTINUE

         IF (FIRSTTYP) THEN
            WRITE (LUNCRT, 1001) ' Interpolation/smoothing options:'
            ITYPE = 18
            CHOICE = MENU (ITYPE) (7 : 13)
         ELSE
            CHOICE = BLANK
         END IF

         CALL SELECT ('Select a curve fit by number or name.',
     >                MXMENU, MENU, .NOT. FIRSTTYP, LUNCRT, LUNKBD,
     >                ITYPE, CHOICE, QUIT)
         IF (QUIT) GO TO 900

         NUMPLT = NUMPLT + 1
         NTPLOT = NPLOT
         LEASTSQRS = FALSE  ! More often than not

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C        Begin loop over (unknown) number of datasets in file
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

C ...    Skip file title if present (identified outside this loop over methods):

         REWIND LUNDAT
         IF (YESTITLE) THEN
  290       CALL GETLINE (LUNDAT, COMMENT, TTLDAT, LAST, IOS)
            IF (LAST == 0) GO TO 290
         END IF

         IDATASET = 0
         FIRSTSET = TRUE
         NOFIT = FALSE   ! Reset below if method just adjusts ordinates
         NOABS = FALSE   !   "    "    "    "    is parametric (Y is not Y(X))

  300    CONTINUE

            IDATASET = IDATASET + 1

C ...       Look for a dataset.  EOF may be normal except first time through.
C           MXPTS - 1 is used here to allow for closing curve for some methods.

            CALL RDXYZ2 (2, LUNDAT, LUNCRT, COMMENT, ICOLUMNS, MXPTS - 1,
     >                   NX, X, Y, Y, FINIS, IER)

            IF (IER == 1) GO TO 700  ! EOF on first read; FINIS used at end
            IF (IER /= 0) GO TO 820  ! RDXYZ explains other errors
            IF (NX  == 0) GO TO 300  ! Has been needed for PROFILE datasets

C           Keep a copy of the raw data because some methods change it in place:

            DO I = 1, NX
               X0 (I) = X (I)
               Y0 (I) = Y (I)
            END DO

C ...       Transform the data?

            IF (NMLIZ) THEN

               CALL GETSCALE (INDEPEND, NDIM, NX, X, Y, Y, SCALE, SHIFT,
     >                        IER)
               CALL USESCALE (NRMLIZ,   NDIM, NX, X, Y, Y, SCALE, SHIFT,
     >                        IER)

C ...          Transform the user-supplied abscissas?

               IF (PLTABS) THEN

                  XEVAL (1 : NEVAL) = XEVAL0 (1 : NEVAL)

                  CALL USESCALE (NRMLIZ, 1, NEVAL, XEVAL, XEVAL, XEVAL,
     >                           SCALE, SHIFT, IER)
               END IF

            END IF

C           Some methods can use the data range for defaults, so always find it:

            XMIN = X (1)
            XMAX = XMIN

            CALL BOUNDS (NX, 1, MXPTS, X, XMIN, XMAX)

            IF (PLTCRV) THEN

C ...          Generate uniformly abscissas in data range for plotting.
C              Some curve fits do not require monotonic abscissas.

               CALL XGRID (NPLOT, 0, XMIN, XMAX, XPLOT)

            END IF

C ...       Set flags for error checking below.  (Not all methods will
C           require all tests, but it's handy to do them all at once.)

            CALL PROTECT (NX, X, Y, ARROW, DISTINCT)

C ...       Compute and evaluate selected fit type.  Some fits require
C           additional prompts, suppressed after the first curve.

            SELECT CASE (ITYPE)

            CASE (1)  ! *** SMOOTHXYZ ***

C              Smooth a noisy curve using a Gaussian kernel.
C              (Written for a surface, but the method is 1-D.)

CCCC           LEASTSQRS = TRUE  ! No; points are smoothed in-place

               IF (FIRSTSET) THEN

                  LEGEND = 'SMOOTHXYZ: FWHM = xxxxxxxxx, NITER = n'

                  FWHM = 0.001 * (XMAX - XMIN)
                  WRITE (LUNCRT, '(4(/, A), ES9.3)')
     >      ' The full width at half maximum of the Gaussian kernel',
     >      ' is ~2.35 sigma, in units of distance from current point.',
     >      ' Use small values for removing high frequency noise.',
     >      ' The default is 0.1% of the X data range: ', FWHM

                  CALL READR (LUNCRT, 'FWHM: ',
     >               LUNKBD, FWHM, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710

                  WRITE (LEGEND(19:27), '(ES9.3)') FWHM

                  NIT = 1
                  CALL READI (LUNCRT, '# iterations (1:9) [<CR> = 1]: ',
     >               LUNKBD, NIT, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710

                  WRITE (LEGEND(38:38), '(I1)') NIT

                  WRITE (LUNLOG, 1100) TRIM (LEGEND)

               END IF

C ...          A temporary "NPLOT" is needed - no fitted curve to evaluate:

               NOFIT  = TRUE
               NOABS  = TRUE
               NTPLOT = NX

               DO I = 1, NX
                  XPLOT (I) = X (I) ! X/YPLOT(*) will be smoothed in place
                  YPLOT (I) = Y (I)
                  WORK (I)  = ZERO  ! For the missing Z argument
               END DO

               CALL SECOND (TIMEA)

               CALL SMOOTHXYZ (NX, XPLOT, YPLOT, WORK, FWHM, NIT)

               CALL SECOND (TIMEB)

            CASE (2)  ! *** ICSSCU ***

C              Spline with smoothing parameter >= 0.

  420          IF (ARROW /= ONE) GO TO 830

               LEASTSQRS = TRUE

               IF (FIRSTSET) THEN
                  SM = 0.001E+0
                  CALL READR (LUNCRT,
     >    'Smoothing parameter SM? (<CR> gives 0.001; SM*NX is used): ',
     >               LUNKBD, SM, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710

                  IF (SM < ZERO .OR. SM > 10.E+0) GO TO 420

                  LEGEND = 'ICSSCU: SM = z.zzzzzzz'
                  WRITE (LEGEND (14:22), '(ES9.3)') SM
                  WRITE (LUNLOG, 1100) TRIM (LEGEND)
               END IF

               CALL SECOND (TIMEA)

               CALL ICSSCU (X, Y, WGHTS, NX, SM*NX, YSM, COEFS, MXPTS,
     >                      WORK, IER)
               IF (IER /= 0) GO TO 860

               CALL ICSEVU (X, YSM, NX, COEFS, MXPTS, X, YFIT,
     >                      NX, IER)
               IF (IER /= 0) GO TO 860

               IF (PLTCRV) THEN
                  CALL ICSEVU (X, YSM, NX, COEFS, MXPTS, XPLOT, YPLOT,
     >                         NPLOT, IER)
                  IF (IER /= 0) GO TO 860
               END IF

               IF (PLTABS) THEN
                  CALL ICSEVU (X, YSM, NX, COEFS, MXPTS, XEVAL, YEVAL,
     >                         NEVAL, IER)
                  IF (IER /= 0) GO TO 860
               END IF

               CALL SECOND (TIMEB)

               CALL PRCOEF (1, NX-1, COEFS (1,1), IDATASET,
     >                      'COEFS (1:NX-1,1):',  IOMIT)
               CALL PRCOEF (1, NX-1, COEFS (1,2), IOMIT,
     >                      'COEFS (1:NX-1,2):',  IOMIT)
               CALL PRCOEF (1, NX-1, COEFS (1,3), IOMIT,
     >                      'COEFS (1:NX-1,3):',  IOMIT)

            CASE (3)  ! *** ICSFKU ***

C              Least squares spline, fixed knots.

               IF (ARROW /= ONE) GO TO 830

               LEASTSQRS = TRUE

C              ICSFKU accepts a maximum of 28 knots.  Pick every "nth"
C              data point as knots:  (Choosing EVERY X seems to cause
C              ICSFKU to fail for reasons unknown at time of writing.)

               IF (FIRSTSET) THEN
                  INCR0 = 2
                  CALL READI (LUNCRT,
     > 'Knot-choosing increment? <CR> gives every 2nd point (28 max): ',
     >               LUNKBD, INCR0, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710

                  INCR = INCR0
  435             NXK  = (NX + (INCR-1)) / INCR

                  IF (NXK > MXNXK) THEN
                     INCR = INCR + 1
                     GO TO 435
                  END IF

                  XK (1) = X (1)
                  J = 1 + INCR

                  DO I = 2, NXK - 1
                     XK (I) = X (J)
                     J = J + INCR
                  END DO
                  XK (NXK) = X (NX)

                  LEGEND = 'ICSFKU:  INCR0 = nnn'
                  WRITE (LEGEND (18:20), '(I3)') INCR0
                  WRITE (LUNLOG, 1100) TRIM (LEGEND)
               END IF

               CALL SECOND (TIMEA)

               CALL ICSFKU (X, Y, NX, IMODE, XK, NXK, YSM, COEFS, MXPTS,
     >                      ERROR, WORK, IER)
               IF (IER /= 0) GO TO 860

               CALL ICSEVU (XK, YSM, NXK, COEFS, MXPTS, X, YFIT,
     >                      NX, IER)
               IF (IER /= 0) GO TO 860

               IF (PLTCRV) THEN
                  CALL ICSEVU (XK, YSM, NXK, COEFS, MXPTS, XPLOT, YPLOT,
     >                         NPLOT, IER)
                  IF (IER /= 0) GO TO 860
               END IF

               IF (PLTABS) THEN
                  CALL ICSEVU (XK, YSM, NXK, COEFS, MXPTS, XEVAL, YEVAL,
     >                         NEVAL, IER)
                  IF (IER /= 0) GO TO 860
               END IF

               CALL SECOND (TIMEB)

               CALL PRCOEF (1, NXK, COEFS (1,1), IDATASET,
     >                      'COEFS (1:NXK,1):',  IOMIT)
               CALL PRCOEF (1, NXK, COEFS (1,2), IOMIT,
     >                      'COEFS (1:NXK,2):',  IOMIT)
               CALL PRCOEF (1, NXK, COEFS (1,3), IOMIT,
     >                      'COEFS (1:NXK,3):',  IOMIT)

            CASE (4)  ! *** ICSVKU ***

C              Least squares spline, variable knots.

               IF (ARROW /= ONE) GO TO 830

               LEASTSQRS = TRUE

C              ICSVKU accepts a maximum of 28 knots.  Pick every "nth"
C              data point as initial guesses.  ICSVKU uses ICSFKU
C              iteratively to optimize the location of the knots.

               IF (FIRSTSET) THEN
                  INCR0 = 2
                  CALL READI (LUNCRT,
     >    'Initial knot-choosing inc.? <CR> = every 2nd pt. (28 max): ',
     >               LUNKBD, INCR0, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710

                  INCR = INCR0
  445             NXK  = (NX + (INCR-1)) / INCR

                  IF (NXK > MXNXK) THEN
                     INCR = INCR + 1
                     GO TO 445
                  END IF

                  XK (1) = X (1)
                  J = 1 + INCR

                  DO I = 2, NXK - 1
                     XK (I) = X (J)
                     J = J + INCR
                  END DO
                  XK (NXK) = X (NX)

                  LEGEND = 'ICSVKU: Initial INC = nn'
                  WRITE (LEGEND (23:24), '(I2)') INCR0
                  WRITE (LUNLOG, 1100) TRIM (LEGEND)
               END IF

               CALL SECOND (TIMEA)

               CALL ICSVKU (X, Y, NX, XK, NXK, YSM, COEFS, MXPTS, ERROR,
     >                      WORK, IER)
               IF (IER /= 0) GO TO 860

               CALL ICSEVU (XK, YSM, NXK, COEFS, MXPTS, X, YFIT,
     >                      NX, IER)
               IF (IER /= 0) GO TO 860

               IF (PLTCRV) THEN
                  CALL ICSEVU (XK, YSM, NXK, COEFS, MXPTS, XPLOT, YPLOT,
     >                         NPLOT, IER)
                  IF (IER /= 0) GO TO 860
               END IF

               IF (PLTABS) THEN
                  CALL ICSEVU (XK, YSM, NXK, COEFS, MXPTS, XEVAL, YEVAL,
     >                         NEVAL, IER)
                  IF (IER /= 0) GO TO 860
               END IF

               CALL SECOND (TIMEB)

               CALL PRCOEF (1, NXK, COEFS (1,1), IDATASET,
     >                      'COEFS (1:NXK,1):',  IOMIT)
               CALL PRCOEF (1, NXK, COEFS (1,2), IOMIT,
     >                      'COEFS (1:NXK,2):',  IOMIT)
               CALL PRCOEF (1, NXK, COEFS (1,3), IOMIT,
     >                      'COEFS (1:NXK,3):',  IOMIT)

            CASE (5)  ! *** PNFIT ***

C              Polynomial fit by linear least squares (QR factorization)

  450          IF (FIRSTSET) THEN
                  LEASTSQRS = TRUE
                  NDEG = -1
                  CALL READI (LUNCRT,
     >               'Enter degree of polynomial to be fitted: ',
     >               LUNKBD, NDEG, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710

                  IF (NDEG < 0) THEN
                     WRITE (LUNCRT, 1001) ' Must be at least 0.'
                     GO TO 450
                  END IF

                  IF (NDEG > 7) THEN
                     WRITE (LUNCRT, 1003)
     >                  ' CAUTION:  Fitting high-order polynomials',
     >                  ' leads to ill-conditioned equations.'
                     IF (NDEG > MXNDEG) THEN
                        WRITE (LUNCRT, '(A, I3, A)')
     >                     ' Degree is being limited to', MXNDEG,
     >                     ' for work-space reasons.'
                        NDEG = MXNDEG
                     END IF
                     WRITE (LUNCRT, 1001) ' Proceeding...'
                  END IF

                  THRU00 = FALSE
                  IF (.NOT. NMLIZ .AND. NDEG > 0) THEN
                     CALL READY (LUNCRT,
     >                  'Must the curve go through (0,0)? (<CR> = No) ',
     >                  LUNKBD, THRU00, DEFAULT, QUIT)
                     IF (QUIT) GO TO 710
                  END IF ! Else it isn't easy to arrange unscaled C (0) to be 0.

                  LEGEND = 'PNFIT: Degree = nn'
                  WRITE (LEGEND (17:18), '(I2)') NDEG
                  IF (THRU00) WRITE (LEGEND (19:24),' (''; C0=0'')')
                  WRITE (LUNLOG, 1100) TRIM (LEGEND)
               END IF

               CALL SECOND (TIMEA)

               CALL PNFIT (NX, X, Y, NDEG, THRU00, MXWORK, WORK, C (0),
     >                     RMSDEV, IER)
               IF (IER /= 0) GO TO 860

               CALL PNEVAL (NDEG, C, NX, X, YFIT)

               IF (PLTCRV) CALL PNEVAL (NDEG, C, NPLOT, XPLOT, YPLOT)

               IF (PLTABS) CALL PNEVAL (NDEG, C, NEVAL, XEVAL, YEVAL)

               CALL SECOND (TIMEB)

               CALL PRCOEF (0, NDEG, C, IDATASET, 'COEFS (0:NDEG):',
     >                      IOMIT)

               IF (NMLIZ) THEN
                  WRITE (LUNLOG, 1003)
     >               ' These coefficients apply to X and Y',
     >               ' transformed to [0,1].'
                  WRITE (LUNLOG, 1010)
     >               ' Shift and scale factors for X:  ',
     >               SHIFT (1), SCALE (1),
     >               ' Shift and scale factors for Y:  ',
     >               SHIFT (2), SCALE (2)
               END IF

            CASE (6)  ! *** FSFIT ***

C              Simple-minded linear-least squares fit of finite sum of
C              sin/cos terms.  Result is really just a pseudo-Fourier series,
C              but it applies to irregular data as well as regular data.

CCCCCC         IF (ARROW /= ONE) GO TO 830  ! Allow duplicate absscissas

  460          IF (FIRSTSET) THEN
                  LEASTSQRS = TRUE
                  WRITE (LUNCRT, 1001)
     >              ' CAUTION:  FSFIT is inefficient for big problems.',
     >              ' But it handles nonuniform data and may be handy.',
     >              BLANK
                  NFCREQ = -1
                  CALL READI (LUNCRT,
     >               'Enter N to estimate Fourier coefs. a, b (0:N): ',
     >               LUNKBD, NFCREQ, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710

                  IF (NFCREQ < 0 .OR. NFCREQ > MXFSFIT) GO TO 460
                  LEGEND = 'FSFIT:  N = nn'
                  WRITE (LEGEND (13:14), '(I2)') NFCREQ
                  WRITE (LUNLOG, 1100) TRIM (LEGEND)
               END IF

               NFCOF = MIN (NFCREQ, (NX - 2) / 2)
               IF (NFCOF < NFCREQ) WRITE (LUNCRT, 1004)
     >            'Warning: No. of coef. pairs is being reduced to ',
     >            NFCOF
                  
C ...          Fourier methods expect abscissas in the range [0,2*pi].
C              USESCALE was already in use with GETSCALE; it can also
C              be used in conjunction with GETXFORM.  (Replaces XFORMX.)

               CALL GETXFORM (XMIN, XMAX, ZERO, TWOPI, SCALE, SHIFT)
               CALL COPY (NX, X, XNORM)
               CALL USESCALE (DENRMLIZ, 1, NX, XNORM, XNORM, XNORM,
     >                        SCALE, SHIFT, IER)

               CALL SECOND (TIMEA)

               CALL FSFIT (NX, XNORM, Y, NFCOF, MXWORK, WORK, A, B,
     >                     RMSDEV, IER)
               IF (IER /= 0) GO TO 860

               CALL FSEVAL (3, NX, XNORM, NFCOF, A, B, YFIT)

               IF (PLTCRV) THEN
                  CALL COPY (NPLOT, XPLOT, XNORM)
                  CALL USESCALE (DENRMLIZ, 1, NPLOT, XNORM, XNORM,
     >                           XNORM, SCALE, SHIFT, IER)
                  CALL FSEVAL (3, NPLOT, XNORM, NFCOF, A, B, YPLOT)
               END IF

               IF (PLTABS) THEN
                  CALL COPY (NEVAL, XEVAL, XNORM)
                  CALL USESCALE (DENRMLIZ, 1, NEVAL, XNORM, XNORM,
     >                           XNORM, SCALE, SHIFT, IER)
                  CALL FSEVAL (3, NEVAL, XNORM, NFCOF, A, B, YEVAL)
               END IF

               CALL SECOND (TIMEB)

               CALL PRCOEF (0, NFCOF, A, IDATASET,'A (0:NFCOF):', IOMIT)
               CALL PRCOEF (0, NFCOF, B, IOMIT,   'B (0:NFCOF):', IOMIT)

            CASE (7)  ! *** FSERIES ***

C              Simple-minded Fourier-type analysis via trapezoidal-rule
C              integrations.  Allows for differing F(1) and F(NX), and
C              for two types of periodic continuation.  Applies to
C              irregular data as well as regular data.  (Written with
C              Macintosh version in mind.)

               IF (ARROW /= ONE) GO TO 830

               IF (FIRSTSET) THEN
                  LEASTSQRS = TRUE
                  WRITE (LUNCRT, 1001)
     >            ' CAUTION:  FSERIES is inefficient for big problems.',
     >            ' But it handles nonuniform data and may be handy.',
     >            BLANK
                  NFCOF = -1
                  DO WHILE (NFCOF < 0 .OR. NFCOF > MXFCOF)
                     CALL READI (LUNCRT,
     >                  'Enter N to apply Fourier coefs. a, b (0:N): ',
     >                  LUNKBD, NFCOF, DEFAULT, QUIT)
                     IF (QUIT) GO TO 710
                  END DO

                  YES = TRUE
                  CALL READY (LUNCRT,
     >               'Are the Ys arbitrary? (<CR>=Yes; No=periodic): ',
     >               LUNKBD, YES, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710
                  IF (YES) THEN
                     MODE = 1
                  ELSE
                     MODE = 3
                  END IF

                  LEGEND = 'FSERIES, mode m:  N = nn'
                  WRITE (LEGEND (15:15), '(I1)') MODE
                  WRITE (LEGEND (23:24), '(I2)') NFCOF
                  WRITE (LUNLOG, 1100) TRIM (LEGEND)
               END IF

               CALL SECOND (TIMEA)

               CALL FSERIES (MODE, NX, X, Y, WORK, NFCOF, A, B, IER)
               IF (IER /= 0) GO TO 860

               CALL FSSUM (MODE, NX, X, NX, X, Y, NFCOF, A, B, YFIT)

               IF (PLTCRV) THEN
                  CALL FSSUM (MODE, NPLOT, XPLOT, NX, X, Y, NFCOF, A, B,
     >                        YPLOT)
               END IF

               IF (PLTABS) THEN
                  CALL FSSUM (MODE, NEVAL, XEVAL, NX, X, Y, NFCOF, A, B,
     >                        YEVAL)
               END IF

               CALL SECOND (TIMEB)

               CALL PRCOEF (0, NFCOF, A, IDATASET,'A (0:NFCOF):', IOMIT)
               CALL PRCOEF (0, NFCOF, B, IOMIT,   'B (0:NFCOF):', IOMIT)

            CASE (8)  ! *** FSCOEF ***

C              One or more pairs of Fourier coefficients by efficient
C              recurrence relation:

  470          IF (ARROW /= ONE) GO TO 830

               LEASTSQRS = TRUE
               RMSDEV = ZERO

               IF (FIRSTSET) THEN
                  WRITE (LUNCRT, 1001)
     >               ' Note: First & last pts. are assumed to "match."',
     >               ' Also: coef. pairs are 0:N where N <= (NPTS-1)/2.'
                  NFCREQ = -1
                  CALL READI (LUNCRT,
     >               'Enter N to estimate Fourier coefs. a, b (0:N): ',
     >               LUNKBD, NFCREQ, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710

                  IF (NFCREQ < 0 .OR. NFCREQ > MXFCOF) GO TO 470
                  LEGEND = 'FSCOEF:  N = nnn'
                  WRITE (LEGEND (14:16), '(I3)') NFCREQ
                  WRITE (LUNLOG, 1100) TRIM (LEGEND)
               END IF

               M = NX - 1
               NFCOF = MIN (NFCREQ, (M - 1) / 2)
               IF (NFCOF < NFCREQ) WRITE (LUNCRT, 1004)
     >            'Warning: # coef. pairs is being reduced to ', NFCOF
                  
               CALL SECOND (TIMEA)

               DO J = 0, NFCOF
                  CALL FSCOEF (3, M, Y, J, A (J), B (J), IER)
                  IF (IER /= 0) GO TO 860
               END DO

C ...          Fourier methods expect abscissas in the range [0, 2*pi].
C              USESCALE was already in use with GETSCALE; it can also be
C              used in conjunction with GETXFORM.  (Replaces XFORMX.)

               CALL GETXFORM (X (1), X (NX), ZERO, TWOPI, SCALE, SHIFT)
               CALL COPY (NX, X, XNORM)
               CALL USESCALE (DENRMLIZ, 1, NX, XNORM, XNORM, XNORM,
     >                        SCALE, SHIFT, IER)

               CALL FSEVAL (3, NX, XNORM, NFCOF, A, B, YFIT)

               IF (PLTCRV) THEN
                  CALL COPY (NPLOT, XPLOT, XNORM)
                  CALL USESCALE (DENRMLIZ, 1, NPLOT, XNORM, XNORM,
     >                           XNORM, SCALE, SHIFT, IER)
                  CALL FSEVAL (3, NPLOT, XNORM, NFCOF, A, B, YPLOT)
               END IF

               IF (PLTABS) THEN
                  CALL COPY (NEVAL, XEVAL, XNORM)
                  CALL USESCALE (DENRMLIZ, 1, NEVAL, XNORM, XNORM,
     >                           XNORM, SCALE, SHIFT, IER)
                  CALL FSEVAL (3, NEVAL, XNORM, NFCOF, A, B, YEVAL)
               END IF

               CALL SECOND (TIMEB)

               CALL PRCOEF (0, NFCOF, A, IDATASET,'A (0:NFCOF):', IOMIT)
               CALL PRCOEF (0, NFCOF, B, IOMIT,   'B (0:NFCOF):', IOMIT)

            CASE (9)  ! *** QHSFIT ***

C              Interpolation by quasi-Hermite spline:

  480          IF (ARROW == ZERO) GO TO 840

               CYCLIC = FALSE
               IF (Y (NX) == Y (1)) THEN
                  CALL READY (LUNCRT,
     >             'Impose periodic boundary conditions? (<CR> = No): ',
     >             LUNKBD, CYCLIC, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710
               END IF

               IF (CYCLIC) THEN
                  METHOD = 'C'
                  LEGEND = 'QHSFIT (periodic)'
               ELSE
                  IF (FIRSTSET) WRITE (LUNCRT, 1003)
     >               ' Akima''s method uses parabolic extrapolation',
     >               ' to add two points off each end.',
     >               ' Normally this suffices, but wild gradients',
     >               ' are possible.',
     >               ' The available monotonic option for the end',
     >               ' intervals may be preferable.'
                  REPLY = 'P'
                  CALL READC (LUNCRT,
     >        'Parabolic end conditions, or Monotonic? (P|M; <CR>=P): ',
     >               LUNKBD, REPLY, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710

                  IF (REPLY == 'P') THEN
                     LEGEND = 'QHSFIT'
                  ELSE IF (REPLY == 'M') THEN
                     LEGEND = 'QHSFIT (monotonic end conditions)'
                  ELSE
                     GO TO 480
                  END IF
                  METHOD = REPLY
               END IF
 
               WRITE (LUNLOG, 1100) TRIM (LEGEND)

               CALL SECOND (TIMEA)

               CALL QHSFIT (NX, X, Y, METHOD,
     >                      COEFS (1,1), COEFS (1,2), COEFS (1,3), IER)
               IF (IER /= 0) GO TO 860

               IF (PLTCRV)
     >            CALL CSEVAL (NX, X, Y, NPLOT, XPLOT, COEFS (1,1),
     >                         COEFS (1,2), COEFS (1,3), YPLOT)
               IF (PLTABS)
     >            CALL CSEVAL (NX, X, Y, NEVAL, XEVAL, COEFS (1,1),
     >                         COEFS (1,2), COEFS (1,3), YEVAL)

               CALL SECOND (TIMEB)

               CALL PRCOEF (1, NX-1, COEFS (1,1), IDATASET,
     >                      'COEFS (1:NX-1,1):',  IOMIT)
               CALL PRCOEF (1, NX-1, COEFS (1,2), IOMIT,
     >                      'COEFS (1:NX-1,2):',  IOMIT)
               CALL PRCOEF (1, NX-1, COEFS (1,3), IOMIT,
     >                      'COEFS (1:NX-1,3):',  IOMIT)

            CASE (10)  ! *** CSFIT ***

C              Cubic interpolating spline with variety of end conditions
C              and proper handling of degenerate cases:

  490          IF (ARROW == ZERO) GO TO 840

               IF (FIRSTSET) THEN
                  WRITE (LUNCRT, 1000) BLANK,
     >            ' Description of spline end conditions:',
     >            '   0: 3rd derivative of spline at endpoint',
     >            '      is to match 3rd derivative of cubic',
     >            '      passing through first/last four data points.',
     >            '   1: 1st derivative of spline at endpoint',
     >            '      is to be user-supplied.',
     >            '   2: 2nd derivative of spline at endpoint',
     >            '      is to be user-supplied.',
     >            '   3: 3rd derivative of spline at endpoint',
     >            '      is to be user-supplied.',
     >            '   4: Cylic, or periodic, case in which spline',
     >            '      and its first three derivatives at X(N)',
     >            '      are to match values at X(1).', BLANK

                  IENDL = 0
                  CALL READI (LUNCRT,
     >               'Choose the  left end condition (0-4; <CR>=0): ',
     >               LUNKBD, IENDL, DEFAULT, QUIT)
                     IF (QUIT) GO TO 710

                  PERIODIC = IENDL == 4
                  IF (IENDL > 0 .AND. .NOT. PERIODIC) THEN
                     DERIVL = ZERO
                     CALL READR (LUNCRT,
     >                  'Enter the value of the derivative there: ',
     >                  LUNKBD, DERIVL, DEFAULT, QUIT)
                     IF (QUIT) GO TO 710
                  END IF

                  IF (PERIODIC) THEN
                     LEGEND = 'CSFIT:  Periodic Case'
                  ELSE
                     IENDR = 0
                     CALL READI (LUNCRT,
     >                 'Choose the right end condition (0-3; <CR>=0): ',
     >                  LUNKBD, IENDR, DEFAULT, QUIT)
                     IF (QUIT) GO TO 710

                     IF (IENDR > 0) THEN
                        DERIVR = ZERO
                        CALL READR (LUNCRT,
     >                     'Enter the value of the derivative there: ',
     >                     LUNKBD, DERIVR, DEFAULT, QUIT)
                        IF (QUIT) GO TO 710
                     END IF

                     LEGEND = 'CSFIT: IENDL=n, IENDR=n'
                     WRITE (LEGEND (14:14), ' (I1)') IENDL
                     WRITE (LEGEND (23:23), ' (I1)') IENDR
                  END IF

                  WRITE (LUNLOG, 1100) TRIM (LEGEND)
               END IF

               IF (PERIODIC) THEN

C ...             Ordinate of last point should match that of first point.
C                 (If it doesn't, we can't insert such a point because we
C                 don't know its abscissa.)

                  IF (Y (NX) /= Y (1)) THEN
                     WRITE (LUNCRT, 1001)
     >               ' Illegal request: End-pt. ordinates do not match.'
                     FIRSTSET = TRUE
                     GO TO 490
                  END IF

                  NDATA = -NX  ! Flag for CSEVAL
               ELSE
                  NDATA =  NX
               END IF

               CALL SECOND (TIMEA)

               CALL CSFIT (NX, X, Y, IENDL, DERIVL, IENDR, DERIVR,
     >                     COEFS (1,1), COEFS (1,2), COEFS (1,3), IER)
               IF (IER /= 0) GO TO 860

               IF (PLTCRV)
     >            CALL CSEVAL (NDATA, X, Y, NPLOT, XPLOT, COEFS (1,1),
     >                         COEFS (1,2), COEFS (1,3), YPLOT)
               IF (PLTABS)
     >            CALL CSEVAL (NDATA, X, Y, NEVAL, XEVAL, COEFS (1,1),
     >                         COEFS (1,2), COEFS (1,3), YEVAL)

               CALL SECOND (TIMEB)

               CALL PRCOEF (1, NX, COEFS (1,1), IDATASET,
     >                      'COEFS (1:NX,1):',  IOMIT)
               CALL PRCOEF (1, NX, COEFS (1,2), IOMIT,
     >                      'COEFS (1:NX,2):',  IOMIT)
               CALL PRCOEF (1, NX, COEFS (1,3), IOMIT,
     >                      'COEFS (1:NX,3):',  IOMIT)

            CASE (11)  ! *** TABLE1 or LINE1D ***

C              Linear interpolation/extrapolation by 1-D table look-up.
C              We have FUNCTION and SUBROUTINE forms - allow for both,
C              if for no other reason than to test them.

               IF (ARROW == ZERO) GO TO 840
               
               IF (NX < 20 .OR. NEVAL < 20) THEN
                  LEGEND = 'TABLE1'
               ELSE
                  LEGEND = 'LINE1D'
               END IF

               CALL SECOND (TIMEA)

               IF (PLTCRV) THEN
                  IF (LEGEND (1:6) == 'TABLE1') THEN
                     INDEX = 1
                     DO I = 1, NPLOT
                        YPLOT (I) =
     >                     TABLE1 (NX, X, Y, INDEX, XPLOT (I), IER)
                        IF (IER /= 0) GO TO 860
                     END DO
                  ELSE
                     CALL LINE1D (NX, X, Y, NPLOT, XPLOT, YPLOT)
                  END IF
               END IF

               IF (PLTABS) THEN
                  IF (LEGEND (1:6) == 'TABLE1') THEN
                     INDEX = 1
                     DO I = 1, NEVAL
                        YEVAL (I) =
     >                     TABLE1 (NX, X, Y, INDEX, XEVAL (I), IER)
                        IF (IER /= 0) GO TO 860
                     END DO
                  ELSE
                     CALL LINE1D (NX, X, Y, NEVAL, XEVAL, YEVAL)
                  END IF
               END IF

               CALL SECOND (TIMEB)

            CASE (12)  ! *** PSFIT ***

C              Parametric interpolating spline (open or closed curves).

               IF (.NOT. PLTCRV) THEN
                  WRITE (LUNCRT, 1001)
     >               ' No plotted points have been requested!'
                  GO TO 200
               END IF

               IF (FIRSTSET) THEN
                  MONO = TRUE
                  WRITE (LUNCRT, 1003) ' Splines for X vs T, Y vs T',
     >               ' may be conventional or monotonic.'
                  CALL READY (LUNCRT,
     >               'Do you want conventional? (Y|N; <CR>=Yes) ',
     >               LUNKBD, MONO, DEFAULT, QUIT)
                  MONO = .NOT. MONO
                  IF (QUIT) GO TO 710

                  CLOSED = FALSE
                  CALL READY (LUNCRT,
     >               'Is curve to be closed smoothly? (Y|N; <CR>=No) ',
     >               LUNKBD, CLOSED, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710
               END IF

               AUGMENTED = FALSE
               IF (CLOSED .AND.
     >            (X (1) /= X (NX) .OR. Y (1) /= Y (NX))) THEN

C ...             Temporarily "close" the data:

                  NX = NX + 1
                  X (NX) = X (1)
                  Y (NX) = Y (1)
                  AUGMENTED = TRUE
               END IF

               IF (.NOT. DISTINCT) GO TO 850

               CALL SECOND (TIMEA)

               IF (MONO) THEN
                  LEGEND = 'PSFIT/MSFIT'
                  CALL PSFIT (NX, X, Y, 'M', CLOSED, IER)
               ELSE
                  LEGEND = 'PSFIT/CSFIT'
                  CALL PSFIT (NX, X, Y, 'C', CLOSED, IER)
               END IF
               IF (IER /= 0) GO TO 860

               IF (PLTCRV) THEN

C ...             Generate uniform points along the arc length:

                  CALL XGRID (NPLOT, 0, TSUBJ (1), TSUBJ (NX), XEVAL)

                  CALL PSTVAL (NPLOT, XEVAL, XPLOT, YPLOT, XPLOT, YPLOT,
     >                         XPLOT, YPLOT, NX, X, Y)
               END IF

               NOABS = TRUE

               IF (AUGMENTED) NX = NX - 1

               CALL SECOND (TIMEB)

            CASE (13)  ! *** WAGFIT ***

C              Wagner functions are appropriate for airfoil geometries,
C              but may find other applications too.
C              The first N are fitted in the least squares sense.

C              The Xs must be on [0,1], either originally, or by normalization:

               IF (ARROW == ZERO) GO TO 840

               IF (ABS (X (1)) > EPS .OR.
     >             ABS (X (NX) - ONE) > EPS) GO TO 870

               IF (FIRSTSET) THEN
                  LEASTSQRS = TRUE
                  NDEG = 0
                  DO WHILE (NDEG < 1 .OR. NDEG > MXWAG)
                     NDEG = 6
                     CALL READI (LUNCRT,
     >               'To fit Wagner functions 1:N, enter N. <CR> = 6: ',
     >                  LUNKBD, NDEG, DEFAULT, QUIT)
                     IF (QUIT) GO TO 710
                  END DO
                  LEGEND = 'WAGFIT:  N = nn'
                  WRITE (LEGEND (14:15), '(I2)') NDEG
                  WRITE (LUNLOG, 1100) TRIM (LEGEND)
               END IF

C ...          Wagner functions are zero at X = 0. and 1.

               Y1 = Y (1)
               CALL SHIFTX (NX, Y, -Y1)

C ...          1st pt. is now [0, 0], but last pt. may not have zero ordinate.
C              To preserve a thick "trailing edge" exactly, a "ramp"
C              function must be subtracted first then added back in:

               YNNORM = -Y (NX)
               RAMP = YNNORM /= ZERO

               IF (RAMP) THEN

C ...             Subtract "ramp" function (in-place):

                  ADD = TRUE

                  CALL BEVAL ('RAMP', 1, YNNORM, ADD, NX, X, Y)

                  YNNORM = -YNNORM
               END IF
                              
               CALL SECOND (TIMEA)

C ...          Set up and solve the linear least squares problem:

               CALL WAGFIT (NX, X, Y, NDEG, MXWORK, WORK, C (1),
     >                      RMSDEV, IER)
               IF (IER /= 0) GO TO 860
 
               IF (RAMP) THEN

C ...             Add back the "ramp" function (in-place):

                  CALL BEVAL ('RAMP', 1, YNNORM, ADD, NX, X, Y)

               END IF
                              
C ...          Restore the vertical shift:

               CALL SHIFTX (NX, Y, Y1)

               ADD = FALSE
               DO J = 1, NDEG
                  B (1) = J
                  B (2) = C (J)

                  CALL BEVAL ('WAGNER', 2, B (1), ADD, NX, X, YFIT)

                  ADD = TRUE
               END DO

               IF (RAMP) THEN
                  CALL BEVAL ('RAMP', 1, YNNORM, ADD, NX, X, YFIT)
               END IF

               CALL SHIFTX (NX, YFIT, Y1)

               IF (PLTCRV) THEN

C                 Abscissas can be assumed to be normalized...

                  ADD = FALSE
                  DO J = 1, NDEG
                     B (1) = J
                     B (2) = C (J)

                     CALL BEVAL ('WAGNER',
     >                           2, B (1), ADD, NPLOT, XPLOT, YPLOT)
                     ADD = TRUE
                  END DO

                  IF (RAMP) THEN
                     CALL BEVAL ('RAMP',
     >                           1, YNNORM, ADD, NPLOT, XPLOT, YPLOT)
                  END IF

                  CALL SHIFTX (NPLOT, YPLOT, Y1)

               END IF

               IF (PLTABS) THEN

                  ADD = FALSE
                  DO J = 1, NDEG
                     B (1) = J
                     B (2) = C (J)

                     CALL BEVAL ('WAGNER',
     >                           2, B (1), ADD, NEVAL, XEVAL, YEVAL)
                     ADD = TRUE
                  END DO

                  IF (RAMP) THEN
                     CALL BEVAL ('RAMP',
     >                           1, YNNORM, ADD, NEVAL, XEVAL, YEVAL)
                  END IF

                  CALL SHIFTX (NEVAL, YEVAL, Y1)

               END IF

               CALL SECOND (TIMEB)
 
               CALL PRCOEF (1, NDEG, C (1), IDATASET, 'COEFS (1:N):',
     >                      IOMIT)

            CASE (14)  ! *** VECFIT ***

C              Linear combination of discrete vectors, which are read
C              from files prompted for interactively.  No point in
C              trying to deal with more than one dataset in one run.
C              Evaluating at other than the original points is not
C              provided for.

               IF (.NOT. FIRSTSET) THEN
                  WRITE (LUNCRT, 1001)
     >               ' VECFIT deals with just one dataset in a file.'
                  GO TO 200
               END IF

               LEASTSQRS = TRUE
               LEGEND = 'VECFIT: NC ='

               CALL SECOND (TIMEA)

               CALL VECFIT (NX, X, Y, LUNCRT, LUNKBD, LUNXGR, LUNLOG,
     >                      WORK, MXVEC, NDEG, C (1), YFIT, RMSDEV,
     >                      IER)
               IF (IER /= 0) GO TO 860

               IF (PLTCRV) THEN
C ...             Not done. Use WORK (*) and C (*) as in VECFIT if needed.
               END IF

               IF (PLTABS) THEN
C ...             Already evaluated at original Xs in VECFIT.
               END IF
 
               CALL SECOND (TIMEB)

               WRITE (LEGEND (13:15), '(I3)') NDEG
               WRITE (LUNLOG, 1100) TRIM (LEGEND)

               CALL PRCOEF (1, NDEG, C(1), IDATASET,
     >                      'Linear coefficients:', IOMIT)

            CASE (15)  ! *** CURV1 ***

C              Cubic interpolating spline under tension, with optional
C              first derivative end conditions (NCAR, 1972):

               IF (ARROW /= ONE) GO TO 830

               IF (FIRSTSET) THEN

                  TENSE = ONE
                  CALL READR (LUNCRT,
     >               'Enter tension factor (0. to 50.; <CR>=1.): ',
     >               LUNKBD, TENSE, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710

                  LEGEND = 'CURV1:  Tension = xx.xxx'
                  WRITE (LEGEND (19:24), '(F6.3)') TENSE

                  TENSE = -TENSE
                  YES = FALSE
                  CALL READY (LUNCRT,
     >  'Do you want to enter end point 1st derivatives? (<CR> = No): ',
     >               LUNKBD, YES, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710

                  IF (YES) THEN
                     TENSE = -TENSE
                     DERIVL = ZERO
                     CALL READR (LUNCRT,
     >                  'Enter 1st derivative at X (1). (<CR> = 0.): ',
     >                  LUNKBD, DERIVL, DEFAULT, QUIT)
                     IF (QUIT) GO TO 710

                     DERIVR = ZERO
                     CALL READR (LUNCRT,
     >                  'Enter 1st derivative at X (N). (<CR> = 0.): ',
     >                  LUNKBD, DERIVR, DEFAULT, QUIT)
                     IF (QUIT) GO TO 710
                  END IF

                  WRITE (LUNLOG, 1100) TRIM (LEGEND)
               END IF

               CALL SECOND (TIMEA)

               CALL CURV1 (NX, X, Y, DERIVL, DERIVR,
     >                     COEFS (1, 1), COEFS (1, 2), TENSE)

               IF (PLTCRV) THEN
                  INDEX = 1
                  DO I = 1, NPLOT
                     YPLOT (I) = CURV2 (XPLOT (I), NX, X, Y,
     >                                  COEFS (1,1), TENSE, INDEX)
                  END DO
               END IF

               IF (PLTABS) THEN
                  INDEX = 1
                  DO I = 1, NEVAL
                     YEVAL (I) = CURV2 (XEVAL (I), NX, X, Y,
     >                                  COEFS (1,1), TENSE, INDEX)
                  END DO
               END IF

               CALL SECOND (TIMEB)

               CALL PRCOEF (1, NX, COEFS (1,1), IDATASET,
     >                      'YP (1:NX):', IOMIT)

            CASE (16)  ! *** MSFIT ***

C              Monotonic cubic interpolating spline:

  550          IF (ARROW == ZERO) GO TO 840

               IF (FIRSTSET) THEN
                  LEGEND = 'MSFIT'
                  WRITE (LUNLOG, 1100) TRIM (LEGEND)
               END IF

               CYCLIC = FALSE
               CALL READY (LUNCRT,
     >            'Impose periodic boundary conditions? (<CR> = No): ',
     >            LUNKBD, CYCLIC, DEFAULT, QUIT)
               IF (QUIT) GO TO 710

               IF (CYCLIC) THEN
                  IF (Y (NX) /= Y (1)) THEN
                     WRITE (LUNCRT, 1001)
     >               ' Illegal request: End-pt. ordinates do not match.'
                     GO TO 550
                  END IF
               END IF

               CALL SECOND (TIMEA)

               CALL MSFIT (NX, X, Y, CYCLIC, COEFS (1,1), COEFS (1,2),
     >                     COEFS (1,3), IER)
               IF (IER /= 0) GO TO 860

               IF (PLTCRV)
     >            CALL CSEVAL (NX, X, Y, NPLOT, XPLOT, COEFS (1,1),
     >                         COEFS (1,2), COEFS (1,3), YPLOT)
               IF (PLTABS)
     >            CALL CSEVAL (NX, X, Y, NEVAL, XEVAL, COEFS (1,1),
     >                         COEFS (1,2), COEFS (1,3), YEVAL)

               CALL SECOND (TIMEB)

               CALL PRCOEF (1, NX, COEFS (1,1), IDATASET,
     >                      'COEFS (1:NX,1):',  IOMIT)
               CALL PRCOEF (1, NX, COEFS (1,2), IOMIT,
     >                      'COEFS (1:NX,2):',  IOMIT)
               CALL PRCOEF (1, NX, COEFS (1,3), IOMIT,
     >                      'COEFS (1:NX,3):',  IOMIT)

            CASE (17)  ! *** PLSFIT ***

C              Parametric local interpolating spline (open or closed curves).

               IF (.NOT. PLTCRV) THEN
                  WRITE (LUNCRT, 1001)
     >               ' No plotted points have been requested!'
                  GO TO 200
               END IF

  560          IF (FIRSTSET) THEN
                  DSTRIB = 'UNIFORM'
                  REPLY = 'T'
                  CALL READC (LUNCRT,
     >               'Tight or loose fit? (T|L; <CR>=Tight) ',
     >               LUNKBD, REPLY, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710

                  IF (REPLY == 'T') THEN
                     LEGEND = 'PLSFIT (tight)'
                     METHOD = 'MONOTONE'
                  ELSE
                     LEGEND = 'PLSFIT (loose)'
                     METHOD = 'BESSEL'
                  END IF

                  CLOSED = FALSE
                  CALL READY (LUNCRT,
     >               'Is curve to be closed smoothly? (Y|N; <CR>=No): ',
     >               LUNKBD, CLOSED, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710

                  WRITE (LUNLOG, 1100) TRIM (LEGEND)

               END IF

               AUGMENTED = FALSE
               IF (CLOSED .AND.
     >            (X (1) /= X (NX) .OR. Y (1) /= Y (NX))) THEN

C ...             Temporarily "close" the data:

                  IF (NX < MXPTS) THEN
                     NX = NX + 1
                     X (NX) = X (1)
                     Y (NX) = Y (1)
                     AUGMENTED = TRUE
                  ELSE
                     WRITE (LUNCRT, 1001)
     >                  ' Not enough space for periodic case.'
                     IF (FIRSTSET) GO TO 560
                     GO TO 300     ! Else skip this dataset
                  END IF
               END IF

               IF (.NOT. DISTINCT) GO TO 850

               CALL SECOND (TIMEA)

               TBEGIN = -ONE
               TEND = -ONE

               CALL PLSFIT (NX, X, Y, TBEGIN, TEND, NPLOT, XPLOT,
     >                      YPLOT, TRUE, CLOSED, METHOD, DSTRIB, IER)
               IF (IER /= 0) GO TO 860

               NOABS = TRUE

               IF (AUGMENTED) NX = NX - 1

               CALL SECOND (TIMEB)

            CASE (18)  ! *** LCSFIT ***

C              Nonparametric form of PLSFIT (local cubic spline method).

               IF (ARROW == ZERO) GO TO 840

               CYCLIC = FALSE
               IF (Y (NX) == Y (1)) THEN
                  CALL READY (LUNCRT,
     >             'Impose periodic boundary conditions? (<CR> = No): ',
     >             LUNKBD, CYCLIC, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710
               END IF

               IF (CYCLIC) THEN
                  METHOD = 'CYCLIC'
                  LEGEND = 'LCSFIT (periodic)'
               ELSE
                  IF (FIRSTSET) THEN
                     REPLY = 'T'
                     CALL READC (LUNCRT,
     >'Loose fit? Tight (monotonic)? Straight lines? (L|T|S; <CR>=T): ',
     >                  LUNKBD, REPLY, DEFAULT, QUIT)
                     IF (QUIT) GO TO 710

                     IF (REPLY == 'T') THEN
                        LEGEND = 'LCSFIT (tight)'
                        METHOD = 'MONOTONE'
                     ELSE IF (REPLY == 'L') THEN
                        LEGEND = 'LCSFIT (loose)'
                        METHOD = 'BESSEL'
                     ELSE
                        LEGEND = 'LCSFIT (linear)'
                        METHOD = 'LINEAR'
                     END IF
                  END IF
               END IF

               CALL SECOND (TIMEA)

               IF (PLTCRV)
     >            CALL LCSFIT (NX, X, Y, TRUE, METHOD, NPLOT, XPLOT,
     >                         YPLOT, YPLOT)
               IF (PLTABS)
     >            CALL LCSFIT (NX, X, Y, TRUE, METHOD, NEVAL, XEVAL,
     >                         YEVAL, YEVAL)

               CALL SECOND (TIMEB)

            CASE (19)  ! *** LSFIT1 ***

C              Local least squares / spline hybrid method.

               IF (ARROW == ZERO) GO TO 840

               IF (FIRSTSET) THEN
                  LEASTSQRS = TRUE
                  CLOSED = FALSE
                  REPLY = 'H'
                  WRITE (LUNCRT, 1003)
     >               ' LSFIT1 first smooths the ordinates F(X) to',
     >               ' give G(X) at the original abscissas.',
     >               ' It can then spline G in two ways, for',
     >               ' interpolation, or suppress the spline',
     >               ' and return G only (no interpolation).', BLANK
                  CALL READC (LUNCRT,
     > 'Hermite or Monotone spline?  Or G only? (H|M|G; <CR>=Hermite) ',
     >               LUNKBD, REPLY, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710

                  IF (REPLY == 'M') THEN
                     LEGEND = 'LSFIT1 (m)'
                     METHOD = 'MONOTONE'
                  ELSE IF (REPLY == 'H') THEN
                     LEGEND = 'LSFIT1'
                     METHOD = 'HERMITE'
                  ELSE
                     LEGEND = 'LSFIT1'
                     METHOD = 'Gonly'  ! Only first character is used
                  END IF

  582             NDEG = 2
                  CALL READI (LUNCRT,
     >  'Degree of local least squares polynomial fits? (1:4; <CR>=2) ',
     >                        LUNKBD, NDEG, DEFAULT, QUIT)
                  IF (NDEG < 1 .OR. NDEG > 4) GO TO 582

  584             NEARPTS = NDEG + 4
                  WRITE (LUNCRT, 1003)
     >               ' The number of pts. in each local fit should be',
     >               ' at least DEGREE + 2',
     >               ' to achieve any smoothing at all.'
                  CALL READI (LUNCRT,
     >'Number of pts. in each local fit? (DEG+2 : 100; <CR>=DEG+4) ',
     >               LUNKBD, NEARPTS, DEFAULT, QUIT)
                  IF (NEARPTS < NDEG + 2 .OR. NEARPTS > MXLSFIT1)
     >               GO TO 584                  ! LSFIT1 has 100 hard-coded

                  WRITE (LEGEND (12:17), 1011) 'Deg: ', NDEG
                  IF (NEARPTS < 10) THEN
                     WRITE (LEGEND (19:24), 1011) 'Pts: ', NEARPTS
                  ELSE
                     WRITE (LEGEND (19:25), 1012) 'Pts: ', NEARPTS
                  END IF

                  CYCLIC = FALSE
                  CALL READY (LUNCRT,
     >               'Impose periodic boundary conditions? (<CR>=No) ',
     >               LUNKBD, CYCLIC, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710

                  WRITE (LUNLOG, 1100) TRIM (LEGEND)

               END IF

               CALL SECOND (TIMEA)

               CALL LSFIT1 (NX, X, Y, NDEG, NEARPTS, METHOD, CYCLIC,
     >                      YFIT, COEFS (1,1), COEFS (1,2), COEFS (1,3),
     >                      IER)
               IF (IER /= 0) GO TO 860

               IF (METHOD (1:1) /= 'G') THEN

                  IF (PLTCRV)
     >               CALL CSEVAL (NX, X, YFIT, NPLOT, XPLOT, COEFS,
     >                            COEFS (1,2), COEFS (1,3), YPLOT)
                  IF (PLTABS)
     >               CALL CSEVAL (NX, X, YFIT, NEVAL, XEVAL, COEFS,
     >                            COEFS (1,2), COEFS (1,3), YEVAL)

                  CALL SECOND (TIMEB)

                  CALL PRCOEF (1, NX, COEFS (1,1), IDATASET,
     >                         'COEFS (1:NX,1):',  IOMIT)
                  CALL PRCOEF (1, NX, COEFS (1,2), IOMIT,
     >                         'COEFS (1:NX,2):',  IOMIT)
                  CALL PRCOEF (1, NX, COEFS (1,3), IOMIT,
     >                         'COEFS (1:NX,3):',  IOMIT)

               ELSE ! There is no fitted curve to evaluate.

                  CALL SECOND (TIMEB)

                  NOFIT  = TRUE
                  NOABS  = TRUE
                  NTPLOT = NX

                  DO I = 1, NX
                     XPLOT (I) = X (I)
                     YPLOT (I) = YFIT (I)
                  END DO

               END IF

            CASE (20)  ! *** QINTERP ***

C              Quadratic interpolation (weighted average where possible).

               IF (ARROW == ZERO) GO TO 840

               IF (FIRSTSET) THEN
                  LEGEND = 'QINTERP'
                  WRITE (LUNLOG, 1100) TRIM (LEGEND)
               END IF

               CALL SECOND (TIMEA)

               IF (PLTCRV) THEN

                  LEFT = 1

                  CALL QINTERP (NX, X, Y, NPLOT, XPLOT, YPLOT, LEFT)

               END IF

               IF (PLTABS) THEN

                  LEFT = 1

                  CALL QINTERP (NX, X, Y, NEVAL, XEVAL, YEVAL, LEFT)

               END IF

               CALL SECOND (TIMEB)

            CASE (21)  ! *** AEBXPLUSC and AXBECX ***

C           Exponential curve fits as for chemical reaction rates.

               IF (FIRSTSET) THEN
                  LEASTSQRS = TRUE
                  CLOSED    = FALSE
                  ARRHENIUS = TRUE
                  CALL READY (LUNCRT,
     >               'Suppress the x^B term? (Y|N; <CR>=Yes) ',
     >               LUNKBD, ARRHENIUS, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710

                  MODE = 0
                  DO WHILE (MODE < 1 .OR. MODE > 2)
                     MODE = 2
                     CALL READI (LUNCRT,
     >            'Mode 1 [e^(Bx)] or mode 2 [e^(B/x)]? (1|2; <CR>=2) ',
     >                  LUNKBD, MODE, DEFAULT, QUIT)
                  END DO
                  IF (QUIT) GO TO 710

                  IF (ARRHENIUS) THEN
                     CE  = ZERO
                     YES = FALSE
                     CALL READY (LUNCRT,
     >                  'Include the + C term? (y|n; <CR>=no) ',
     >                  LUNKBD, YES, DEFAULT, QUIT)
                     IF (QUIT) GO TO 710
                     IF (YES) CE = ONE

                     IF (MODE == 1) THEN
                        LEGEND = 'Ae^(Bx) + C'
                        IF (.NOT. YES) LEGEND(8:11) = BLANK
                     ELSE
                        LEGEND = 'Ae^(B/x) + C'
                        IF (.NOT. YES) LEGEND(9:12) = BLANK
                     END IF
                  ELSE
                     IF (MODE == 1) THEN
                        LEGEND = 'Ax^Be(Cx)'
                     ELSE
                        LEGEND = 'Ax^Be(C/x)'
                     END IF
                  END IF

                  SHOWITERS = TRUE
                  CALL READY (LUNCRT,
     >               'Show iterations? (y|n; <CR>=yes) ',
     >               LUNKBD, SHOWITERS, DEFAULT, QUIT)
               END IF

               IER = 1
               IF (.NOT. SHOWITERS) IER = 0

               CALL SECOND (TIMEA)

               IF (ARRHENIUS) THEN

                  CALL AeBxplusC (MODE, NX, X, Y, AE, BE, CE, SSQMIN,
     >                            IER)

                  WRITE (LUNCRT, '(A, I3, 4ES17.8)')
     >              'IER, A, B, C, SSQMIN:', IER, AE, BE, CE, SSQMIN

                  IF (IER /= 0) GO TO 860

                  DO I = 1, NX
                     IF (MODE == 1) THEN
                        YFIT (I) = AE * EXP (BE * X (I)) + CE
                     ELSE
                        YFIT (I) = AE * EXP (BE / X (I)) + CE
                     END IF
                  END DO

                  IF (PLTCRV) THEN
                     DO I = 1, NPLOT
                        IF (MODE == 1) THEN
                           YPLOT (I) = AE * EXP (BE * XPLOT (I)) + CE
                        ELSE
                           YPLOT (I) = AE * EXP (BE / XPLOT (I)) + CE
                        END IF
                     END DO
                  END IF

                  IF (PLTABS) THEN
                     DO I = 1, NEVAL
                        IF (MODE == 1) THEN
                           YEVAL (I) = AE * EXP (BE * XEVAL (I)) + CE
                        ELSE
                           YEVAL (I) = AE * EXP (BE / XEVAL (I)) + CE
                        END IF
                     END DO
                  END IF

               ELSE

                  CALL AxBeCx (MODE, NX, X, Y, AE, BE, CE, SSQMIN, IER)

                  WRITE (LUNCRT, '(A, I3, 4ES17.8)')
     >               'IER, A, B, C, SSQMIN:', IER, AE, BE, CE, SSQMIN

                  IF (IER /= 0) GO TO 860

                  DO I = 1, NX
                     IF (MODE == 1) THEN
                        YFIT (I) = AE * X (I)**BE * EXP (CE * X (I))
                     ELSE
                        YFIT (I) = AE * X (I)**BE * EXP (CE / X (I))
                     END IF
                  END DO

                  IF (PLTCRV) THEN
                     DO I = 1, NPLOT
                        IF (MODE == 1) THEN
                           YPLOT (I) =
     >                        AE * XPLOT (I)**BE * EXP (CE * XPLOT (I))
                        ELSE
                           YPLOT (I) =
     >                        AE * XPLOT (I)**BE * EXP (CE / XPLOT (I))
                        END IF
                     END DO
                  END IF

                  IF (PLTABS) THEN
                     DO I = 1, NEVAL
                        IF (MODE == 1) THEN
                           YEVAL (I) =
     >                        AE * XEVAL (I)**BE * EXP (CE * XEVAL (I))
                        ELSE
                           YEVAL (I) =
     >                        AE * XEVAL (I)**BE * EXP (CE / XEVAL (I))
                        END IF
                     END DO
                  END IF

               END IF

               CALL SECOND (TIMEB)

               C (1) = AE
               C (2) = BE
               C (3) = CE

               CALL PRCOEF (1, 3, C (1), IDATASET, 'A, B, C:', IOMIT)

            CASE (22)  ! *** AXBPLUSC ***

C              Power law curve fit as for certain viscous flow properties.

               IF (FIRSTSET) THEN
                  LEASTSQRS = TRUE
                  CLOSED    = FALSE
                  CE  = ZERO
                  YES = FALSE
                  CALL READY (LUNCRT,
     >               'Include the + C term? (y|n; <CR>=no) ',
     >               LUNKBD, YES, DEFAULT, QUIT)
                  IF (QUIT) GO TO 710
                  IF (YES) CE = ONE

                  SHOWITERS = TRUE
                  CALL READY (LUNCRT,
     >               'Show iterations? (y|n; <CR>=yes) ',
     >               LUNKBD, SHOWITERS, DEFAULT, QUIT)
               END IF

               IER = 1
               IF (.NOT. SHOWITERS) IER = 0

               CALL SECOND (TIMEA)

               CALL AxBplusC (NX, X, Y, AE, BE, CE, SSQMIN, IER)

               WRITE (LUNCRT, '(A, I3, 4ES17.8)')
     >            'IER, A, B, C, SSQMIN:', IER, AE, BE, CE, SSQMIN

               IF (IER /= 0) GO TO 860

               DO I = 1, NX
                  YFIT (I) = AE * X (I)**BE + CE
               END DO

               IF (PLTCRV) THEN
                  DO I = 1, NPLOT
                     YPLOT (I) = AE * XPLOT (I)**BE + CE
                  END DO
               END IF

               IF (PLTABS) THEN
                  DO I = 1, NEVAL
                     YEVAL (I) = AE * XEVAL (I)**BE + CE
                  END DO
               END IF

               CALL SECOND (TIMEB)

               C (1) = AE
               C (2) = BE
               C (3) = CE

               CALL PRCOEF (1, 3, C (1), IDATASET, 'A, B, C:', IOMIT)

            CASE (23)  ! *** LEGENDRE ***

C              Linear combination of Legendre polynomials 0:n.

               IF (FIRSTSET) THEN
                  NDEG = 0
                  DO WHILE (NDEG < 1 .OR. NDEG > MXLEGENDRE)
                     NDEG = 6
                     CALL READI (LUNCRT,
     >           'To fit Legendre polynomials 0:N, enter N. <CR> = 6: ',
     >                           LUNKBD, NDEG, DEFAULT, QUIT)
                     IF (QUIT) GO TO 710
                  END DO
                  LEGEND = 'LEGENDRE:  N = nn'
                  WRITE (LEGEND (16:17), '(I2)') NDEG
                  WRITE (LUNLOG, 1100) TRIM (LEGEND)

                  LEASTSQRS = TRUE
                  CLOSED    = FALSE
               END IF

               CALL SECOND (TIMEA)

               CALL LEGENDRE (TRUE, NX, NDEG, X, Y, LEGENDRE_COEFS,
     >                        SSQMIN, TRUE, NEVAL, XEVAL, YEVAL, IER)

               WRITE (LUNCRT, '(A, I3, ES17.8)') 'IER, SSQMIN:', IER,
     >            SSQMIN

               IF (IER /= 0) GO TO 860

CC             DO I = 1, NX
CC                YFIT (I) = AE * X (I)**BE + CE
CC             END DO

               CALL LEGENDRE (FALSE, NX, NDEG, X, Y, LEGENDRE_COEFS,
     >                        SSQMIN, TRUE, NX, X, YFIT, IER)

               IF (PLTCRV) THEN
CC                DO I = 1, NPLOT
CC                   YPLOT (I) = AE * XPLOT (I)**BE + CE
CC                END DO
                  CALL LEGENDRE (FALSE, NX, NDEG, X, Y, LEGENDRE_COEFS,
     >                           SSQMIN, TRUE, NPLOT, XPLOT, YPLOT, IER)
               END IF

               CALL SECOND (TIMEB)

               CALL PRCOEF (0, NDEG, LEGENDRE_COEFS, IDATASET, 'COEFS:',
     >                      IOMIT)

            CASE (24)  ! *** Radial basis functions ***

C              Interpolation of smooth data (not noisy data; no smoothing):

               IF (FIRSTSET) THEN
                  SHAPE = ZERO
                  DO WHILE (SHAPE <= ZERO)
                     SHAPE = 3.
                     CALL READR (LUNCRT,
     >                           'Shape parameter > 0.; <CR> = 3: ',
     >                           LUNKBD, SHAPE, DEFAULT, QUIT)
                     IF (QUIT) GO TO 710
                  END DO
                  LEGEND = 'Radial basis functions:  SHAPE = rr.rr'
                  WRITE (LEGEND (34:38), '(F5.2)') SHAPE
                  WRITE (LUNLOG, 1100) TRIM (LEGEND)

                  LEASTSQRS = TRUE
                  CLOSED    = FALSE
               END IF

               CALL SECOND (TIMEA)

               CALL RADIAL_BASIS_1D_WEIGHTS (NX, X, Y, SHAPE, WGHTS,
     >                                       IER)
               IF (IER /= 0) GO TO 860

               IF (PLTCRV) THEN
CC   >            CALL CSEVAL (NDATA, X, Y, NPLOT, XPLOT, COEFS (1,1),
CC   >                         COEFS (1,2), COEFS (1,3), YPLOT)

                  CALL RADIAL_BASIS_1D_EVAL (NX, X, WGHTS, SHAPE, NPLOT,
     >                                       XPLOT, YPLOT)
               END IF

               IF (PLTABS) THEN
CC   >            CALL CSEVAL (NDATA, X, Y, NEVAL, XEVAL, COEFS (1,1),
CC   >                         COEFS (1,2), COEFS (1,3), YEVAL)

                  CALL RADIAL_BASIS_1D_EVAL (NX, X, WGHTS, SHAPE, NEVAL,
     >                                       XEVAL, YEVAL)
               END IF

               CALL SECOND (TIMEB)

            END SELECT

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C           End of case statement for available curve fits/evaluations.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

            DELTAT = TIMEB - TIMEA

            IF (DELTAT < SMALL) THEN
               IF (.NOT. LEASTSQRS) THEN  ! Avoid 2 lines on screen per dataset
                  WRITE (LUNCRT, 1004)
     >            'CPU time negligible for processing dataset', IDATASET
               END IF
            ELSE
               WRITE (LUNCRT, 1006)
     >            ' CPU seconds for processing dataset ', IDATASET,
     >            ': ', DELTAT
            END IF

C ...       Save results.

            IF (FIRSTTYP .AND. YESPLT) THEN  ! Raw data is to be plotted

               SHORTLEG = 'Dataset nn'
               WRITE (SHORTLEG (9:10), '(I2)') IDATASET
               IF (NMLIZ) SHORTLEG (11:) = '  (normalization applied)'

               CALL QPLCRV (LUNPLT, NX, X0, Y0, 'END CURVE', OMIT,
     >                      BLANK, OMIT, OMIT, OMIT, OMIT, OMIT,
     >                      OMIT, SHORTLEG, 'SYMBOLS')
            END IF

C           Evaluations of curve fit at uniform intervals for plotting:

            IF (PLTCRV .OR. NOFIT) THEN  ! NOFIT means plot adjusted data

               IF (NMLIZ) CALL USESCALE (DENRMLIZ, NDIM, NTPLOT, XPLOT,
     >                                  YPLOT, YPLOT, SCALE, SHIFT, IER)

               IF (YESPLT) THEN

                  IF (FIRSTSET) THEN
                     LINE = LTYPES (NUMPLT - MXLTYP*(NUMPLT/(MXLTYP+1)))
                     CALL QPLCRV (LUNPLT, NTPLOT, XPLOT, YPLOT,
     >                            'END CURVE', OMIT, BLANK, OMIT, OMIT,
     >                            OMIT, OMIT, OMIT, OMIT, LEGEND, LINE)

                  ELSE  ! No need for a new legend entry for the same fit type

                     CALL QPLCRV (LUNPLT, NTPLOT, XPLOT, YPLOT,
     >                            'END CURVE', OMIT, BLANK, OMIT, OMIT,
     >                            OMIT, OMIT, OMIT, OMIT, BLANK, LINE)
                  END IF

               END IF

               IF (YESOUT) THEN  ! Save results in SMOOTH format

                  WRITE (LUNOUT, 1013) NTPLOT, TRIM (LEGEND)
                  WRITE (LUNOUT, 1014)
     >               (XPLOT (I), YPLOT (I), I = 1, NTPLOT)
               END IF

            END IF

C ...       Evaluations at specified abscissas:

            IF (PLTABS .AND. .NOT. NOABS) THEN

               IF (NMLIZ) THEN

                  XEVAL (1 : NEVAL) = XEVAL0 (1 : NEVAL)

                  CALL USESCALE (DENRMLIZ, 1, NEVAL, YEVAL, YEVAL,
     >                           YEVAL, SCALE (2), SHIFT (2), IER)
               END IF

               IF (YESPLT)
     >            CALL QPLCRV (LUNPLT, NEVAL, XEVAL, YEVAL, 'END CURVE',
     >                OMIT, BLANK, OMIT, OMIT, OMIT, OMIT, OMIT, OMIT,
     >                'Fit at specified abscissas', 'SYMBOLS')

               IF (YESOUT) THEN
                  WRITE (LUNOUT, 1015) NEVAL, TRIM (LEGEND)
                  WRITE (LUNOUT, 1014)
     >               (XEVAL (I), YEVAL (I), I = 1, NEVAL)
               END IF

            END IF

C ...       Deviations at data abscissas for least-squares methods:

            IF (LEASTSQRS) THEN

               IF (NMLIZ) THEN
                  CALL USESCALE (DENRMLIZ, 1, NX, YFIT, YFIT,
     >                           YFIT, SCALE (2), SHIFT (2), IER)
               END IF

               SSQ    = ZERO
               DEVMAX = ZERO

               DO I = 1, NX
                  YDEV (I) = YFIT (I) - Y0 (I)
                  DEVMAX   = MAX (ABS (YDEV (I)), DEVMAX)
                  SSQ      = YDEV (I) ** 2  +  SSQ
               END DO

               RMSDEV = SQRT (SSQ / NX)

C ...          Log the fit statistics and show them on the screen as well:
    
               WRITE (LUNLOG, 1020) SSQ, RMSDEV, DEVMAX           ! Log file
               IF (IDATASET == 1) WRITE (LUNCRT, 1000)            ! Screen
               WRITE (LUNCRT, 1030) IDATASET, SSQ, RMSDEV, DEVMAX !  "  "

               IF (PLTDEV) THEN

C ...             The following relies on above legend assignments being
C                 within LEGEND (1:25).

                  WRITE (LEGEND (27:91), 1009)
     >               'Sum of Squares:', SSQ, '  RMS Dev:', RMSDEV,
     >               '  Max Dev:', DEVMAX 

                  WRITE (LUNDEV, 1000) ' $OPTIONS'
                  WRITE (LUNDEV, 1040) ' LEGEND=', LEGEND
                  LEGEND (26:91) = BLANK
                  WRITE (LUNDEV, 1000) ' XCOL=1, YCOL=3', ' $END',
     >               '!            X           Fit    Fit - Data'
                  WRITE (LUNDEV, 1050)
     >               (X0 (I), YFIT (I), YDEV (I), I = 1, NX)

               END IF

            END IF

C ...       Look for another dataset (curve) in current file:

            FIRSTSET = FALSE

            IF (.NOT. FINIS)
     >   GO TO 300

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C                End of loop over datasets in input file
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  700    CONTINUE

C ...    Come here upon EOF found reading dataset.  (Will have dropped
C        through for a normal, indefinite NX case, because FINIS is TRUE -
C        hence the test on NX.)

         IF (NX == 0) IDATASET = IDATASET - 1
         IF (IDATASET == 0) GO TO 805

  710    CONTINUE

C ...    Come here after cancelling out of a curve fit option.
C        Go back for another option:

         FIRSTTYP = FALSE

      GO TO 200

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C                End of loop over curve fit selections
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C                           Error handling:
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  800 WRITE (LUNCRT, 1003) ' Error on initial read of file ', DATAFILE
      WRITE (LUNCRT, 1004) 'Status code: ', IOS
      GO TO 999

  805 WRITE (LUNCRT, 1003) ' No datasets found in ', DATAFILE
      GO TO 100

  810 WRITE (LUNCRT, 1003) ' Error reading abscissas from ', XFILE
      WRITE (LUNCRT, 1004) 'Status code: ', IOS
      GO TO 999

  820 WRITE (LUNCRT, 1004) 'Error reading dataset ', IDATASET
      GO TO 999

  830 WRITE (LUNCRT, 1004)
     >   'Abscissas must be monotone increasing. Skipping dataset no.',
     >   IDATASET
      GO TO 300

  840 WRITE (LUNCRT, 1004)
     >   'Abscissas must be monotonic. Skipping dataset no.', IDATASET
      GO TO 300

  850 WRITE (LUNCRT, 1004)
     >   'Successive points must be distinct. Skipping dataset no.',
     >   IDATASET
      GO TO 300

  860 WRITE (LUNCRT, 1003) ' Bad return: ', LEGEND (1:7)
      WRITE (LUNCRT, 1004) 'IER =', IER
      IF (LEGEND (1:1)=='I') WRITE (LUNCRT, 1000) BLANK,
     >   ' IER   IMSL Error Description',
     >   '  33   Extrapolation required at minimum abscissa',
     >   '  34   Extrapolation required at maximum abscissa',
     >   ' 129   Improperly dimensioned workspace array',
     >   ' 130   Too few data points (less than 4)',
     >   ' 131   Abscissas not monotonic increasing'
      GO TO 999

  870 WRITE (LUNCRT, 1000) BLANK,
     >   ' Wagner functions require the data to be on [0,1].',
     >   ' Repeat and invoke the normalization option this time.'
      GO TO 999

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C               Normal termination for current data file
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  900 CONTINUE

C ... Allow for repeating with a new data file (or the same one).

      YES = FALSE
      WRITE (LUNCRT, 1001)
     > ' You may repeat with a new, or the same, file (new plot frame).'
      CALL READY (LUNCRT,
     >   'Do you want to process more data? (Y/N; <CR>=No) ',
     >   LUNKBD, YES, DEFAULT, QUIT)

      CLOSE (LUNDAT)

      IF (YES) THEN
         IF (YESPLT) CALL QPLNXY (LUNPLT, 0, X, Y, 'END FRAME')
         IF (PLTDEV) CALL QPLNXY (LUNDEV, 0, X, Y, 'END FRAME')
         FIRSTFIL = FALSE
         GO TO 100
      END IF

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C                End of loop over multiple data files
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      CALL SECOND (TIME2)
      TOTALT = TIME2 - TIME1

      DO LUN = LUNCRT, LUNLOG, LUNLOG - LUNCRT

         WRITE (LUN, 1005)    ' Total CPU secs:', TOTALT

         IF (YESOUT) THEN
            WRITE (LUN, 1002) ' Reusable output:  ', OUTFILE
            CLOSE (LUNOUT)
         END IF

         IF (YESPLT) THEN
            WRITE (LUN, 1000) ' Plottable data:   smooth.plt'
            CLOSE (LUNPLT)
         END IF

         IF (PLTDEV) THEN
            WRITE (LUN, 1000) ' Deviations plot:  smooth.dev'
            CLOSE (LUNDEV)
         END IF

      END DO

      WRITE (LUNCRT, 1000)    ' Log file:         smooth.log'
      CLOSE (LUNLOG)

  999 WRITE (LUNCRT, 1000)

C*****STOP ! Avoid system dependencies.

C     Formats:
C     --------

 1000 FORMAT (A)
 1001 FORMAT (/, A)
 1002 FORMAT (A, A)
 1003 FORMAT (/, (A, A))
 1004 FORMAT (/, 1X, A, I4)
 1005 FORMAT (/, A, F7.2)
 1006 FORMAT (/, A, I2, A, F6.2)
 1009 FORMAT (3(A, ES10.3))
 1010 FORMAT (A, 2ES24.16)
 1011 FORMAT (A, I1)
 1012 FORMAT (A, I2)
!1013 FORMAT (I3, ' ! Uniform Xs: ', A)
 1013 FORMAT ('# ', I3, ' ! Uniform Xs: ', A)
 1014 FORMAT (2ES24.16)
!1015 FORMAT (I3, ' ! Specified Xs: ', A)
 1015 FORMAT ('# ', I3, ' ! Specified Xs: ', A)
 1020 FORMAT (/, ' Sum of squares =', ES10.3, '  RMS dev =', ES10.3,
     >        '  Max |dev| =', ES10.3)
 1030 FORMAT ('    Curve ', I2, ':   SUMSQRS =', ES10.3,
     >        '  RMS dev =', ES10.3, '  Max |dev| =', ES10.3)
 1040 FORMAT (A, '''', A, ''',')
 1050 FORMAT (3ES14.6)
 1060 FORMAT (3A)
 1100 FORMAT (/, 35('-'), //, 1X, A, /, ' Computed coefficients:')

C     Internal procedure for program SMOOTH:

      CONTAINS

!        ------------------------------------------------------------
         SUBROUTINE PRCOEF (IBGN, IEND, COEFS, IDATASET, TEXT, IOMIT)
!        ------------------------------------------------------------
!
!        Print coefficients from the various smoothing/interpolating methods.
!
!        Arguments:
!        NAME      DIM    TYPE DESCRIPTION
!        IBGN       -      I   Beginning index of COEFS.
!        IEND       -      I   Ending      "   "    "  .
!        COEFS  IBGN:IEND  R   Coefficients.
!        IDATASET    -     I   Number of dataset.
!        TEXT      (*)     C   Descriptive text.
!        IOMIT      -      I   Arguments equal to IOMIT will be suppressed.
!
!        AUTHOR:  Ronald Langhi, Sterling Software, Palo Alto, CA  (Feb. 85)
!
!        -----------------------------------------------------------------------

!        Arguments:

         INTEGER    IDATASET, IBGN, IEND, IOMIT
         REAL       COEFS (IBGN : IEND)
         CHARACTER  TEXT * (*)

!        Local constants:

         INTEGER, PARAMETER :: MXRNGE = 200

!        Local variables:

         INTEGER    I, IRANGE, NEWEND

!        Execution:

         IRANGE = IEND - IBGN + 1

         IF (IRANGE > MXRNGE) THEN
            NEWEND = IBGN + MXRNGE - 1
         ELSE
            NEWEND = IEND
         END IF

         WRITE (LUNLOG, 1040)

         IF (IDATASET /= IOMIT) WRITE (LUNLOG, 1020) IDATASET

         WRITE (LUNLOG, 1040) TEXT, COEFS (IBGN : NEWEND)

         IF (IRANGE > MXRNGE) WRITE (LUNLOG, 1060)
     >      ' Et cetera.  Last coefficient index: ', IRANGE

 1020    FORMAT (' Dataset ', I2, /)
 1040    FORMAT (1X, A, /, (1X, 3ES24.16))
 1060    FORMAT (/, A, I3)

         END SUBROUTINE PRCOEF

      END PROGRAM SMOOTH

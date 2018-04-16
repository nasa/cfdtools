C+---------------------------------------------------------------------
C
      SUBROUTINE PSFIT (NGEOM, XGEOM, YGEOM, METHOD, CLOSED, IER)
C
C  ACRONYM: Parametric Spline FIT (for evaluating Y vs. X or X/Y vs. T)
C           -          -      ---
C  PURPOSE:
C                PSFIT fits an interpolating parametric cubic spline to
C             the given open or closed curve represented  by  abscissas
C             which are not necessarily monotonic, and by corresponding
C             ordinates.    (The parametric spline actually consists of
C             two splines, one for the abscissas and one for the ordin-
C             ates, each fitted against the same values of a parametric
C             variable - cumulative chord length - generated here.)
C
C                Evaluation for given abscissas of a specified subcurve
C             of the fitted parametric spline can be done using PSEVAL,
C             and probably TSUBJ, which returns the value of the  para-
C             metric variable T at a given data point J.  (PSEVAL needs
C             an interval in which to estimate the T which  corresponds
C             to the given X, then Y is computed from this T.)
C
C                Evaluation for a given value or values of T (a simpler
C             problem) is also provided for by subroutine PSTVAL.
C
C                This package was originally written for  working  with
C             airfoil sections with either sharp or rounded leading and
C             trailing edges.
C
C  METHOD:
C                The parametric variable is the usual approximation  to
C             arc length: cumulative chord length from the first point.
C
C                The splines for X vs. T and Y vs. T are interpolatory,
C             with a choice of conventional or monotonic.  (A monotonic
C             spline follows the data better in difficult cases such as
C             adjacent points with equal Y where flatness is  required.
C             However, monotonic splines are not appropriate for curva-
C             ture calculations, because the second derivatives are not
C             necessarily continuous.)
C
C                In either case,  the  curve  may be closed smoothly by
C             invoking periodic boundary conditions in the spline fits.
C
C                This implementation is intended to permit  evaluations
C             of Y for given values of either X or T.   It  necessarily
C             uses an internal COMMON area for the spline coefficients,
C             since this is the only way that the zero-finder  employed
C             by PSEVAL can get at the information it needs.  
C
C                One consequence is that PSFIT and PSEVAL or PSTVAL can
C             be used on only one curve at a time (unless copies of the
C             internal COMMON data  are kept - not recommended - appli-
C             cation programs should normally NOT reference /PSWORK/).
C
C  ARGUMENTS:
C   ARG       TYPE I/O/S  DIM    DESCRIPTION
C  NGEOM        I    I     -     Number of data points.
C  XGEOM,YGEOM  R    I   NGEOM   Coordinates of the data points.
C  METHOD     C*(*)  I     -     'C' or 'CC': conventional spline fits
C                                             using IENDL=IENDR=0;
C                                'M' or 'MM': monotonic spline fits.
C                                'MC': monotonic for X vs. T, but
C                                      conventional for Y vs. T;
C                                'CM': converse of 'MC'.
C  CLOSED       L    I     -     .TRUE. means the data points form a
C                                closed curve. In this case, the calling
C                                program must duplicate the first data
C                                point (XGEOM (1), YGEOM (1)) as
C                                (XGEOM (NGEOM), YGEOM (NGEOM)), and
C                                smoothness in the spline will be sought
C                                here as at all the other data points;
C                                .FALSE. means the end-points are treated
C                                as distinct, and are not joined (though
C                                in fact they MAY be the same point, as
C                                in the case of an airfoil with a sharp
C                                trailing edge).
C  IER          I    O     -     Error return code:
C                                IER = 0: No errors were detected;
C                                IER = 1: Too few data points:
C                                         (NGEOM < 2 for open curve or
C                                          NGEOM < 3 for closed curve);
C                                IER = 4: Closed curve case only:
C                                         First and last data points do
C                                         not match;
C                                IER = 5: Too many data points: NGEOM
C                                         exceeds MAXPS used in the
C                                         internal COMMON block.
C
C  INTERNAL COMMON BLOCK:
C   /PSWORK/ Work space required for parametric spline interpolation.
C            Row dimensions are set to MAXPS by PARAMETER statement.
C            (/PSWORK/ is normally NOT referenced by the calling program.)
C
C   VAR    TYPE I/O/S DIM  DESCRIPTION
C   PSCURV  R     O MAXPS  Copy of XGEOM for FUNCTION GETT, which is
C                          an argument of a zero-finder and is itself
C                          limited to a single argument.
C   TGEOM   R     O MAXPS  Values of parametric variable for the data
C                          points, generated here; available to the
C                          user via FUNCTION TSUBJ (to avoid COMMON).
C   PSCOFS  R     O MAXPS  Cubic spline coefficients from CSFIT, one
C                   *3*2   set for XGEOM vs. TGEOM, one for YGEOM vs.
C                          TGEOM.
C   NCURV   I     O   -    Copy of NGEOM.
C
C  PROCEDURES:
C   CHORD   Safeguarded calculation of arc lengths
C   CSFIT   Fits conventional interpolating spline
C   MSFIT   Fits monotonic interpolating spline
C
C  RELATED MODULES:
C   PSEVAL  Evaluates parametric spline set up by PSFIT for 
C           ordinates and derivatives at arbitrary abscissas (Xs)
C           on a specified monotonic sub-curve.  Uses FUNCTION
C           GETT internally, which should not concern the user.
C   TSUBJ   FUNCTION routine which may be used to determine intervals
C           defining subcurves for evaluation purposes.  Both GETT
C           and TSUBJ are in the same source module as PSEVAL.
C   PSTVAL  Evaluates parametric spline set up by PSFIT for given
C           value(s) of T (as opposed to X).  Returns X, Y, Y', Y''.
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77
C
C  HISTORY:
C     Jan. 83  R.Lefkowitz  Original implementation using conventional
C                           splines from IMSL, which did not permit
C                           proper handling of closed curves.
C     Sep. 84    "    "     Handled closed curves properly via CSFIT.
C     Sep. 85  D.Saunders   Changed from index number to cumulative chord
C                           length for the parametric variable.
C     Jan. 87    "    "     Introduced monotonic spline option (prompted
C                           by a mesh generation requirement).
C     Sep. 87    "    "     Made it clear that /PSWORK/ is internal.
C     Oct. 91    "    "     Introduced CHORD for more careful arc-length
C                           calculations; added METHOD argument in place
C                           of NGEOM < 0 for indicating monotonic option;
C                           introduced IMPLICIT NONE; other cosmetics.
C     Aug. 93    "    "     Provided METHOD='MC' (and 'CM') option to
C                           accommodate redistribution of blunt-nosed
C                           airfoils where 'C' can allow the curve fit
C                           to exceed the data range at the leading edge.
C
C  AUTHOR: Rosalie C. Lefkowitz, Sterling Software/NASA Ames, Moffett Field, CA
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   NGEOM, IER
      REAL
     >   XGEOM (NGEOM), YGEOM (NGEOM)
      LOGICAL
     >   CLOSED
      CHARACTER
     >   METHOD * (*)

C     Internal COMMON.

      INTEGER
     >   MAXPS
      PARAMETER
     >  (MAXPS = 1001)
      INTEGER
     >   NCURV
      REAL
     >   PSCURV, TGEOM, PSCOFS, XLOC
      COMMON /PSWORK/
     >   PSCURV (MAXPS), TGEOM (MAXPS), PSCOFS (MAXPS,3,2), XLOC, NCURV

C     Local constants.

      REAL
     >   ZERO
      PARAMETER
     >  (ZERO = 0.E+0)

C     Local variables.

      INTEGER
     >   IENDL, J
      CHARACTER
     >   FIT * 2

C     Procedures.

      REAL
     >   CHORD
      EXTERNAL
     >   CHORD, CSFIT, MSFIT


C     Execution.

C     Too few data points?

      IF (NGEOM .LT. 2 .OR. (CLOSED .AND. NGEOM .LT. 3)) THEN
         IER = 1
         GO TO 999
      END IF

C     Too many data points?

      IF (NGEOM .GT. MAXPS) THEN
         IER = 5
         GO TO 999
      END IF


C     Before METHOD = 'MC' or 'CM' were allowed, applications are assumed
C     to have passed single-character strings:

      FIT = METHOD
      IF (LEN (METHOD) .EQ. 1) FIT (2 : 2) = FIT (1 : 1)


C     Set up parametric variable values, and transfer original abscissas
C     to internal COMMON for later use by zero-finder called by PSEVAL.

      NCURV = NGEOM
      TGEOM (1) = ZERO
      PSCURV (1) = XGEOM (1)

      DO 20, J = 2, NGEOM
         TGEOM (J) = CHORD (XGEOM, YGEOM, J-1, J) + TGEOM (J-1) 
         PSCURV (J) = XGEOM (J)
   20 CONTINUE


C     For conventional splines, smooth closure means periodic end conditions.

      IF (CLOSED) THEN
         IENDL = 4
      ELSE
         IENDL = 0
      END IF


C     Fit spline for XGEOM vs. TGEOM.

      IF (FIT (1 : 1) .EQ. 'C') THEN

C        Conventional spline:

         CALL CSFIT (NGEOM, TGEOM, XGEOM, IENDL, ZERO, 0, ZERO,
     >      PSCOFS (1,1,1), PSCOFS (1,2,1), PSCOFS (1,3,1), IER)

      ELSE

C        Monotonic spline: CLOSED used for "cyclic" is misleading, but valid.

         CALL MSFIT (NGEOM, TGEOM, XGEOM, CLOSED,
     >      PSCOFS (1,1,1), PSCOFS (1,2,1), PSCOFS (1,3,1), IER)
      END IF

      IF (IER .NE. 0) GO TO 999


C     Likewise for YGEOM vs. TGEOM.

      IF (FIT (2 : 2) .EQ. 'C') THEN

         CALL CSFIT (NGEOM, TGEOM, YGEOM, IENDL, ZERO, 0, ZERO,
     >      PSCOFS (1,1,2), PSCOFS (1,2,2), PSCOFS (1,3,2), IER)

      ELSE

         CALL MSFIT (NGEOM, TGEOM, YGEOM, CLOSED,
     >      PSCOFS (1,1,2), PSCOFS (1,2,2), PSCOFS (1,3,2), IER)

      END IF


  999 RETURN
      END

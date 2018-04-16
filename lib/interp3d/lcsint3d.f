C+----------------------------------------------------------------------
C
      FUNCTION LCSINT3D (NX, NY, NZ, IDIMX, IDIMY, XTBL, YTBL, ZTBL,
     >                   FTBL, X, Y, Z, IER)
C
C     ONE-LINER: Local Cubic Spline INTerpolation in 3 Dimensions
C                -     -     -      ---              - -
C
C     DESCRIPTION:
C
C        LCSINT3D applies 1-dimensional local cubic spline techniques
C     to a 3-dimensional rectangular table in place of the linear
C     interpolations of module TABLE3, from which it is adapted.
C     
C        The hope is to provide a utility which may suffice for applications
C     requiring reasonably smooth functions, such as those involving gradient
C     methods of optimization.  In the absence of reliable 3-D methods,
C     the present approach should be preferable to TABLE3.  It is understood
C     that a given 3-D interpolation using nonlinear 1-D methods typically
C     depends on the order in which the 1-D interpolations are applied (in
C     this case the order X, Y, Z).  But the assumption here is that smooth
C     variation of the result with respect to small changes in the target
C     point is more important than maximum precision, which is probably
C     lacking in the data table anyway.
C
C        In view of these assumptions, only the "monotonic" choice of local
C     method is offered here, to be sure of avoiding wild excursions in
C     the 1-D fits.  Thus the calling sequence matches that of TABLE3.
C
C        Note the rectangularity requirement:  generalizing LCSINT2D (which
C     handles a form of pseudorectangularity) was considered intractable.
C
C        "Local cubic spline" methods are essentially 4-point methods
C     which avoid calculating and storing all spline coefficients when
C     any are needed by generating them on the fly from local data.
C     They provide first derivative continuity across data points but
C     cannot guarantee second derivative continuity.
C
C        As for TABLE3, LCSINT3D returns a single result per call.
C     The rectangular table coordinates need not be uniform, but they
C     must be monotonic increasing or monotonic decreasing in each
C     dimension.
C
C     ARGUMENTS:
C
C        ARG    DIM     TYPE I/O/S DESCRIPTION
C        NX      -        I    I   Number of "columns" in the table, >= 2
C        NY      -        I    I   Number of "rows", >= 2
C        NZ      -        I    I   Number of "layers", >= 2
C        IDIMX,  -        I    I   Declared first two dimensions of FTBL (*,*,*)
C        IDIMY
C        XTBL   NX        R    I   Rectangular table coordinates in the
C        YTBL   NY                 1st, 2nd, and 3rd dimensions
C        ZTBL   NZ
C        FTBL  IDIMX,     R   I    The table function values should be in
C             IDIMY,NZ             FTBL (1:NX, 1:NY, 1:NZ)
C        X,Y,Z   -        R    I   Point at which interpolation is required
C        IER     -        I    O   Error return code - see below
C        LCSINT3D         R    O   The interpolated value at (X, Y, Z) is
C                                  returned as the FUNCTION value.
C                                  LCSINT3D = 0. if an error is detected
C                                  (rather than being left undefined).
C     PROCEDURES:
C
C        INTERVAL   1-D search utility
C        LCSFIT     1-D local cubic spline utility
C
C     ERROR HANDLING:
C
C        IER = 0 if no error was detected;
C              = 1 if NX < 2;
C              = 2 if NY < 2;
C              = 3 if NZ < 2.
C        There is no check for monotonicity in X(*), Y(*), and Z(*).
C
C     ENVIRONMENT:  FORTRAN 77 + IMPLICIT NONE + names <= 8 characters + END DO
C
C     HISTORY:
C
C        04/17/93   DAS   Adaptation of TABLE3.
C
C     AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   IDIMX, IDIMY, NX, NY, NZ, IER
      REAL
     >   LCSINT3D, X, Y, Z, XTBL (NX), YTBL (NY), ZTBL (NZ),
     >   FTBL (IDIMX, IDIMY, NZ)

C     Local constants:

      REAL
     >   ONE, ZERO
      LOGICAL
     >   NEW
      CHARACTER
     >   METHOD * 1
      PARAMETER
     >  (ONE    = 1.E+0,
     >   ZERO   = 0.E+0,
     >   NEW    = .TRUE.,  ! New data for each call LCSFIT
     >   METHOD = 'M')     ! Monotonic option

C     Local variables:

      INTEGER
     >   J, K, I1, J1, K1, I2, J2, K2, NI, NJ, NK
      REAL
     >   FYATX (4), FZATXY (4), FXYZ, SLOPE

C     Procedures:

      EXTERNAL
     >   INTERVAL, LCSFIT

C     Execution:

      LCSINT3D = ZERO
      IER = 0

      IF (NX .LT. 2) THEN
         IER = 1
         GO TO 99
      END IF

      IF (NY .LT. 2) THEN
         IER = 2
         GO TO 99
      END IF

      IF (NZ .LT. 2) THEN
         IER = 3
         GO TO 99
      END IF

C     Locate the cell containing (X, Y, Z), and isolate the relevant
C     table entries.  LCSFIT uses at most 4 points per interpolation
C     (3 at a boundary; 2 are also possible; 1 is not allowed).

      I1 = 1

      CALL INTERVAL (NX, XTBL, X, SIGN (ONE, XTBL (2) - XTBL (1)), I1)

      I2 = MIN (I1 + 2, NX)
      I1 = MAX (I1 - 1, 1)
      NI = I2 - I1 + 1

      J1 = 1

      CALL INTERVAL (NY, YTBL, Y, SIGN (ONE, YTBL (2) - YTBL (1)), J1)

      J2 = MIN (J1 + 2, NY)
      J1 = MAX (J1 - 1, 1)
      NJ = J2 - J1 + 1

      K1 = 1
      CALL INTERVAL (NZ, ZTBL, Z, SIGN (ONE, ZTBL (2) - ZTBL (1)), K1)

      K2 = MIN (K1 + 2, NZ)
      K1 = MAX (K1 - 1, 1)
      NK = K2 - K1 + 1

C     Perform 1-D interpolations along the X direction, then Y, then Z:

      DO K = K1, K2

         DO J = J1, J2

            CALL LCSFIT (NI, XTBL (I1), FTBL (I1, J, K), NEW, METHOD,
     >                   1, X, FYATX (J - J1 + 1), SLOPE)
         END DO

         CALL LCSFIT (NJ, YTBL (J1), FYATX, NEW, METHOD,
     >                1, Y, FZATXY (K - K1 + 1), SLOPE)
      END DO

      CALL LCSFIT (NK, ZTBL (K1), FZATXY, NEW, METHOD,
     >             1, Z, FXYZ, SLOPE)

      LCSINT3D = FXYZ

   99 RETURN
      END

      SUBROUTINE D1PJIS (CMEM, IMEM, DMEM, ITS, U, V, IPJ, NWORK, IER)

C     PURPOSE:
C        Point on Joined Initialize to Surface.  Initialize the Point
C        on Joined IPJ to be the point determined by parameter values
C        U and V on Trimmed Surface ITS.  This Trimmed Surface
C        pointer is copied into the IPJ point entity and the U,V
C        and X,Y,Z values are calculated and stored into the IPJ
C        point entity.
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        ITS      Pointer to Trimmed Surface entity.  If ITS = 0, then
C                 all values in the PJ entity are reset to zero.
C        U, V     Parameter values on Trimmed Surface
C        IPJ      Pointer to Point on Joined entity
C
C     OUTPUT:
C        NWORK    Amount of working storage used/needed.
C        IER      Error flag.  If negative, the entity IPJ is not
C                 updated.
C                 0     No errors detected.
C
C                 -2    IPJ does not point to a Point on Joined
C                       entity.
C
C                 -6  = ITS does not point to a Trimmed Surface entity
C
C                 -7  = ITS->IBFS does not point to a B-Spline
C                       function entity
C
C                 -9  = Insufficient DMEM for working storage or
C                       insufficient IMEM to allocate working storage
C                       entity.
C
C                 -19 = The Surface B-spline function entity
C                       has K(i) .LE. 0 for some i
C                 -20 = The Surface B-spline function entity
C                       has a number of B-spline coefficients
C                       with respect to the ith independent
C                       variable is less than Ki for some i
C                 -21 = The Surface B-spline function entity
C                       has invalid knot set
C                 -22 = The Surface B-spline function entity
C                       has denominator = 0
C                 -23 = U,V is inside too small an interval
C                 -24 = U or V is out of range.
C                 -25 = The Surface B-spline function entity
C                       has NDOM .NE. 2
C                 -26 = The Surface B-spline function entity
C                       has NRNG .LT. 1
C                 -27 = The Surface B-spline function entity
C                       does not have 3 dependent variables
C
C                 -999= Unexpected error occurred in a called routine
C
C     CALLS:
C        D0TSFP
C        D0PJI2
C        D2STAI
C        D2STAD
C
C
C     HISTORY:
C        16Jul92  D. Parsons   Created.
C
C     ------

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER           ITS, IPJ, NWORK, IER
      DOUBLE PRECISION  U, V

      INTEGER           IBFS, IEG, IDMY, IERX, IA(2)
      DOUBLE PRECISION  T, X, Y, Z, DA(6)
      CHARACTER         LABELX

      IER = 0

      IF (ITS .NE. 0) THEN

C        Get the pointer to the Surface B-Spline from the ITS entity

         CALL D0TSFP (CMEM, IMEM, DMEM, ITS, LABELX, IBFS, IDMY, IDMY,
     +         IDMY, IDMY, IDMY, IERX)
         IF (IERX .NE. 0) THEN
            IER = -6
            GOTO 9000
         ENDIF

C        Calculate the X, Y, and Z parameter values from the Surface
C        Spline

         CALL D0PJI2 (CMEM, IMEM, DMEM, IBFS, U, V,
     +         NWORK, X, Y, Z, IER)
         IF (IER .NE. 0) THEN
            GOTO 9000
         ENDIF

      ELSE
         U = 0.0D0
         V = 0.0D0
         X = 0.0D0
         Y = 0.0D0
         Z = 0.0D0
      ENDIF

C     Write the Trimmed Surface (ITS), U, V, X, Y, and Z values into
C     the Point entity and zero the Edge and T values.

      IEG = 0
      T   = 0.0D0

      IA(1) = IEG
      IA(2) = ITS

      CALL D2STAI (CMEM, IMEM, DMEM, IA, 1, 2, IPJ, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9000
      ENDIF

      DA(1) = T
      DA(2) = U
      DA(3) = V
      DA(4) = X
      DA(5) = Y
      DA(6) = Z

      CALL D2STAD (CMEM, IMEM, DMEM, DA, 1, 6, IPJ, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9000
      ENDIF

 9000 CONTINUE

      RETURN
      END

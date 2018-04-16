      SUBROUTINE D1PJIC (CMEM, IMEM, DMEM, IEG, T, IPJ, NWORK, IER)

C     PURPOSE:
C        Point on Joined Initialize to Curve.  Initialize the Point
C        on Joined IPJ to be the point determined by parameter value
C        T on Edge IEG.  Edge IEG points to a Loop entity, which
C        points to a Trimmed Surface entity.  This Trimmed Surface
C        pointer is copied into the IPJ point entity and the U,V
C        and X,Y,Z values are calculated.
C
C        If ITS is already defined, but does not match the
C        Edge->Loop->ITS, an error will be returned.
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IEG      Pointer to Edge entity.  If IEG = 0, then
C                 the Edge and T values in the PJ entity are
C                 reset to zero.
C        T        Parameter value on Edge curve
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
C                 -3  = IEG does not point to an Edge entity.
C                 -4  = IEG->IBFC does not point to a B-spline
C                       Function entity.
C
C                 -5  = IEG->ILP does not point to a Loop entity
C
C                 -6  = IEG->ILP->ITS does not point to a Trimmed
C                       Surface entity
C
C                 -7  = IEG->ILP->ITS->IBFS does not point to a
C                       B-Spline function entity
C
C                 -8  = ITS is not zero and IEG->LOOP->ITS are
C                       different.
C
C                 -9  = Insufficient DMEM for working storage or
C                       insufficient IMEM to allocate working storage
C                       entity.
C
C                 -10 = The Curve B-spline function entity
C                       has C(3) .LE. 0
C                 -11 = The Curve B-spline function entity
C                       has C(4) < C(3)
C                 -12 = The Curve B-spline function entity
C                       has invalid knot set
C                 -13 = The Curve B-spline function entity
C                       has denominator = 0
C                 -14 = T is inside too small an interval
C                 -15 = T out of range.
C                 -16 = The Curve B-spline function entity
C                       has NDOM .NE. 1
C                 -17 = The Curve B-spline function entity
C                       has NRNG .LT. 1
C                 -18 = The Curve B-spline function entity
C                       does not have 2 dependent variables
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
C        D0EGFP
C        D0LPFP
C        D0PJFP
C        D0PJI1
C        D0PJIS
C        D2STEI
C        D2STED
C
C     HISTORY:
C        16Jul92  D. Parsons   Created.
C
C     ------

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IEG, IPJ, NWORK, IER
      DOUBLE PRECISION  T

      CHARACTER         LABELX
      INTEGER           IDMY, ILP, ITSE, ITSP, IERX
      INTEGER           ITS, IBFC
      DOUBLE PRECISION  DDMY, U, V

C     ------

      IER = 0

      IF (IEG .NE. 0) THEN

C        Find the IEG->LOOP->ITS Trimmed Surface

         CALL D0EGFP (CMEM, IMEM, DMEM, IEG, LABELX, IBFC, ILP, IDMY,
     +         IDMY, IDMY, IERX)
         IF (IERX .NE. 0) THEN
            IER = -3
            GOTO 9000
         ENDIF

         CALL D0LPFP (CMEM, IMEM, DMEM, ILP, LABELX, ITSE, IDMY,
     +      IDMY, IDMY, IERX)
         IF (IERX .NE. 0) THEN
            IER = -5
            GOTO 9000
         ENDIF

         IF (ITSE .EQ. 0) THEN
            IER = -6
            GOTO 9000
         ENDIF

C        Check to see if ITS is already defined for this Point

         CALL D0PJFP (CMEM, IMEM, DMEM, IPJ, LABELX, IDMY, ITSP,
     +         DDMY, DDMY, DDMY, DDMY, DDMY, DDMY, IDMY, IERX)
         IF (IERX .NE. 0) THEN
C           (This shouldn't happen!)
            IER = -2
            GOTO 9000
         ENDIF

         IF (ITSP .NE. 0) THEN
            IF (ITSP .NE. ITSE) THEN
               IER = -8
               GOTO 9000
            ENDIF
         ENDIF

         ITS = ITSE

C        Calculate the U and V parameter values from the Edge Curve
C        Spline

         CALL D0PJI1 (CMEM, IMEM, DMEM, IBFC, T, NWORK, U, V, IER)
         IF (IER .NE. 0) THEN
            GOTO 9000
         ENDIF

C        Calculate X,Y,Z values for the point and store
C        ITS, U, V, X, Y, and Z values to the Point entity

         CALL D1PJIS (CMEM, IMEM, DMEM, ITS, U, V, IPJ,
     +         NWORK, IER)
         IF (IER .NE. 0) THEN
            GOTO 9000
         ENDIF

      ENDIF

C     Write the Edge and T values into the Point entity.

      CALL D2STEI (CMEM, IMEM, DMEM, IEG, 1, IPJ, IER)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9000
      ENDIF

      CALL D2STED (CMEM, IMEM, DMEM, T, 1, IPJ, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9000
      ENDIF

 9000 CONTINUE

      RETURN
      END

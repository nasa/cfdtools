      SUBROUTINE D2PJDC (CMEM, IMEM, DMEM, DT, IPJ, IER)

C     PURPOSE:
C        Point on Joined Delta Move on Curve.  Change the Point
C        on Joined IPJ by adding DT to the current parameter value
C        of the same Edge on which it is currently located.  Error
C        if not on an Edge.
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        DT       Change to be made to the parameter value on
C                 Edge curve
C        IPJ      Pointer to Point on Joined entity
C
C     OUTPUT:
C        IER      Error flag.  If negative, the entity IPJ is not
C                 updated.
C                 0     No errors detected.
C                 -1    Dynamic Memory is corrupt or uninitialized.
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
C                 -9  = Insufficient DMEM for working storage.
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
C        D1MBAD
C        D0PJFP
C        D1PJIC
C        DTERR
C
C     HISTORY:
C        16Jul92  D. Parsons   Created.
C
C     ------

C     Long name alias:
C        ENTRY D2_DELTA_MOVE_POINT_CURVE (CMEM, IMEM, DMEM,
C    +         DT, IPJ, IER)

      EXTERNAL          D1MBAD
      LOGICAL           D1MBAD

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IPJ, IER
      DOUBLE PRECISION  DT

      INTEGER           IEG, IDMY, IERX, NWORK
      DOUBLE PRECISION  T, DDMY
      CHARACTER         LABELX

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2PJDC'/

C     ------

      IER = 0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

C     Check the validity of the IPJ entity and get the old values

      CALL D0PJFP (CMEM, IMEM, DMEM, IPJ, LABELX, IEG, IDMY, T,
     +      DDMY, DDMY, DDMY, DDMY, DDMY, IDMY, IERX)
      IF (IERX .NE. 0) THEN
         IER = -2
         GOTO 9000
      ENDIF

      IF (IEG .EQ. 0) THEN
         IER = -3
         GOTO 9000
      ENDIF

C     Update the Point on Joined entity

      T = T+DT

      CALL D1PJIC (CMEM, IMEM, DMEM, IEG, T, IPJ, NWORK, IER)

      IF (IER .NE. 0) GOTO 9000

      GOTO 9999

C     Error Handling

 9000 CONTINUE

      IF (IER .EQ. -999) THEN
         CALL DTERR (5, SUBNAM, IER, 0)
      ELSE IF (IER .EQ. -9) THEN
         CALL DTERR (2, SUBNAM, IER, NWORK)
      ELSE
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

 9999 CONTINUE

      RETURN
      END

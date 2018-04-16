      SUBROUTINE D2PJIS (CMEM, IMEM, DMEM, ITS, U, V, IPJ, IER)

C     PURPOSE:
C        Point on Joined Initialize on Surface.  Initialize the Point
C        on Joined IPJ to be the point determined by parameter values
C        (U,V) on Trimmed Surface ITS.
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        ITS      Pointer to Trimmed Surface entity.  If ITS = 0, then
C                 the Trimmed Surface and (U,V) values in the PJ entity
C                 are reset to zero.
C        (U,V)    Parameter values on Trimmed Surface
C        IPJ      Pointer to Point on Joined entity
C
C     OUTPUT:
C        IER      Error flag.  If negative, the entity IPJ is not
C                 updated.
C                 0     No errors detected.
C                 -1  = Dynamic memory is corrupt or uninitialized.
C                 -2  = IPJ does not point to a Point on Joined
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
C        D1MBAD
C
C     HISTORY:
C        09Jul92  D. Parsons   Created.
C
C     ------

C     Long name alias:
C        ENTRY D2_INITIALIZE_POINT_SURFACE (CMEM, IMEM, DMEM,
C    +         ITS, U, V, IPJ, IER)
      EXTERNAL          D1MBAD
      LOGICAL           D1MBAD

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER           ITS, IPJ, IER
      DOUBLE PRECISION  U, V

      INTEGER           IDMY, IERX, NWORK
      DOUBLE PRECISION  DDMY
      CHARACTER         LABELX

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2PJIS'/

C     ------

      IER = 0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

C     Check the validity of the IPJ entity

      CALL D0PJFP (CMEM, IMEM, DMEM, IPJ, LABELX, IDMY, IDMY, DDMY,
     +      DDMY, DDMY, DDMY, DDMY, DDMY, IDMY, IERX)
      IF (IERX .NE. 0) THEN
         IER = -2
         GOTO 9000
      ENDIF

C     Update the Point on Joined entity

      CALL D1PJIS (CMEM, IMEM, DMEM, ITS, U, V, IPJ, NWORK, IER)

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

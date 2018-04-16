      SUBROUTINE D2PJDS (CMEM, IMEM, DMEM, DU, DV, IPJ, IER)

C     PURPOSE:
C        Point on Joined Delta Move on Surface.  Change the Point
C        on Joined IPJ by adding Du and DV to the current surface
C        parameter values on the same Surface on which it is currently
C        located.  If moving away from an Edge, the Edge pointer is
C        nulled.  Error if not on a Surface.
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        (DU,DV)  Change to be made to the parameter values on
C                 Trimmed Surface
C        IPJ      Pointer to Point on Joined entity
C
C     OUTPUT:
C        IER      Error flag.  If negative, the entity IPJ is not
C                 updated.
C                 0     No errors detected.
C
C                 -1    Dynamic memory corrupt or uninitialized.
C
C                 -2    IPJ does not point to a Point on Joined
C                       entity.
C
C                 -6  = ITS does not point to a Trimmed Surface entity
C
C                 -7  = ITS->IBFS does not point to a B-Spline
C                       function entity
C
C                 -9  = Insufficient DMEM for working storage.
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
C        D1PJIS
C        DTERR
C
C     HISTORY:
C        16Jul92  D. Parsons   Created.
C
C     ------

C     Long name alias:
C        ENTRY D2_DELTA_MOVE_POINT_SURFACE (CMEM, IMEM, DMEM,
C    +      DU, DV, IPJ, IER)

      EXTERNAL          D1MBAD
      LOGICAL           D1MBAD

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IPJ, IER
      DOUBLE PRECISION  DU, DV

      INTEGER           ITS, IERX, IDMY, NWORK
      DOUBLE PRECISION  U, V, DDMY
      CHARACTER         LABELX

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2PJDS'/

C     ------

      IER = 0

C     Set the tolerance for determining whether (xc,yc,zc) = (xs,ys,zs)

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

C     Check the validity of the IPJ entity and get ITS, U, and V

      CALL D0PJFP (CMEM, IMEM, DMEM, IPJ, LABELX, IDMY, ITS, DDMY, 
     +      U, V, DDMY, DDMY, DDMY, IDMY, IERX)
      IF (IERX .NE. 0) THEN
         IER = -2
         GOTO 9000
      ENDIF

      IF (ITS .EQ. 0) THEN
         IER = -6
         GOTO 9000
      ENDIF

      U = U+DU
      V = V+DV

C     Update the Point on Joined entity

      CALL D1PJIS (CMEM, IMEM, DMEM, ITS, U, V, IPJ, NWORK, IER)

      IF (IER .NE. 0) GOTO 9000

      GOTO 9999

 9000 CONTINUE

C     Error Handling

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


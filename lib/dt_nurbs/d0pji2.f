      SUBROUTINE D0PJI2 (CMEM, IMEM, DMEM, IBFS, U, V,
     +      NWORK, X, Y, Z, IER)

C     PURPOSE:
C        Evaluate a B-spline Surface at (U,V).
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IBFS     Pointer to B-Spline Surface entity.
C        U,V      Parameter values on Trimmed Surface
C
C     OUTPUT:
C        NWORK    Amount of working storage used/needed
C        X,Y,Z    3-Space coordinates of the point
C        IER      Error flag.
C
C                 0     No errors detected.
C
C                 -7  = IBFS does not point to a B-Spline
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
C
C     CALLS:
C        D0BFFP
C        D2FEBD
C        DTERPT
C        D1DEFE
C        D2ERR
C        DTNPVL
C        D2UNLD
C        D2ERAS
C        D2TSLD
C
C     HISTORY:
C        16Jul92  D. Parsons   Created.
C
C     ------

      EXTERNAL          D2TSLD
      LOGICAL           D2TSLD

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

C =====================================================================
C
C     INCLUDE FILE DITYID.INC
C
C     Reserved Entity Type ID Numbers
C
C     HISTORY
C        12May92     D. Parsons     created
C
C                          Double Precision Array
      INTEGER     ENTDA
      PARAMETER  (ENTDA = 253)

C =====================================================================

      INTEGER           IBFS, NWORK, IER
      DOUBLE PRECISION  U, V, X, Y, Z

      INTEGER           IDMY, IERX
      INTEGER           CSTART, CEND, IWKST, IWKEND, IDWORK
      INTEGER           NRNG, NDOM, KMAX, NZERO, I, K
      DOUBLE PRECISION  EV(4), EX(2)
      CHARACTER         LABELX
      LOGICAL           INITIZ

C     ------

      IER = 0

      IDWORK = 0
      X = 0.0D0
      Y = 0.0D0
      Z = 0.0D0

C     Fetch bounds and lock the B-spline function

      CALL D0BFFP (CMEM, IMEM, DMEM, IBFS, LABELX, IDMY, IDMY, IERX)
      IF (IERX .NE. 0) THEN
         IER = -7
         GOTO 9000
      ENDIF

      CALL D2FEBD (CMEM, IMEM, DMEM, IBFS, CSTART, CEND, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9000
      ENDIF

C     Set NDOM = C(1)

      NDOM = DMEM(CSTART)

      IF (NDOM .NE. 2) THEN
         IER = -25
         GOTO 9000
      ENDIF

C     Set NRNG = C(2)

      NRNG = DMEM(CSTART+1)

      IF (NRNG .EQ. 0) THEN
         IER = -26
         GOTO 9000
      ELSEIF (NRNG .LT. 0) THEN
C        Compensate for rational spline
         NRNG = -1-NRNG
      ENDIF

      IF (NRNG .NE. 3) THEN
         IER = -27
         GOTO 9000
      ENDIF

C     Define and lock some working storage

      KMAX = 0
      NZERO = 1

      DO 10 I = 1, NDOM
         K = DMEM(CSTART+I+1)
         IF (K .LE. 0) THEN
            IER = -19
            GOTO 9000
         ENDIF
         KMAX = MAX(KMAX,K)
         IF (NDOM .GE. 2) NZERO = NZERO * K
   10 CONTINUE

      NRNG =  ABS(DMEM(CSTART+1))
      NWORK = NZERO * (NRNG+1) + 3 * KMAX + NDOM

      INITIZ = .FALSE.

      CALL D1DEFE (CMEM, IMEM, DMEM, ENTDA, 0, 0, NWORK, INITIZ,
     +      IDWORK, IERX)
      IF (IERX .NE. 0) THEN
         IER = -9
         GOTO 9000
      ENDIF

      CALL D2FEBD (CMEM, IMEM, DMEM, IDWORK, IWKST, IWKEND, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9000
      ENDIF

C     Evaluate C(U,V)

      EV(1) = 0.0D0
      EV(2) = 0.0D0
      EV(3) = 0.0D0
      EV(4) = 0.0D0

      EX(1) = U
      EX(2) = V

      CALL DTERPT(0)
      CALL DTNPVL (EX, 1, DMEM(CSTART), DMEM(IWKST), NWORK, EV, IERX)
      CALL DTERPT(1)

      IF (IERX .NE. 0) THEN
         IF       (IERX .EQ. -1) THEN
C           C(3) .LE. 0)
            IER = -19
         ELSE IF  (IERX .EQ. -6) THEN
C           C(4) < C(3)
            IER = -20
         ELSE IF  (IERX .EQ. -8) THEN
C           Invalid knot set
            IER = -21
         ELSE IF  (IERX .EQ. -10) THEN
C           Denominator = 0
            IER = -22
         ELSE IF  (IERX .EQ. -38) THEN
C           Interval too small
            IER = -23
         ELSE IF  (IERX .EQ. -50) THEN
C           U,V out of range
            IER = -24
         ELSE IF  (IERX .EQ. -51) THEN
C           C(1) .NE. 2
            IER = -25
         ELSE IF  (IERX .EQ. -52) THEN
C           C(2) = 0
            IER = -26
         ELSE
            CALL DTERR (1, 'DTNPVL', IERX, 0)
            IER = -999
         ENDIF

         GOTO 9000
      ENDIF

      X = EV(1)
      Y = EV(2)
      Z = EV(3)

C     Unlock the B-spline Function entity

      CALL D2UNLD (CMEM, IMEM, DMEM, IBFS, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9000
      ENDIF

C     Erase the Working storage entity

      CALL D2ERAS (CMEM, IMEM, DMEM, IDWORK, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9000
      ENDIF

      GOTO 9999

 9000 CONTINUE

C     Error recovery:

C     Try to unlock the B-spline Surface Function entities

      CALL DTERPT(0)

      IF (D2TSLD (CMEM, IMEM, DMEM, IBFS))
     +      CALL D2UNLD (CMEM, IMEM, DMEM, IBFS, IERX)

C     Try to delete the Working storage entity

      CALL D2ERAS (CMEM, IMEM, DMEM, IDWORK, IERX)

      CALL DTERPT(1)

 9999 CONTINUE

      RETURN
      END

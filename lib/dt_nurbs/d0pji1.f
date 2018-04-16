      SUBROUTINE D0PJI1 (CMEM, IMEM, DMEM, IBFC, T, NWORK, U, V, 
     +      IER)

C     PURPOSE:
C        Evaluate a B-spline Curve at T.
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IBFC     Pointer to B-Spline Function Curve entity
C        T        Parameter value on Edge curve
C
C     OUTPUT:
C        NWORK    Amount of working storage used/needed
C        U, V     Paramter values from Curve IBFC(T)
C        IER      Error flag.  If negative, the entity IPJ is not
C                 updated.
C                 0     No errors detected.
C
C                 -4  = IBFC does not point to a valid 
C                       B-spline
C                       Function entity.
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
C                 -999= Unexpected error occurred in a called routine
C
C     CALLS:
C        D0BFFP
C        D2FEBD
C        DTERPT
C        D1DEFE
C        DTERR
C        DTSPVL
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

      INTEGER           IBFC, NWORK, IER
      DOUBLE PRECISION  T, U, V

      INTEGER           IDMY, IERX
      INTEGER           CSTART, CEND, IWKST, IWKEND, IDWORK
      INTEGER           KORD, NRNG
      DOUBLE PRECISION  EV(3)
      CHARACTER         LABELX
      LOGICAL           INITIZ

C     ------

      IER = 0

      IDWORK = 0
      U = 0.0D0
      V = 0.0D0

C     Fetch bounds and lock the B-spline function

      CALL D0BFFP (CMEM, IMEM, DMEM, IBFC, LABELX, IDMY, IDMY, IERX)
      IF (IERX .NE. 0) THEN
         IER = -4
         GOTO 9000
      ENDIF

      CALL D2FEBD (CMEM, IMEM, DMEM, IBFC, CSTART, CEND, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9000
      ENDIF

C     Set NRNG = C(2)

      NRNG = DMEM(CSTART+1)

      IF (NRNG .EQ. 0) THEN
         IER = -17
         GOTO 9000
      ELSEIF (NRNG .LT. 0) THEN
C        Compensate for rational spline
         NRNG = -1-NRNG
      ENDIF

      IF (NRNG .NE. 2) THEN
         IER = -18
         GOTO 9000
      ENDIF

C     Define and lock some working storage

C     Set KORD = C(3)

      KORD = DMEM(CSTART+2)
      IF (KORD .LE. 0) THEN
         IER = -10
         GOTO 9000
      ENDIF

      INITIZ = .FALSE.
      NWORK  = 5 * KORD - 2

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

C     Evaluate C(T)

      EV(1) = 0.0D0
      EV(2) = 0.0D0
      EV(3) = 0.0D0

      CALL DTERPT(0)
      CALL DTSPVL (T, DMEM(CSTART), DMEM(IWKST), NWORK, EV, IERX)
      CALL DTERPT(1)

      IF (IERX .NE. 0) THEN
         IF       (IERX .EQ. -1) THEN
C           C(3) .LE. 0)
            IER = -10
         ELSE IF  (IERX .EQ. -6) THEN
C           C(4) < C(3)
            IER = -11
         ELSE IF  (IERX .EQ. -8) THEN
C           Invalid knot set
            IER = -12
         ELSE IF  (IERX .EQ. -10) THEN
C           Denominator = 0
            IER = -13
         ELSE IF  (IERX .EQ. -38) THEN
C           Interval too small
            IER = -14
         ELSE IF  (IERX .EQ. -50) THEN
C           T out of range
            IER = -15
         ELSE IF  (IERX .EQ. -51) THEN
C           C(1) .NE. 1
            IER = -16
         ELSE IF  (IERX .EQ. -52) THEN
C           C(2) = 0
            IER = -17
         ELSE
            CALL DTERR (1, 'DTSPVL', IERX, 0)
            IER = -999
         ENDIF

         GOTO 9000
      ENDIF

      U = EV(1)
      V = EV(2)

C     Unlock the B-spline Function entity

      CALL D2UNLD (CMEM, IMEM, DMEM, IBFC, IERX)
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

C     Try to unlock the B-spline Curve Function entities

      CALL DTERPT(0)

      IF (D2TSLD (CMEM, IMEM, DMEM, IBFC))
     +      CALL D2UNLD (CMEM, IMEM, DMEM, IBFC, IERX)

C     Try to delete the Working storage entity

      CALL D2ERAS (CMEM, IMEM, DMEM, IDWORK, IERX)

      CALL DTERPT(1)

 9999 CONTINUE

      RETURN
      END


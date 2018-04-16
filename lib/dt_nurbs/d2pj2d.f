      SUBROUTINE D2PJ2D (CMEM, IMEM, DMEM, IPJ, U, V, IER)

C     PURPOSE:
C        Point on Joined 2-Dimensional parametric coordinates.
C        Retrieve the current surface parameter coordinates of
C        Point on Joined IPJ into U and V
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IPJ      Point on Joined Surface entity ID
C
C     OUTPUT:
C        U,V      The parameter values from the entity
C        IER      Error flag.
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = IPJ does not point to a Point on Joined
C                       entity.
C                 -3  = IPJ->ITS = 0.
C
C     CALLS:
C        D1MBAD
C        D0PJFP
C        DTERR
C
C     HISTORY:
C        09Jul92  D. Parsons   Created.
C
C     ------

C     Long name alias:
C        ENTRY D2_POINT_JOINED_2D (CMEM, IMEM, DMEM, IPJ,
C    +      U, V, IER)

      EXTERNAL D1MBAD
      LOGICAL  D1MBAD

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IPJ, ITS, IER
      DOUBLE PRECISION  U, V, DDMY
      CHARACTER         LABELX
      INTEGER           IDMY, IERX

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2PJ2D'/

C     ------

      IER = 0
      U   = 0.0D0
      V   = 0.0D0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

C     Fetch the data from the entity

      CALL D0PJFP (CMEM, IMEM, DMEM, IPJ, LABELX, IDMY, ITS,
     +      DDMY, U, V, DDMY, DDMY, DDMY, IDMY, IERX)
      IF (IERX .NE. 0) THEN
         IER = -2
         GOTO 9000
      ENDIF

      IF (ITS .EQ. 0) THEN
         U = 0.0D0
         V = 0.0D0
         IER = -3
         GOTO 9000
      ENDIF

      GOTO 9999

 9000 CONTINUE

C     Error Handling

      CALL DTERR (1, SUBNAM, IER, 0)

 9999 CONTINUE

      RETURN
      END

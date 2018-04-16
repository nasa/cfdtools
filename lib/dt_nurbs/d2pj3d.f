      SUBROUTINE D2PJ3D (CMEM, IMEM, DMEM, IPJ, X, Y, Z, IER)

C     PURPOSE:
C        Point on Joined 3-Dimensional rectangular coordinates.
C        Retrieve the current coordinates of Point on Joined IPJ
C        into X, Y, and Z.
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
C        X,Y,Z    The rectangular coordinates of the point
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
C        ENTRY D2_POINT_JOINED_3D (CMEM, IMEM, DMEM, IPJ,
C    +      X, Y, Z, IER)

      EXTERNAL D1MBAD
      LOGICAL  D1MBAD

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IPJ, ITS, IER
      DOUBLE PRECISION  X, Y, Z, DDMY
      CHARACTER         LABELX
      INTEGER           IDMY, IERX

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2PJ3D'/

C     ------

      IER = 0
      X = 0.0D0
      Y = 0.0D0
      Z = 0.0D0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

C     Fetch the data from the entity

      CALL D0PJFP (CMEM, IMEM, DMEM, IPJ, LABELX, IDMY, ITS,
     +      DDMY, DDMY, DDMY, X, Y, Z, IDMY, IERX)
      IF (IERX .NE. 0) THEN
         IER = -2
         GOTO 9000
      ENDIF

      IF (ITS .EQ. 0) THEN
         X = 0.0D0
         Y = 0.0D0
         Z = 0.0D0
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

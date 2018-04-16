      SUBROUTINE D2LPIX (CMEM, IMEM, DMEM, ILP, IXLP, IER)

C     PURPOSE:
C        Find the index IXLP of Loop ILP in the list of Loops of
C        the associated Trimmed Surface.
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        ILP      Loop Entity ID
C
C     OUTPUT:
C        IXLP     Loop index into Trimmed Surface entity
C        IER      Error flag.
C                 -1  = Dynamic memory is corrupt or uninitialized.
C                 -2  = ILP does not point to a Loop entity
C
C     CALLS:
C        D1MBAD
C        D0LPFP
C
C     HISTORY:
C        02Jul92  D. Parsons     Created
C
C     ------

C     Long name alias:
C        ENTRY D2_LOOP_INDEX (CMEM, IMEM, DMEM, ILP, IXLP, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      EXTERNAL          D1MBAD
      LOGICAL           D1MBAD

      INTEGER     ILP, IXLP, IER
      INTEGER     ITS, MAXEG, IERX, IHDR
      CHARACTER   LABEL

      IER  = 0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

C     Extract the information

      CALL D0LPFP (CMEM, IMEM, DMEM, ILP, LABEL, ITS, IXLP,
     +      MAXEG, IHDR, IERX)

      IF (IERX .NE. 0) THEN
         IER  = -2
         IXLP = 0
      ENDIF

 9000 CONTINUE

      IF (IER .NE. 0) THEN
         CALL DTERR (1, 'D2LPIX', IER, 0)
      ENDIF

      RETURN
      END

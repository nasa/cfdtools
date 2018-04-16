      SUBROUTINE D2EGIX (CMEM, IMEM, DMEM, IEG, IXEG, IER)

C     PURPOSE:
C        Find the index IXEG of Edge IEG in the list of segments of
C        the associated loop.
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IEG      Edge Entity ID
C
C     OUTPUT:
C        IXEG     Edge index into Loop entity
C        IER      Error flag.
C                 -1  = Dynamic memory is corrupt or uninitialized.
C                 -2  = IEG does not point to a edge entity
C
C     CALLS:
C        D1MBAD
C        D0EGFP
C
C     HISTORY:
C        02Jul92  D. Parsons     Created
C
C     ------

C     Long name alias:
C        ENTRY D2_EDGE_INDEX (CMEM, IMEM, DMEM, IEG, IXEG, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      EXTERNAL          D1MBAD
      LOGICAL           D1MBAD

      INTEGER     IEG, IXEG, IER
      INTEGER     IBFC, ILP, JEG, IERX, IHDR
      CHARACTER   LABEL

      IER  = 0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

C     Extract the information

      CALL D0EGFP (CMEM, IMEM, DMEM, IEG, LABEL, IBFC, ILP, IXEG,
     +      JEG, IHDR, IERX)

      IF (IERX .NE. 0) THEN
         IER  = -2
         IXEG = 0
      ENDIF

 9000 CONTINUE

      IF (IER .NE. 0) THEN
         CALL DTERR (1, 'D2EGIX', IER, 0)
      ENDIF

      RETURN
      END

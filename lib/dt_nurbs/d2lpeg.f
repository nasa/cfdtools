      SUBROUTINE D2LPEG (CMEM, IMEM, DMEM, ILP, IXEG, IEG, IER)

C     PURPOSE:
C        Find the IXEGth Edge of Loop ILP and return its pointer in
C        IEG.
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        ILP      Loop Entity ID
C        IXEG     Edge Index
C
C     OUTPUT:
C        IEG      The IXEGth Edge
C        IER      Error flag.
C                 -1    = Dynamic memory is corrupt or uninitialized.
C                 -2    = ILP does not point to a Loop entity
C                 -3    = IXEG < 1 or IXEG > MAXEG
C                 -999  = Unexpected error in a lower-level routine
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
C        ENTRY D2_LOOP_EDGE (CMEM, IMEM, DMEM, ILP, IXEG, IEG, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      EXTERNAL          D1MBAD
      LOGICAL           D1MBAD

      INTEGER     ILP, IXEG, IEG, IER
      INTEGER     ITS, IXLP, MAXEG, IERX
      INTEGER     JSUB, IELM, IHDR
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
         GOTO 9000
      ENDIF

      IF ((IXEG .LT. 1) .OR. (IXEG .GT. MAXEG)) THEN
         IER = -3
         GOTO 9000
      ENDIF

      JSUB = 3+IXEG

      CALL D2FEEI (CMEM, IMEM, DMEM, ILP, JSUB, IELM, IERX)

      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9000
      ENDIF

      IEG = IELM

 9000 CONTINUE

      IF (IER .NE. 0) THEN
         IF (IER .EQ. -999) THEN
            CALL DTERR (5, 'D2LPEG', IER, 0)
         ELSE
            CALL DTERR (1, 'D2LPEG', IER, 0)
         ENDIF
         IEG = 0
      ENDIF

      RETURN
      END

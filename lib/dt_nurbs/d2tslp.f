      SUBROUTINE D2TSLP (CMEM, IMEM, DMEM, ITS, IXLP, ILP, IER)

C     PURPOSE:
C        Find the IXLPth Loop of Trimmed Surface ITS and return its
C        pointer in ILP.
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        ITS      Trimmed Surface Entity ID
C        IXLP     Edge Index
C
C     OUTPUT:
C        ILP      The IXLPth Loop
C        IER      Error flag.
C                 -1    = Dynamic memory is corrupt or uninitialized.
C                 -2    = ITS does not point to a Trimmed Surface
C                         entity
C                 -3    = IXLP < 1 or IXLP > MAXLP
C                 -999  = Unexpected error in a lower-level routine
C
C     CALLS:
C        D1MBAD
C        D0TSFP
C
C     HISTORY:
C        02Jul92  D. Parsons     Created
C
C     ------

C     Long name alias:
C        ENTRY D2_TRIMMED_SURFACE_LOOP (CMEM, IMEM, DMEM, ITS,
C    +      IXLP, ILP, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      EXTERNAL          D1MBAD
      LOGICAL           D1MBAD

      INTEGER     ITS, IXLP, ILP, IER
      INTEGER     IBFS, IJS, IXTS, IORIEN, MAXLP, IERX
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

      CALL D0TSFP (CMEM, IMEM, DMEM, ITS, LABEL, IBFS, IJS, IXTS,
     +      IORIEN, MAXLP, IHDR, IERX)

      IF (IERX .NE. 0) THEN
         IER  = -2
         GOTO 9000
      ENDIF

      IF ((IXLP .LT. 1) .OR. (IXLP .GT. MAXLP)) THEN
         IER = -3
         GOTO 9000
      ENDIF

      JSUB = 5+IXLP

      CALL D2FEEI (CMEM, IMEM, DMEM, ITS, JSUB, IELM, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9000
      ENDIF

      ILP = IELM

 9000 CONTINUE

      IF (IER .NE. 0) THEN
         IF (IER .EQ. -999) THEN
            CALL DTERR (5, 'D2TSLP', IER, 0)
         ELSE
            CALL DTERR (1, 'D2TSLP', IER, 0)
         ENDIF
         ILP = 0
      ENDIF

      RETURN
      END

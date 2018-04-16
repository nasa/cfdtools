      SUBROUTINE D2JSTS (CMEM, IMEM, DMEM, IJS, IXTS, ITS, IER)

C     PURPOSE:
C        Find the IXTSth Trimmed Surface of Joined Surface IJS
C        and return its pointer in ITS.
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IJS      Joined Surface Entity ID
C        IXTS     Trimmed Surface Index
C
C     OUTPUT:
C        ITS      The IXTSth Trimmed Surface
C        IER      Error flag.
C                 -1    = Dynamic memory is corrupt or uninitialized.
C                 -2    = IJS does not point to a Joined Surface
C                         entity
C                 -3    = IXTS < 1 or IXTS > MAXTS
C                 -999  = Unexpected error in a lower-level routine
C
C     CALLS:
C        D1MBAD
C        D0JSFP
C
C     HISTORY:
C        02Jul92  D. Parsons     Created
C
C     ------

C     Long name alias:
C        ENTRY D2_JOINED_SURFACE_TRIMMED_SURFACE (CMEM, IMEM, DMEM,
C    +      IJS, IXTS, ITS, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      EXTERNAL          D1MBAD
      LOGICAL           D1MBAD

      INTEGER     IJS, IXTS, ITS, IER
      INTEGER     ICLS, MAXTS, IERX
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

      CALL D0JSFP (CMEM, IMEM, DMEM, IJS, LABEL, ICLS, MAXTS,
     +      IHDR, IERX)

      IF (IERX .NE. 0) THEN
         IER  = -2
         GOTO 9000
      ENDIF

      IF ((IXTS .LT. 1) .OR. (IXTS .GT. MAXTS)) THEN
         IER = -3
         GOTO 9000
      ENDIF

      JSUB = 2+IXTS

      CALL D2FEEI (CMEM, IMEM, DMEM, IJS, JSUB, IELM, IERX)

      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9000
      ENDIF

      ITS = IELM

 9000 CONTINUE

      IF (IER .NE. 0) THEN
         IF (IER .EQ. -999) THEN
            CALL DTERR (5, 'D2JSTS', IER, 0)
         ELSE
            CALL DTERR (1, 'D2JSTS', IER, 0)
         ENDIF
         ITS = 0
      ENDIF

      RETURN
      END


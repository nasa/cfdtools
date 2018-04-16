      SUBROUTINE D2TSIX (CMEM, IMEM, DMEM, ITS, IXTS, IER)

C     PURPOSE:
C        Find the index IXTS of Trimmed Surface ITS in the list
C        of Trimmed Surfaces of the associated Joined Surface.
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        ITS      Trimmed Surface Entity ID
C
C     OUTPUT:
C        IXTS     Trimmed Surface index into Loop entity
C        IER      Error flag.
C                 -1  = Dynamic memory is corrupt or uninitialized.
C                 -2  = ITS does not point to a Trimmed Surface entity
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
C        ENTRY D2_TRIMMED_SURFACE_INDEX (CMEM, IMEM, DMEM, ITS,
C    +      IXTS, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      EXTERNAL          D1MBAD
      LOGICAL           D1MBAD

      INTEGER     ITS, IXTS, IER
      INTEGER     IBFS, IJS, IORIEN, MAXLP, IERX, IHDR
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
         IXTS = 0
      ENDIF

 9000 CONTINUE

      IF (IER .NE. 0) THEN
         CALL DTERR (1, 'D2TSIX', IER, 0)
      ENDIF

      RETURN
      END


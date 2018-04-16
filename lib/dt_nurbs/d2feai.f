      SUBROUTINE D2FEAI (CMEM, IMEM, DMEM, IDE, JLO, JHI, IA, IER)

C     PURPOSE:
C        Fetch subarray of integer data from the data space of
C        entity IDE
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IDE      Pointer to entity to fetch
C        JLO      Relative index into the integer space
C        JHI      Relative index into the integer space
C
C     OUTPUT:
C        IA       Integer array to put data into
C        IER      Error flag.  If negative, the subarray is unchanged
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = Null pointer (IDE = 0)
C                 -3  = Garbage value in pointer (points outside valid
C                       range in IMEM)
C                 -4  = Either a deleted entity or a garbage value in
C                       pointer
C                 -5  = Pointer refers to a deleted entity
C                 -8  = JLO <= 0
C                 -9  = JHI > LENI
C                 -10 = JLO > JHI
C
C     CALLS:
C        D0PTR, DTERR
C
C     HISTORY:
C        10Mar92  D. Parsons   Created.
C
C     ------

C     Long name alias:
C        ENTRY D2_FETCH_IMEM_SUBARRAY (CMEM, IMEM, DMEM, IDE, JLO,
C    +         JHI, IA, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*), IDE, JLO, JHI, IA(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IHDR, ITYP, LENI, JBLO, JBHI, IER, I

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2FEAI'/

C     ------

      IER = 0

C     Call the general pointer-check utility routine

      CALL D0PTR (CMEM, IMEM, DMEM, IDE, 0, 0, IHDR, ITYP, IER)
      IF (IER .NE. 0) GOTO 9000


      IF (JLO .LE. 0) THEN
         IER = -8
         GOTO 9000
      ENDIF

      IF (JHI .LT. JLO) THEN
         IER = -10
         GOTO 9000
      ENDIF

C     Get the actual integer array bounds

      JBLO = IMEM(IHDR+2)
      LENI = ABS(IMEM(IHDR+5))
      JBHI = JBLO + LENI - 1

      IF (JHI .GT. LENI) THEN
         IER = -9
         GOTO 9000
      ENDIF

C     Get the Subarray from IMEM

      DO 100 I = JLO, JHI
         IA(I-JLO+1) = IMEM(JBLO+I-1)
  100 CONTINUE

      GOTO 9999

C     Error return

 9000 CONTINUE

      CALL DTERR (1, SUBNAM, IER, 0)

 9999 CONTINUE

      RETURN
      END

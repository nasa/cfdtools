      SUBROUTINE D2STAD (CMEM, IMEM, DMEM, DA, JLO, JHI, IDE, IER)

C     PURPOSE:
C        Store subarray of double precision data into the data space of
C        entity IDE
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        DA       Double precision array to put into DMEM
C        JLO      Relative index into the double precision space
C        JHI      Relative index into the double precision space
C        IDE      Pointer to entity to to store
C
C     OUTPUT:
C        IER      Error flag.  If negative, the space remains unchanged
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = Null pointer (IDE = 0)
C                 -3  = Garbage value in pointer (points outside valid
C                       range in IMEM)
C                 -4  = Either a deleted entity or a garbage value in
C                       pointer
C                 -5  = Pointer refers to a deleted entity
C                 -8  = JLO <= 0
C                 -9  = JHI > LEND
C                 -10 = JLO > JHI
C
C     CALLS:
C        D0PTR, DTERR
C
C     HISTORY:
C        11Mar92  D. Parsons   Created.
C
C     ------

C     Long name alias:
C        ENTRY D2_STORE_DMEM_SUBARRAY (CMEM, IMEM, DMEM, DA, JLO,
C    +         JHI, IDE, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*), IDE, JLO, JHI
      DOUBLE PRECISION  DMEM(*), DA(*)

      INTEGER           IHDR, ITYP, LEND, JBLO, JBHI, IER, I

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2STAD'/

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

C     Get the actual double precision array bounds

      JBLO = IMEM(IHDR+3)
      LEND = ABS(IMEM(IHDR+6))
      JBHI = JBLO + LEND - 1

      IF (JHI .GT. LEND) THEN
         IER = -9
         GOTO 9000
      ENDIF

C     Store the array into DMEM

      DO 100 I = JLO, JHI
         DMEM(JBLO+I-1) = DA(I-JLO+1)
  100 CONTINUE

      GOTO 9999

C     Error return

 9000 CONTINUE

      CALL DTERR (1, SUBNAM, IER, 0)

 9999 CONTINUE

      RETURN
      END

      SUBROUTINE D1STAC (CMEM, IMEM, DMEM, CA, JLO, JHI, IHDR, IER)

C     PURPOSE:
C        Store subarray of character data into the data space of
C        entity IDE
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        CA       Character array to put into CMEM
C        JLO      Relative index into the character space
C        JHI      Relative index into the character space
C        IHDR     Header index of entity to store to
C
C     OUTPUT:
C        IER      Error flag.  If negative, the space remains unchanged
C                 -1 = JLO <= 0
C                 -2 = JHI > LENC
C                 -3 = JLO > JHI
C
C     CALLS:
C        none
C
C     HISTORY:
C        28May93  D. Parsons  Extracted from D2STAC
C        23Mun93  D. Parsons  Fixed possibility of STRLEN > LEN(STR)
C
C     ------

      CHARACTER         CMEM*(*), CA*(*)
      INTEGER           IMEM(*), JLO, JHI
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IHDR, LENC, JBLO, JBHI, IER     
      INTEGER           MINLEN

C     ------

      IER = 0

      IF (JLO .LE. 0) THEN
         IER = -1
         GOTO 9000
      ENDIF

      IF (JHI .LT. JLO) THEN
         IER = -3
         GOTO 9000
      ENDIF

C     Get the actual character array bounds

      JBLO = IMEM(IHDR+1)
      LENC = ABS(IMEM(IHDR+4))
      JBHI = JBLO + LENC - 1

      IF (JHI .GT. LENC) THEN
         IER = -2
         GOTO 9000
      ENDIF

C     Store the array into CMEM

      MINLEN = MIN((JHI-JLO+1),LEN(CA))
      CMEM(JBLO+JLO-1:JBLO+JHI-1) = CA(1:MINLEN)

 9000 CONTINUE

      RETURN
      END

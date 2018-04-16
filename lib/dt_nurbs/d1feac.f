      SUBROUTINE D1FEAC (CMEM, IMEM, DMEM, IHDR, JLO, JHI, CA, IER)

C     PURPOSE:
C        Fetch subarray of character data from the data space of
C        entity IDE
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IHDR     Header index of entity to fetch
C        JLO      Relative index into the character space
C        JHI      Relative index into the character space
C
C     OUTPUT:
C        CA       Character array to put data into
C        IER      Error flag.  If negative, the subarray is unchanged
C                 -8  = JLO <= 0
C                 -9  = JHI > LENC
C                 -10 = JLO > JHI
C
C     CALLS:
C        (none)
C
C     HISTORY:
C        29Jun93  D. Parsons  Extracted from D2FEAC.
C
C     ------
      CHARACTER         CMEM*(*), CA*(*)
      INTEGER           IMEM(*), JLO, JHI
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IHDR, LENC, JBLO, JBHI, IER
      INTEGER           LENSTR

C     ------

      IER = 0 
      CA  = ' '
      
      LENSTR = LEN(CA)

      IF (JLO .LE. 0) THEN
         IER = -8
         GOTO 9000
      ENDIF

      IF (JHI .LT. JLO) THEN
         IER = -10
         GOTO 9000
      ENDIF

C     Get the actual character array bounds

      JBLO = IMEM(IHDR+1)
      LENC = ABS(IMEM(IHDR+4))
      JBHI = JBLO + LENC - 1

      IF (JHI .GT. LENC) THEN
         IER = -9
         GOTO 9000
      ENDIF

C     Get the Subarray from CMEM

      CA(1:MIN(LENSTR,JHI-JLO+1)) = 
     +      CMEM(JBLO+JLO-1:JBLO+MIN(JLO+LENSTR-1,JHI)-1)

 9000 CONTINUE
 
      IF (IER .LT. 0) CA = ' '

      RETURN
      END

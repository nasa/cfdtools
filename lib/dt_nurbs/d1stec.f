      SUBROUTINE D1STEC (CMEM, IMEM, DMEM, CELM, JSUB, IHDR, IER)

C     PURPOSE:
C        Store element of character data into the data space of
C        entity IDE
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        CELM     Character to put into CMEM
C        JSUB     Relative index into the character space
C        IHDR     Header index of entity to store
C
C     OUTPUT:
C        IER      Error flag.  If negative, the space remains unchanged
C                 -8 =  JSUB <= 0
C                 -9 =  JSUB > LENC
C
C     CALLS:
C        none
C
C     HISTORY:
C        12Mar92  D. Parsons   Created.
C        09Jun93  D. Parsons   Extracted from D2STEC
C
C     ------

      CHARACTER         CMEM*(*), CELM
      INTEGER           IMEM(*), JSUB
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IHDR, LENC, JBLO, JBHI, IER

C     ------

      IER = 0

      IF (JSUB .LE. 0) THEN
         IER = -8
         GOTO 9000
      ENDIF

C     Get the actual character array bounds

      JBLO = IMEM(IHDR+1)
      LENC = ABS(IMEM(IHDR+4))
      JBHI = JBLO + LENC - 1

      IF (JSUB .GT. LENC) THEN
         IER = -9
         GOTO 9000
      ENDIF

C     Store the array into CMEM

      CMEM(JBLO+JSUB-1:JBLO+JSUB-1) = CELM

 9000 CONTINUE

      RETURN
      END

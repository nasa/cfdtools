      SUBROUTINE D1FEEC (CMEM, IMEM, DMEM, IHDR, JSUB, CELM, IER)

C     PURPOSE:
C        Fetch element of character data from the data space of
C        entity IDE
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IHDR     Header index of entity to fetch
C        JSUB     Relative index into the character space
C
C     OUTPUT:
C        CELM     Character to put data into
C        IER      Error flag.  If negative, CELM is unchanged
C                 -8 =  JSUB <= 0
C                 -9 =  JSUB > LENC
C
C     CALLS:
C        (none)
C
C     HISTORY:
C        11Mar92  D. Parsons   Created.
C        29Jun93  D. Parsons   Extracted from D2FEEC
C
C     ------

      CHARACTER         CMEM*(*), CELM
      INTEGER           IMEM(*), JSUB
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IHDR, LENC, JBLO, IER

C     ------

      IER = 0

      IF (JSUB .LE. 0) THEN
         IER = -8
         GOTO 9000
      ENDIF

C     Get the actual character array bounds

      JBLO = IMEM(IHDR+1)
      LENC = ABS(IMEM(IHDR+4))

      IF (JSUB .GT. LENC) THEN
         IER = -9
         GOTO 9000
      ENDIF

C     Get the element from CMEM

      CELM = CMEM(JBLO+JSUB-1:JBLO+JSUB-1)

 9000 CONTINUE

      RETURN
      END

      SUBROUTINE D1STED (CMEM, IMEM, DMEM, DELM, JSUB, IHDR, IER)

C     PURPOSE:
C        Store element of double precision data into the data space of
C        entity IDE
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        DELM     Double precision to put into DMEM
C        JSUB     Relative index into the double precision space
C        IHDR     Pointer index of entity to store
C
C     OUTPUT:
C        IER      Error flag.  If negative, the space remains unchanged
C                 -8 =  JSUB <= 0
C                 -9 =  JSUB > LEND
C
C     CALLS:
C        (none)
C
C     HISTORY:
C        12Mar92  D. Parsons   Created.
C        29Jun93  D. Parsons   Extracted from D2STED
C
C     ------
      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*), JSUB
      DOUBLE PRECISION  DMEM(*), DELM

      INTEGER           IHDR, LEND, JBLO, JBHI, IER

C     ------

      IER = 0

      IF (JSUB .LE. 0) THEN
         IER = -8
         GOTO 9000
      ENDIF

C     Get the actual double precision array bounds

      JBLO = IMEM(IHDR+3)
      LEND = ABS(IMEM(IHDR+6))
      JBHI = JBLO + LEND - 1

      IF (JSUB .GT. LEND) THEN
         IER = -9
         GOTO 9000
      ENDIF

C     Store the element into DMEM

      DMEM(JBLO+JSUB-1) = DELM

 9000 CONTINUE

      RETURN
      END

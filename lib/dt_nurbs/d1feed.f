      SUBROUTINE D1FEED (CMEM, IMEM, DMEM, IHDR, JSUB, DELM, IER)

C     PURPOSE:
C        Fetch element of double precision data from the data
C        space of entity IDE
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IHDR     Header index of entity to fetch
C        JSUB     Relative index into the double precision space
C
C     OUTPUT:
C        DELM     Double precision to put data into
C        IER      Error flag.  If negative, DELM is unchanged
C                 -8 =  JSUB <= 0
C                 -9 =  JSUB > LEND
C
C     CALLS:
C        (none)
C
C     HISTORY:
C        11Mar92  D. Parsons   Created.
C        29Jun93  D. Parsons   Extracted guts from D2FEED.
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

C     Get the element from DMEM

      DELM = DMEM(JBLO+JSUB-1)

C     Error return

 9000 CONTINUE

      RETURN
      END

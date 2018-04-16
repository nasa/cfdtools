      SUBROUTINE D1FEEI (CMEM, IMEM, DMEM, IHDR, JSUB, IELM, IER)

C     PURPOSE:
C        Fetch element of integer data from the data space of
C        entity with header index IHDR
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IHDR     Header index on entity
C        JSUB     Relative index into the integer space
C
C     OUTPUT:
C        IELM     Integer to put data into
C        IER      Error flag.  If negative, IELM is unchanged
C                 -1 =  JSUB <= 0
C                 -2 =  JSUB > LENI
C
C     CALLS:
C        none
C
C     HISTORY:
C        28May93  D. Parsons  Extraced from D2FEEI
C
C     ------

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*), JSUB, IELM
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IHDR, LENI, JBLO, JBHI, IER

C     ------

      IER = 0

      IF (JSUB .LE. 0) THEN
         IER = -1
         GOTO 9000
      ENDIF

C     Get the actual integer array bounds

      JBLO = IMEM(IHDR+2)
      LENI = ABS(IMEM(IHDR+5))
      JBHI = JBLO + LENI - 1

      IF (JSUB .GT. LENI) THEN
         IER = -2
         GOTO 9000
      ENDIF

C     Get the element from IMEM

      IELM = IMEM(JBLO+JSUB-1)

 9000 CONTINUE

      RETURN
      END

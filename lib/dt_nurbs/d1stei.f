      SUBROUTINE D1STEI (CMEM, IMEM, DMEM, IELM, JSUB, IHDR, IER)

C     PURPOSE:
C        Store element of integer data into the data space of
C        entity IDE
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IELM     Integer to put into IMEM
C        JSUB     Relative index into the integer space
C        IHDR     Header index to entity to store
C
C     OUTPUT:
C        IER      Error flag.  If negative, the space remains unchanged
C                 -1 =  JSUB <= 0
C                 -2 =  JSUB > LENI
C
C     CALLS:
C        none
C
C     HISTORY:
C        01Jun93  D. Parsons   Extracted from D2STEI
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

C     Store the element into IMEM

      IMEM(JBLO+JSUB-1) = IELM

 9000 CONTINUE
      RETURN
      END

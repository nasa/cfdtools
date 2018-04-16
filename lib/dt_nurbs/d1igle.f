      SUBROUTINE D1IGLE (CMEM, IMEM, DMEM, IGEHDR, NPAR, NSTR, NSTRCH, 
     +                   NREAL, IER)

C     PURPOSE:
C        Fetch IGES Entity Length  (shadow routine)
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IGEHDR   Header index of IGES entity
C
C     OUTPUT:
C        NPAR     Number of Parameters
C        NSTR     Number of string parameters to allocate
C        NSTRCH   Total number of characters in strings
C        NREAL    Number of double precision parameters
C        IER      Error flag.
C                 -3  = IGE entity has no associated character sequence 
C                       entity
C
C     CALLS:
C        D1FEEI   Fetch Integer Element
C
C
C     HISTORY:
C        17Jun93  D. Parsons   Created.
C
C     ------

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

C     Get the Reserved Entity Type ID Numbers

C =====================================================================
C
C     INCLUDE FILE DITYID.INC
C
C     Reserved Entity Type ID Numbers
C
C     HISTORY
C        12May92     D. Parsons     created
C
C                          Character String Sequence
      INTEGER     ENTCQ
      PARAMETER  (ENTCQ = 249)

C =====================================================================

      INTEGER         IGEHDR, NPAR, NSTR, NSTRCH, NREAL, IER
      INTEGER         ICQHDR, ITYP, ICQ, IERX
      CHARACTER*6     SUBNAM
      DATA            SUBNAM /'D2IGLE'/ 
      
      IER    = 0
      IERX   = 0
      NPAR   = 0
      NSTR   = 0
      NSTRCH = 0
      NREAL  = 0

C     Fetch the NPAR from the IMEM space.
      
      CALL D1FEEI (CMEM, IMEM, DMEM, IGEHDR, 3, NPAR, IERX) 
      IF (IERX .NE. 0) NPAR = 0

C     Compute the NSTR from the length of the ICQ IMEM
      
      CALL D1FEEI (CMEM, IMEM, DMEM, IGEHDR, 2, ICQ, IERX)
      IF (IERX .NE. 0) THEN
          IER = -3
          GOTO 9000
      ENDIF
      
      CALL D0PTR (CMEM, IMEM, DMEM, ICQ, ENTCQ, 0, ICQHDR,
     +    ITYP, IERX)
      IF (IERX .NE. 0) THEN
          IER = -3
          GOTO 9000
      ENDIF

      NSTR = ABS(IMEM(ICQHDR+5)) - 1
      
C     Compute the NSTRCH from the length of the ICQ CMEM

      NSTRCH = ABS(IMEM(ICQHDR+4))
      
C     Compute the NREAL from the length of the DMEM

      NREAL = ABS(IMEM(IGEHDR+6))
      

 9000 CONTINUE

      IF (IER .NE. 0) THEN
         NPAR   = 0
         NSTR   = 0
         NSTRCH = 0
         NREAL  = 0
      ENDIF

      RETURN
      END

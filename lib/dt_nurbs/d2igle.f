      SUBROUTINE D2IGLE (CMEM, IMEM, DMEM, IGE, NPAR, NSTR, NSTRCH, 
     +                   NREAL, IER)

C     PURPOSE:
C        Fetch IGES Entity Length
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IGE      Pointer to IGES entity
C
C     OUTPUT:
C        NPAR     Number of Parameters
C        NSTR     Number of string parameters to allocate
C        NSTRCH   Total number of characters in strings
C        NREAL    Number of double precision parameters
C        IER      Error flag.  If negative, the entity is not defined.
C                 -1  = Dynamic memory is corrupt or uninitialized.
C                 -2  = IGE does not point to a valid IGES Entity  
C                 -3  = IGE entity has no associated character sequence 
C                       entity
C
C     CALLS:
C        D0PTR    Pointer check
C        D1FEEI   Fetch Integer Element
C        DTERR    Error handler
C
C
C     HISTORY:
C        17Jun93  D. Parsons   Created.
C
C     ------

C     Long name alias:
C        ENTRY D2_IGES_ENTITY_LENGTH (CMEM, IMEM, DMEM, IGE, NPAR, 
C    +                   NSTR, NSTRCH, NREAL, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      EXTERNAL D1MBAD
      LOGICAL  D1MBAD

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
C                          IGES Entity Data
      INTEGER     ENTGE
      PARAMETER  (ENTGE = 238)

C =====================================================================

      INTEGER         IGE, NPAR, NSTR, NSTRCH, NREAL, IER
      INTEGER         IGEHDR, ITYP, IERX
      CHARACTER*6     SUBNAM
      DATA            SUBNAM /'D2IGLE'/ 
      
      IER    = 0
      IERX   = 0
      NPAR   = 0
      NSTR   = 0
      NSTRCH = 0
      NREAL  = 0

C     Error checking

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF
      
      CALL D0PTR (CMEM, IMEM, DMEM, IGE, ENTGE, 0, IGEHDR,
     +    ITYP, IERX)
      IF (IERX .NE. 0) THEN
          IER = -2
          GOTO 9000
      ENDIF 
      
      CALL D1IGLE (CMEM, IMEM, DMEM, IGEHDR, NPAR, NSTR, NSTRCH, 
     +                   NREAL, IER)      

 9000 CONTINUE

      IF (IER .NE. 0) THEN
         CALL DTERR (1, SUBNAM, IER, 0)
         NPAR   = 0
         NSTR   = 0
         NSTRCH = 0
         NREAL  = 0
      ENDIF

      RETURN
      END

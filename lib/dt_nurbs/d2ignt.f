      SUBROUTINE D2IGNT (CMEM, IMEM, DMEM, IGI, ITYPE, JDE, IGE, IER)
      
C   Fetch the next IGE entity pointer to an entity of type ITYPE from
C   IGES Index Entity IGI, beginning with position JDE.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM     Dynamically managed character data
C     IMEM     Dynamically managed integer data
C     DMEM     Dynamically managed double precision data
C
C   INPUT:
C     IGI      Pointer to the IGES Index Entity
C     ITYPE    IGES Entity type number to search for
C     JDE      Starting Directry Entry pointer in IGI.  If JDE is
C              not an ODD integer, the search will begin with the
C              next ODD integer (JDE+1).
C
C   OUTPUT: 
C     IGE      Pointer to the next IGES Entity of type ITYPE.
C              (0 if none found)
C     JDE      Index position of IGE in IGI, JDE MUST BE ODD.
C     IER      Returned error value.  Zero implies no errors.  
C              -1  = Dynamic memory is corrupt or uninitialized
C              -2  = IGI is a null pointer
C              -3  = Garbage value found in IGI
C              -4  = Ambiguous--garbage pointer or deleted entity
C              -5  = IGI points to a deleted entity
C              -6  = IGI does not point to an IGES Index Entity
C                    
C   CALLS:
C     D0PTR
C     DTERR
C
C   HISTORY:
C     9/15/93  D. Parsons  Created
C
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_NEXT_TYPE (CMEM, IMEM, DMEM, IGI, ITYPE, JDE, 
C    +   IGE, IER)

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      
      INTEGER     ENTGI
      PARAMETER  (ENTGI = 237)

      INTEGER  IGI, ITYPE, JDE, IGE, IER
      INTEGER  IGIHDR, ITYP, IERX, INDEX 
      INTEGER  NSL, NDE
      
      CHARACTER*6    SUBNAM
      DATA SUBNAM   /'D2IGNT'/
      
C****

      IER  = 0
      IERX = 0
      IGE  = 0
      
C     Check pointer and type
      CALL D0PTR (CMEM, IMEM, DMEM, IGI, ENTGI, 0, IGIHDR, ITYP, IER)
      
      IF (IER .NE. 0) GOTO 9900
      
C     Find out the number of Directory Entry Entities (NDE)

      CALL D1IGLI (CMEM, IMEM, DMEM, IGIHDR, NSL, NDE, IERX)
      
      JDE = MAX(1,JDE)
      IF (MOD(JDE,2) .EQ. 0) JDE = JDE+1
      
      DO 100 INDEX = JDE, NDE, 2
         CALL D2IGFI (CMEM, IMEM, DMEM, IGI, 'I', INDEX, IGE, IERX)
         IF (IGE .NE. 0) THEN
            CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'D', 1, ITYP, IERX)
         ELSE
            ITYP = 0
         ENDIF
         IF (ITYP .EQ. ITYPE) GOTO 200
  100 CONTINUE
  
C        Not found

         IGE = 0
         
  200 CONTINUE
  
      JDE = INDEX

C     Error reporting section

 9900 CONTINUE

      IF (IER .NE. 0) THEN
          CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      RETURN
      END

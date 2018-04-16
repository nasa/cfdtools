      LOGICAL FUNCTION D0IGDO (ITYPE, IFORM, IGDO, IGHI)
C
C     Stupid little function to decide if this is an entity to "DO"
C
C     INPUT
C        ITYPE    The entity type
C        IGDO     The level we want to "DO"
C        IGHI     The maximum level encountered so far
C
C     OUTPUT
C        D0IGDO   .TRUE. means do it, .FALSE. means don't
C        IGHI     The MAX of IGHI and the level of this entity.
C
C     HISTORY
C        8/24/93  D. Parsons     Created

      INTEGER  ITYPE, IFORM, IGDO, IGHI
         
      D0IGDO = .FALSE.
      
C     NASA-IGES-NURBS Entity Types

      IF ((ITYPE .EQ.   0) .OR.
     +    (ITYPE .EQ. 102) .OR.
     +    (ITYPE .EQ. 124 .AND. IFORM .GE. 0 .AND. IFORM .LE. 1) .OR.
     +    (ITYPE .EQ. 126) .OR.
     +    (ITYPE .EQ. 128) .OR.
     +    (ITYPE .EQ. 141) .OR.
     +    (ITYPE .EQ. 142) .OR.
     +    (ITYPE .EQ. 143) .OR.
     +    (ITYPE .EQ. 212 .AND. IFORM .EQ. 0) .OR.
     +    (ITYPE .EQ. 314) .OR.
     +    (ITYPE .EQ. 402 .AND. (IFORM .EQ. 1 .OR. IFORM .EQ. 7 .OR.
     +                           IFORM .EQ. 14 .OR. IFORM .EQ. 15)) .OR.
     +    (ITYPE .EQ. 406 .AND. IFORM .EQ. 15)) THEN
         IGHI = MAX(IGHI,1)
         IF (IGDO .GE. 1) D0IGDO = .TRUE.
         
C        NASA-IGES Entity Types

      ELSE IF (
     +    (ITYPE .EQ. 100) .OR.
     +    (ITYPE .EQ. 104) .OR.
     +    (ITYPE .EQ. 106 .AND. IFORM .GE. 1 .AND. IFORM .LE. 3) .OR.
     +    (ITYPE .EQ. 110) .OR.
     +    (ITYPE .EQ. 116) .OR.
     +    (ITYPE .EQ. 308) .OR.
     +    (ITYPE .EQ. 408)) THEN
         IGHI = MAX(IGHI,2)
         IF (IGDO .GE. 2) D0IGDO = .TRUE.
         
C     Unrestricted Entity Types  
      
      ELSE
         IGHI = MAX(IGHI,3)
         IF (IGDO .GE. 3) D0IGDO = .TRUE.
      ENDIF   
         
      RETURN
      END

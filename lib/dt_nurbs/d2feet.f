      SUBROUTINE D2FEET (CMEM, IMEM, DMEM, IDE, IDTYP, IER)

C     PURPOSE:
C        Fetch Entity Type ID number
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IDE      Pointer to entity to fetch
C
C     OUTPUT:
C        IDTYP    Requested entity type ID number
C        IER      Error flag.  If negative, the entity pointer
C                 is invalid.
C                 -1 = Dynamic memory is corrupt or uninitialized.
C                 -2 = Null pointer
C                 -3 = IDE does not hold a pointer
C                    = IDE is either not a pointer or points to a deleted
C                 -4   entity (ambiguous case)
C                 -5 = IDE points to a deleted entity
C
C     CALLS:
C        D0PTR
C
C     HISTORY:
C        12Feb92  D. Parsons   Created.
C
C     ------

C     Long name alias:
C        ENTRY D2_FETCH_ENTITY_TYPE (CMEM, IMEM, DMEM, IDE, IDTYP, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*), IDE
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IHDR, ITYP, IER

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2FEET'/

C     ------

C     Decode pointer and fetch the entity type

      CALL D0PTR (CMEM, IMEM, DMEM, IDE, 0, 0, IHDR, ITYP, IER)

 9900 CONTINUE

      IF (IER .EQ. 0) THEN
         IDTYP = ITYP
      ELSE
         IDTYP = 0
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      RETURN
      END

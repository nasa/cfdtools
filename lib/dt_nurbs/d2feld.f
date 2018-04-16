      SUBROUTINE D2FELD (CMEM, IMEM, DMEM, IDE, LEND, IER)

C     PURPOSE:
C        Fetch Length of Entity's double precision space
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
C        LEND     Length of the entity's double precision space
C        IER      Error flag.  If negative, the length is set to 0.
C                 -1 =  Dynamic memory is corrupt or uninitialized.
C                 -2 =  Null pointer (IDE = 0)
C                 -3 =  Garbage value in pointer (points outside valid
C                       range in IMEM)
C                 -4 =  Either a deleted entity or a garbage value in
C                       pointer
C                 -5 =  Pointer refers to a deleted entity
C
C     CALLS:
C        D0PTR
C
C     HISTORY:
C        27Feb92  D. Parsons   Created.
C        09Apr92  D. Parsons  Add check of dynamic memory.
C
C     ------

C     Long name alias:
C        ENTRY D2_FETCH_DMEM_LENGTH (CMEM, IMEM, DMEM, IDE, LEND, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*), IDE
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IHDR, ITYP, LEND, IER

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2FELD'/

C     ------

      IER = 0

C     Call the general pointer-check utility routine

      CALL D0PTR (CMEM, IMEM, DMEM, IDE, 0, 0, IHDR, ITYP, IER)

 9900 CONTINUE

      IF (IER .EQ. 0) THEN
         LEND = ABS(IMEM(IHDR+6))
      ELSE
         LEND = 0
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      RETURN
      END

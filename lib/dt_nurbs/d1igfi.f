      SUBROUTINE  D1IGFI (CMEM, IMEM, DMEM, IHDR, IGLBL, IPAR,
     +       IVAL, IER)

C     Fetch an integer element from an IGES entity or an IGES Index
C     entity.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IHDR    Pointer to an IGES entity or IGES global entity.
C     IGLBL   Global entity flag.
C             .TRUE.    IHDR points to an IGES global entity.
C             .FALSE.   IHDR points to a normal IGES entity.
C     IPAR    Parameter number to fetch.
C   OUTPUT:
C     IVAL    Integer element fetched.
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             -99 = Entity IHDR is internally inconsistent.
C
C *** WARNING:  NO INPUT ARGUMENT ERROR CHECKING IS DONE HERE ***
C
C   CALLS:
C     D1FEEI   Fetch an integer element
C     D1FEEC   Fetch a character element
C
C   HISTORY:
C     7/26/93  D. Parsons   Created
C     8/26/93  D. Parsons   Adjusted PAR values to being with 0 (Entity
C                           Type ID)
C*****
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)

      INTEGER        IHDR, IPAR, IVAL, IER
      INTEGER        IMEMX, ITCDX, IERX
      CHARACTER*(1)  TYPCOD
      LOGICAL        IGLBL

C ***
      IER = 0

C     Set index into IMEM (IMEMX) and index into TYPCODs (ITCDX)

      IF (IGLBL) THEN
         IMEMX = IPAR+3
         ITCDX = IPAR-2
      ELSE
         IMEMX = IPAR+19
         ITCDX = IPAR+24
      ENDIF

      IF (IPAR .EQ. 0) THEN

C        Entity type Number

         CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 5, IVAL, IERX)
         IF (IERX .NE. 0) GOTO 9900

      ELSE

C        Get the type-code.

         CALL D1FEEC (CMEM, IMEM, DMEM, IHDR, ITCDX, TYPCOD,
     +      IERX)
         IF (IERX .NE. 0) GOTO 9900

         IF (TYPCOD .EQ. 'O') THEN
C           Nothing to fetch
            IVAL = 0
            GOTO 9999
         ENDIF

C        Fetch the value

         CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, IMEMX, IVAL, IERX)
         IF (IERX .NE. 0) GOTO 9900

      ENDIF

      GOTO 9999

 9900 CONTINUE
      IER = -99

 9999 CONTINUE

      RETURN
      END

      SUBROUTINE  D1IGSI (CMEM, IMEM, DMEM, IHDR, IGLBL, TYPCOD, IPAR,
     +       IVAL, IER)

C     Store an integer element into an IGES entity or an IGES Index
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
C     TYPCOD  IGES-context data type code.
C     IPAR    Parameter number to replace.
C     IVAL    Integer element to store.
C   OUTPUT:
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             -99 = Entity IHDR is internally inconsistent.
C
C *** WARNING:  NO INPUT ARGUMENT ERROR CHECKING IS DONE HERE ***
C
C   CALLS:
C     D1STEI   Store an integer element
C     D1FEEC   Fetch a character element
C     D1STEC   Store a character element
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
      INTEGER        IMEMX, ITCDX, IERX, IOVAL, ICQ, ICQHDR, ITYP
      CHARACTER*(1)  TYPCOD, OTCOD
      LOGICAL        IGLBL

C =====================================================================
C
C     INCLUDE FILE DITYID.INC
C
C     Reserved Entity Type ID Numbers
C
C     HISTORY
C        12May92     D. Parsons     created
C
C                          Character Sequence Entity
      INTEGER     ENTCQ
      PARAMETER  (ENTCQ = 249)
C =====================================================================
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
         CALL D1STEI (CMEM, IMEM, DMEM, IVAL, 5, IHDR, IERX)
         IF (IERX .NE. 0) GOTO 9900

      ELSE

C        Get the old type-code.
         CALL D1FEEC (CMEM, IMEM, DMEM, IHDR, ITCDX, OTCOD,
     +      IERX)
         IF (IERX .NE. 0) GOTO 9900

C        Get the old value.
         CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, IMEMX, IOVAL,
     +      IERX)
         IF (IERX .NE. 0) GOTO 9900

C        If we are changing TYPCOD, free the old
         IF  (OTCOD .EQ. 'C') THEN
C           This was a Character string--free it, ignore any errors
            CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 2, ICQ, IERX)
            IF (IERX .EQ. 0)
     +         CALL D0PTR (CMEM, IMEM, DMEM, ICQ, ENTCQ, 0,
     +               ICQHDR, ITYP, IERX)
            IF (IERX .EQ. 0)
     +         CALL D1CQDL (CMEM, IMEM, DMEM, ICQHDR, IOVAL, IERX)

         ELSE IF ((OTCOD .EQ. 'E') .OR.
     +       (OTCOD .EQ. 'D')) THEN
C           This was an index into DMEM--free it, ignore any errors
            CALL D0IGFD (CMEM, IMEM, DMEM, IHDR, IOVAL, IERX)

         ELSE IF ((OTCOD .EQ. 'T') .OR.
     +       (OTCOD .EQ. 'R')) THEN
C           This was an ICQ pointer--free it, ignore any errors
            CALL DTERPT(0)
            CALL D2ERAS (CMEM, IMEM, DMEM, IOVAL, IERX)
            CALL DTERPT(1)

         ENDIF

C        Store the value

         IF (TYPCOD .NE. 'O') THEN
            CALL D1STEI (CMEM, IMEM, DMEM, IVAL, IMEMX, IHDR, IERX)
            IF (IERX .NE. 0) GOTO 9900
         ENDIF

C        Update the type-code

         IF (TYPCOD .NE. OTCOD) THEN
            CALL D1STEC (CMEM, IMEM, DMEM, TYPCOD, ITCDX, IHDR,
     +         IERX)
            IF (IERX .NE. 0) GOTO 9900
         ENDIF
      ENDIF

      GOTO 9999

 9900 CONTINUE
      IER = -99

 9999 CONTINUE

      RETURN
      END

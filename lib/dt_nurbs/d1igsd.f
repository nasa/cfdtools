      SUBROUTINE  D1IGSD (CMEM, IMEM, DMEM, IHDR, IGLBL, TYPCOD, IPAR,
     +       DVAL, NEED, IER)

C     Store a double precision element into an IGES entity or an IGES Index
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
C     DVAL    Double precision element to store.
C   OUTPUT:
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             -12 = Insufficient storage for this item.
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

      INTEGER        IHDR, IPAR, IER
      DOUBLE PRECISION DVAL
      INTEGER        IMEMX, ITCDX, DMEMX, IERX, IVAL, ICQHDR, ITYP
      INTEGER        MPAR, MREAL, MSTR, MSTRCH, NEED
      CHARACTER*(1)  TYPCOD, OTCOD
      LOGICAL        IGLBL

      INTEGER ENTCQ
      PARAMETER (ENTCQ=249)

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

C     Get the old type-code.

      CALL D1FEEC (CMEM, IMEM, DMEM, IHDR, ITCDX, OTCOD,
     +   IERX)
      IF (IERX .NE. 0) GOTO 9900

C     Get the old index or value
      CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, IMEMX, IVAL, IERX)
      IF (IERX .NE. 0) GOTO 9900

C     Get the index into DMEM
      IF ((OTCOD .NE. 'E') .AND.
     +     (OTCOD .NE. 'D')) THEN
C        We are changing TYPCOD, free the old
         IF  (OTCOD .EQ. 'C') THEN
C           This was a Character string--free it, ignore any errors
            CALL D0PTR (CMEM, IMEM, DMEM, IVAL, ENTCQ, 0, ICQHDR, ITYP,
     +         IERX)
            IF (IERX .EQ. 0)
     +         CALL D1CQDL (CMEM, IMEM, DMEM, ICQHDR, IVAL, IERX)

         ELSE IF ((OTCOD .EQ. 'T') .OR.
     +       (OTCOD .EQ. 'R')) THEN
C           This was an ICQ pointer--free it, ignore any errors
            CALL DTERPT(0)
            CALL D2ERAS (CMEM, IMEM, DMEM, IVAL, IERX)
            CALL DTERPT(1)
         ENDIF
         CALL D0IGND (CMEM, IMEM, DMEM, IHDR, DMEMX, IERX)
         IF (IERX .NE. 0) THEN
            IF (IERX .EQ. -2) THEN
C              All allocated DMEM has been assigned
C              Try to expand entity
               MPAR   = 0
               MSTR   = 0
               MSTRCH = 0
               MREAL  = 1
               CALL D1IGXE (CMEM, IMEM, DMEM, IHDR, MPAR, MSTR, MSTRCH,
     +                   MREAL, NEED, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -12
                  GOTO 9999
               ENDIF
               CALL D0IGND (CMEM, IMEM, DMEM, IHDR, DMEMX, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -12
                  GOTO 9999
               ENDIF
            ELSE
               GOTO 9900
            ENDIF
         ENDIF
C        Update the IMEM index
         CALL D1STEI (CMEM, IMEM, DMEM, DMEMX, IMEMX, IHDR, IERX)
         IF (IERX .NE. 0) GOTO 9900
      ELSE
         DMEMX = IVAL
      ENDIF

C     Store the value
      IF (TYPCOD .NE. 'O') THEN
         CALL D1STED (CMEM, IMEM, DMEM, DVAL, DMEMX, IHDR, IERX)
         IF (IERX .NE. 0) GOTO 9900
      ENDIF

C     Update the type-code
      IF (TYPCOD .NE. OTCOD) THEN
         CALL D1STEC (CMEM, IMEM, DMEM, TYPCOD, ITCDX, IHDR,
     +      IERX)
         IF (IERX .NE. 0) GOTO 9900
      ENDIF

      GOTO 9999

 9900 CONTINUE
      IER = -99

 9999 CONTINUE

      RETURN
      END

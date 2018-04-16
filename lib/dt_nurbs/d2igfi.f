      SUBROUTINE D2IGFI (CMEM, IMEM, DMEM, IGIE, SECT, ITEM, IVAL, IER)

C     Fetch an integer element from an IGES entity or an IGES Index
C     entity.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IGIE    Pointer to an IGES entity or IGES Index entity.
C     SECT    IGES "section" to which the request is directed.
C     ITEM    Index into the section to which the request is directed.
C   OUTPUT:
C     IVAL    Fetched integer element.
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             -1  = Dynamic memory is corrupt or uninitialized
C             -2  = IGIE is a null pointer
C             -3  = Garbage value found in IGIE
C             -4  = Ambiguous--garbage pointer or deleted entity
C             -5  = IGIE points to a deleted entity
C             -6  = IGIE does not point to an IGES entity or an
C                   IGES Index entity.
C             -7  = SECT is invalid for this type of entity (IGES
C                   or Index).
C             -8  = ITEM is out of range with respect to this IGES
C                   data structure and store request.
C             -99 = Entity IGIE is internally inconsistent.
C
C
C   CALLS:
C     D0PTR    Pointer check
C     D1IGLI   Fetch length parameters of IGES Index Entity
C     D1FEEI   Fetch an integer element
C     D1IGFI   Fetch an integer from parameter storage
C     D1IGLE   Fetch length parameters of IGES Entity
C     DTERR    Error handler
C
C   HISTORY:
C     7/26/93  D. Parsons   Created
C     8/26/93  D. Parsons   Adjusted PAR values to being with 0 (Entity
C                           Type ID)
C     9/28/93  P. Kraushar  Reorganize to permit access to Index info
C                           from Entity
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_FETCH_INTEGER (CMEM, IMEM, DMEM, TYPCOD, IVAL,
C     +                          IGIE, SECT, ITEM, IER)

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)

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

C                          IGES Index
      INTEGER     ENTGI
      PARAMETER  (ENTGI = 237)

C =====================================================================

      CHARACTER*1    SECT
      INTEGER        IGIE, ITEM, IVAL, IER

      INTEGER        NSL, NDE, NPAR, NSTR, NSTRCH, NREAL, IERX
      INTEGER        IHDR, IGHDR, IGLBL, ITYP, INDEX, IGI
      LOGICAL        IGLBF


      CHARACTER*6    SUBNAM
      DATA SUBNAM   /'D2IGFI'/

C****
      IER  = 0
      IVAL = 0

C     Check pointer and type, then locate relevant data spaces
      CALL D0PTR (CMEM, IMEM, DMEM, IGIE, 0, 0, IHDR, ITYP, IER)
      IF (IER .LT. 0) GOTO 9900

      IF ((ITYP .NE. ENTGE) .AND. (ITYP .NE. ENTGI)) THEN
         IER = -6
         GOTO 9900
      ENDIF

      IF (SECT .EQ. 'P') THEN

C        *************************
C        *** Parameter section ***
C        *************************

         IF (ITYP .NE. ENTGE) THEN
            IER = -7
            GOTO 9900
         ENDIF

         CALL D1IGLE (CMEM, IMEM, DMEM, IHDR, NPAR, NSTR, NSTRCH,
     +                NREAL, IERX)
         IF (IERX .NE. 0) THEN
            IER = -99
            GOTO 9900
         ENDIF

         IF ((ITEM .LT. 0) .OR. (ITEM .GT. NPAR)) THEN
            IER = -8
            GOTO 9900
         ENDIF
         IGLBF = .FALSE.
         CALL D1IGFI (CMEM, IMEM, DMEM, IHDR, IGLBF, ITEM, IVAL, IER)

      ELSE IF (SECT .EQ. 'D') THEN

C           ***********************
C           *** Directory entry ***
C           ***********************

         IF (ITYP .NE. ENTGE) THEN
            IER = -7
            GOTO 9900
         ENDIF

         IF ((ITEM .EQ. 16) .OR.
     +       (ITEM .EQ. 17) .OR.
     +       (ITEM .EQ. 18) .OR.
     +       (ITEM .GT. 24) .OR.
     +       (ITEM .LT. 1)) THEN
            IER = -8
            GOTO 9900
         ENDIF

         IF (ITEM .GE. 21) THEN

C           Value is a subfield of Directory Entry 9

            CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 13, IVAL, IERX)
            IF (IERX .NE. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF

            IVAL = MOD (IVAL/100**(24-ITEM), 100)

         ELSE

            IF (ITEM .EQ. 20) THEN
               INDEX = 14
               CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, INDEX, IVAL, IERX)
               IVAL = IVAL+1
            ELSE
               IF (ITEM .LE. 10) THEN
                  INDEX = ITEM+4
               ELSE IF (ITEM .EQ. 11) THEN
C                 ITEM 11 = ITEM 1 = ENTITY TYPE
                  INDEX = 5
               ELSE IF (ITEM .LE. 15) THEN
                  INDEX = ITEM+3
               ELSE IF (ITEM .EQ. 19) THEN
                  INDEX = ITEM
               ENDIF
               CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, INDEX, IVAL, IERX)
            ENDIF

            IF (IERX .NE. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF

         ENDIF

      ELSE IF (SECT .EQ. 'I') THEN

C        *********************
C        *** Index section ***
C        *********************

         IF (ITYP .EQ. ENTGE) THEN
            IER = -99
            CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 1, IGI, IERX)
            IF (IERX .NE. 0) GOTO 9900
            CALL D0PTR (CMEM, IMEM, DMEM, IGI, ENTGI, 0, IHDR, ITYP,
     +                  IERX)
            IF (IERX .NE. 0) GOTO 9900
            IER = 0
         ENDIF

         CALL D1IGLI (CMEM, IMEM, DMEM, IHDR, NSL, NDE, IERX)
         IF (IERX .NE. 0) THEN
            IER = -99
            GOTO 9900
         ENDIF

         INDEX = (ITEM+1)/2

         IF ((MOD(ITEM,2) .NE. 1) .OR. 
     +       (ITEM .GT. NDE) .OR.
     +       (ITEM .LT. 1)) THEN
            IER = -8
            GOTO 9900
         ENDIF

C        Fetch the value

         CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, INDEX+2, IVAL, IER)
         IF (IERX .NE. 0) THEN
            IER = -99
            GOTO 9900
         ENDIF

      ELSE IF (SECT .EQ. 'G') THEN

C        **********************
C        *** Global section ***
C        **********************

         IF (ITYP .EQ. ENTGE) THEN
            IER = -99
            CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 1, IGI, IERX)
            IF (IERX .NE. 0) GOTO 9900
            CALL D0PTR (CMEM, IMEM, DMEM, IGI, ENTGI, 0, IHDR, ITYP,
     +                  IERX)
            IF (IERX .NE. 0) GOTO 9900
            IER = 0
         ENDIF

         IF ((ITEM .GT. 25) .OR. (ITEM .LT. 3)) THEN
            IER = -8
            GOTO 9900
         ENDIF

C        Get the pointer to the Global IGES entity

         CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 1, IGLBL, IERX)
         IF (IERX .NE. 0) THEN
            IER = -99
            GOTO 9900
         ENDIF

         CALL D0PTR (CMEM, IMEM, DMEM, IGLBL, ENTGE, 0, IGHDR,
     +         ITYP, IERX)
         IF (IERX .NE. 0) THEN
            IER = -99
            GOTO 9900
         ENDIF

C        Fetch the element.
         IGLBF = .TRUE.
         CALL D1IGFI (CMEM, IMEM, DMEM, IGHDR, IGLBF, ITEM,
     +      IVAL, IER)

      ELSE

C        *** Invalid section specified ***

         IER = -7
         GOTO 9900
      ENDIF

C     Error reporting section

 9900 CONTINUE

      IF (IER .NE. 0) THEN
          CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      RETURN
      END

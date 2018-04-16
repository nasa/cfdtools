      SUBROUTINE D2IGFC (CMEM, IMEM, DMEM, IGIE, SECT, ITEM,
     +   CSTR, CLEN, IER)

C     Fetch a string from an IGES entity or an IGES Index entity.
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
C     CSTR    String value fetched.
C     CLEN    Length of used portion of CSTR.
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             >0  = String was truncated.  Original length = IER.
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
C                   data structure and store request
C             -9  = TYPCOD for data is non-character.
C
C             -99 = Entity IGIE is internally inconsistent.
C
C
C   CALLS:
C     D0PTR    Pointer check
C     D1IGLI   Fetch length parameters of IGES Index Entity
C     D1FEAC   Fetch a character array
C     D1FEEC   Fetch an character element
C     D1FEEI   Fetch an integer element
C     D1CQFE   Fetch a string from a Character Sequence Entity
C     D1IGLE   Fetch length parameters of IGES Entity
C
C     D1STEC   Store a character element
C     DTERR    Error handler
C
C   HISTORY:
C     7/26/93   D. Parsons   Created
C     8/26/93   D. Parsons   Adjusted PAR values to being with 0 (Entity
C                            Type ID)
C     9/28/93   P. Kraushar  Reorganized to permit access to Index info
C                            from IGES Entity
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_FETCH_CHARS (CMEM, IMEM, DMEM, IGIE, SECT, ITEM,
C     +   CSTR, CLEN, IER)

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

C                          Character Sequence Entity
      INTEGER     ENTCQ
      PARAMETER  (ENTCQ = 249)
C =====================================================================

      CHARACTER*1    TYPCOD, SECT
      CHARACTER*(*)  CSTR
      INTEGER        IGIE, ITEM, CLEN, IER

      INTEGER        NSL, NDE, NPAR, NSTR, NSTRCH, NREAL, IERX
      INTEGER        IHDR, IGHDR, IGLBL, ITYP, JLO, JHI, IGI
      INTEGER        ICQ, ICQHDR, ICQNDX, LENSTR

      CHARACTER*6    SUBNAM
      DATA SUBNAM   /'D2IGFC'/

C****
      IER = 0
      CSTR = ' '

C     Check pointer and type, then locate relevant data spaces
      CALL D0PTR (CMEM, IMEM, DMEM, IGIE, 0, 0, IHDR, ITYP, IER)
      IF (IER .LT. 0) GOTO 9900

      IF ((ITYP .NE. ENTGE) .AND. (ITYP .NE. ENTGI)) THEN
         IER = -6
         GOTO 9900
      ENDIF

      LENSTR = LEN(CSTR)

      IF (SECT .EQ. 'P') THEN

C        *************************
C        *** Parameter section ***
C        *************************

         IF (ITYP .NE. ENTGE) THEN
            IER = -7
            GOTO 9900
         ENDIF

         CALL D1IGLE (CMEM, IMEM, DMEM, IHDR, NPAR, NSTR, NSTRCH,
     +                   NREAL, IERX)
         IF (IERX .NE. 0) THEN
            IER = -99
            GOTO 9900
         ENDIF

         IF ((ITEM .GT. NPAR) .OR. (ITEM .LT. 1)) THEN
            IER = -8
            GOTO 9900
         ENDIF

C        Get the type-code.

         CALL D1FEEC (CMEM, IMEM, DMEM, IHDR, ITEM+24, TYPCOD, IERX)
         IF (IERX .NE. 0) THEN
            IER = -99
            GOTO 9900
         ENDIF

         IF (TYPCOD .EQ. 'O') THEN
            CLEN = 0
            CSTR = ' '
            GOTO 9900
         ENDIF

         IF ((TYPCOD .NE. 'C'))  THEN
            IER = -9
            GOTO 9900
         ENDIF

C        Fetch the index

         CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, ITEM+19, ICQNDX,
     +      IERX)
         IF (IERX .NE. 0) THEN
            IER = -99
            GOTO 9900
         ENDIF

C        Fetch the ICQ

         CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 2, ICQ, IERX)
         IF (IERX .NE. 0) THEN
            IER = -99
            GOTO 9900
         ENDIF

         CALL D0PTR (CMEM, IMEM, DMEM, ICQ, ENTCQ, 0, ICQHDR,
     +         ITYP, IERX)
         IF (IERX .NE. 0) THEN
            IER = -99
            GOTO 9900
         ENDIF

C        Fetch the string

         CALL D1CQFE (CMEM, IMEM, DMEM, ICQHDR, ICQNDX, CSTR,
     +      CLEN, IERX)
         IF (IERX .GT. 0) THEN
            IER = IERX
         ELSE IF (IERX .LT. 0) THEN
            IER = -99
            GOTO 9900
         ENDIF

      ELSE IF (SECT .EQ. 'D') THEN

C        ***********************
C        *** Directory entry ***
C        ***********************

         IF (ITYP .NE. ENTGE) THEN
            IER = -7
            GOTO 9900
         ENDIF

         IF ((ITEM .EQ. 16) .OR.
     +       (ITEM .EQ. 17) .OR.
     +       (ITEM .EQ. 18)) THEN
            JLO = (ITEM-16)*8 + 1
            JHI = JLO+MIN(8,LENSTR)-1
            CALL D1FEAC (CMEM, IMEM, DMEM, IHDR, JLO, JHI, CSTR,
     +         IERX)
            IF (LENSTR .LT. 8) THEN
               IER = 8
               CLEN = LENSTR
               GOTO 9900
            ELSE IF (IERX .LT. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF
         ELSE
            IER = -8
            GOTO 9900
         ENDIF
         CLEN = 8

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

         IF ((ITEM .LT. 1) .OR. (ITEM .GT. 25)) THEN
            IER = -8
            GOTO 9900
         ENDIF

         IF ((ITEM .EQ. 1) .OR. (ITEM .EQ. 2)) THEN

C           Delimiter characters are stored in IGI, not ICQ

            CALL D1FEEC (CMEM, IMEM, DMEM, IHDR, ITEM,
     +         CSTR(1:1), IERX)
            CLEN = 1
            IF (IERX .NE. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF

         ELSE

C           Get the pointer to the Global IGES entity

            CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 1, IGLBL, IERX)
            IF (IERX .NE. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF

            CALL D0PTR (CMEM, IMEM, DMEM, IGLBL, ENTGE, 0, IGHDR,
     +            ITYP, IERX)
            IF (IERX .NE. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF

C           Get the type-code.

            CALL D1FEEC (CMEM, IMEM, DMEM, IGHDR, ITEM-2, TYPCOD,
     +            IERX)
            IF (IERX .NE. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF

            IF (TYPCOD .EQ. 'O') THEN
               CLEN = 0
               CSTR = ' '
               GOTO 9900
            ENDIF

            IF ((TYPCOD .NE. 'C'))  THEN
               IER = -9
               GOTO 9900
            ENDIF

C           Fetch the index

            CALL D1FEEI (CMEM, IMEM, DMEM, IGHDR, ITEM+3,
     +            ICQNDX, IERX)
            IF (IERX .NE. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF

C           Fetch the ICQ

            CALL D1FEEI (CMEM, IMEM, DMEM, IGHDR, 2, ICQ, IERX)
            IF (IERX .NE. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF

            CALL D0PTR (CMEM, IMEM, DMEM, ICQ, ENTCQ, 0, ICQHDR,
     +         ITYP, IERX)
            IF (IERX .NE. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF

C           Fetch the string

            CALL D1CQFE (CMEM, IMEM, DMEM, ICQHDR, ICQNDX, CSTR,
     +         CLEN, IERX)
            IF (IERX .GT. 0) THEN
               IER = IERX
            ELSE IF (IERX .LT. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF


         ENDIF

      ELSE IF (SECT .EQ. 'S') THEN

C        *********************
C        *** Start section ***
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

         IF ((ITEM .LT. 0) .OR. (ITEM .GT. NSL)) THEN
            IER = -8
            GOTO 9900
         ENDIF

         IF (ITEM .EQ. 0) THEN
            CALL D1FEEC (CMEM, IMEM, DMEM, IHDR, 3, CSTR, IERX)
            IF (IERX .LT. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF
            CLEN = 1
         ELSEIF (ITEM .LE. NSL) THEN
            JLO = 4+(ITEM-1)*72
            CLEN = MIN(72,LENSTR)
            JHI = JLO+CLEN-1
            CALL D1FEAC (CMEM, IMEM, DMEM, IHDR, JLO, JHI, CSTR,
     +         IERX)
            IF (LENSTR .LT. 72) THEN
               IER = 72
               GOTO 9900
            ELSE IF (IERX .LT. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF
         ENDIF

      ELSE

C        *** Invalid section specified ***

         IER = -7
         GOTO 9900

      ENDIF


C     Error reporting section

 9900 CONTINUE

      IF (IER .NE. 0) THEN
          CALL DTERR (1, SUBNAM, IER, 0)
          IF (IER .LT. 0) THEN
            CLEN = 0
            CSTR = ' '
          ENDIF
      ENDIF

      RETURN
      END

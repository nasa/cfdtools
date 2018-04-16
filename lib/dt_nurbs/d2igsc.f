      SUBROUTINE D2IGSC (CMEM, IMEM, DMEM, TYPCOD, CSTR, IGIE, SECT,
     +      ITEM, IER)

C     Store a string into an IGES entity or an IGES Index entity.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     TYPCOD  IGES-context data type code.  Only 'C' and 'O' are
C             valid here.
C     CSTR    String value to store.
C     IGIE    Pointer to an IGES entity or IGES Index entity.
C     SECT    IGES "section" to which the request is directed.
C     ITEM    Index into the section to which the request is directed.
C   OUTPUT:
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             >0  = String was truncated to length = IER to fit.
C             -1  = Dynamic memory is corrupt or uninitialized
C             -2  = IGIE is a null pointer
C             -3  = Garbage value found in IGIE
C             -4  = Ambiguous--garbage pointer or deleted entity
C             -5  = IGIE points to a deleted entity
C             -6  = IGIE does not point to an IGES entity or an
C                   IGES Index entity.
C             -7  = TYPCOD not = 'O' or 'C'.
C             -8  = SECT is invalid for this type of entity (IGES
C                   or Index).
C             -9  = ITEM is out of range with respect to this IGES
C                   data structure and store request.
C             -10 = Attempt to change TYPCOD is invalid.
C             -11 = TYPCOD = 'O' only permitted in SECT='P'
C             -12 = Insufficient storage for this item.
C
C             -99 = Entity IGIE is internally inconsistent.
C
C
C   CALLS:
C     D0PTR    Pointer check
C     D1IGLI   Fetch length parameters of IGES Index Entity
C     D1IGLE   Fetch length parameters of IGES Entity
C     D1STAC   Store a character array
C     D1STEC   Store a character element
C     D1FEEI   Fetch an integer element
C     D1FEEC   Fetch an character element
C     D1CQST   Store a string to a Character Sequence Entity
C     DTERR    Error handler
C
C   HISTORY:
C     7/26/93   D. Parsons   Created
C     8/6/93    D. Parsons   Add NPAR expansion capability.
C     8/26/93   D. Parsons   Adjusted PAR values to being with 0 (Entity
C                            Type ID)
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_STORE_CHARS (CMEM, IMEM, DMEM, TYPCOD, CSTR,
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

C                          Character Sequence Entity
      INTEGER     ENTCQ
      PARAMETER  (ENTCQ = 249)
C =====================================================================

      CHARACTER*1    TYPCOD, SECT
      CHARACTER*(*)  CSTR
      INTEGER        IGIE, ITEM, IER

      INTEGER        NSL, NDE, NPAR, NSTR, NSTRCH, NREAL, IERX
      INTEGER        MPAR, MSTR, MSTRCH, MREAL, NEED
      INTEGER        IHDR, IGHDR, IGLBL, ITYP, JLO, JHI
      INTEGER        ICQ, ICQHDR, ICQNDX, LENSTR, IVAL
      CHARACTER*(1)  OTCOD


      CHARACTER*6    SUBNAM
      DATA SUBNAM   /'D2IGSC'/

C****
      IER = 0


C     Check pointer and type, then locate relevant data spaces
      CALL D0PTR (CMEM, IMEM, DMEM, IGIE, 0, 0, IHDR, ITYP, IER)
      IF (IER .LT. 0) GOTO 9900

      IF ((ITYP .NE. ENTGE) .AND. (ITYP .NE. ENTGI)) THEN
         IER = -6
         GOTO 9900
      ENDIF

      IF ((TYPCOD .NE. 'O') .AND. (TYPCOD .NE. 'C')) THEN
         IER = -7
         GOTO 9900
      ENDIF

      IF ((TYPCOD .EQ. 'O') .AND. (SECT .NE. 'P')) THEN
         IER = -11
         GOTO 9900
      ENDIF

      LENSTR = LEN(CSTR)

      IF (ITYP .EQ. ENTGI) THEN

C        *************************
C        *** IGES INDEX ENTITY ***
C        *************************

         CALL D1IGLI (CMEM, IMEM, DMEM, IHDR, NSL, NDE, IERX)
         IF (IERX .NE. 0) THEN
            IER = -99
            GOTO 9900
         ENDIF

         IF (SECT .EQ. 'S') THEN

C           *** Start section ***

            IF (ITEM .EQ. 0) THEN
               CALL D1STEC (CMEM, IMEM, DMEM, CSTR(1:1), 3, IHDR, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF
            ELSE IF (ITEM .LE. NSL) THEN
               JLO = 4+(ITEM-1)*72
               JHI = JLO+MIN(72,LENSTR)-1
               CALL D1STAC (CMEM, IMEM, DMEM, CSTR, JLO, JHI, IHDR,
     +            IERX)
               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF

               IF (LENSTR .GT. 72) THEN
                  IER = 72
               ENDIF
            ELSE
               IER = -9
               GOTO 9900
            ENDIF

         ELSE IF (SECT .EQ. 'G') THEN

C           *** Global section ***

            IF ((ITEM .LT. 1) .OR. (ITEM .GT. 25)) THEN
               IER = -9
               GOTO 9900
            ENDIF

            IF ((ITEM .EQ. 1) .OR. (ITEM .EQ. 2)) THEN

C              Delimiter characters are stored in IGI, not ICQ

               CALL D1STEC (CMEM, IMEM, DMEM, CSTR(1:1), ITEM, IHDR,
     +               IERX)
               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF

            ELSE

C              Get the pointer to the Global IGES entity

               CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 1, IGLBL, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF

               CALL D0PTR (CMEM, IMEM, DMEM, IGLBL, ENTGE, 0, IGHDR,
     +               ITYP, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF

C              Get the old type-code.

               CALL D1FEEC (CMEM, IMEM, DMEM, IGHDR, ITEM-2, OTCOD,
     +               IERX)
               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF

               IF (TYPCOD .NE. OTCOD) THEN
                  IER = -10
                  GOTO 9900
               ENDIF

C              Fetch the index

               CALL D1FEEI (CMEM, IMEM, DMEM, IGHDR, ITEM+3, ICQNDX,
     +               IERX)
               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF

C              Fetch the ICQ

               CALL D1FEEI (CMEM, IMEM, DMEM, IGHDR, 2, ICQ, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF


               CALL D0PTR (CMEM, IMEM, DMEM, ICQ, ENTCQ, 0, ICQHDR,
     +               ITYP, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF

C              Update the ICQ entry

               CALL D1CQST (CMEM, IMEM, DMEM, ICQHDR, ICQNDX,
     +               LENSTR, CSTR, NEED, IERX)
               IF (IERX .NE. 0) THEN
                  IF ((IERX .EQ. -13) .OR. (IERX .EQ. -14)) THEN
                     IER = -12
                  ELSE
                     IER = -99
                  ENDIF
                  GOTO 9900
               ENDIF

            ENDIF


         ELSE

C           *** Invalid section specified ***

            IER = -8
            GOTO 9900
         ENDIF
      ELSE

C        *******************
C        *** IGES ENTITY ***
C        *******************

         CALL D1IGLE (CMEM, IMEM, DMEM, IHDR, NPAR, NSTR, NSTRCH,
     +                   NREAL, IERX)
         IF (IERX .NE. 0) THEN
            IER = -99
            GOTO 9900
         ENDIF


         IF (SECT .EQ. 'D') THEN

C           *** Directory entry ***

            IF ((ITEM .EQ. 16) .OR.
     +          (ITEM .EQ. 17) .OR.
     +          (ITEM .EQ. 18)) THEN
               JLO = (ITEM-16)*8 + 1
               JHI = JLO+7
               CALL D1STAC (CMEM, IMEM, DMEM, CSTR, JLO, JHI,
     +               IHDR, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF

               IF (LENSTR .GT. 8) THEN
                  IER = 8
               ENDIF

            ELSE
               IER = -9
               GOTO 9900
            ENDIF
         ELSE IF (SECT .EQ. 'P') THEN

C           *** Parameter section ***

            IF (ITEM .LT. 1) THEN
               IER = -9
               GOTO 9900
            ENDIF

            IF (ITEM .GT. NPAR) THEN

C              Try to expand entity

               MPAR   = ITEM-NPAR
               MSTR   = 0
               MSTRCH = 0
               MREAL  = 0
               CALL D1IGXE (CMEM, IMEM, DMEM, IHDR, MPAR, MSTR, MSTRCH,
     +                   MREAL, NEED, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -12
                  GOTO 9900
               ENDIF
            ENDIF

C           Get the old type-code.
            CALL D1FEEC (CMEM, IMEM, DMEM, IHDR, ITEM+24, OTCOD, IERX)
            IF (IERX .NE. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF

C           Fetch the ICQ and it's Header pointer
            CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 2, ICQ, IERX)
            IF (IERX .NE. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF

            CALL D0PTR (CMEM, IMEM, DMEM, ICQ, ENTCQ, 0, ICQHDR,
     +            ITYP, IERX)
            IF (IERX .NE. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF

C           Get the old index or IVAL
            CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, ITEM+19, IVAL,
     +          IERX)
            IF (IERX .NE. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF

C           Get the index into the ICQ space
            IF (OTCOD .NE. 'C') THEN

               IF ((OTCOD .EQ. 'D') .OR.
     +             (OTCOD .EQ. 'E')) THEN
C                 This was a DMEM pointer--free it, ignore any errors
                  CALL D0IGFD (CMEM, IMEM, DMEM, IHDR, IVAL, IERX)
               ELSE IF ((OTCOD .EQ. 'T') .OR.
     +                  (OTCOD .EQ. 'R')) THEN
C                 This was an ICQ pointer--free it, ignore any errors
                  CALL DTERPT(0)
                  CALL D2ERAS (CMEM, IMEM, DMEM, IVAL, IERX)
                  CALL DTERPT(1)
               ENDIF

C              Find the next free position in the ICQ
               CALL D0CQNF (CMEM, IMEM, DMEM, ICQHDR, ICQNDX, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF
            ELSE
               ICQNDX = IVAL
            ENDIF

C           Update the ICQ entry
            IF (TYPCOD .EQ. 'C') THEN
               CALL D1CQST (CMEM, IMEM, DMEM, ICQHDR, ICQNDX,
     +               LENSTR, CSTR, NEED, IERX)
               IF (IERX .NE. 0) THEN
                  IF ((IERX .EQ. -13) .OR. (IERX .EQ. -14)) THEN
                     IER = -11
                  ELSE
                     IER = -99
                  ENDIF
                  GOTO 9900
               ENDIF
C              Update the ICQ index
               CALL D1STEI (CMEM, IMEM, DMEM, ICQNDX, ITEM+19, IHDR,
     +            IERX)
               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF
            ENDIF

C           Update the type-code
            IF (TYPCOD .NE. OTCOD) THEN
               CALL D1STEC (CMEM, IMEM, DMEM, TYPCOD, ITEM+24, IHDR,
     +               IERX)
               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF
            ENDIF

         ELSE

C           *** Invalid section specified ***

            IER = -8
            GOTO 9900
         ENDIF
      ENDIF


C     Error reporting section

 9900 CONTINUE

      IF (IER .EQ. -12) THEN
          CALL DTERR (2, SUBNAM, IER, NEED)
      ELSE IF (IER .NE. 0) THEN
          CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      RETURN
      END

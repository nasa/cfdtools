      SUBROUTINE D2IGSI (CMEM, IMEM, DMEM, TYPCOD, IVAL, IGIE, SECT,
     +      ITEM, IER)

C     Store an integer element into an IGES entity or an IGES Index
C     entity.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     TYPCOD  IGES-context data type code.  Only 'I', 'P', 'L',
C             'T', 'R', '*', and 'O' are valid here.
C     IVAL    Integer element to store.
C     IGIE    Pointer to an IGES entity or IGES Index entity.
C     SECT    IGES "section" to which the request is directed.
C     ITEM    Index into the section to which the request is directed.
C   OUTPUT:
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             -1  = Dynamic memory is corrupt or uninitialized
C             -2  = IGIE is a null pointer
C             -3  = Garbage value found in IGIE
C             -4  = Ambiguous--garbage pointer or deleted entity
C             -5  = IGIE points to a deleted entity
C             -6  = IGIE does not point to an IGES entity or an
C                   IGES Index entity.
C             -7  = TYPCOD not = 'I', 'P', 'L', 'T', 'R', '*', or 'O'.
C             -8  = SECT is invalid for this type of entity (IGES
C                   or Index).
C             -9  = ITEM is out of range with respect to this IGES
C                   data structure and store request.
C             -10 = Attempt to change Global TYPCOD is invalid.
C             -12 = SECT 'I' specified, TYPCOD not '*'.
C             -13 = TYPCOD = '*', but IVAL is not a valid MEM pointer to
C                   an IGES Entity.
C             -14 = TYPCOD = 'L', but IVAL is not zero or one
C             -15 = TYPCOD = 'T' or 'R', but IVAL is not a valid
C                   Character Seqence entity.
C             -16 = TYPCOD = '*', but SECT does not equal 'I'.
C             -17 = TYPCOD = 'R', 'T', or 'O' only applies to parameter
C                   data, but SECT is not 'P'.
C             -18 = Insufficient storage to perform needed expansion.
C             -99 = Entity IGIE is internally inconsistent.
C
C
C   CALLS:
C     D0PTR    Pointer check
C     D1IGLI   Fetch length parameters of IGES Index Entity
C     D1FEEI   Fetch an integer element
C     D1FEEC   Fetch a character element
C     D1IGSI   Store an integer to parameter storage
C     D1STEI   Store an integer element
C     D1IGLE   Fetch length parameters of IGES Entity
C     DTERR    Error handler
C
C   HISTORY:
C     7/26/93  D. Parsons   Created.
C     8/6/93   D. Parsons   Add parameter-section expansion capability.
C     8/26/93  D. Parsons   Adjusted PAR values to being with 0 (Entity
C                           Type ID)
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_STORE_INTEGER (CMEM, IMEM, DMEM, TYPCOD, IVAL,
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
      INTEGER        IGIE, ITEM, IVAL, IER

      INTEGER        NSL, NDE, NPAR, NSTR, NSTRCH, NREAL, IERX
      INTEGER        MPAR, MSTR, MSTRCH, MREAL, NEED
      INTEGER        IHDR, IGHDR, IGLBL, ITYP
      CHARACTER*(1)  OTCOD
      INTEGER        IGEHDR, ICQHDR
      LOGICAL        IGLBF
      INTEGER        IVAL9, INDEX, IVALB


      CHARACTER*6    SUBNAM
      DATA SUBNAM   /'D2IGSI'/

C****
      IER = 0


C     Check pointer and type, then locate relevant data spaces
      CALL D0PTR (CMEM, IMEM, DMEM, IGIE, 0, 0, IHDR, ITYP, IER)
      IF (IER .LT. 0) GOTO 9900

      IF ((ITYP .NE. ENTGE) .AND. (ITYP .NE. ENTGI)) THEN
         IER = -6
         GOTO 9900
      ENDIF

C     Verify that this value matches the TYPCOD
      IF (TYPCOD .EQ. '*') THEN
         IF (IVAL .NE. 0) THEN
            CALL D0PTR (CMEM, IMEM, DMEM, IVAL, ENTGE, 0, IGEHDR, ITYP,
     +                  IERX)
            IF (IERX .NE. 0) THEN
               IER = -13
               GOTO 9900
            ENDIF
         ELSE
            IGEHDR = 0
         ENDIF
         IF (SECT .NE. 'I') THEN
            IER = -16
            GOTO 9900
         ENDIF
      ELSE IF (TYPCOD .EQ. 'L') THEN
         IF ((IVAL .NE. 0) .AND. (IVAL .NE. 1)) THEN
            IER = -14
            GOTO 9900
         ENDIF
      ELSE IF ((TYPCOD .EQ. 'T') .OR. (TYPCOD .EQ. 'R')) THEN
         CALL D0PTR (CMEM, IMEM, DMEM, IVAL, ENTCQ, 0, ICQHDR, ITYP,
     +      IERX)
         IF (IERX .NE. 0) THEN
            IER = -15
            GOTO 9900
         ENDIF
         IF (SECT .NE. 'P') THEN
            IER = -17
            GOTO 9900
         ENDIF
      ELSE IF (TYPCOD .EQ. 'O') THEN
         IF (SECT .NE. 'P') THEN
            IER = -17
            GOTO 9900
         ENDIF
      ELSE IF ((TYPCOD .NE. 'I') .AND. (TYPCOD .NE. 'P')) THEN
         IER = -7
         GOTO 9900
      ENDIF

      IF (ITYP .EQ. ENTGI) THEN

C        *************************
C        *** IGES INDEX ENTITY ***
C        *************************

         CALL D1IGLI (CMEM, IMEM, DMEM, IHDR, NSL, NDE, IERX)
         IF (IERX .NE. 0) THEN
            IER = -99
            GOTO 9900
         ENDIF

         IF (SECT .EQ. 'G') THEN

C           *** Global section ***

            IF ((ITEM .GT. 25) .OR. (ITEM .LT. 3)) THEN
               IER = -9
               GOTO 9900
            ENDIF

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

C           Get the old type-code.

            CALL D1FEEC (CMEM, IMEM, DMEM, IGHDR, ITEM-2, OTCOD, IERX)
            IF (IERX .NE. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF

            IF (TYPCOD .NE. OTCOD) THEN
               IER = -10
               GOTO 9900
            ENDIF

C           Update the entry

            IGLBF = .TRUE.
            CALL D1IGSI (CMEM, IMEM, DMEM, IGHDR, IGLBF, TYPCOD, ITEM,
     +         IVAL, IER)

         ELSE IF (SECT .EQ. 'I') THEN

C           *** Index to IGES Entity ***

            INDEX = (ITEM+1)/2

            IF ((MOD(ITEM,2) .NE. 1) .OR. 
     +          (ITEM .GT. NDE) .OR.
     +          (ITEM .LT. 1)) THEN
               IER = -9
               GOTO 9900
            ENDIF

            IF (TYPCOD .NE. '*') THEN
               IER = -12
               GOTO 9900
            ENDIF

C           Store the value

            CALL D1STEI (CMEM, IMEM, DMEM, IVAL, INDEX+2, IHDR, IERX)
            IF (IERX .NE. 0) THEN
               IER = -99
               GOTO 9900
            ENDIF

C           Point this entity back to the IGES Index Entity

            IF (IGEHDR .NE. 0) THEN
               CALL D1STEI (CMEM, IMEM, DMEM, IGIE, 1, IGEHDR, IERX)
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
     +          (ITEM .EQ. 18) .OR.
     +          (ITEM .GT. 24) .OR.
     +          (ITEM .LT. 1)) THEN
               IER = -9
               GOTO 9900
            ENDIF

            IF (ITEM .GE. 21) THEN

C              Value is a subfield of Directory Entry 9

               CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 13, IVAL9, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF

               IVALB = MOD (IVAL9/100**(24-ITEM), 100)
               IVAL9 = IVAL9 + (IVAL-IVALB)*100**(24-ITEM)

               CALL D1STEI (CMEM, IMEM, DMEM, IVAL9, 13, IHDR, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF

            ELSE

               IF (ITEM .EQ. 20) THEN
                  INDEX = 14
                  CALL D1STEI (CMEM, IMEM, DMEM, IVAL-1, INDEX, IHDR, 
     +                         IERX)
               ELSE
                  IF (ITEM .LE. 10) THEN
                     INDEX = ITEM+4
                  ELSE IF (ITEM .EQ. 11) THEN
C                    ITEM 11 = ITEM 1 = ENTITY TYPE
                     INDEX = 5
                  ELSE IF (ITEM .LE. 15) THEN
                     INDEX = ITEM+3
                  ELSE IF (ITEM .EQ. 19) THEN
                     INDEX = ITEM
                  ENDIF
                  CALL D1STEI (CMEM, IMEM, DMEM, IVAL, INDEX, IHDR, 
     +                         IERX)
               ENDIF

               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF

            ENDIF


         ELSE IF (SECT .EQ. 'P') THEN

C           *** Parameter section ***

            IF (ITEM .LT. 0) THEN
               IER = -9
               GOTO 9900
            ENDIF
            
            IF (ITEM .GT. NPAR) THEN

C              Try to expand entity

               MPAR   = ITEM-NPAR
               MSTR   = 0
               MSTRCH = 0
               MREAL  = 0
               CALL D1IGXE (CMEM, IMEM, DMEM, IHDR, MPAR, MSTR,
     +                   MSTRCH, MREAL, NEED, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -18
                  NEED = MPAR
                  GOTO 9900
               ENDIF
            ENDIF

            IGLBF = .FALSE.
            CALL D1IGSI (CMEM, IMEM, DMEM, IHDR, IGLBF, TYPCOD, ITEM,
     +         IVAL, IER)

         ELSE

C           *** Invalid section specified ***

            IER = -8
            GOTO 9900
         ENDIF
      ENDIF


C     Error reporting section

 9900 CONTINUE

      IF (IER .EQ. -18) THEN
          CALL DTERR (2, SUBNAM, IER, NEED)
      ELSE IF (IER .NE. 0) THEN
          CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      RETURN
      END

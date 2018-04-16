      SUBROUTINE D2IGSD (CMEM, IMEM, DMEM, TYPCOD, DVAL, IGIE, SECT,
     +      ITEM, IER)

C     Store a double precision element into an IGES entity or an
C     IGES Index entity.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     TYPCOD  IGES-context data type code.  Only 'E', 'D' and 'O' are
C             valid here.
C     DVAL    Double Precision element to store.
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
C             -7  = TYPCOD not = 'E', 'D', or 'O'.
C             -8  = SECT is invalid for this type of entity (IGES
C                   or Index).
C             -9  = ITEM is out of range with respect to this IGES
C                   data structure and store request.
C             -10 = Attempt to change Global TYPCOD is invalid.
C             -12 = Insufficient storage for this item.
C
C             -99 = Entity IGIE is internally inconsistent.
C
C
C   CALLS:
C     D0PTR    Pointer check
C     D1IGLI   Fetch length parameters of IGES Index Entity
C     D1FEEI   Fetch an integer element
C     D1FEEC   Fetch a character element
C     D1IGSD   Store a double precision element to parameter storage
C     D1IGLE   Fetch length parameters of IGES Entity
C     DTERR    Error handler
C
C   HISTORY:
C     7/26/93  D. Parsons   Created
C     8/6/93   D. Parsons   Add parameter-section expansion capability.
C     8/26/93  D. Parsons   Adjusted PAR values to being with 0 (Entity
C                           Type ID)
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_STORE_DOUBLE (CMEM, IMEM, DMEM, TYPCOD, DVAL,
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

      CHARACTER*1    TYPCOD, SECT
      INTEGER        IGIE, ITEM, IER
      DOUBLE PRECISION DVAL

      INTEGER        NSL, NDE, NPAR, NSTR, NSTRCH, NREAL, IERX
      INTEGER        MPAR, MSTR, MSTRCH, MREAL, NEED
      INTEGER        IHDR, IGHDR, IGLBL, ITYP
      LOGICAL        IGLBF
      CHARACTER*(1)  OTCOD

      CHARACTER*6    SUBNAM
      DATA SUBNAM   /'D2IGSD'/

C****
      IER = 0

C     Check pointer and type, then locate relevant data spaces
      CALL D0PTR (CMEM, IMEM, DMEM, IGIE, 0, 0, IHDR, ITYP, IER)
      IF (IER .LT. 0) GOTO 9900

      IF ((ITYP .NE. ENTGE) .AND. (ITYP .NE. ENTGI)) THEN
         IER = -6
         GOTO 9900
      ENDIF

      IF ((TYPCOD .NE. 'E') .AND.
     +    (TYPCOD .NE. 'D') .AND.
     +    (TYPCOD .NE. 'O'))     THEN
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
C           Do not allow changing of Global TYPCODs
               IER = -10
               GOTO 9900
            ENDIF

C           Update the entry

            IGLBF = .TRUE.
            CALL 	D1IGSD (CMEM, IMEM, DMEM, IGHDR, IGLBF, TYPCOD, ITEM,
     +         DVAL, NEED, IER)
            IF (IER .NE. 0) GOTO 9900
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


         IF (SECT .EQ. 'P') THEN

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

            IGLBF = .FALSE.
            CALL 	D1IGSD (CMEM, IMEM, DMEM, IHDR, IGLBF, TYPCOD, ITEM,
     +         DVAL, NEED, IER)
            IF (IER .NE. 0) THEN
               NEED = 1
               GOTO 9900
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

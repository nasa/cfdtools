      SUBROUTINE D2IGFT (CMEM, IMEM, DMEM, IGIE, SECT, ITEM,
     +   TYPCOD, IER)

C    Fetch parameter TYPCOD from an IGES entity or IGES Index entity.
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
C     TYPCOD  IGES-context data type code.
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             -1  = Dynamic memory is corrupt or uninitialized
C             -2  = IGIE is a null pointer
C             -3  = Garbage value found in IGIE
C             -4  = Ambiguous--garbage pointer or deleted entity
C             -5  = IGIE points to a deleted entity
C             -6  = IGIE does not point to an IGES entity or an
C                   IGES Index entity.
C             -7  = SECT is invalid for this type of entity;
C                   must be 'G' if IGIE points to IGES Index entity,
C                   or 'D' or 'P' if IGIE points to an IGES entity.
C             -8  = ITEM is out of range with respect to this IGES
C                   data structure.
C             -99 = Entity IGIE is internally inconsistent.
C
C
C   CALLS:
C     D0PTR    Pointer check
C     D1IGLI   Fetch length parameters of IGES Index Entity
C     D1FEEI   Fetch an integer element
C     D1FEEC   Fetch a character element
C     D1IGLE   Fetch length parameters of IGES Entity
C     DTERR    Error handler
C
C   HISTORY:
C     7/26/93  D. Parsons   Created
C     8/26/93  D. Parsons   Adjusted PAR values to being with 0 (Entity
C                           Type ID)
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_FETCH_TYPE (CMEM, IMEM, DMEM, IGIE, SECT, ITEM,
C     +   TYPCOD, IER)

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

      INTEGER        NSL, NDE, NPAR, NSTR, NSTRCH, NREAL, IERX
      INTEGER        IHDR, IGHDR, IGLBL, ITYP

      CHARACTER*6    SUBNAM
      DATA SUBNAM   /'D2IGFT'/

C****
      IER    = 0
      TYPCOD = ' '


C     Check pointer and type, then locate relevant data spaces
      CALL D0PTR (CMEM, IMEM, DMEM, IGIE, 0, 0, IHDR, ITYP, IER)
      IF (IER .LT. 0) GOTO 9900

      IF ((ITYP .NE. ENTGE) .AND. (ITYP .NE. ENTGI)) THEN
         IER = -6
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

            IF ((ITEM .GT. 25) .OR. (ITEM .LT. 1)) THEN
               IER = -8
               GOTO 9900
            ENDIF    

C           Get the pointer to the Global IGES entity

            IF ((ITEM .EQ. 1) .OR. (ITEM .EQ. 2)) THEN
               TYPCOD = 'C'
            ELSE
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

C              Get the type-code.

               CALL D1FEEC (CMEM, IMEM, DMEM, IGHDR, ITEM-2, TYPCOD, 
     +                      IERX)
               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF
            ENDIF

         ELSE

C           *** Invalid section specified ***

            IER = -7
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

            IF ((ITEM .LT. 1) .OR. (ITEM .GT. 24)) THEN
               IER = -8
               GOTO 9900
            ELSE IF ((ITEM .EQ. 16) .OR.
     +               (ITEM .EQ. 17) .OR.
     +               (ITEM .EQ. 18) .OR.
     +               (ITEM .GE. 21)) THEN
               TYPCOD = 'C'
            ELSE
               TYPCOD = 'I'
            ENDIF

         ELSE IF (SECT .EQ. 'P') THEN

C           *** Parameter section ***

            IF ((ITEM .GT. NPAR) .OR. (ITEM .LT. 0)) THEN
               IER = -8
               GOTO 9900
            ENDIF

            IF (ITEM .EQ. 0) THEN
               TYPCOD = 'I'
            ELSE

C              Get the old type-code.

               CALL D1FEEC (CMEM, IMEM, DMEM, IHDR, ITEM+24, TYPCOD,
     +            IERX)
               IF (IERX .NE. 0) THEN
                  IER = -99
                  GOTO 9900
               ENDIF
            ENDIF

         ELSE

C           *** Invalid section specified ***

            IER = -7
            GOTO 9900
         ENDIF
      ENDIF


C     Error reporting section

 9900 CONTINUE

      IF (IER .NE. 0) THEN
          CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      RETURN
      END

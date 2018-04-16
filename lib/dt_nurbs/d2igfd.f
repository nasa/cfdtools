      SUBROUTINE D2IGFD (CMEM, IMEM, DMEM, IGIE, SECT, ITEM, DVAL, IER)

C     Fetch a double precision element from an IGES entity or an
C     IGES Index entity.
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
C     DVAL    Fetched double precision element.
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
C             -9 =  TYPCOD associated with this data is not 'D', 'E',
C                   or 'O'.
C
C             -99 = Entity IGIE is internally inconsistent.
C
C
C   CALLS:
C     D0PTR    Pointer check
C     D1IGLI   Fetch length parameters of IGES Index Entity
C     D1FEEI   Fetch an integer element
C     D1IGDI   Fetch a double precision element from parameter storage
C     D1FEED   Fetch a double precision element
C     D1IGLE   Fetch length parameters of IGES Entity
C     DTERR    Error handler
C
C   HISTORY:
C     7/26/93 D. Parsons   Created
C     8/26/93 D. Parsons   Adjusted PAR values to being with 0 (Entity
C                          Type ID)
C     9/28/93 P. Kraushar  Reorganized to permit access to Index info
C                          from IGES Entity
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_FETCH_DOUBLE (CMEM, IMEM, DMEM, IGIE, SECT,
C    +      ITEM, DVAL, IER)

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
      INTEGER        IGIE, ITEM, IER
      DOUBLE PRECISION DVAL

      INTEGER        NPAR, NSTR, NSTRCH, NREAL, IERX
      INTEGER        IHDR, IGHDR, IGLBL, ITYP, IGI
      LOGICAL        IGLBF

      CHARACTER*1    TYPCOD

      CHARACTER*6    SUBNAM
      DATA SUBNAM   /'D2IGFD'/

C****
      IER  = 0
      DVAL = 0.0D0

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
     +                   NREAL, IERX)
         IF (IERX .NE. 0) THEN
            IER = -99
            GOTO 9900
         ENDIF

         IF ((ITEM .GT. NPAR) .OR. (ITEM .LT. 1)) THEN
            IER = -8
            GOTO 9900
         ENDIF

         IGLBF = .FALSE.
         CALL D1IGFD (CMEM, IMEM, DMEM, IHDR, IGLBF, ITEM,
     +      DVAL, IER)

      ELSE IF (SECT .EQ. 'G') THEN

C           **********************
C           *** Global section ***
C           **********************

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
     +      ITYP, IERX)
         IF (IERX .NE. 0) THEN
            IER = -99
            GOTO 9900
         ENDIF

C        Get the type-code.
         CALL D1FEEC (CMEM, IMEM, DMEM, IGHDR, ITEM-2, TYPCOD,
     +         IERX)
         IF (IERX .NE. 0) THEN
            IER = -99
            GOTO 9900
         ENDIF

         IF (TYPCOD .EQ. 'O') THEN
            DVAL = 0.0D0
            GOTO 9900
         ENDIF

C        Fetch the element.
         IGLBF = .TRUE.
         CALL D1IGFD (CMEM, IMEM, DMEM, IGHDR, IGLBF, ITEM,
     +      DVAL, IER)

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

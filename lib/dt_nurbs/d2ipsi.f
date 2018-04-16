      SUBROUTINE D2IPSI (CMEM, IMEM, DMEM, NUM, IA, INCIA, IGE, TYPCOD,
     +      JPAR, INCPAR, IER)

C     Store an integer array into an IGES entity parameter space.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     NUM     Number of integers to store
C     IA      Array of integers to store
C     INCIA   Increment to IA, elements (stride).  INCIA = 0
C             results in a "fill" of the same value.
C     IGE     Pointer to an IGES Entity
C     TYPCOD  IGES-context data type code.  Only 'I', 'L',
C             'P', 'A', 'T', and 'O' are valid here.
C     JPAR    Parameter number of 1st item location.
C     INCPAR  Increment to Parameter numbers (stride)
C   OUTPUT:
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             -1  = Dynamic memory is corrupt or uninitialized
C             -2  = IGE is a null pointer
C             -3  = Garbage value found in IGE
C             -4  = Ambiguous--garbage pointer or deleted entity
C             -5  = IGE points to a deleted entity
C             -6  = IGE does not point to an IGES entity.
C             -7  = NUM < 0.
C             -8  = INCIA < 0.
C             -9  = TYPCOD not = 'I', 'L', 'P', 'A', 'T', or 'O'.
C             -10 = JPAR is out of range.
C             -11 = Insufficient storage to perform needed expansion.
C             -12 = INCPAR < 0.
C             -13 = TYPCOD = 'T' but one or more values is not a MEM pointer
C                   to a Character Sequence Entity.
C             -14 = TYPCOD = 'L' but one or more values is not zero or one.
C             -99 = Entity IGE is internally inconsistent.
C
C
C   CALLS:
C     D0PTR    Pointer check
C     D1IGLE   Fetch length parameters of IGES Entity
C     D1IGSI   Store an integer element to an IGES Entity
C     DTERR    Error handler
C
C   HISTORY:
C     7/26/93 D. Parsons   Created
C     8/6/93  D. Parsons   Add parameter-section expansion capability.
C     8/26/93 D. Parsons   Adjusted PAR values to being with 0 (Entity
C                          Type ID)
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_PD_STORE_INTS (CMEM, IMEM, DMEM, NUM, IA, INCIA,
C    +   IGE, TYPCOD, JPAR, INCPAR, IER)

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

C                          Character String Sequence
      INTEGER     ENTCQ
      PARAMETER  (ENTCQ = 249)

C =====================================================================

      CHARACTER*1    TYPCOD
      INTEGER        NUM, IA(*), INCIA, IGE, JPAR, INCPAR, IER

      INTEGER        NPAR, NSTR, NSTRCH, NREAL, IERX
      INTEGER        MPAR, MSTR, MSTRCH, MREAL
      INTEGER        IHDR, ITYP, I, IPAR, IVAL, NEED
      LOGICAL        IGLBL

      CHARACTER*6    SUBNAM
      DATA SUBNAM   /'D2IPSI'/

C****
      IER = 0

C     Check pointer and type, then locate relevant data spaces
      CALL D0PTR (CMEM, IMEM, DMEM, IGE, ENTGE, 0, IHDR, ITYP, IER)
      IF (IER .LT. 0) GOTO 9900

      IF (NUM .LT. 0) THEN
         IER = -7
         GOTO 9900
      ENDIF

      IF (INCIA .LT. 0) THEN
         IER = -8
         GOTO 9900
      ENDIF

      IF ((TYPCOD .NE. 'O') .AND.
     +    (TYPCOD .NE. 'I') .AND.
     +    (TYPCOD .NE. 'P') .AND.
     +    (TYPCOD .NE. 'L') .AND.
     +    (TYPCOD .NE. 'A'))     THEN
         IER = -9
         GOTO 9900
      ENDIF

      IF (INCPAR .LT. 0) THEN
         IER = -12
         GOTO 9900
      ENDIF

      IF (TYPCOD .EQ. 'L') THEN
         DO 80 I=1,1+(NUM-1)*INCIA,INCIA
            IF (IA(I) .LT. 0 .OR. IA(I) .GT. 1) THEN
               IER = -14
               GOTO 9900
            ENDIF
   80    CONTINUE
      ELSE IF (TYPCOD .EQ. 'T') THEN
         DO 90 I=1,1+(NUM-1)*INCIA,INCIA
            IF (IA(I) .NE. 0) THEN
               CALL D0PTR (CMEM, IMEM, DMEM, IA(I), ENTCQ, 0,
     +                     IVAL, ITYP, IERX)
               IF (IERX .NE. 0) THEN
                  IER = -13
                  GOTO 9900
               ENDIF
            ENDIF
   90    CONTINUE
      ENDIF

C     Fetch the defined lengths

      CALL D1IGLE (CMEM, IMEM, DMEM, IHDR, NPAR, NSTR, NSTRCH,
     +                NREAL, IERX)
      IF (IERX .NE. 0) THEN
         IER = -99
         GOTO 9900
      ENDIF

      IF (JPAR .LT. 0) THEN
         IER = -10
         GOTO 9900
      ENDIF

      NEED = JPAR + (NUM-1)*INCPAR

      IF (NEED .GT. NPAR) THEN

C        Try to expand entity

         MPAR   = NEED-NPAR
         MSTR   = 0
         MSTRCH = 0
         MREAL  = 0
         CALL D1IGXE (CMEM, IMEM, DMEM, IHDR, MPAR, MSTR, MSTRCH,
     +                MREAL, NEED, IERX)
         IF (IERX .NE. 0) THEN
            IER = -11
            GOTO 9900
         ENDIF
      ENDIF

      IGLBL = .FALSE.

C     Store the data

      IF (TYPCOD .EQ. 'O') THEN

         DO 100 I = 1, NUM

            IPAR = JPAR+(I-1)*INCPAR

            CALL D1IGSI (CMEM, IMEM, DMEM, IHDR, IGLBL, TYPCOD,
     +         IPAR, 0, IER)
            IF (IER .NE. 0) GOTO 9900

  100    CONTINUE

      ELSE

         DO 200 I = 1, NUM

            IPAR = JPAR+(I-1)*INCPAR
            IVAL = IA(1+(I-1)*INCIA)

            CALL D1IGSI (CMEM, IMEM, DMEM, IHDR, IGLBL, TYPCOD,
     +         IPAR, IVAL, IER)
            IF (IER .NE. 0) GOTO 9900

  200    CONTINUE

      ENDIF

C     Error reporting section

 9900 CONTINUE

      IF (IER .EQ. -11) THEN
          CALL DTERR (2, SUBNAM, IER, NEED)
      ELSE IF (IER .NE. 0) THEN
          CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      RETURN
      END

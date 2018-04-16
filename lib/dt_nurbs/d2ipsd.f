      SUBROUTINE D2IPSD (CMEM, IMEM, DMEM, NUM, DA, INCDA, IGE, TYPCOD,
     +      JPAR, INCPAR, IER)

C     Store a double precision array into an IGES entity parameter space.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     NUM     Number of double precision elements to store
C     DA      Array of double precision elements to store
C     INCDA   Increment to DA, elements (stride).  INCDA = 0
C             results in a "fill" of the same value.
C     IGE     Pointer to an IGES Entity
C     TYPCOD  IGES-context data type code.  Only 'E', 'D', and 'O'
C             are valid here.
C     JPAR    Parameter number of 1st item location.
C     INCPAR  Increment to Parameter numbers (stride)
C   OUTPUT:
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             -1  = Dynamic memory is corrupt or uninitDAlized
C             -2  = IGE is a null pointer
C             -3  = Garbage value found in IGE
C             -4  = Ambiguous--garbage pointer or deleted entity
C             -5  = IGE points to a deleted entity
C             -6  = IGE does not point to an IGES entity.
C             -7  = NUM < 0.
C             -8  = INCDA < 0.
C             -9  = TYPCOD not = 'E', 'D', or 'O'.
C             -10 = JPAR < 2.
C             -11 = INCPAR < 0.
C             -12 = Insufficient storage reserved for DP.
C             -99 = Entity IGE is internally inconsistent.
C
C
C   CALLS:
C     D0PTR    Pointer check
C     D1IGLE   Fetch length parameters of IGES Entity
C     D1IGSD   Store a double precision element to an IGES Entity
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
C     ENTRY D2_IGES_PD_STORE_DBLS (CMEM, IMEM, DMEM, NUM, DA, INCDA,
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

C =====================================================================

      CHARACTER*1    TYPCOD
      INTEGER        NUM, INCDA, IGE, JPAR, INCPAR, IER
      DOUBLE PRECISION DA(*), DVAL

      INTEGER        NPAR, NSTR, NSTRCH, NREAL, IERX
      INTEGER        MPAR, MSTR, MSTRCH, MREAL, NEED
      INTEGER        IHDR, ITYP, I, IPAR
      LOGICAL        IGLBL

      CHARACTER*6    SUBNAM
      DATA SUBNAM   /'D2IPSD'/

C****
      IER = 0

C     Check pointer and type, then locate relevant data spaces
      CALL D0PTR (CMEM, IMEM, DMEM, IGE, ENTGE, 0, IHDR, ITYP, IER)
      IF (IER .LT. 0) GOTO 9900

      IF (NUM .LT. 0) THEN
         IER = -7
         GOTO 9900
      ENDIF

      IF (INCDA .LT. 0) THEN
         IER = -8
         GOTO 9900
      ENDIF

      IF ((TYPCOD .NE. 'O') .AND.
     +    (TYPCOD .NE. 'E') .AND.
     +    (TYPCOD .NE. 'D'))     THEN
         IER = -9
         GOTO 9900
      ENDIF

      IF (INCPAR .LT. 0) THEN
         IER = -11
         GOTO 9900
      ENDIF

C     Fetch the defined lengths

      CALL D1IGLE (CMEM, IMEM, DMEM, IHDR, NPAR, NSTR, NSTRCH,
     +                NREAL, IERX)
      IF (IERX .NE. 0) THEN
         IER = -99
         GOTO 9900
      ENDIF

      IF (JPAR .LT. 1) THEN
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
            IER = -12
            GOTO 9900
         ENDIF
      ENDIF

      IGLBL = .FALSE.

C     Store the data

      IF (TYPCOD .EQ. 'O') THEN

         DO 100 I = 1, NUM

            IPAR = JPAR+(I-1)*INCPAR

            CALL D1IGSD (CMEM, IMEM, DMEM, IHDR, IGLBL, TYPCOD,
     +         IPAR, 0.0D0, NEED, IERX)
            IF (IERX .NE. 0) THEN
               IER = -13
               NEED = NUM-I+1
               GOTO 9900
            ENDIF

  100    CONTINUE

      ELSE

         DO 200 I = 1, NUM

            IPAR = JPAR+(I-1)*INCPAR
            DVAL = DA(1+(I-1)*INCDA)

            CALL D1IGSD (CMEM, IMEM, DMEM, IHDR, IGLBL, TYPCOD,
     +         IPAR, DVAL, NEED, IERX)
            IF (IERX .NE. 0) THEN
               IER = -13
               NEED = NUM-I+1
               GOTO 9900
            ENDIF

  200    CONTINUE

      ENDIF


C     Error reporting section

 9900 CONTINUE

      IF (IER .EQ. -13) THEN
         CALL DTERR (2, SUBNAM, IER, NEED)
      ELSE IF (IER .NE. 0) THEN
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      RETURN
      END

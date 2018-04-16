      SUBROUTINE D2IPFI (CMEM, IMEM, DMEM, NUM, IGE, JPAR, INCPAR,
     +      INCIA, IA, IER)

C     Fetch an integer array from an IGES entity parameter space.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     NUM     Number of integers to fetch
C     IGE     Pointer to an IGES Entity
C     JPAR    Parameter number of 1st item location.
C     INCPAR  Increment to Parameter numbers (stride)
C     INCIA   Increment to IA, elements (stride).  INCIA = 0
C             results in a "fill" of the same value.
C   OUTPUT:
C     IA      Array of integers fetched
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             -1  = Dynamic memory is corrupt or uninitialized
C             -2  = IGE is a null pointer
C             -3  = Garbage value found in IGE
C             -4  = Ambiguous--garbage pointer or deleted entity
C             -5  = IGE points to a deleted entity
C             -6  = IGE does not point to an IGES entity.
C             -7  = NUM < 0.
C             -8  = INCIA < 0.
C             -10 = JPAR is out of range.
C             -11 = JPAR + (NUM-1) * INCPAR is out of range.
C             -12 = INCPAR < 0.
C             -99 = Entity IGE is internally inconsistent.
C
C
C   CALLS:
C     D0PTR    Pointer check
C     D1IGLE   Fetch length parameters of IGES Entity
C     D1IGFI   Fetch an integer element from an IGES Entity
C     DTERR    Error handler
C
C   HISTORY:
C     7/26/93  D. Parsons   Created
C     8/26/93  D. Parsons   Adjusted PAR values to being with 0 (Entity
C                           Type ID)
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_PD_FETCH_INTS (CMEM, IMEM, DMEM, NUM, IGE, JPAR, INCPAR,
C     +      INCIA, IA, IER)

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

      INTEGER        NUM, IA(*), INCIA, IGE, JPAR, INCPAR, IER

      INTEGER        NPAR, NSTR, NSTRCH, NREAL, IERX
      INTEGER        IHDR, ITYP, I, IPAR, IVAL
      LOGICAL        IGLBL

      CHARACTER*6    SUBNAM
      DATA SUBNAM   /'D2IPFI'/

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

      IF (INCPAR .LT. 0) THEN
         IER = -12
         GOTO 9900
      ENDIF

C     Fetch the defined lengths

      CALL D1IGLE (CMEM, IMEM, DMEM, IHDR, NPAR, NSTR, NSTRCH,
     +                NREAL, IERX)
      IF (IERX .NE. 0) THEN
         IER = -99
         GOTO 9900
      ENDIF

      IF ((JPAR .LT. 0) .OR. (JPAR .GT. NPAR)) THEN
         IER = -10
         GOTO 9900
      ENDIF

      IF ((JPAR+(NUM-1)*INCPAR) .GT. NPAR) THEN
         IER = -11
         GOTO 9900
      ENDIF

      IGLBL = .FALSE.

C     Fetch the data

      DO 100 I = 1, NUM

         IPAR = JPAR+(I-1)*INCPAR

         CALL D1IGFI (CMEM, IMEM, DMEM, IHDR, IGLBL, IPAR,
     +       IVAL, IER)
         IF (IER .NE. 0) GOTO 9900

         IA(1+(I-1)*INCIA) = IVAL

  100 CONTINUE


C     Error reporting section

 9900 CONTINUE

      IF (IER .NE. 0) THEN
          CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      RETURN
      END

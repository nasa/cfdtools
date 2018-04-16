      SUBROUTINE D2BFDF (CMEM, IMEM, DMEM, LEND, LABEL, IBF, IER)

C     PURPOSE:
C        B-spline Function DeFine.  Create a B-spline Function data
C        entity with room for a C-array of length LEND and with label
C        LABEL (which can be a null string).
C
C        B-Spline Function Entity:
C
C           CMEM  [1..len()]  LABEL
C
C           DMEM  [1..LEND]   C-Array
C
C           IMEM  (not used)
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        LEND     Length of the C-Array vector
C        LABEL    Label for the entity
C
C     OUTPUT:
C        IBF      New entity ID
C        IER      Error flag.  If negative, the entity is not defined.
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = Insufficient space available for character data
C                       (LABEL)
C                 -3  = Insufficient space available for entity header
C                 -4  = Insufficient space available for double
C                       precision data
C                 -5  = Negative double precision data length requested
C                 -999= Unexpected error in a called subroutine
C
C     CALLS:
C        D1MBAD
C        D0TRMC
C        D1DEFE
C        D2STAC
C        DTERR
C
C
C     HISTORY:
C        12May92  D. Parsons   Created.
C        18May92  D. Parsons   Make error returns sequential.
C        11Aug92  D. Parsons   Change D2DEFE call to D1DEFE.
C
C     ------

C     Long name alias:
C        ENTRY D2_BSPLINE_FUNCTION_DEFINE (CMEM, IMEM, DMEM, LEND
C    +         LABEL, IBF, IER)

      EXTERNAL D0TRMC
      INTEGER  D0TRMC

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      EXTERNAL D1MBAD
      LOGICAL  D1MBAD

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
C                          B-spline Function
      INTEGER     ENTBF
      PARAMETER  (ENTBF = 246)

C =====================================================================

      INTEGER        LEND, IBF, IER, IERX
      CHARACTER*(*)  LABEL
      INTEGER        IDTYP, LENC, LENI, IDE
      LOGICAL        INITIZ

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2BFDF'/

C     ------

      IER = 0
      IBF = 0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

      IF (LEND .LT. 0) THEN
         IER = -5
         GOTO 9000
      ENDIF

      IDTYP  = ENTBF
      LENC   = D0TRMC(LABEL)
      LENI   = 0
      INITIZ = .TRUE.

C     Call D1DEFE to define a new entity

      CALL D1DEFE (CMEM, IMEM, DMEM, IDTYP, LENC, LENI, LEND, INITIZ,
     +      IDE, IER)
      IF (IER .NE. 0) GOTO 9000

C     Store the LABEL in the Character memory

      IF (LENC .GT. 0) THEN
         CALL D2STAC (CMEM, IMEM, DMEM, LABEL, 1, LENC, IDE, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
      ENDIF

      IBF = IDE

      GOTO 9999

C     Error Handling

 9000 CONTINUE

      IF (IER .EQ. -999) THEN
         CALL DTERR (5, SUBNAM, IER, 0)
      ELSE
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      IBF = 0

 9999 CONTINUE

      RETURN
      END

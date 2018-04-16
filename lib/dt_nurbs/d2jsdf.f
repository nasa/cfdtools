      SUBROUTINE D2JSDF (CMEM, IMEM, DMEM, MAXTS, LABEL, IJS, IER)

C     PURPOSE:
C        Joined Surface DeFine.  Create a Joined Surface data entity
C        with room for MAXTS Trimmed surfaces.
C
C        Output a pointer to  the newly allocated Joined Surface in IJS.
C
C        Joined Surface Entity:
C
C           CMEM  [1..len()]     LABEL
C
C           DMEM  (not used)
C
C           IMEM  [1]            ICLS  (closed flag -1, 0, +1)
C                 [2]            MAXTS
C                 [3..2+MAXTS]   ITS(IXTS)
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        MAXTS    Maximum Trimmed Surface pointers
C        LABEL    Label for the entity
C
C     OUTPUT:
C        IJS      New entity ID
C        IER      Error flag.  If negative, the entity is not defined.
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = Insufficient space available for character data
C                       (LABEL)
C                 -3  = Insufficient space available for integer data
C                 -4  = MAXTS < 0
C                 -999= Unexpected error in a called subroutine
C
C     CALLS:
C        D1MBAD
C        D0TRMC
C        D1DEFE
C        D2STAC
C        D2STAI
C        DTERR
C
C     HISTORY:
C        22Jun92  D. Parsons   Created.
C        11Aug92  D. Parsons   Change D2DEFE call to D1DEFE.
C
C     ------

C     Long name alias:
C        ENTRY D2_JOINED_SURFACE_DEFINE (CMEM, IMEM, DMEM, MAXTS,
C    +      LABEL, IJS, IER)

      EXTERNAL D0TRMC, D1MBAD
      INTEGER  D0TRMC
      LOGICAL  D1MBAD

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

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
C                          Joined Surface
      INTEGER     ENTJS
      PARAMETER  (ENTJS = 240)

C =====================================================================

      INTEGER        MAXTS, IJS, ICLS, IER
      CHARACTER*(*)  LABEL
      INTEGER        IDTYP, LENC, LENI, LEND, IDE
      INTEGER        IERX, IA(2)
      LOGICAL        INITIZ

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2JSDF'/

C     ------

      IER = 0
      IJS = 0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

      IDTYP  = ENTJS
      LENC   = D0TRMC(LABEL)
      LENI   = MAXTS+2
      LEND   = 0
      INITIZ = .TRUE.

C     Check the validity of MAXTS

      IF (MAXTS .LT. 0) THEN
         IER = -4
         GOTO 9000
      ENDIF

C     Call D1DEFE to define a new entity

      CALL D1DEFE (CMEM, IMEM, DMEM, IDTYP, LENC, LENI, LEND, INITIZ,
     +      IDE, IERX)
      IF (IERX .NE. 0) THEN
         IF       (IERX .EQ. -2) THEN
            IER = -2
         ELSE IF  (IERX .EQ. -3) THEN
            IER = -3
         ENDIF
         GOTO 9000
      ENDIF

C     Store the LABEL in the Character memory

      IF (LENC .GT. 0) THEN
         CALL D2STAC (CMEM, IMEM, DMEM, LABEL, 1, LENC, IDE, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
      ENDIF

C     Store the Integer data into the Integer memory

      ICLS = 0

      IA(1) = ICLS
      IA(2) = MAXTS

      CALL D2STAI (CMEM, IMEM, DMEM, IA, 1, 2, IDE, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9000
      ENDIF

      IJS = IDE

      GOTO 9999

C     Error Handling

 9000 CONTINUE

      IF (IER .EQ. -999) THEN
         CALL DTERR (5, SUBNAM, IER, 0)
      ELSE
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      IJS = 0

 9999 CONTINUE

      RETURN
      END

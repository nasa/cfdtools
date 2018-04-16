      SUBROUTINE D2LPDF (CMEM, IMEM, DMEM, MAXEG, ITS, IXLP, LABEL,
     +      ILP, IER)

C     PURPOSE:
C        LooP DeFine.  Create a Loop of a Trimmed Surface data entity,
C        with space for MAXEG Edge pointers, for the IXLPth Loop in
C        Trimmed Surface ITS, and with label LABEL.  Any of the pointers
C        may be null, represented by a value of zero, during this call.
C
C        Output a pointer to the newly allocated Loop in ILP.
C
C        Loop Entity:
C
C           CMEM  [1..len()]     LABEL
C
C           DMEM  (not used)
C
C           IMEM  [1]            ITS
C                 [2]            IXLP
C                 [3]            MAXEG
C                 [4..3+MAXEG]   IEG(IXEG)
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        MAXEG    Maximum edge pointers
C        ITS      Pointer to Trimmed Surface Entity
C        IXLP     Loop index in ITS Trimmed Surface Entity
C        LABEL    Label for the entity
C
C     OUTPUT:
C        ILP      New entity ID
C        IER      Error flag.  If negative, the entity is not defined.
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = Insufficient space available for character data
C                       (LABEL)
C                 -3  = Insufficient space available for integer data
C                 -4  = MAXEG < 0
C                 -5  = ITS does not point to a Trimmed Surface
C                       entity
C                 -6  = IXLP is not in the range [1..MAXLP] as defined
C                       in entity ITS
C                 -7  = Index location IXLP in Trimmed Surface ITS
C                       is already assigned (non-zero)
C                 -999= Unexpected error in a called subroutine
C
C     CALLS:
C        D1MBAD
C        D0TRMC
C        D0TSFP
C        D1DEFE
C        D2FEEI
C        D2STAC
C        D2STAI
C        D2STEI
C        DTERR
C
C     HISTORY:
C        20May92  D. Parsons   Created.
C        11Aug92  D. Parsons   Change D2DEFE call to D1DEFE.
C
C     ------

C     Long name alias:
C        ENTRY D2_LOOP_DEFINE (CMEM, IMEM, DMEM, MAXEG, ITS, IXLP,
C     LABEL, ILP, IER)

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
C                          Loop of a Trimmed Surface
      INTEGER     ENTLP
      PARAMETER  (ENTLP = 242)

C =====================================================================

      INTEGER        MAXEG, ITS, IXLP, ILP, IER
      CHARACTER*(*)  LABEL
      CHARACTER      LABELX
      INTEGER        IDTYP, LENC, LENI, LEND, IDE, IA(3), IHDR
      INTEGER        IBFS, IJS, IXTS, IORIEN, MAXLP, IERX, LPPTR
      LOGICAL        INITIZ

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2LPDF'/

C     ------

      IER = 0
      ILP = 0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

      IDTYP  = ENTLP
      LENC   = D0TRMC(LABEL)
      LENI   = MAXEG+3
      LEND   = 0
      INITIZ = .TRUE.

C     Check the validity of MAXEG

      IF (MAXEG .LT. 0) THEN
         IER = -4
         GOTO 9000
      ENDIF

C     Check the validity of ITS

      IF (ITS .NE. 0) THEN
         CALL D0TSFP (CMEM, IMEM, DMEM, ITS, LABELX, IBFS, IJS,
     +         IXTS, IORIEN, MAXLP, IHDR, IERX)
         IF (IERX .NE. 0) THEN
C           ITS does not point to a Trimmed Surface Entity
            IER = -5
            GOTO 9000
         ENDIF
      ENDIF

C     Check the validity of IXLP

      IF (IXLP .NE. 0) THEN
         IF ((ITS .EQ. 0)
     +         .OR. (IXLP .GT. MAXLP)
     +         .OR. (IXLP .LT. 0)) THEN
            IER = -6
            GOTO 9000
         ENDIF
         CALL D2FEEI (CMEM, IMEM, DMEM, ITS, 5+IXLP, LPPTR, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
         IF (LPPTR .NE. 0) THEN
            IER = -7
            GOTO 9000
         ENDIF
      ELSE
         IF (ITS .NE. 0) THEN
            IER = -6
            GOTO 9000
         ENDIF
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

C     Store the pointer to this new entity in the corresponding
C        Trimmed Surface list

      IF (ITS .NE. 0) THEN
         CALL D2STEI (CMEM, IMEM, DMEM, IDE, 5+IXLP, ITS, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
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

      IA(1) = ITS
      IA(2) = IXLP
      IA(3) = MAXEG

      CALL D2STAI (CMEM, IMEM, DMEM, IA, 1, 3, IDE, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9000
      ENDIF

      ILP = IDE

      GOTO 9999

C     Error Handling

 9000 CONTINUE

      IF (IER .EQ. -999) THEN
         CALL DTERR (5, SUBNAM, IER, 0)
      ELSE
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      ILP = 0

 9999 CONTINUE

      RETURN
      END

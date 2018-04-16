      SUBROUTINE D2EGDF (CMEM, IMEM, DMEM, IBFC, ILP, IXEG, JEG, LABEL,
     +      IEG, IER)

C     PURPOSE:
C        EdGe DeFine.  Create an Edge of a Trimmed Surface data entity,
C        with B-spline curve Function pointer IBFC, which is the IXEGth
C        segment of the Loop pointer ILP, is joined to Edge JEG, and
C        has label LABEL.  Any of the pointers may be null, represented
C        by a value of zero, during this call.
C
C        Output pointer to the newly allocated Edge in IEG.
C
C        Edge Entity:
C
C           CMEM  [1..len()]  LABEL
C
C           DMEM  (not used)
C
C           IMEM  [1]         IBFC
C                 [2]         ILP
C                 [3]         IXEG
C                 [4]         JEG
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        IBFC     Pointer to B-spline curve function
C        ILP      Pointer to Loop Entity
C        IXEG     Edge index in ILP Loop Entity
C        JEG      Pointer to adjoining Edge Entity
C        LABEL    Label for the entity
C
C     OUTPUT:
C        IEG      New entity ID
C        IER      Error flag.  If negative, the entity is not defined.
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = Insufficient space available for character data
C                       (LABEL)
C                 -3  = Insufficient space available for integer data
C                 -4  = IBFC does not point to a B-spline curve entity
C                 -5  = ILP does not point to a Loop entity
C                 -6  = IXEG is not in the range [1..MAXEG] as defined
C                       in entity ILP
C                 -7  = Index location IXEG in Loop ILP is already
C                       assigned (non-zero)
C                 -8  = JEG does not point to an Edge entity
C                 -999= Unexpected error in a called subroutine
C
C     CALLS:
C        D0BFFP
C        D0EGFP
C        D0LPFP
C        D1MBAD
C        D0TRMC
C        D1DEFE
C        D2FEEI
C        D2STAC
C        D2STAI
C        D2STEI
C        DTERR
C
C
C     HISTORY:
C        21May92  D. Parsons   Created.
C        11Aug92  D. Parsons   Change D2DEFE call to D1DEFE.
C
C     ------

C     Long name alias:
C        ENTRY D2_EDGE_DEFINE (CMEM, IMEM, DMEM, IBFC, ILP, IXEG, JEG, LABEL,
C    +      IEG, IER)

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
C                          Edge of a Trimmed Surface
      INTEGER     ENTEG
      PARAMETER  (ENTEG = 243)

C =====================================================================

      INTEGER        IBFC, IXEG, ILP, JEG, IER
      INTEGER        IBFCX, IXEGX, ILPX, JEGX, LENDX
      CHARACTER*(*)  LABEL
      CHARACTER      LABELX
      INTEGER        IDTYP, LENC, LENI, LEND, IDE
      INTEGER        ITS, IXLP, MAXEG, IHDR
      INTEGER        IERX, EGPTR, IA(4)
      LOGICAL        INITIZ

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2EGDF'/

C     ------

      IER = 0
      IEG = 0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

      IDTYP  = ENTEG
      LENC   = D0TRMC(LABEL)
      LENI   = 4
      LEND   = 0
      INITIZ = .TRUE.

C     Check the validity of IFBC

      IF (IBFC .NE. 0) THEN
         CALL D0BFFP (CMEM, IMEM, DMEM, IBFC, LABELX, LENDX, IHDR,
     +         IERX)
         IF (IERX .NE. 0) THEN
C           IFBC does not point to a B-spline Entity
            IER = -4
            GOTO 9000
         ENDIF
      ENDIF

C     Check the validity of ILP

      IF (ILP .NE. 0) THEN
         CALL D0LPFP (CMEM, IMEM, DMEM, ILP, LABELX, ITS, IXLP,
     +         MAXEG, IHDR, IERX)
         IF (IERX .NE. 0) THEN
C           ILP does not point to a Loop entity
            IER = -5
            GOTO 9000
         ENDIF
      ENDIF

C     Check the validity of IXEG

      IF (IXEG .NE. 0) THEN
         IF ((ILP .EQ. 0)
     +         .OR. (IXEG .GT. MAXEG)
     +         .OR. (IXEG .LT. 0)) THEN
            IER = -6
            GOTO 9000
         ENDIF
         CALL D2FEEI (CMEM, IMEM, DMEM, ILP, 3+IXEG, EGPTR, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
         IF (EGPTR .NE. 0) THEN
            IER = -7
            GOTO 9000
         ENDIF
      ELSE
         IF (ILP .NE. 0) THEN
            IER = -6
            GOTO 9000
         ENDIF
      ENDIF

C     Check the validity of JEG

      IF (JEG .NE. 0) THEN
         CALL D0EGFP (CMEM, IMEM, DMEM, JEG, LABELX, IBFCX, ILPX,
     +         IXEGX, JEGX, IHDR, IERX)
         IF (IERX .NE. 0) THEN
C           JEG does not point to an Edge entity
            IER = -8
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
C        Loop list

      IF (ILP .NE. 0) THEN
         CALL D2STEI (CMEM, IMEM, DMEM, IDE, 3+IXEG, ILP, IERX)
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

      IA(1) = IBFC
      IA(2) = ILP
      IA(3) = IXEG
      IA(4) = JEG

      CALL D2STAI (CMEM, IMEM, DMEM, IA, 1, 4, IDE, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9000
      ENDIF

      IEG = IDE

      GOTO 9999

C     Error Handling

 9000 CONTINUE

      IF (IER .EQ. -999) THEN
         CALL DTERR (5, SUBNAM, IER, 0)
      ELSE
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      IEG = 0

 9999 CONTINUE

      RETURN
      END


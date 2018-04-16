      SUBROUTINE D2TSDF (CMEM, IMEM, DMEM, MAXLP, IBFS, IJS, IXTS,
     +      IORIEN, LABEL, ITS, IER)

C     PURPOSE:
C        Trimmed Surface DeFine.  Create a Trimmed Surface of a
C        Joined Surface data entity with room for MAXLP Loops, with
C        surface B-spline function IBFS, fro the IXTSth Trimmed
C        Surface in Joined Surface IJS, with orientation flag IORIEN,
C        and with labe LABEL.  Any of the pointers may be null,
C        represented by a value of zero, during this call.
C
C        Output a pointer to the newly allocated Trimmed Surface in
C        ITS.
C
C        Loop Entity:
C
C           CMEM  [1..len()]     LABEL
C
C           DMEM  (not used)
C
C           IMEM  [1]            IBFS
C                 [2]            IJS
C                 [3]            IXTS
C                 [4]            IORIEN
C                 [5]            MAXLP
C                 [6..5+MAXLP]   ITS(IXLP)
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        MAXLP    Maximum number of Loop entities
C        IBFS     Pointer to B-Spline surface function
C        IJS      Pointer to Joined Surface entity
C        IXTS     Trimmed Surface index in IJS Joined Surface
C                 entity
C        IORIEN   Orientation flag:
C                 +1 The normal for the trimmed surface agrees
C                    with the normal of the underlying surface
C                 -1 The normal for the trimmed surface is opposite
C                    to the normal of the underlying surface
C        LABEL    Label for the entity
C
C     OUTPUT:
C        ITS      New entity ID
C        IER      Error flag.  If negative, the entity is not defined.
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = Insufficient space available for character data
C                       (LABEL)
C                 -3  = Insufficient space available for integer data
C                 -4  = MAXLP < 0.
C                 -5  = IBFS does not point to a B-spline surface entity
C                 -6  = IJS does not point to a Joined Surface entity
C                 -7  = IXTS is not in the range [1..MAXTS] as defined
C                       in entity IJS
C                 -8  = Index location IXTS in Joined Surface IJS is
C                       already assigned (non-zero)
C                 -9  = IORIEN not equal to 1 or -1
C                 -999= Unexpected error in a called subroutine
C
C     CALLS:
C        D0JSFP
C        D1MBAD
C        D0TRMC
C        D1DEFE
C        D2FEED
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
C        ENTRY D2_TRIMMED_SURFACE_DEFINE (CMEM, IMEM, DMEM, MAXLP,
C           IBFS, IJS, IXTS, IORIEN, LABEL, ITS, IER)

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
C                          Trimmed Surface of a Joined Surface
      INTEGER     ENTTS
      PARAMETER  (ENTTS = 241)

C =====================================================================

      INTEGER           MAXLP, IBFS, IJS, IXTS, IORIEN, ITS, IER
      CHARACTER*(*)     LABEL
      CHARACTER         LABELX
      INTEGER           IDTYP, LENC, LENI, LEND, IDE
      INTEGER           ICLS, MAXTS, IA(5), IHDR
      INTEGER           IERX, LENDX, TSPTR
      DOUBLE PRECISION  INDEP
      LOGICAL           INITIZ

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2TSDF'/

C     ------

      IER = 0
      ITS = 0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

C     Check MAXLP .GE. 0

      IF (MAXLP .LT. 0) THEN
         IER = -4
         GOTO 9000
      ENDIF

      IDTYP  = ENTTS
      LENC   = D0TRMC(LABEL)
      LENI   = MAXLP+5
      LEND   = 0
      INITIZ = .TRUE.

C     Check the validity of IFBS

      IF (IBFS .NE. 0) THEN
         CALL D0BFFP (CMEM, IMEM, DMEM, IBFS, LABELX, LENDX,
     +         IHDR, IERX)
         IF (IERX .NE. 0) THEN
C           IFBS does not point to a B-spline Entity
            IER = -5
            GOTO 9000
         ELSE
C           Make a cursory check... is this a Surface?
            CALL D2FEED (CMEM, IMEM, DMEM, IBFS, 1, INDEP, IERX)
            IF ((IERX .NE. 0) .OR.
     +            ((INDEP .NE. 0) .AND. (INDEP .NE. 2))) THEN
               IER = -5
               GOTO 9000
            ENDIF
         ENDIF
      ENDIF

C     Check the validity of IJS

      IF (IJS .NE. 0) THEN
         CALL D0JSFP (CMEM, IMEM, DMEM, IJS, LABELX, ICLS, MAXTS,
     +         IHDR, IERX)
         IF (IERX .NE. 0) THEN
C           IJS does not point to a Joined Surface entity
            IER = -6
            GOTO 9000
         ENDIF
      ENDIF

C     Check the validity of IXTS

      IF (IXTS .NE. 0) THEN
         IF ((IJS .EQ. 0)
     +         .OR. (IXTS .GT. MAXTS)
     +         .OR. (IXTS .LT. 0)) THEN
            IER = -7
            GOTO 9000
         ENDIF
         CALL D2FEEI (CMEM, IMEM, DMEM, IJS, 2+IXTS, TSPTR, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
         IF (TSPTR .NE. 0) THEN
            IER = -8
            GOTO 9000
         ENDIF
      ELSE
         IF (IJS .NE. 0) THEN
            IER = -7
            GOTO 9000
         ENDIF
      ENDIF

      IF ((IORIEN .NE. 1) .AND. (IORIEN .NE. -1)) THEN
         IER = -9
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

C     Store the pointer to this new entity in the corresponding
C        Trimmed Surface list

      IF (IJS .NE. 0) THEN
         CALL D2STEI (CMEM, IMEM, DMEM, IDE, 2+IXTS, IJS, IERX)
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

      IA(1) = IBFS
      IA(2) = IJS
      IA(3) = IXTS
      IA(4) = IORIEN
      IA(5) = MAXLP

      CALL D2STAI (CMEM, IMEM, DMEM, IA, 1, 5, IDE, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9000
      ENDIF

      ITS = IDE

      GOTO 9999

C     Error Handling

 9000 CONTINUE

      IF (IER .EQ. -999) THEN
         CALL DTERR (5, SUBNAM, IER, 0)
      ELSE
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      ITS = 0

 9999 CONTINUE

      RETURN
      END

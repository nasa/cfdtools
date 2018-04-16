      SUBROUTINE D2PJDF (CMEM, IMEM, DMEM, LABEL, IPJ, IER)

C     PURPOSE:
C        Point on a Joined surface DeFine.  Create a Point on a Joined
C        Surface data entity.
C
C        Output a pointer to the newly allocated Point on a Joined
C        Surface in IPJ.
C
C        Point on a Joined Surface Entity:
C
C           CMEM  [1..len()]     LABEL
C
C           DMEM  [1]            T (Parameter value for Edge)
C                 [2,3]          U,V (Parameter values for Surface)
C                 [4,5,6]        X,Y,Z (3-space coordinates)
C
C           IMEM  [1]            IEG
C                 [2]            ITS
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        LABEL    Label for the entity
C
C     OUTPUT:
C        IPJ      New entity ID
C        IER      Error flag.  If negative, the entity is not defined.
C                 -1  = Dynamic memory is corrupt or uninitialized
C                 -2  = Insufficient space available for character data
C                       (LABEL)
C                 -3  = Insufficient space available for integer data
C                 -4  = Insufficient space available for double
C                       precision data
C                 -999= Unexpected error in a called subroutine
C
C     CALLS:
C        D1MBAD
C        D0TRMC
C        D1DEFE
C        D2STAC
C        DTERR
C
C     HISTORY:
C        22Jun92  D. Parsons   Created.
C        11Aug92  D. Parsons   Change D2DEFE call to D1DEFE.
C
C     ------

C     Long name alias:
C        ENTRY D2_POINT_JOINED_DEFINE (CMEM, IMEM, DMEM, LABEL, IPJ,
C    +      IER)

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
C                          Point on a Joined Surface
      INTEGER     ENTPJ
      PARAMETER  (ENTPJ = 239)

C =====================================================================

      INTEGER        IPJ, IER
      CHARACTER*(*)  LABEL
      INTEGER        IDTYP, LENC, LENI, LEND, IDE
      INTEGER        IERX
      LOGICAL        INITIZ

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2PJDF'/

C     ------

      IER = 0
      IPJ = 0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

      IDTYP  = ENTPJ
      LENC   = D0TRMC(LABEL)
      LENI   = 2
      LEND   = 6
      INITIZ = .TRUE.

C     Call D1DEFE to define a new entity

      CALL D1DEFE (CMEM, IMEM, DMEM, IDTYP, LENC, LENI, LEND, INITIZ,
     +      IDE, IERX)
      IF (IERX .NE. 0) THEN
         IF       (IERX .EQ. -2) THEN
            IER = -2
         ELSE IF  (IERX .EQ. -3) THEN
            IER = -3
         ELSE IF  (IERX .EQ. -4) THEN
            IER = -4
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

      IPJ = IDE

      GOTO 9999

C     Error Handling

 9000 CONTINUE

      IF (IER .EQ. -999) THEN
         CALL DTERR (5, SUBNAM, IER, 0)
      ELSE
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      IPJ = 0

 9999 CONTINUE

      RETURN
      END

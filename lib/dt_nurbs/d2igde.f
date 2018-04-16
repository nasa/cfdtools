      SUBROUTINE D2IGDE (CMEM, IMEM, DMEM, NPAR, NSTR, NSTRCH, NREAL,
     +                   IGI, JDE, IGE, IER)

C     PURPOSE:
C        IGES Define Entity.  Create an IGES entity.
C
C        Output pointer to the newly allocated entity in IGE.
C
C        IGES Entity: If not Global (IGLOBL = .FALSE.)
C
C           CMEM  [1..8]      Directory Entry 16.
C                 [9..16]     Directory Entry 17.
C                 [17..24]    Directory Entry 18.
C                 [25]        TYPCOD for Parameter 2.
C                 [26]        TYPCOD for Parameter 3.
C                 [27]        TYPCOD for Parameter 4.
C                  .                   .
C                  .                   .
C
C           IMEM  [1]         Pointer to the IGES Index entity.
C                 [2]         Pointer to the Character Sequence entity.
C                 [3]         Total number of parameter data values
C                 [4]         DE 1, Entity Type Number
C                 [5]         DE 2, Parameter Data initial line number
C                 [6]         DE 3, Structure
C                 [7]         DE 4, Line Font Pattern
C                 [8]         DE 5, Level
C                 [9]         DE 6, View
C                 [10]        DE 7, Transformation Matrix
C                 [11]        DE 8, Label Display Associativity
C                 [12]        DE 9, Status Number
C                 [13]        DE 10, Sequence Number
C                 [14]        DE 12, Line Weight Number
C                 [15]        DE 13, Color Number
C                 [16]        DE 14, Parameter Line Count Number
C                 [17]        DE 15, Form Number
C                 [18]        DE 19, Entity Subscript Number
C                 [19]        Parameter 2 Value or reference
C                 [20]        Parameter 3 Value or reference
C                  .                     .
C                  .                     .
C
C           DMEM  [1]         1st Real-valued parameter
C                 [2]         2nd Real-valued parameter
C                  .                     .
C                  .                     .
C
C
C        IGES Entity: If Global (IGLOBL = .TRUE.)
C
C           CMEM  [1]
C                  .          TYPCOD's for Global parameters.
C                  .          See table in documentation.
C                 [23]
C
C
C           IMEM  [1]         Pointer to the IGES Index entity.
C                 [2]         Pointer to the Character Sequence entity.
C                 [3]         = 23; Total number of parameter data values
C                 [4]         = -11 Global entity type number
C                 [5]         Pointer into CQ entity for product ID
C                  .                     .
C                  .          pointers or values for Global parameters.
C                  .          See table in documentation.
C                  .                     .
C                 [27]        Pointer into CQ entity for mod Date
C
C           DMEM  [1]         Model Space Scale
C                 [2]         Width of maximum line weight
C                 [3]         Minimum user-intended resolution
C                 [4]         Approximate maximum coordinate
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        NPAR     Number of Parameters
C        NSTR     Number of string parameters to allocate
C        NSTRCH   Total number of characters in strings
C        NREAL    Number of double precision parameters
C        IGLOBL   Global flag (=1 indicates special Global entity)
C        IGI      Pointer to Index entity
C        JDE      IGES DE (directory entry) sequence number of this 
C                 new entity.  Always odd.
C
C     OUTPUT:
C        IGE      New IGES entity ID
C        ICQ      New Character sequence entity ID
C        IER      Error flag.  If negative, the entity is not defined.
C                 -1  = Dynamic memory is corrupt or uninitialized.
C                 -2  = Insufficient space available for character data.
C                 -3  = Insufficient space available for integer data.
C                 -4  = Insufficient space available for double
C                       precision data.
C                 -5  = NPAR < 0.
C                 -6  = NSTR < 0.
C                 -7  = NSTRCH < 0.
C                 -8  = NREAL < 0.
C                 -9  = IGI not a valid IGES Index Entity.
C                 -10 = JDE < 1.
C                 -11 = JDE > NDE as defined when the Index was created.
C                 -12 = JDE even.
C                 -999= Unexpected error in a called subroutine.
C
C     CALLS:
C        D0PTR    Pointer check
C        D1IGDE   Define IGES entity (shadow)
C        D1STEI   Store integer element
C        DTERR    Error handler
C
C
C     HISTORY:
C        16Jun93  D. Parsons   Created.
C
C     ------

C     Long name alias:
C        ENTRY D2_IGES_DEFINE_ENTITY (CMEM, IMEM, DMEM, NPAR, NSTR,
C    +            NSTRCH, NREAL, IGI, JDE, IGE, IER)

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
C                          IGES Entity Data
      INTEGER     ENTGE
      PARAMETER  (ENTGE = 238)

C                          IGES Index
      INTEGER     ENTGI
      PARAMETER  (ENTGI = 237)

C =====================================================================

      INTEGER         NPAR, NSTR, NSTRCH, NREAL, IGI, JDE, IGE, IER
      INTEGER         IHDR, ITYP, IERX, NDE, ICQ, NEED
      LOGICAL         IGLOBL
      CHARACTER*6     SUBNAM
      DATA            SUBNAM /'D2IGDE'/
      IER  = 0
      IERX = 0
      IGE  = 0
      ICQ  = 0

C     Error checking

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

      IF (NPAR .LT. 0) THEN
          IER = -5
          GOTO 9000
      ENDIF

      IF (NSTR .LT. 0) THEN
          IER = -6
          GOTO 9000
      ENDIF

      IF (NSTRCH .LT. 0) THEN
          IER = -7
          GOTO 9000
      ENDIF

      IF (NREAL .LT. 0) THEN
          IER = -8
          GOTO 9000
      ENDIF

      CALL D0PTR (CMEM, IMEM, DMEM, IGI, ENTGI, 0, IHDR,
     +    ITYP, IERX)
      IF (IERX .NE. 0) THEN
          IER = -9
          GOTO 9000
      ENDIF

      IF (JDE .LT. 1) THEN
          IER = -10
          GOTO 9000
      ENDIF
      
      CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 2, NDE, IERX) 
      IF (JDE .GT. NDE) THEN
          IER = -11
          GOTO 9000
      ENDIF 
      
      IF (MOD(JDE,2) .NE. 1) THEN
         IER = -12
         GOTO 9000
      ENDIF

C     Define the new IGES Entity

      IGLOBL = .FALSE.
      CALL D1IGDE (CMEM, IMEM, DMEM, NPAR, NSTR, NSTRCH, NREAL,
     +      IGLOBL, IGI, IGE, ICQ, NEED, IER)
      IF (IER .NE. 0) GOTO 9000

C     Update the IGES Index (JDE)

      CALL D1STEI (CMEM, IMEM, DMEM, IGE, 2+(JDE+1)/2, IHDR, IERX)

 9000 CONTINUE

      IF (IER .NE. 0) THEN
         IF (IER .EQ. -999) THEN
            CALL DTERR (5, SUBNAM, IER, 0)
         ELSE IF ((IER .EQ. -2) .OR. 
     +            (IER .EQ. -3) .OR. 
     +            (IER .EQ. -4)) THEN
            CALL DTERR (2, SUBNAM, IER, NEED)
         ELSE
            CALL DTERR (1, SUBNAM, IER, 0)
         ENDIF
         IGE = 0
      ENDIF

      RETURN
      END

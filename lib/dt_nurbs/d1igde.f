      SUBROUTINE D1IGDE (CMEM, IMEM, DMEM, NPAR, NSTR, NSTRCH, NREAL,
     +      IGLOBL, IGI, IGE, ICQ, NEED, IER)

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
C                 [4]         Pointer to first free location in DMEM
C                 [5]         DE 1, Entity Type Number
C                 [6]         DE 2, Parameter Data initial line number
C                 [7]         DE 3, Structure
C                 [8]         DE 4, Line Font Pattern
C                 [9]         DE 5, Level
C                 [10]        DE 6, View
C                 [11]        DE 7, Transformation Matrix
C                 [12]        DE 8, Label Display Associativity
C                 [13]        DE 9, Status Number
C                 [14]        DE 10, Sequence Number
C                 [15]        DE 12, Line Weight Number
C                 [16]        DE 13, Color Number
C                 [17]        DE 14, Parameter Line Count Number
C                 [18]        DE 15, Form Number
C                 [19]        DE 19, Entity Subscript Number
C                 [20]        Parameter 2 Value or reference
C                 [21]        Parameter 3 Value or reference
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
C                 [4]         Pointer to first free location in DMEM
C                 [5]         = -11 Global entity type number
C                 [6]         Pointer into CQ entity for product ID
C                  .                     .
C                  .          pointers or values for Global parameters.
C                  .          See table in documentation.
C                  .                     .
C                 [28]        Pointer into CQ entity for mod Date
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
C
C     OUTPUT:
C        IGE      New IGES entity ID
C        ICQ      New Character sequence entity ID    
C        NEED     Amount of space needed, if insufficient space available.
C        IER      Error flag.  If negative, the entity is not defined.
C                 -2  = Insufficient space available for character data.
C                 -3  = Insufficient space available for integer data.
C                 -4  = Insufficient space available for double
C                       precision data.
C                 -999= Unexpected error in a called subroutine.
C
C     CALLS:
C        D1DEFE   Define a new entity
C        D1CQDF   Define a character sequence entity
C        D2STAC   Store a character array
C        D2STAI   Store an integer array
C        D2ERAS   Erase an entity
C
C
C     HISTORY:
C        26May93  D. Parsons   Created.
C        29Jun93  D. Parsons   Added DMEM free-space stack.
C
C     ------


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
C                          IGES Entity Data
      INTEGER     ENTGE
      PARAMETER  (ENTGE = 238)

C                          Character String Sequence
      INTEGER     ENTCQ
      PARAMETER  (ENTCQ = 249)

C =====================================================================

      INTEGER     NPAR, NSTR, NSTRCH, NREAL, IGI, IGE, ICQ, IER
      INTEGER     IDTYP, LENC, LENI, LEND, IERX, NEED
      CHARACTER   INSTR
      LOGICAL     IGLOBL, INITIZ
      INTEGER     IA(28), GLBLIA(24), JLO, JHI, I
      CHARACTER   CA*23,  GLBLCA*23
      DOUBLE PRECISION  DELM

      DATA        GLBLCA /'CCCCIIIIICEICIECEECCIIC'/
      DATA        GLBLIA /-11,
     +               1,2,3,4,0,0,0,0,0,5,1,0,6,0,2,7,3,4,8,9,0,0,10/

C     ------

C     Call D1DEFE to define a new entity

      IDTYP  = ENTGE
      INITIZ = .TRUE.
      NEED   = 0
      IER    = 0

      IF (IGLOBL) THEN
         LENC   = 23
         LENI   = 28
      ELSE
         LENC   = 24 + NPAR
         LENI   = 19 + NPAR
      ENDIF

      LEND   = NREAL

      CALL D1DEFE (CMEM, IMEM, DMEM, IDTYP, LENC, LENI, LEND, INITIZ,
     +      IGE, IER)
      IF (IER .NE. 0) THEN
         IF (IER .EQ. -2) THEN
            NEED = LENC
         ELSE IF (IER .EQ. -3) THEN
            NEED = LENI
         ELSE IF (IER .EQ. -4) THEN
            NEED = LEND
         ENDIF
         GOTO 9000
      ENDIF


C     Call D1CQDF to define a Character Sequence

      INSTR = ';'

      CALL D1CQDF (CMEM, IMEM, DMEM, NSTR, NSTRCH, INSTR, ICQ, IER)
      IF (IER .NE. 0) THEN
         GOTO 9000
      ENDIF


C     Fill in standard data

      IA(1) = IGI
      IA(2) = ICQ
      IA(3) = NPAR
      IA(4) = 1

      JLO = 1

      IF (IGLOBL) THEN

C        Load the TYPCOD's for Global data

         CA = GLBLCA
         JHI = 23

         CALL D2STAC (CMEM, IMEM, DMEM, CA, JLO, JHI, IGE, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF

         DO 10, I = 5, 28
            IA(I) = GLBLIA(I-4)
   10    CONTINUE

         JHI = 28
      ELSE
         JHI = 4
      ENDIF

      CALL D2STAI (CMEM, IMEM, DMEM, IA, JLO, JHI, IGE, IERX)
      IF (IERX .NE. 0) THEN
         IER = -999
         GOTO 9000
      ENDIF
      
C     Initialize the DMEM free-space stack
      IF (NREAL .GT. 0) THEN
         DO 20 I = 1, NREAL-1 
            DELM = I+1
            CALL D2STED (CMEM, IMEM, DMEM, DELM, I, IGE, IERX)
            IF (IERX .NE. 0) THEN
               IER = -999
               GOTO 9000
            ENDIF
   20    CONTINUE
         DELM = 0
         CALL D2STED (CMEM, IMEM, DMEM, DELM, NREAL, IGE, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
      ELSE
         CALL D2STEI (CMEM, IMEM, DMEM, 0, 4, IGE, IERX)
         IF (IERX .NE. 0) THEN
            IER = -999
            GOTO 9000
         ENDIF
      ENDIF

 9000 CONTINUE

      IF (IER .NE. 0) THEN
C        An error occurred.  Try to erase any created entities.

         IF (IGE .NE. 0) THEN
C           Try to release the defined IGES entity
            CALL D2ERAS (CMEM, IMEM, DMEM, IGE, IERX)
            IGE = 0
         ENDIF

         IF (ICQ .NE. 0) THEN
C           Try to release the defined character sequence entity
            CALL D2ERAS (CMEM, IMEM, DMEM, ICQ, IERX)
            ICQ = 0
         ENDIF
      ENDIF

      RETURN
      END

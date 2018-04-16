      SUBROUTINE D2IGDI (CMEM, IMEM, DMEM, NSL, NDE, INIT, IGI, IER)

C     PURPOSE:
C        IGES Define Index.  Create an IGES Index entity.
C
C        Output pointer to the newly allocated Index entity in IGI.
C
C        IGES Index Entity:
C
C
C           CMEM  [1]         Parameter delimiter
C                 [2]         Record delimiter
C                 [3]         Compression flag
C                             'C' = compressed
C                             otherwise  = fixed
C                 [4..75]     Start section line 1
C                  .
C                  .
C                 [...NSL*72+3]  Start section line NSL
C
C
C           IMEM  [1]         Pointer to Special IGES Entity (GLOBAL)
C                 [2]         NDE = Number of Directory section lines
C                 [3]         Pointer to IGES DE 1
C                  .
C                  .
C                 [NDE/2+2]   Pointer to IGES DE (NDE/2)
C
C
C           DMEM  (not used)
C
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C
C     INPUT:
C        NSL      Number of Start Section lines allocated
C        NDE      Number of Directory Section lines allocated, which
C                 is twice the number of IGES Entities that can be
C                 placed in the index.
C        INIT     Global parameter nitialization flag.
C                 .FALSE. = Do not initialize
C                 .TRUE.  = Initialize Global parameters to library
C                           default values.
C
C     OUTPUT:
C        IGI      New Index entity ID
C        IER      Error flag.  If negative, the entity is not defined.
C                 -1  = Dynamic memory is corrupt or uninitialized.
C                 -2  = Insufficient space available for character data.
C                 -3  = Insufficient space available for integer data.
C                 -4  = Insufficient space available for double
C                       precision data.
C                 -5  = NSL < 0.
C                 -6  = NDE < 0.
C                 -7  = NDE not even.
C                 -999= Unexpected error in a called subroutine.
C
C     CALLS:
C        D1DEFE   Define entity
C        D1IGDE   Define IGES Entity (for GLOBAL data)
C        D2STAC   Store character array
C        D2STAI   Store integer array
C        D2CQST   Store character sequence string
C        D2STAD   Store double precision array
C        DTERR    error handler
C        DTERAS   Erase entity
C
C
C     HISTORY:
C        20May93  D. Parsons   Created.
C        29Jun93  D. Parsons   Added DMEM free-space stack.
C
C     ------

C     Long name alias:
C        ENTRY D2_IGES_DEFINE_INDEX (CMEM, IMEM, DMEM, NSL, NDE, INIT,
C       +    IGI, IER)

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      EXTERNAL D1MBAD, DTJCON
      LOGICAL  D1MBAD
      INTEGER  DTJCON

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

      INTEGER     NSL, NDE, IGI, JSEQ, LENSTR, IER, NEED
      INTEGER     IGE, ICQ, IERX, IDTYP, LENC, LENI, LEND
      INTEGER     NPAR, NSTR, NSTRCH, NREAL, JLO, JHI, IA(9)
      DOUBLE PRECISION DA(4)
      LOGICAL     INIT, INITIZ, IGLOBL
      CHARACTER*6 SUBNAM
      CHARACTER   CA*3, STR*13

      DATA SUBNAM /'D2IGDI'/

C     ------

      IER  = 0
      IERX = 0
      IGI  = 0
      IGE  = 0
      ICQ  = 0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9000
         ENDIF
      ENDIF

      IDTYP  = ENTGI
      LENC   = NSL*72 + 3
      LENI   = NDE/2  + 2
      LEND   = 0
      INITIZ = .TRUE.

C     Check the validity of input arguments

      IF (NSL .LT. 0) THEN
         IER = -5
         GOTO 9000
      ENDIF

      IF (NDE .LT. 0) THEN
         IER = -6
         GOTO 9000
      ENDIF

      IF (MOD(NDE,2) .NE. 0) THEN
          IER = -7
          GOTO 9000
      ENDIF


C     Call D1DEFE to define a new entity

      CALL D1DEFE (CMEM, IMEM, DMEM, IDTYP, LENC, LENI, LEND, INITIZ,
     +      IGI, IERX)
      IF (IERX .NE. 0) THEN
         IF       (IERX .EQ. -2) THEN
            IER = -2
         ELSE IF  (IERX .EQ. -3) THEN
            IER = -3
         ENDIF
         IGI = 0
         GOTO 9000
      ENDIF

C     Call D1IGDE to define a new IGES Entity for the GLOBAL data

      NPAR   = 23
      NSTR   = 10
      NSTRCH = 28
      NREAL  = 4
      IGLOBL = .TRUE.

      CALL D1IGDE (CMEM, IMEM, DMEM, NPAR, NSTR, NSTRCH, NREAL,
     +      IGLOBL, IGI, IGE, ICQ, NEED, IER)

C     IER codes are the same for "insufficient memory"

      IF (IER .NE. 0) GOTO 9000

C     Store the data into the entities
C     Write the entity pointers into the IMEM

      CA = ',;F'
      JLO = 1
      JHI = 3
      CALL D2STAC (CMEM, IMEM, DMEM, CA, JLO, JHI, IGI, IERX)
      IF (IERX .NE. 0) THEN
          IER = -999
          GOTO 9000
      ENDIF

      IA(1) = IGE
      IA(2) = NDE
      JLO = 1
      JHI = 2
      CALL D2STAI (CMEM, IMEM, DMEM, IA, JLO, JHI, IGI, IERX)
      IF (IERX .NE. 0) THEN
          IER = -999
          GOTO 9000
      ENDIF

C     Write the global defaults, if requested

      IF (INIT) THEN
C         CHARACTER DEFAULTS:
C           Default Product ID  (NULL)
          JSEQ   = 1
          LENSTR = 0
          CALL D2CQST (CMEM, IMEM, DMEM, ICQ, JSEQ, LENSTR, STR, IERX)
          IF (IERX .NE. 0) THEN
              IER = -999
              GOTO 9000
          ENDIF

C           Default File Name  (NULL)
          JSEQ   = 2
          LENSTR = 0
          CALL D2CQST (CMEM, IMEM, DMEM, ICQ, JSEQ, LENSTR, STR, IERX)
          IF (IERX .NE. 0) THEN
              IER = -999
              GOTO 9000
          ENDIF

C           Default System ID  (NULL)
          JSEQ   = 3
          LENSTR = 0
          CALL D2CQST (CMEM, IMEM, DMEM, ICQ, JSEQ, LENSTR, STR, IERX)
          IF (IERX .NE. 0) THEN
              IER = -999
              GOTO 9000
          ENDIF

C           Default Preprocessor Version  (NULL)
          JSEQ   = 4
          LENSTR = 0
          CALL D2CQST (CMEM, IMEM, DMEM, ICQ, JSEQ, LENSTR, STR, IERX)
          IF (IERX .NE. 0) THEN
              IER = -999
              GOTO 9000
          ENDIF

C           Default Product IGEntification  (NULL)
          JSEQ   = 5
          LENSTR = 0
          CALL D2CQST (CMEM, IMEM, DMEM, ICQ, JSEQ, LENSTR, STR, IERX)
          IF (IERX .NE. 0) THEN
              IER = -999
              GOTO 9000
          ENDIF

C           Default Units  (INCHES)
          JSEQ   = 6
          LENSTR = 2
          STR = 'IN'
          CALL D2CQST (CMEM, IMEM, DMEM, ICQ, JSEQ, LENSTR, STR, IERX)
          IF (IERX .NE. 0) THEN
              IER = -999
              GOTO 9000
          ENDIF

C           Default Date & Time of Generation  (NULL)
          JSEQ   = 7
          LENSTR = 13
          STR = '000000.000000'
          CALL D2CQST (CMEM, IMEM, DMEM, ICQ, JSEQ, LENSTR, STR, IERX)
          IF (IERX .NE. 0) THEN
              IER = -999
              GOTO 9000
          ENDIF

C           Default Author  (NULL)
          JSEQ   = 8
          LENSTR = 0
          CALL D2CQST (CMEM, IMEM, DMEM, ICQ, JSEQ, LENSTR, STR, IERX)
          IF (IERX .NE. 0) THEN
              IER = -999
              GOTO 9000
          ENDIF

C           Default Organization  (NULL)
          JSEQ   = 9
          LENSTR = 0
          CALL D2CQST (CMEM, IMEM, DMEM, ICQ, JSEQ, LENSTR, STR, IERX)
          IF (IERX .NE. 0) THEN
              IER = -999
              GOTO 9000
          ENDIF

C           Default Date & Time of Modification  (NULL)
          JSEQ   = 10
          LENSTR = 13
          STR = '000000.000000'
          CALL D2CQST (CMEM, IMEM, DMEM, ICQ, JSEQ, LENSTR, STR, IERX)
          IF (IERX .NE. 0) THEN
              IER = -999
              GOTO 9000
          ENDIF


C         INTEGER DEFAULTS:
C           Default Number of Binary Bits for Integer
          IA(1) = 32
C           Default Maximum Power of Ten for Real
          IA(2) = 0
C         Default Number of Significant Digits for Real
          IA(3) = DTJCON(10)
C           Default Maximum Power of Ten for Double Precision
          IA(4) = 0
C           Default Number of Significant Digits for Double Precision
          IA(5) = DTJCON(11)

          JLO = 10
          JHI = 14
          CALL D2STAI (CMEM, IMEM, DMEM, IA(1), JLO, JHI, IGE, IERX)
          IF (IERX .NE. 0) THEN
              IER = -999
              GOTO 9000
          ENDIF

C           Default Unit Flag (1 = INCHES)
          IA(6) = 1

          JLO = 17
          JHI = 17
          CALL D2STAI (CMEM, IMEM, DMEM, IA(6), JLO, JHI, IGE, IERX)
          IF (IERX .NE. 0) THEN
              IER = -999
              GOTO 9000
          ENDIF

C           Default Maximum Number of Line Weight Gradations
          IA(7) = 1

          JLO = 19
          JHI = 19
          CALL D2STAI (CMEM, IMEM, DMEM, IA(7), JLO, JHI, IGE, IERX)
          IF (IERX .NE. 0) THEN
              IER = -999
              GOTO 9000
          ENDIF

C           Default Specification Version (9 = Ver. 5.1)
          IA(8) = 9
C           Default Drafting Standard
          IA(9) = 0

          JLO = 26
          JHI = 27
          CALL D2STAI (CMEM, IMEM, DMEM, IA(8), JLO, JHI, IGE, IERX)
          IF (IERX .NE. 0) THEN
              IER = -999
              GOTO 9000
          ENDIF

C         DOUBLE PRECISION DEFAULTS:
C           Model Space Scale
          DA(1) = 1.0D0
C           Width of Maximum Line Weight
          DA(2) = 0.0D0
C           Minimum User-Intended Resolution
          DA(3) = 0.0D0
C           Approximate Maximum Coordinate
          DA(4) = 0.0D0

          JLO = 1
          JHI = 4
C           Zero the free-space pointer          
          CALL D2STEI (CMEM, IMEM, DMEM, 0, 4, IGE, IERX)
          IF (IERX .NE. 0) THEN
             IER = -999
              GOTO 9000
          ENDIF
          CALL D2STAD (CMEM, IMEM, DMEM, DA(1), JLO, JHI, IGE, IERX)
          IF (IERX .NE. 0) THEN
              IER = -999
              GOTO 9000
          ENDIF

      ENDIF

C     Done

      GOTO 9999

C     Error Handling

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
      ENDIF

C     Try to erase the IGES Index Entity if started

      IF (IGI .NE. 0) THEN
         CALL D2ERAS (CMEM, IMEM, DMEM, IGI, IERX)
      ENDIF

C     Try to erase the IGES (Global) Entity if started

      IF (IGE .NE. 0) THEN
         CALL D2ERAS (CMEM, IMEM, DMEM, IGE, IERX)
      ENDIF

C     Try to erase the Character Sequence Entity if started

      IF (ICQ .NE. 0) THEN
         CALL D2ERAS (CMEM, IMEM, DMEM, ICQ, IERX)
      ENDIF

      IGI = 0

 9999 CONTINUE

      RETURN
      END

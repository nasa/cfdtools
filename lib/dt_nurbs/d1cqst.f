      SUBROUTINE D1CQST (CMEM, IMEM, DMEM, IHDR, JSEQ, LENSTR, STR,
     +                   NEED, IER)

C     Shadow routine for D2CQST
C     Store a string into a Character Sequence (CQ) entity
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IHDR    Header Index pointer to a Character String
C             Sequence entity
C     JSEQ    Ordinal number of string to store to ICQ
C     LENSTR  Length of STR.  LENSTR = 0 to remove the field.
C     STR     String to store in JSEQth position in ICQ
C   OUTPUT:
C     NEED    If IER = -13 or -14, NEED is the needed amount
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             -8  = JSEQ is less than or equal to zero
C             -9  = LENSTR is less than zero
C             -10 = Entity ICQ is internally inconsistent
C             -11 = Character space is locked
C             -12 = Integer space is locked
C             -13 = Insufficient Character space for needed expansion
C             -14 = Insufficient Integer space for needed expansion
C
C   CALLS:    D1FEEI  Fetch Integer Element
C             D1STEI  Store Integer Element
C             D1DEFX  Expand Entity Definition
C             D1STAC  Store Character Array
C             D1STEC  Store Character Element
C
C   HISTORY:
C     14Jun93 D. Parsons   Created
C     28Jun93 D. Parsons   Add support for deleted & undefined
C                          strings.  These will be designated with
C                          -(endpos+1) in their associated IMEM
C                          positions.
C*****
C

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER IHDR, JSEQ, LENSTR, NEED, IER
      CHARACTER STR*(*)

      INTEGER I, NCFREE, JLO, JHI, IXPTR, NMORE, NLESS
      INTEGER LENI, LENC, MXSEQ, CURLEN, LUSED, LBE4
      INTEGER STPOS, ENDPOS
      INTEGER MOREC, MOREI, MORED, IERX
      LOGICAL INITIZ

C****
      IER   = 0
      NEED = 0

      IF (JSEQ .LE. 0) THEN
          IER = -8
          GOTO 9900
      ENDIF

      IF (LENSTR .LT. 0) THEN
          IER = -9
          GOTO 9900
      ENDIF

      MOREC = 0
      MOREI = 0
      MORED = 0

      MXSEQ = 0

      STPOS  = 1
      ENDPOS = 1
      LBE4   = 0

C     Get LENI and LENC

      LENI = ABS(IMEM(IHDR+5))
      LENC = ABS(IMEM(IHDR+4))

C     Get the current number of CQ fields

      IF (LENI .GE. 1) THEN
          CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 1, MXSEQ, IERX)
      ENDIF
      MXSEQ = MAX(0,MXSEQ)
C
      IF (JSEQ .GT. MXSEQ) THEN
C         We need to add more CQ fields.
          IF (JSEQ+1 .GT. LENI) THEN
C             We need to add more integer space to do that
              MOREI = JSEQ+1 - LENI
          ENDIF
          CURLEN = 0
          IF (MXSEQ .GT. 0) THEN
              CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, MXSEQ+1, LBE4, IERX)
              IF (IERX .NE. 0) THEN
                  IER = -10
                  GOTO 9900
              ENDIF
              IF (LBE4 .LT. 0)   LBE4   = -(LBE4+1)
          ENDIF
      ELSE
C         Get the current length of this field
          CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, JSEQ+1, ENDPOS, IERX)
          IF (IERX .NE. 0) THEN
              IER = -10
              GOTO 9900
          ENDIF
          IF (JSEQ .GT. 1) THEN
              CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, JSEQ, LBE4, IERX)
              IF (IERX .NE. 0) THEN
                  IER = -10
                   GOTO 9900
              ENDIF
          ENDIF
          IF (ENDPOS .LT. 0) ENDPOS = -(ENDPOS+1)
          IF (LBE4 .LT. 0)   LBE4   = -(LBE4  +1)
          CURLEN = ENDPOS-LBE4
      ENDIF

C     Get the current used length of all fields

      IF (MXSEQ .GT. 0) THEN
          CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, MXSEQ+1, LUSED, IERX)
          IF (IERX .NE. 0) THEN
              IER = -10
              GOTO 9900
          ENDIF
          IF (LUSED .LT. 0) LUSED = -(LUSED+1)
      ELSE
          LUSED = 0
      ENDIF

      IF (LENSTR .GT. CURLEN) THEN
C         We need to expand this field and maybe the entire entity

          NCFREE = LENC-LUSED
          IF (NCFREE .LT. 0) THEN
              IER = -10
              GOTO 9900
          ENDIF

          NMORE = LENSTR - CURLEN
          MOREC = MAX(NMORE-NCFREE,0)
      ELSE
          NLESS = CURLEN - LENSTR
      ENDIF

      IF ((MOREC .GT. 0) .OR. (MOREI .GT. 0)) THEN
C         Expand the entity
          INITIZ = .TRUE.
          CALL D1DEFX (CMEM, IMEM, DMEM, IHDR, MOREC, MOREI, MORED,
     +        INITIZ, NEED, IERX)
          IF (IERX .NE. 0) THEN
              IF ((IERX .EQ. -8) .OR. (IERX .EQ. -9)) THEN
C                 Locked Character or Integer space
                  IER = IERX-3
              ELSE IF ((IERX .EQ. -11) .OR. (IERX .EQ. -12)) THEN
C                 Insufficient Character or Integer space
                  IER = IERX-2
              ELSE
                  IER = -999
              ENDIF
              GOTO 9900
          ENDIF
      ENDIF

      IF (LENSTR .GT. CURLEN) THEN
C         We need to shift the characters right to make room for
C             the expanded field

          IXPTR = IMEM(IHDR+1)
          DO 20 I = IXPTR+LUSED+NMORE-1, IXPTR+LBE4+CURLEN+NMORE, -1
              CMEM(I:I) = CMEM(I-NMORE:I-NMORE)
   20     CONTINUE

C         Update the ending position of this string in CQ
          CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, JSEQ+1, ENDPOS, IERX)
          IF (IERX .NE. 0) THEN
              IER = -10
              GOTO 9900
          ENDIF
          IF (ENDPOS .LT. 0) ENDPOS = -(ENDPOS+1)
          CALL D1STEI (CMEM, IMEM, DMEM, ENDPOS+NMORE, JSEQ+1, IHDR, 
     +        IERX)
          IF (IERX .NE. 0) THEN
              IER = -10
              GOTO 9900
          ENDIF
          
C         Update the subsequent ending positions in CQ
          DO 30 I = JSEQ+2, MXSEQ+1
              CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, I, ENDPOS, IERX)
              IF (IERX .NE. 0) THEN
                      IER = -10
              GOTO 9900
              ENDIF
              IF (ENDPOS .LT. 0) THEN
                 CALL D1STEI (CMEM, IMEM, DMEM, ENDPOS-NMORE, I, IHDR,
     +                 IERX)
              ELSE
                 CALL D1STEI (CMEM, IMEM, DMEM, ENDPOS+NMORE, I, IHDR,
     +                 IERX)
              ENDIF
              IF (IERX .NE. 0) THEN
                  IER = -10
                  GOTO 9900
              ENDIF
   30     CONTINUE
      ELSE IF (LENSTR .LT. CURLEN) THEN
C        We need to shift the characters left to remove excess space
C           from the field.
         
          IXPTR = IMEM(IHDR+1)
          DO 40 I = IXPTR+LBE4+LENSTR, IXPTR+LUSED-NLESS-1
              CMEM(I:I) = CMEM(I+NLESS:I+NLESS)
   40     CONTINUE
   
C        Asterisk out the end of the line

          DO 45 I = IXPTR+LUSED-NLESS, IXPTR+LUSED-1
              CMEM(I:I) = '*'
   45     CONTINUE

C         Update the ending position of this string in CQ
          CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, JSEQ+1, ENDPOS, IERX)
          IF (IERX .NE. 0) THEN
              IER = -10
              GOTO 9900
          ENDIF
          IF (ENDPOS .LT. 0) ENDPOS = -(ENDPOS+1)
          CALL D1STEI (CMEM, IMEM, DMEM, ENDPOS-NLESS, JSEQ+1, IHDR, 
     +        IERX)
          IF (IERX .NE. 0) THEN
              IER = -10
              GOTO 9900
          ENDIF
          
C         Update the ending positions in CQ
          DO 50 I = JSEQ+2, MXSEQ+1
              CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, I, ENDPOS, IERX)
              IF (IERX .NE. 0) THEN
                      IER = -10
              GOTO 9900
              ENDIF
              IF (ENDPOS .LT. 0) THEN
                 CALL D1STEI (CMEM, IMEM, DMEM, ENDPOS+NLESS, I, IHDR,
     +                 IERX)
              ELSE
                 CALL D1STEI (CMEM, IMEM, DMEM, ENDPOS-NLESS, I, IHDR,
     +                 IERX)
              ENDIF
              IF (IERX .NE. 0) THEN
                  IER = -10
                  GOTO 9900
              ENDIF
   50     CONTINUE
      ENDIF

C     Update the "number of strings", if necessary
      IF (JSEQ .GT. MXSEQ) THEN
          CALL D1STEI (CMEM, IMEM, DMEM, JSEQ, 1, IHDR, IERX)
          IF (IERX .NE. 0) THEN
              IER = -10
              GOTO 9900
          ENDIF

C         Mark subsequent strings as uninitialized
          DO 60 I = MAX(MXSEQ+2,3), JSEQ
              CALL D1STEI (CMEM, IMEM, DMEM, -(LUSED+1), I, IHDR, IERX)
              IF (IERX .NE. 0) THEN
                  IER = -10
                  GOTO 9900
              ENDIF
   60     CONTINUE

C         Store the new ending position in CQ
          CALL D1STEI (CMEM, IMEM, DMEM, LUSED+LENSTR, JSEQ+1, IHDR,
     +                 IERX)
          IF (IERX .NE. 0) THEN
              IER = -10
              GOTO 9900
          ENDIF
      ENDIF

C     Store the string

      JLO = LBE4+1
      JHI = JLO+LENSTR-1

C     JHI < JLO if either LENSTR and CURLEN = 0. That is, no space
C         is to be assigned for this field.
      IF (JHI .GE. JLO) THEN
          CALL D1STAC (CMEM, IMEM, DMEM, STR, JLO, JHI, IHDR, IERX)
          IF (IERX .NE. 0) THEN
              IER = -10
              GOTO 9900
          ENDIF
      ENDIF

 9900 CONTINUE

      RETURN
      END


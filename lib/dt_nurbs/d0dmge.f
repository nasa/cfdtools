      SUBROUTINE D0DMGE (CMEM, IMEM, DMEM, IDE, LUNIT, IER)

C     IGES Entity Check


      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER           IDE, LUNIT, IER, IERX, I1, I2
      INTEGER           LENC, LENI, LEND, ITYPE, INDEX
      INTEGER           NPAR, NSTR, NSTRCH, NREAL, IVAL, IREAL, ICHAR
      INTEGER           IGI, IHDR, ITYP, ICQ, IPAR, INEXT, I, THIS
      DOUBLE PRECISION  NEXT
      LOGICAL           GLOBLF


C                          IGES Entity Data
      INTEGER     ENTGE
      PARAMETER  (ENTGE = 238)

C                          IGES Index
      INTEGER     ENTGI
      PARAMETER  (ENTGI = 237)

C                          Character Sequence Entity
      INTEGER     ENTCQ
      PARAMETER  (ENTCQ = 249)

      CHARACTER*(23)    TYPCOD, XTYPE
      DATA XTYPE        /'CCCCIIIIICEICIECEECCIIC'/


      IER   = 0
      IERX  = 0
      I1   = IDE/256
      I2   = MOD(ABS(IDE),256)

      CALL D2IGLE (CMEM, IMEM, DMEM, IDE, NPAR, NSTR, NSTRCH,
     +                   NREAL, IERX)
      IF (IERX .NE. 0) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 100) I1, I2, IERX
            IER = -999
         ELSE
            IER = -1
            GOTO 9999
         ENDIF
      ENDIF

C     Is this a GLOBAL IGES entity, or not?
C        If entity type = -11 then this is a GLOBAL IGES entity

      CALL D2FEEI (CMEM, IMEM, DMEM, IDE, 5, ITYPE, IERX)
      IF (ITYPE .EQ. -11) THEN
         GLOBLF = .TRUE.
      ELSE
         GLOBLF = .FALSE.
      ENDIF
      
      IF (GLOBLF) THEN
         IF (NPAR .NE. 23) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 120) I1, I2, NPAR
               IER = -999
            ELSE
               IER = -3
               GOTO 9999
            ENDIF
         ENDIF
         IF (NSTR .NE. 10) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 130) I1, I2, NSTR
               IER = -999
            ELSE
               IER = -4
               GOTO 9999
            ENDIF
         ENDIF
C         IF (NSTRCH .NE. ??) THEN
C            IF (LUNIT .GE. 0) THEN
C               WRITE (LUNIT, 140) I1, I2, NSTRCH
C               IER = -999
C            ELSE
C               IER = -5
C               GOTO 9999
C            ENDIF
C         ENDIF
         IF (NREAL .NE. 4) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 150) I1, I2, NREAL
               IER = -999
            ELSE
               IER = -6
               GOTO 9999
            ENDIF
         ENDIF
      ENDIF

      CALL D2FELC (CMEM, IMEM, DMEM, IDE, LENC, IERX)
      IF (GLOBLF .AND. (LENC .NE. 23)) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 110) I1, I2, LENC, 23
            IER = -999
         ELSE
            IER = -2
            GOTO 9999
         ENDIF
      ELSE IF (.NOT. GLOBLF .AND. (LENC .NE. NPAR+24)) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 110) I1, I2, LENC, NPAR+24
            IER = -999
         ELSE
            IER = -2
            GOTO 9999
         ENDIF
      ENDIF


      CALL D2FELI (CMEM, IMEM, DMEM, IDE, LENI, IERX)
      IF (GLOBLF .AND. (LENI .NE. 28)) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 160) I1, I2, LENI
            IER = -999
         ELSE
            IER = -7
            GOTO 9999
         ENDIF
      ELSE IF (.NOT. GLOBLF .AND. (LENI .NE. NPAR+19)) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 170) I1, I2, LENI, NPAR+19
            IER = -999
         ELSE
            IER = -8
            GOTO 9999
         ENDIF
      ENDIF

      CALL D2FELD (CMEM, IMEM, DMEM, IDE, LEND, IERX)
      IF (LEND .NE. NREAL) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 180) I1, I2, LEND, NREAL
            IER = -999
         ELSE
            IER = -9
            GOTO 9999
         ENDIF
      ENDIF

C     Check pointers
      CALL D2FEEI (CMEM, IMEM, DMEM, IDE, 1, IGI, IERX)
      CALL D0PTR (CMEM, IMEM, DMEM, IGI, ENTGI, 0, IHDR,
     +    ITYP, IERX)
      IF (IERX .NE. 0) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 190) I1, I2, IGI, IERX
            IER = -999
         ELSE
            IER = -10
            GOTO 9999
         ENDIF
      ENDIF

      CALL D2FEEI (CMEM, IMEM, DMEM, IDE, 2, ICQ, IERX)
      CALL D0PTR (CMEM, IMEM, DMEM, ICQ, ENTCQ, 0, IHDR,
     +    ITYP, IERX)
      IF (IERX .NE. 0) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 200) I1, I2, ICQ, IERX
            IER = -999
         ELSE
            IER = -11
            GOTO 9999
         ENDIF
      ENDIF

C     Check TYPCODs
      IF (GLOBLF) THEN
         CALL D2FEAC (CMEM, IMEM, DMEM, IDE, 1, 23, TYPCOD, IERX)
         IF (TYPCOD .NE. XTYPE) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 210) I1, I2, TYPCOD, XTYPE
               IER = -999
            ELSE
               IER = -12
               GOTO 9999
            ENDIF
         ENDIF
      ELSE
         DO 10 IPAR = 1, NPAR
            CALL D2FEEC (CMEM, IMEM, DMEM, IDE, IPAR+24, TYPCOD(1:1),
     +            IERX)
            CALL D2FEEI (CMEM, IMEM, DMEM, IDE, IPAR+19, IVAL, IERX)
            IF (TYPCOD(1:1) .EQ. 'C') THEN
               IF (IVAL .GT. LENC) THEN
                  IF (LUNIT .GE. 0) THEN
                     WRITE (LUNIT, 230) I1, I2, TYPCOD(1:1), IVAL, LENC
                     IER = -999
                  ELSE
                     IER = -14
                     GOTO 9999
                  ENDIF
               ENDIF
            ELSE IF (TYPCOD(1:1) .EQ. 'I') THEN
C              No test for plain integer
            ELSE IF (TYPCOD(1:1) .EQ. 'P' .OR.
     +               TYPCOD(1:1) .EQ. '*') THEN
               CALL D0PTR (CMEM, IMEM, DMEM, IVAL, ENTGE, 0, IHDR,
     +            ITYP, IERX)
               IF (IERX .NE. 0) THEN
                  IF (LUNIT .GE. 0) THEN
                     WRITE (LUNIT, 240) I1, I2, TYPCOD(1:1), IVAL, IERX
                     IER = -999
                  ELSE
                     IER = -15
                     GOTO 9999
                  ENDIF
               ENDIF
            ELSE IF (TYPCOD(1:1) .EQ. 'L') THEN
               IF (IVAL .NE. 1 .AND. IVAL .NE. 0) THEN
                  IF (LUNIT .GE. 0) THEN
                     WRITE (LUNIT, 250) I1, I2, TYPCOD(1:1), IVAL
                     IER = -999
                  ELSE
                     IER = -16
                     GOTO 9999
                  ENDIF
               ENDIF
            ELSE IF (TYPCOD(1:1) .EQ. 'A') THEN
C              No test for ambiguous integer
            ELSE IF (TYPCOD(1:1) .EQ. 'E' .OR.
     +               TYPCOD(1:1) .EQ. 'D') THEN
               IF (IVAL .GT. LEND) THEN
                  IF (LUNIT .GE. 0) THEN
                     WRITE (LUNIT, 260) I1, I2, TYPCOD(1:1), IVAL, LEND
                     IER = -999
                  ELSE
                     IER = -17
                     GOTO 9999
                  ENDIF
               ENDIF
            ELSE IF (TYPCOD(1:1) .EQ. 'T' .OR.
     +               TYPCOD(1:1) .EQ. 'R') THEN
               CALL D0PTR (CMEM, IMEM, DMEM, IVAL, ENTCQ, 0, IHDR,
     +            ITYP, IERX)
               IF (IERX .NE. 0) THEN
                  IF (LUNIT .GE. 0) THEN
                     WRITE (LUNIT, 270) I1, I2, TYPCOD(1:1), IVAL, IERX
                     IER = -999
                  ELSE
                     IER = -18
                     GOTO 9999
                  ENDIF
               ENDIF
            ELSE IF (TYPCOD(1:1) .EQ. 'O' .OR.
     +               TYPCOD(1:1) .EQ. ' ') THEN
C              No test for omitted value
            ELSE
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 220) I1, I2, IPAR, TYPCOD(1:1)
                  IER = -999
               ELSE
                  IER = -13
                  GOTO 9999
               ENDIF
            ENDIF
   10    CONTINUE
      ENDIF

C     Check the DMEM stack
      CALL D2FEEI (CMEM, IMEM, DMEM, IDE, 4, INEXT, IERX)
      NEXT = INEXT
      DO 20 I = 1, NREAL
         IF (INT(NEXT) .EQ. 0) GOTO 30
         IF (INT(NEXT) .GT. LEND) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 280) I1, I2, I, NEXT, LEND
               IER = -999
            ELSE
               IER = -19
               GOTO 9999
            ENDIF
         ENDIF
         THIS = NEXT
         CALL D2FEED (CMEM, IMEM, DMEM, IDE, THIS, NEXT, IERX)
   20 CONTINUE
         IF (INT(NEXT) .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 290) I1, I2, NEXT
               IER = -999
            ELSE
               IER = -20
               GOTO 9999
            ENDIF
         ENDIF
   30 CONTINUE

C     Count TYPCODs
      IREAL = 0
      ICHAR = 0        
      DO 40 I = 1, NPAR
         IF (GLOBLF) THEN
            INDEX = I
         ELSE
            INDEX = I+24
         ENDIF
         CALL D2FEEC (CMEM, IMEM, DMEM, IDE, INDEX, TYPCOD(1:1), IERX)
         IF (TYPCOD(1:1) .EQ. 'C') THEN
            ICHAR = ICHAR+1
         ELSE IF ((TYPCOD(1:1) .EQ. 'E') .OR.
     +            (TYPCOD(1:1) .EQ. 'D')) THEN
            IREAL = IREAL+1
         ENDIF
   40 CONTINUE

      IF (IREAL .GT. NREAL) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 300) I1, I2, IREAL, NREAL
            IER = -999
         ELSE
            IER = -21
            GOTO 9999
         ENDIF
      ENDIF

C     EOP and EOR are considered Global parameters!!
      IF (GLOBLF) NSTR = NSTR+2
      IF (ICHAR .GT. NSTR) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 310) I1, I2, ICHAR, NSTR
            IER = -999
         ELSE
            IER = -22
            GOTO 9999
         ENDIF
      ENDIF


 9999 CONTINUE

      RETURN

  100 FORMAT (' IDE: ',I3,'#',I3.3,
     +        ' GE: Unable to get parameter lengths; IER = ',I5,'.')
  110 FORMAT (' IDE: ',I3,'#',I3.3,' GE: (Global) LENC = ',I5,
     +        ' != ',I5,'.')
  120 FORMAT (' IDE: ',I3,'#',I3.3,' GE: (Global) NPAR = ',I5,
     +        ' != 23.')
  130 FORMAT (' IDE: ',I3,'#',I3.3,' GE: (Global) NSTR = ',I5,
     +        ' != 10.')
  140 FORMAT (' IDE: ',I3,'#',I3.3,' GE: (Global) NSTRCH = ',I5,
     +        ' != ??.')
  150 FORMAT (' IDE: ',I3,'#',I3.3,' GE: (Global) NREAL = ',I5,
     +        ' != 4.')
  160 FORMAT (' IDE: ',I3,'#',I3.3,' GE: (Global) LENI = ',I5,
     +        ' != 28.')
  170 FORMAT (' IDE: ',I3,'#',I3.3,' GE: LENI = ',I5,
     +        ' != NPAR+18= ',I5,'.')
  180 FORMAT (' IDE: ',I3,'#',I3.3,' GE: LEND = ',I5,
     +        ' != NREAL = ',I5,'.')
  190 FORMAT (' IDE: ',I3,'#',I3.3,' GE: IGI = ',I5,
     +        'Not a valid IGES Index; D0PTR returned ',I5,'.')
  200 FORMAT (' IDE: ',I3,'#',I3.3,' GE: ICQ = ',I5,
     +        'Not a valid Character Sequence; D0PTR returned ',I5,'.')
  210 FORMAT (' IDE: ',I3,'#',I3.3,' GE: (Global) TYPCODs = ',A23,
     +        ' != Expected = ',A23,'.')
  220 FORMAT (' IDE: ',I3,'#',I3.3,' GE: IPAR = ',I5,' TYPCOD = ',A1,
     +        ' Invalid TYPCOD.')
  230 FORMAT (' IDE: ',I3,'#',I3.3,' GE: IPAR = ',I5,' TYPCOD = ',A1,
     +        ' IVAL = ',I5,' > LENC = ',I5,'.')
  240 FORMAT (' IDE: ',I3,'#',I3.3,' GE: IPAR = ',I5,' TYPCOD = ',A1,
     +        ' IVAL = ',I5,
     +        'Not a valid IGES Entity; D0PTR returned ',I5,'.')
  250 FORMAT (' IDE: ',I3,'#',I3.3,' GE: IPAR = ',I5,' TYPCOD = ',A1,
     +        ' IVAL = ',I5,'Not a valid "LOGICAL" value.')
  260 FORMAT (' IDE: ',I3,'#',I3.3,' GE: IPAR = ',I5,' TYPCOD = ',A1,
     +        ' IVAL = ',I5,' > LEND = ',I5,'.')
  270 FORMAT (' IDE: ',I3,'#',I3.3,' GE: IPAR = ',I5,' TYPCOD = ',A1,
     +        ' IVAL = ',I5,
     +        'Not a valid Character Sequence; D0PTR returned ',I5,'.')
  280 FORMAT (' IDE: ',I3,'#',I3.3,' GE: DMEM stack element ',I5,' = ',
     +        F5.1,' > LEND = ',I5,'.')
  290 FORMAT (' IDE: ',I3,'#',I3.3,' GE: Last DMEM stack element ',
     +        ' = ',F5.1,' needs to be zero.')
  300 FORMAT (' IDE: ',I3,'#',I3.3,' GE: Number of real TYPCODs = ',I5,
     +        ' Exceeds NREAL = ',I5,'.')
  310 FORMAT (' IDE: ',I3,'#',I3.3,' GE: Number of char TYPCODs = ',I5,
     +        ' Exceeds NSTR = ',I5,'.')

      END

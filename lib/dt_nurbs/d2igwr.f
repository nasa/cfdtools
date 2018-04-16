      SUBROUTINE D2IGWR (CMEM, IMEM, DMEM, IGI, LVLO, FILNAM, LUNIT,
     +                   LVLE, IER)

C     Read an IGES file into the DTRC Dynamic Memory.
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C     INPUT:
C        IGI      Dynamic memory pointer to the new IGES Index entity.
C        LVLO     The IGES entity types to accept as input:
C                 1  NASA-IGES-NURBS only
C                 2  NASA-IGES
C                 3  Unrestricted
C        FILNAM   The name of the output IGES file to write.
C        LUNIT    The I/O unit to use when working on the IGES file.
C                 If LUNIT < 0 on input, a line will be displayed on
C                 Unit 6, giving the process status for each entity
C                 processed.
C                 IF LUNIT > 0, no status will be displayed.
C     OUTPUT:
C        LVLE     The Level of IGES entities encountered.  See LVLO
C                 for meaning.
C        IER      Returned error value.  Zero implies no errors.
C                 Otherwise,
C                 +1 LVLE > LVLO.  Some non-NASA-IGES or non-NASA-IGES-
C                    NURBS entities were written to the file as Null
C                    entities.  If other NASA-IGES entities point to 
C                    these entities, the output file may be invalid for 
C                    some applications.
C                 -1 Defective structure or data in Flag or Start
C                    section information.
C                 -2 Defective structure or data in Global section
C                    information.
C                 -n (n an odd number, not 1)
C                    Defective structure or data in Directory Entry of
C                    IGES Entity with DE pointer (n-2).
C                 -n (n an even number, not 2)
C                    Defective structure or data in Parameter Entry of
C                    IGES Entity with DE ptr (n-3).
C                 -9999  Unlikely I/O error occurred during copying of
C                    scratch file to final output.  Output file may
C                    exist, but is incomplete.
C
C   CALLS:
C     DTJCON   Integer Constants
C     D0IGDO   Decided whether this is an entity to "DO"
C     D0IGWP   Write an IGES parameter
C     D2IGLI   Get Index lengths
C     D2IGLE   Get Entity lengths
C     D2IGFT   Fetch parameter type
C     D2IGFC   Fetch Character parameter
C     D2IGFD   Fetch Double Precision parameter
C     D2IGFI   Fetch Integer parameter
C     DTERR    Error Handling
C
C   HISTORY:
C     9/5/93      D. Parsons   Created
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_WRITE (CMEM, IMEM, DMEM, IGI, LVLO, FILNAM, LUNIT,
C     +                   LVLE, IER)

      EXTERNAL          D0IGDO, DTJCON
      LOGICAL           D0IGDO
      INTEGER           DTJCON

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      CHARACTER*(80) LINE, FILNAM
      INTEGER        LUNIT, LVLO, IGI, LVLE, IER

      CHARACTER*(72) CVAL
      CHARACTER      EOP, EOR, TYPE, SECT, TYPCOD
      INTEGER        ALUNIT, DLUNIT, IERX, I, CLEN, SEQ, DSEQ
      INTEGER        NSL, NDE, NPAR
      INTEGER        ICURLN, IPCURP, IOCHKA, IOCHKD, MSGUN
      INTEGER        IDLINE, IPLINE, PSEQ, PSSTRT, IDE, IGE
      INTEGER        NSTR, NSTRCH, NREAL
      INTEGER        IPAR, IVAL
      INTEGER        ITYPE, IFORM, IPSTRT, NGLBL
      LOGICAL        EORFLG, STATUS, EXFLAG
      DOUBLE PRECISION RVAL
      CHARACTER*8    CDVAL(3)
      INTEGER        IDVAL(15)

      CHARACTER*6    SUBNAM
      DATA SUBNAM   /'D2IGWR'/

C****

C     LUNIT is the unit number passed in by the user.  May be negative.
C        User must define this because on some systems that file has to
C        be attached to that unit in the Job Control Language.
C     ALUNIT is ABS(LUNIT)
C     DLUNIT is a scratch file for direct access writing.

      IF (LUNIT .LT. 0) THEN
         STATUS = .TRUE.
         ALUNIT = IABS(LUNIT)
         MSGUN  = DTJCON(6)
      ELSE
         STATUS = .FALSE.
         ALUNIT = LUNIT
      ENDIF

      IER  = 0
      IERX = 0
      LVLE = 0
      LINE = ' '

C     Find a free unit number for DLUNIT

      DLUNIT = 100
   10 CONTINUE
         INQUIRE (UNIT=DLUNIT, OPENED=EXFLAG)
         IF (EXFLAG) THEN
            DLUNIT = DLUNIT+1
            GOTO 10
         ENDIF

C     Open a scratch file for DIRECT ACCESS processing.

      OPEN (DLUNIT, ACCESS='DIRECT', FORM='FORMATTED', RECL=80,
     +      STATUS='SCRATCH', IOSTAT=IOCHKD)
      IF (IOCHKD .NE. 0) THEN
C        An error occurred in the OPEN statement
         IF (STATUS) WRITE (MSGUN, '(A,I6,/)')
     +      ' FATAL ERROR: Unable to open scratch file; IOCHKD = ',
     +      IOCHKD
         IER = -1
         GOTO 9010
      ENDIF

C     Open the output IGES file

      OPEN (ALUNIT, FILE=FILNAM, IOSTAT=IOCHKA)
      IF (IOCHKA .NE. 0) THEN
C        An error occurred in the OPEN statement
         IF (STATUS) WRITE (MSGUN, '(A,A,A,I6,/)')
     +      ' FATAL ERROR: Unable to open ',FILNAM,' IOCHKA = ',IOCHKA
         IER = -1
         GOTO 9010
      ENDIF

C     *** INITIALIZE COUNTERS AND POINTERS ***

      CALL D2IGLI (CMEM, IMEM, DMEM, IGI, NSL, NDE, IERX)
      IF (IERX .NE. 0) THEN
         IER = -1
         GOTO 9020
      ENDIF

      ICURLN = 1

C     Write the Start Section

      SECT = 'S'
      DO 100 I = 1, NSL
         CALL D2IGFC (CMEM, IMEM, DMEM, IGI, SECT, I, LINE, CLEN, IERX)
         IF (IERX .NE. 0) THEN
            IF (STATUS) WRITE (MSGUN, '(A,I3,/)')
     +         ' FATAL ERROR: Unable to fetch the Start section line ',
     +         I
            IER = -1
            GOTO 9020
         ENDIF
         WRITE (DLUNIT, '(A72,''S'',BZ,I7)',REC=ICURLN,IOSTAT=IOCHKD)
     +         LINE(1:CLEN), I
         IF (IOCHKD .NE. 0) THEN
            IF (STATUS) WRITE (MSGUN, '(A,I3,A,I6,/)')
     +         ' FATAL ERROR: Unable to write Start section line ',
     +         I,'; IOCHKD = ',IOCHKD
            IER = -1
            GOTO 9010
         ENDIF
         LINE = ' '
         ICURLN = ICURLN+1
  100 CONTINUE

C     Write the Global Section

      SECT   = 'G'
      IPCURP = 1
      SEQ    = 1
      DSEQ   = 0
      EOP    = ','
      EOR    = ';'
      LINE   = ' '
      EORFLG = .FALSE.

      DO 200 IPAR = 1, 25
         CALL D2IGFT (CMEM, IMEM, DMEM, IGI, SECT, IPAR, TYPCOD,
     +         IERX)
         IF (IERX .NE. 0) THEN
            IF (STATUS) WRITE (MSGUN, '(A,A,I2,/)')
     +         ' FATAL ERROR: Unable to fetch Global TYPCOD;',
     +         ' IPAR = ',IPAR
            IER = -2
            GOTO 9020
         ENDIF

         IF (TYPCOD .EQ. 'C') THEN
            CALL D2IGFC (CMEM, IMEM, DMEM, IGI, SECT, IPAR, CVAL, CLEN,
     +         IERX)
            TYPE = 'C'
         ELSE IF (TYPCOD .EQ. 'I') THEN
            CALL D2IGFI (CMEM, IMEM, DMEM, IGI, SECT, IPAR, IVAL, IERX)
            TYPE = 'I'
         ELSE IF (TYPCOD .EQ. 'E') THEN
            CALL D2IGFD (CMEM, IMEM, DMEM, IGI, SECT, IPAR, RVAL, IERX)
            TYPE = 'R'
         ELSE
            IF (STATUS) WRITE (MSGUN, '(A,A,A,I2,/)')
     +         ' FATAL ERROR: Invalid Global TYPCOD = ',TYPCOD,
     +         '; IPAR =', IPAR
            IER = -2
            GOTO 9010
         ENDIF
         IF (IERX .NE. 0) THEN
            IF (STATUS) WRITE (MSGUN, '(A,A,I2,/)')
     +         ' FATAL ERROR: Unable to fetch Global Parameter',
     +         ' IPAR =', IPAR
            IER = -2
            GOTO 9020
         ENDIF

         IF (IPAR .EQ. 1) EOP = CVAL(1:1)
         IF (IPAR .EQ. 2) EOR = CVAL(1:1)

         IF (IPAR .EQ. 25) EORFLG = .TRUE.
         CALL D0IGWP (LINE, DLUNIT, ICURLN, IPCURP, EOP, EOR,
     +                   EORFLG, TYPE, SECT, SEQ, DSEQ,
     +                   IVAL, RVAL, CVAL(1:MAX(1,CLEN)))
  200 CONTINUE

      NGLBL = ICURLN-NSL-1

C     Process Directory and Parameter Sections


C     Process Entities

      IDLINE = ICURLN
      IPLINE = ICURLN + NDE
      PSSTRT = IPLINE
      IPSTRT = IPLINE
      PSEQ   = 1
      DSEQ   = 1

      DO 500 IDE = 1, NDE, 2
C        Get the address of the IGES Entity
         SECT = 'I'
         CALL D2IGFI (CMEM, IMEM, DMEM, IGI, SECT, IDE, IGE, IERX)
         IF (IERX .NE. 0) THEN
            IF (STATUS) WRITE (MSGUN, '(A,/)')
     +         ' FATAL ERROR: Unable to fetch pointer to IGES entity.'
            IER = -IDE-2
            GOTO 9020
         ENDIF

         IF (IGE .NE. 0) THEN
            SECT = 'D'
            DO 310 I = 1, 19
               IF (I .LE. 10) THEN
                  CALL D2IGFI (CMEM, IMEM, DMEM, IGE, SECT, I,
     +               IDVAL(I), IERX)
               ELSE IF (I .EQ. 11) THEN
C                 (Skipped parameter)
               ELSE IF (I .LE. 15) THEN
                  CALL D2IGFI (CMEM, IMEM, DMEM, IGE, SECT, I,
     +               IDVAL(I-1), IERX)
               ELSE IF (I .LE. 18) THEN
                  CALL D2IGFC (CMEM, IMEM, DMEM, IGE, SECT, I,
     +               CDVAL(I-15), CLEN, IERX)
               ELSE
                  CALL D2IGFI (CMEM, IMEM, DMEM, IGE, SECT, I,
     +               IDVAL(I-4), IERX)
               ENDIF
               IF (IERX .NE. 0) THEN
                  IF (STATUS) WRITE (MSGUN, '(A,/)')
     +                  ' FATAL ERROR: Unable to fetch Directory entry'
                  IER = -IDE-2
                  GOTO 9020
               ENDIF
  310       CONTINUE

            ITYPE = IDVAL(1)
            IFORM = IDVAL(14)

         ELSE
            ITYPE = 0
            IFORM = 0
         ENDIF

         SECT = 'P'
         DSEQ = IDE
         IPCURP = 1
         LINE = ' '

         IF (ITYPE .EQ. 0) LVLE = MAX(LVLE,1)

         IF ((IGE .NE. 0) .AND. (ITYPE .NE. 0) .AND.
     +       (D0IGDO(ITYPE, IFORM, LVLO, LVLE))) THEN

            IF (STATUS) WRITE (MSGUN,'(A,I4,A,I3,A,I5)')
     +         ' D2IGWR: Processing Entity (IDE) Number ',IDE,
     +         '; Entity Type = ',ITYPE, ' Form = ',IFORM

            CALL D2IGLE (CMEM, IMEM, DMEM, IGE, NPAR, NSTR, NSTRCH,
     +                   NREAL, IERX)
            IF (IERX .NE. 0) THEN
               IF (STATUS) WRITE (MSGUN, '(A,A,/)')
     +            ' FATAL ERROR: Unable to fetch IGES entity',
     +            ' parameter lengths.'
               IER = -IDE-2
               GOTO 9020
            ENDIF

         ELSE

            IF (STATUS) WRITE (MSGUN,'(A,I4,A,I3,A,I5)')
     +         ' D2IGWR: Skipping Entity (IDE) Number ',IDE,
     +         '; Entity Type = ',ITYPE, ' Form = ',IFORM

            DO 320 I = 1, 15
               IDVAL(I) = 0
  320       CONTINUE
            DO 330 I = 1, 3
               CDVAL(I) = ' '
  330       CONTINUE

            ITYPE  = 0
            IFORM  = 0

            NPAR   = 0
            NSTR   = 0
            NSTRCH = 0
            NREAL  = 0

         ENDIF

C        Store the entity type

         IF (NPAR .EQ. 0) THEN
            EORFLG = .TRUE.
         ELSE
            EORFLG = .FALSE.
         ENDIF

         CALL D0IGWP (LINE, DLUNIT, IPLINE, IPCURP, EOP, EOR,
     +          EORFLG, 'I', SECT, PSEQ, DSEQ,
     +          ITYPE, RVAL, CVAL(1:MAX(1,CLEN)))

C        Process parameters 1 through NPAR

         DO 400 IPAR = 1, NPAR
            CALL D2IGFT (CMEM, IMEM, DMEM, IGE, SECT, IPAR,
     +         TYPCOD, IERX)
            IF (IERX .NE. 0) THEN
               IF (STATUS) WRITE (MSGUN, '(A,A,I3,/)')
     +            ' FATAL ERROR: Unable to fetch Parameter',
     +            ' TYPCOD; IPAR = ',IPAR
               IER = -IDE-3
               GOTO 9020
            ENDIF

            IF (TYPCOD .EQ. 'C') THEN
               CALL D2IGFC (CMEM, IMEM, DMEM, IGE, SECT, IPAR,
     +         CVAL, CLEN, IERX)
               TYPE = 'C'
            ELSE IF (
     +         (TYPCOD .EQ. 'I') .OR.
     +         (TYPCOD .EQ. 'P') .OR.
     +         (TYPCOD .EQ. 'L') .OR.
     +         (TYPCOD .EQ. 'T') .OR.
     +         (TYPCOD .EQ. 'R') .OR.
     +         (TYPCOD .EQ. '*') ) THEN
               CALL D2IGFI (CMEM, IMEM, DMEM, IGE, SECT, IPAR,
     +                      IVAL, IERX)
               TYPE = 'I'
            ELSE IF (
     +         (TYPCOD .EQ. 'E') .OR.
     +         (TYPCOD .EQ. 'D') ) THEN
               CALL D2IGFD (CMEM, IMEM, DMEM, IGE, SECT, IPAR,
     +            RVAL, IERX)
               TYPE = 'R'
            ELSE IF (TYPCOD .EQ. 'O') THEN
               TYPE = 'O'
            ELSE
               IF (STATUS) WRITE (MSGUN, '(A,A,A,I3,/)')
     +            ' FATAL ERROR: Invalid Global TYPCOD = ',TYPCOD,
     +            '; IPAR =', IPAR
               IER = -IDE-3
               GOTO 9010
            ENDIF
            IF (IERX .NE. 0) THEN
               IF (STATUS) WRITE (MSGUN, '(A,A,I3,/)')
     +            ' FATAL ERROR: Unable to fetch Parameter',
     +            ' IPAR =', IPAR
               IER = -IDE-3
               GOTO 9020
            ENDIF

            IF (IPAR .EQ. NPAR) THEN
               EORFLG = .TRUE.
            ELSE
               EORFLG = .FALSE.
            ENDIF

            CALL D0IGWP (LINE, DLUNIT, IPLINE, IPCURP, EOP, EOR,
     +             EORFLG, TYPE, SECT, PSEQ, DSEQ,
     +             IVAL, RVAL, CVAL(1:MAX(1,CLEN)))
  400    CONTINUE

C        Write the Directory lines

C        Write 1st line

         IDVAL(2)  = IPSTRT-PSSTRT+1
         IDVAL(10) = DSEQ
         IDVAL(13) = IPLINE-IPSTRT
         WRITE (LINE(1:80),'(8I8,I8.8,''D'',I7.7)')
     +         (IDVAL(I),I=1,10)

         WRITE (DLUNIT, '(A80)', REC=IDLINE) LINE

C        Write 2nd line

         WRITE (LINE(1:80),'(5I8,3A8,I8,''D'',I7.7)')
     +      IDVAL(1), (IDVAL(I),I=11,14), (CDVAL(I),I=1,3),
     +      IDVAL(15), DSEQ+1

         WRITE (DLUNIT, '(A80)', REC=IDLINE+1) LINE

         DSEQ   = DSEQ+2
         IDLINE = IDLINE+2
         IPSTRT = IPLINE

  500 CONTINUE

C     Write Terminate line

      WRITE (LINE(1:80),'(''S'',I7.7,''G'',I7.7,''D'',I7.7,''P'',I7.7,
     +                T73,''T'',I7.7)')
     +          NSL, NGLBL, NDE, PSEQ-1, 1
      WRITE (DLUNIT, '(A80)', REC=IPLINE) LINE

C     Copy the scratch file into the IGES file.

      REWIND (ALUNIT)
      REWIND (DLUNIT)

      ICURLN = 1

  600 CONTINUE
         READ (DLUNIT, '(A80)', REC=ICURLN, IOSTAT=IOCHKD) LINE
         WRITE (ALUNIT, '(A80)', IOSTAT=IOCHKA) LINE
         ICURLN = ICURLN+1
         IF ((IOCHKD .NE. 0) .OR. (IOCHKA .NE. 0)) THEN
            IER = -9999
            GOTO 9010
         ENDIF
      IF (LINE(73:73) .NE. 'T') GOTO 600

      IF (LVLE .GT. LVLO) THEN
         IER = 1
         CALL DTERR (0, SUBNAM, IER, 0)
      ENDIF
      GOTO 9999


C                 ** ERROR HANDLING **

 9010 CONTINUE

C        A direct error occurred

         CALL DTERR (1, SUBNAM, IER, 0)
         GOTO 9999

 9020 CONTINUE

C        An error occurred in a called subroutine

         CALL DTERR (4, SUBNAM, IER, 0)
         GOTO 9999

 9999 CONTINUE

      CLOSE (ALUNIT)
      CLOSE (DLUNIT)

      RETURN
      END

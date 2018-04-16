      SUBROUTINE D2IGRD (CMEM, IMEM, DMEM, FILNAM, LUNIT, LVLI,
     +                   IGI, LVLE, IER)

C     Read an IGES file into the DTRC Dynamic Memory.
C
C     DYNAMIC MEMORY INPUT-OUTPUT:
C        CMEM     Dynamically managed character data
C        IMEM     Dynamically managed integer data
C        DMEM     Dynamically managed double precision data
C     INPUT:
C        FILNAM   The name of the IGES file.
C        LUNIT    The I/O unit to use when working on the IGES file.
C                 If LUNIT < 0 on input, a line will be displayed on
C                 Unit DTJCON(6), giving the process status for each
C                 entity processed.
C                 IF LUNIT > 0, no status will be displayed.
C        LVLI     The IGES entity types to accept as input:
C                 1  NASA-IGES-NURBS only
C                 2  NASA-IGES
C                 3  Unrestricted
C     OUTPUT:
C        IGI      Dynamic memory pointer to the new IGES Index entity.
C        LVLE     The Level of IGES entities encountered.  See LVLI
C                 for meaning.
C        IER      Returned error value.  Zero implies no errors.
C                 Otherwise,
C                 +1 LVLE > LVLI.  Some non-NASA-IGES or non-NASA-IGES-
C                    NURBS entities were ignored.  If valid NASA-IGES
C                    entities point to these ignored entities, such
C                    otherwise valid entities should be deleted from the
C                    IGES data structure to obtain a perfectly consistent
C                    dataset.
C                 -1 Ambiguous, see printed error message.
C                 -n An error occurred while reading or processing line
C                    'n' of the IGES File.
C                 Note that IER=-1 can occur if any start-up processing
C                 fails, including failure to initialize the Dynamic
C                 Memory (D2INIT), or insufficient memory to initialize
C                 the IGES Index entity, etc.
C
C   CALLS:
C     DTJCON   Integer Constants
C     D0IGDO   Decided whether this is an entity to "DO"
C     D0IGRE   Read an IGES entity into dynamic memory
C     D0IGRP   Read an IGES parameter
C     D2IGDI   Define IGES Index entity
C     D2IGSI   Store integer
C     D2IGSD   Store double precision
C     D2IGSC   Store character
C     DTERR    Error Handling
C
C   HISTORY:
C     9/5/93     D. Parsons   Created
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_READ (CMEM, IMEM, DMEM, FILNAM, LUNIT, LVLI,
C     +                   IGI, LVLE, IER)

      EXTERNAL          D0IGDO, DTJCON
      LOGICAL           D0IGDO
      INTEGER           DTJCON

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      CHARACTER*(80) FILNAM
      INTEGER        LUNIT, LVLI, IGI, LVLE, IER

      CHARACTER*(80) LINE, STRING
      CHARACTER      EOP, EOR, TYPE
      CHARACTER*8    CDVAL(3)
      INTEGER        ALUNIT, DLUNIT, IERX, I, IDVAL(15)
      INTEGER        NSL, NDE, NPAR, NGLBL, NTERM, IPARST
      INTEGER        NLINES, ICURLN, IPCURP, IOCHK
      INTEGER        IPAR, IVAL, ILEN, IENT, ITYPE, IFORM, MSGUN
      LOGICAL        INIT, EORFLG, STATUS, EXFLAG
      DOUBLE PRECISION RVAL

      CHARACTER*6    SUBNAM
      DATA SUBNAM   /'D2IGRD'/

C****

C     LUNIT is the unit number passed in by the user.  May be negative.
C        User must define this because on some systems that file has to
C        be attached to that unit in the Job Control Language.
C     ALUNIT is ABS(LUNIT)
C     DLUNIT is a scratch file for direct access reading.

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

C     Find a free unit number for DLUNIT

      DLUNIT = 100
   10 CONTINUE
         INQUIRE (UNIT=DLUNIT, OPENED=EXFLAG)
         IF (EXFLAG) THEN
            DLUNIT = DLUNIT+1
            GOTO 10
         ENDIF

C     Copy the IGES file to a direct-access file

      OPEN (ALUNIT, FILE=FILNAM, STATUS='OLD', IOSTAT=IOCHK)

      IF (IOCHK .NE. 0) THEN
C        An error occurred in the OPEN statement
         IF (STATUS) WRITE (MSGUN,'(A,A8,A,I6,/)')
     +      ' FATAL ERROR: Unable to open ',FILNAM,' IOCHK = ',IOCHK
         IER = -1
         GOTO 9010
      ENDIF

      OPEN (DLUNIT, ACCESS='DIRECT', FORM='FORMATTED', RECL=80,
     +      STATUS='SCRATCH', IOSTAT=IOCHK)

      IF (IOCHK .NE. 0) THEN
C        An error occurred in the OPEN statement
         IF (STATUS) WRITE (MSGUN,'(A,A,I6,/)')
     +      ' FATAL ERROR: Unable to open direct access scratch file; ',
     +      'IOCHK = ',IOCHK
         IER = -1
         GOTO 9010
      ENDIF

C     *** 1ST PASS; INITIALIZE COUNTERS AND POINTERS ***

      REWIND (ALUNIT)
      REWIND (DLUNIT)

      NSL   = 0
      NDE   = 0
      NPAR  = 0
      NGLBL = 0
      NTERM = 0

      ICURLN = 0

  300 CONTINUE

         ICURLN = ICURLN+1
         READ (ALUNIT, '(A80)', IOSTAT=IOCHK) LINE
         IF (IOCHK .NE. 0) THEN
            IF (STATUS) WRITE (MSGUN,'(A,I6,A,I6,/)')
     +         ' FATAL ERROR: Unable to read line ',ICURLN,
     +         ' IOCHK = ',IOCHK
            IER = -ICURLN
            GOTO 9010
         ENDIF

C        Copy sequential file to direct access scratch file
         WRITE (DLUNIT, '(A80)', REC=ICURLN, IOSTAT=IOCHK) LINE
         IF (IOCHK .NE. 0) THEN
            IF (STATUS) WRITE (MSGUN,'(A,I6,A,I6,/)')
     +         ' FATAL ERROR: Unable to write line ',ICURLN,
     +         ' to scratch file; IOCHK = ',IOCHK
            IER = -ICURLN
            GOTO 9010
         ENDIF

         IF     (LINE(73:73) .EQ. 'S') THEN
            NSL = NSL+1
         ELSEIF (LINE(73:73) .EQ. 'G') THEN
            NGLBL = NGLBL+1
         ELSEIF (LINE(73:73) .EQ. 'D') THEN
            NDE = NDE+1
         ELSEIF (LINE(73:73) .EQ. 'P') THEN
            NPAR = NPAR+1
         ELSEIF (LINE(73:73) .EQ. 'T') THEN
            NTERM = NTERM+1
         ELSE
            IF (STATUS) WRITE (MSGUN,'(A,A,/)')
     +         ' FATAL ERROR: Unrecognized Section ID ',LINE(73:73)
            IER = -ICURLN
            GOTO 9010
         ENDIF

      IF (LINE(73:73) .NE. 'T') GOTO 300

      NLINES = ICURLN

C     Define the IGES Index entity

      INIT = .TRUE.
      CALL D2IGDI (CMEM, IMEM, DMEM, NSL, NDE, INIT, IGI, IERX)
      IF (IERX .NE. 0) THEN
         IF (STATUS) WRITE (MSGUN,'(A,/)')
     +      ' FATAL ERROR: Unable to initialize IGES Index Entity'
         IER = -1
         GOTO 9020
      ENDIF

C     *** 2ND PASS; PROCESS IGES ENTITY ***

      CLOSE (ALUNIT)
      REWIND (DLUNIT)

C     Read and store the "START" section lines.

      DO 400 ICURLN = 1, NSL
         READ (DLUNIT, '(A80)', REC=ICURLN, IOSTAT=IOCHK) LINE
         IF ((LINE(73:73) .NE. 'S') .OR. (IOCHK .NE. 0)) THEN
            IF (STATUS) WRITE (MSGUN,'(A,A,A,I6,/)')
     +         ' FATAL ERROR: Unexpected Section ID ',LINE(73:73),
     +         ' or unable to read line ',ICURLN,
     +         ' IOCHK = ',IOCHK
            IER = -ICURLN
            GOTO 9010
         ENDIF
         CALL D2IGSC (CMEM, IMEM, DMEM, 'C', LINE(1:72), IGI, 'S',
     +                ICURLN, IERX)
         IF (IERX .NE. 0) THEN
            IF (STATUS) WRITE (MSGUN,'(A,/)')
     +         ' FATAL ERROR: Unable to write START section'
            IER = -ICURLN
            GOTO 9020
         ENDIF
  400 CONTINUE

C     Read and store the "GLOBAL" section parameters.

      ICURLN = NSL+1

      READ (DLUNIT, '(A80)', REC=ICURLN, IOSTAT=IOCHK) LINE
      IF ((LINE(73:73) .NE. 'G') .OR. (IOCHK .NE. 0)) THEN
            IF (STATUS) WRITE (MSGUN,'(A,A,A,I6,A,I6,/)')
     +         ' FATAL ERROR: Unexpected Section ID ',LINE(73:73),
     +         ' or unable to read line ',ICURLN,
     +         ' IOCHK = ',IOCHK
         IER = -ICURLN
         GOTO 9010
      ENDIF

      EOP = ','
      EOR = ';'

      IPCURP = 1

      DO 500 IPAR = 1, 25
         CALL D0IGRP (LINE, DLUNIT, ICURLN, IPCURP, EOP, EOR,
     +                EORFLG, TYPE, IVAL, RVAL, STRING, ILEN)
C        Is this parameter supposed to be character?
         IF ((IPAR .GE. 1 .AND. IPAR .LE. 6) .OR.
     +       (IPAR .EQ. 12) .OR.
     +       (IPAR .EQ. 15) .OR.
     +       (IPAR .EQ. 18) .OR.
     +       (IPAR .EQ. 21) .OR.
     +       (IPAR .EQ. 22) .OR.
     +       (IPAR .EQ. 25) )          THEN
            IF (ILEN .GT. 0) THEN
               CALL D2IGSC (CMEM, IMEM, DMEM, 'C', STRING(1:ILEN), IGI,
     +                     'G', IPAR, IERX)
               IF (IERX .NE. 0) THEN
                  IF (STATUS) WRITE (MSGUN,'(A,/)')
     +               ' FATAL ERROR: Unable to write Global Section'
                  IER = -ICURLN
                  GOTO 9020
               ENDIF
            ENDIF
            IF ((IPAR .EQ. 1) .AND. (TYPE .EQ. 'C')) EOP = STRING(1:1)
            IF ((IPAR .EQ. 2) .AND. (TYPE .EQ. 'C')) EOR = STRING(1:1)

         ELSE
C           Integer or Real?
            IF ((IPAR .GE. 7 .AND. IPAR .LE. 11) .OR.
     +          (IPAR .EQ. 14) .OR.
     +          (IPAR .EQ. 16) .OR.
     +          (IPAR .EQ. 23) .OR.
     +          (IPAR .EQ. 24) )       THEN
               IF (TYPE .NE. 'I') THEN
C                 This parameter needs to be an integer
                  IF (STATUS) WRITE (MSGUN,'(A,A,/)')
     +               ' FATAL ERROR: Non-integer found in Global',
     +               ' integer position'
                  IER = -ICURLN
                  GOTO 9010
               ENDIF
               CALL D2IGSI (CMEM, IMEM, DMEM, 'I', IVAL, IGI, 'G',
     +                      IPAR, IERX)
               IF (IERX .NE. 0) THEN
                  IF (STATUS) WRITE (MSGUN,'(A,/)')
     +                ' FATAL ERROR: Unable to write Global integer.'
                  IER = -ICURLN
                  GOTO 9020
               ENDIF
            ELSE
               IF (TYPE .NE. 'R') RVAL = IVAL
               CALL D2IGSD (CMEM, IMEM, DMEM, 'E', RVAL, IGI, 'G',
     +                      IPAR, IERX)
               IF (IERX .NE. 0) THEN
                  IF (STATUS) WRITE (MSGUN,'(A,/)')
     +                ' FATAL ERROR: Unable to write Global real.'
                  IER = -ICURLN
                  GOTO 9020
               ENDIF
            ENDIF
         ENDIF  
         
C        If End-of-record encountered prematurely, leave default values.

         IF (EORFLG) GOTO 510

         IF (LINE(73:73) .NE. 'G') THEN
            IF (STATUS) WRITE (MSGUN,'(A,A,/)')
     +         ' FATAL ERROR: Unexpected Section ID ',LINE(73:73)
            IER = -ICURLN
            GOTO 9010
         ENDIF
  500 CONTINUE
  510 CONTINUE


C     Read and store the entities.

      DO 600 IENT = 1, NDE, 2
         ICURLN  = NSL+NGLBL+IENT
         IPARST = NSL+NGLBL+NDE+1

C        1st Line of directory

         READ (DLUNIT, '(A80)', REC=ICURLN, IOSTAT=IOCHK) LINE
         IF ((LINE(73:73) .NE. 'D') .OR. (IOCHK .NE. 0)) THEN
            IF (STATUS) WRITE (MSGUN,'(A,I6,/)')
     +         ' FATAL ERROR: Unable to read Directory line ',ICURLN
            IER = -ICURLN
            GOTO 9010
         ENDIF

         READ (LINE(1:80),'(BZ,9I8,1X,I7)', IOSTAT=IOCHK)
     +         (IDVAL(I),I=1,10)
         IF (IOCHK .NE. 0) THEN
            IER = -ICURLN
            GOTO 9010
         ENDIF

C        2nd Line of directory

         ICURLN = ICURLN+1
         READ (DLUNIT, '(A80)', REC=ICURLN, IOSTAT=IOCHK) LINE
         IF ((LINE(73:73) .NE. 'D') .OR. (IOCHK .NE. 0)) THEN
            IF (STATUS) WRITE (MSGUN,'(A,I6,/)')
     +         ' FATAL ERROR: Unable to read Directory line ',ICURLN
            IER = -ICURLN
            GOTO 9010
         ENDIF

         READ (LINE(9:72),'(4I8,3A8,I8)', IOSTAT=IOCHK)
     +         (IDVAL(I),I=11,14), (CDVAL(I),I=1,3), IDVAL(15)
         IF (IOCHK .NE. 0) THEN
            IER = -ICURLN
            GOTO 9010
         ENDIF

         ITYPE = IDVAL(1)
         IFORM = IDVAL(14)

         IF (ITYPE .EQ. 0) LVLE = MAX(LVLE,1)

         IF ((ITYPE .NE. 0) .AND.
     +      D0IGDO (ITYPE, IFORM, LVLI, LVLE)) THEN
C           Process this entity
            IF (STATUS) WRITE (MSGUN,'(A,I4,A,I3,A,I5)')
     +         ' D2IGRD: Processing Entity (IDE) Number ',IENT,
     +         '; Entity Type = ',ITYPE, ' Form = ',IFORM
            CALL D0IGRE (CMEM, IMEM, DMEM, DLUNIT, IENT, IGI, ITYPE,
     +                   IDVAL, CDVAL, IPARST, IERX)
            IF (IERX .NE. 0) THEN
               IF (STATUS) WRITE (MSGUN,'(A,/)')
     +            ' FATAL ERROR: Unable to read entity'
               IER = IERX
               GOTO 9020
            ENDIF
         ELSE
C           This entity is NULL or is excluded by the LVLI parameter
            IF (STATUS) WRITE (MSGUN,'(A,I4,A,I3,A,I5)')
     +         ' D2IGRD: Skipping Entity (IDE) Number ',IENT,
     +         '; Entity Type = ',ITYPE, ' Form = ',IFORM
            CALL D2IGSI (CMEM, IMEM, DMEM, '*', 0, IGI, 'I',
     +                   IENT, IERX)
            IF (IERX .NE. 0) THEN
               IF (STATUS) WRITE (MSGUN,'(A,/)')
     +            ' FATAL ERROR: Unable to write NULL entity'
               IER = -ICURLN
               GOTO 9020
            ENDIF
         ENDIF
  590    CONTINUE

  600 CONTINUE

      IF (LVLE .GT. LVLI) THEN
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

      CLOSE (DLUNIT)

      RETURN
      END

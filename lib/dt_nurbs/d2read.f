      SUBROUTINE D2READ (CMEM, IMEM, DMEM, FILNAM, LUNIT, LTYPTR, NUMPA,
     +                   MEMPA, LENPA, IER)

C     Read a D2-File-Format file back into dynamic memory.  Normally such
C     a file is created by D2WRIT.  If successful, all entities are recreated
C     and all old pointers are replaced with the corresponding new pointers.
C     Old data type numbers may be optionally translated into new ones under
C     the control of a type translation array.
C
C     D2READ creates a translation table for old pointers to new pointers.
C     Hence, D2READ temporarily requires more space than that occupied by the
C     the data entities themselves.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     FILNAM  File name to use for the file to be written.
C     LUNIT   FORTRAN logical unit number to use during file operations.
C     LTYPTR  Optional type translation array.  
C             LTYPTR(N) >  0 means translate file user data type N to LTYPTR(N).
C             LTYPTR(N) =  0 means end LTYPTR, but treat all M >= N as if -1.
C             LTYPTR(N) = -1 means search for existing matching type name.
C             LTYPTR(N) = -2 means end LTYPTR, but treat all M >= N as if -3.
C             LTYPTR(N) = -3 means do not translate this type.  It is an error.
C             The maximum size LTYPTR could be is 191.  To omit it entirely,
C             but cause type translation to occur by matching names, supply
C             a single integer argument of value zero.  To omit LTYPTR entirely,
C             and cause all user-defined types to be errors, supply a single
C             integer argument of value minus two.
C     NUMPA   Length of array MEMPA.
C   OUTPUT:
C     MEMPA   Array of new MEM pointers to the entities that were selected to be
C             written to the file and in the same order.  (Additional entities
C             may have been written merely because they were linked to these.)
C     LENPA   Length of list of selected pointers in MEMPA, or, if zero,
C             an indication that all active entities were selected.
C     IER     Returned error code.  Zero implies no errors.
C             -1  = Unable to open given file on given unit.
C             -2  = Insufficient character data space available.
C             -3  = Insufficient integer data space available.
C             -4  = Insufficient double precision space available.
C             -5  = Unable to read first line.
C             -6  = File does not begin with 'D2 File Format'.
C             -7  = D2 File Format version of file is later than D2READ version.
C             -8  = Unable to interpret format info in first line.
C             -9  = Value out of range in LTYPTR.
C             -10 = Reading or formatting error in second line.
C             -11 = Impossible values in second line.
C             -12 = Untranslatable user-defined data type encountered.
C             -13 = Impossible error from D1DEFE
C             -14 = Impossible error from D1DEFX
C             for error code numbers < -30:
C             -10*n-0 = Read error at line n (in list of pointers section).
C             -10*n-1 = Read error at line n (in list of user types section).
C             -10*n-2 = Read error at line n (in blank line preceding entity).
C             -10*n-3 = Read error at line n (in initial line of entity).
C             -10*n-4 = Read error at line n (in character data of entity).
C             -10*n-5 = Read error at line n (in integer data of entity).
C             -10*n-6 = Read error at line n (in double prec. data of entity).
C             -10*n-7 = Unknown pointer error at line n.
C             -10*n-8 = Unknown user data type at line n.
C             -10*n-9 = Unknown character code in line n.
C
C   CALLS:
C     D1DEFE
C     D2ERAS
C     D2FETN
C     DTERPT
C     DTERR
C
C   HISTORY:
C     20Sep93  P. Kraushar   Created.
C*****
C
C   Long Name Alias:
C     ENTRY D2_READ_DT_FILE (CMEM, IMEM, DMEM, FILNAM, LUNIT, LTYPTR,
C    +                       NUMPA, MEMPA, LENPA, IER)
      CHARACTER CMEM*(*), FILNAM*(*)
      INTEGER IMEM(*), LTYPTR(*), MEMPA(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER LUNIT, NUMPA, LENPA, IER

      CHARACTER VERSN*18, FORMI*12, FORMIA*12, FORMD*12, FORMC*5
      CHARACTER CH*40, CHX*40, CCH*33, PFLAG(10)
      CHARACTER SUBNAM*6, LINE1*38, TYPNAM*74, MEMTYP*74
      INTEGER NUMENL, NUMENT, NUMUT, NUMC, NUMI, NUMD, NI, LIW, ND
      INTEGER LDW, LDD, LDE, LINUM, IPTR, IHDR, IH, I, ILIM, J, K, L, M
      INTEGER LPTR(10), DTJCON, LINLEN, NHASH, IHASH0, LUTR(191), IERX
      INTEGER NUMCRE, ITYP
      DOUBLE PRECISION VER, DTMCON
      LOGICAL ISOPEN
      EXTERNAL DTJCON, DTMCON

      DATA VERSN /'D2 File Format 1.0'/
      DATA FORMC /'(A,A)'/
      DATA CCH /'@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_~'/
      DATA NUMENL, NUMENT, NUMUT, NUMC, NUMI, NUMD /6*0/
      DATA SUBNAM /'D2READ'/

C****
      IER = -1
      OPEN (LUNIT, FILE=FILNAM, ERR=9900)

C     Process first line of file      
      IER = -5
      READ (LUNIT, 10, ERR=9900) LINE1
   10 FORMAT(A)
      IF (LINE1(1:14) .NE. VERSN(1:14)) THEN
         IER = -6
         GOTO 9900
      ENDIF
      IF (LINE1(1:18) .NE. VERSN) THEN
         READ (LINE1(15:18), '(F4.0)') VER
         IF (VER .GT. 1.0) THEN
            IER = -7
            GOTO 9900
         ENDIF
      ENDIF
C     Reconstruct formats used in file from info in rest of first line
      IER = -8
      READ (LINE1(19:38), 20, ERR=9900) LINLEN, LIW, LDW, LDD, LDE
   20 FORMAT (5I4)
      NI = MIN (LINLEN/(LIW+1), 10)
      ND = LINLEN/LDW
      WRITE (FORMIA, 30) NI, LIW
   30 FORMAT ('(',I2.2,'(I',I2.2,',A1))')
C     FORMIA ends up looking like '(06(I11,A1))' for 80-char line, 32-bit int
      FORMI = FORMIA
      FORMI(9:10) = '1X'
      WRITE (FORMD, 40) ND, LDW, LDD, LDE
   40 FORMAT ('(',I1,'G',I2.2,'.',I2.2,'E',I2.2,')')
C     FORMD ends up looking like '(3G23.16E3)' for 80-char line, 64-bit dble

C     Process second line of file to get various file counts
      IER = -10
      READ (LUNIT, FORMI, ERR=9900) NUMENL, NUMENT, NUMUT, NUMC, NUMI,
     +                              NUMD
      IF (NUMENL .LT. 0 .OR. NUMENT .LT. 0 .OR. NUMUT .LT. 0
     + .OR. NUMC .LT. 0 .OR. NUMI .LT. 0 .OR. NUMD .LT. 0 
     + .OR. NUMUT .GT. 191) THEN
         IER = -11
         GOTO 9900
      ENDIF
C     Is there enough free space to read in the file?
      IF (NUMC .GT. IMEM(1)-IMEM(4)) THEN
         IER = -2
         GOTO 9900
      ENDIF
      NHASH = NUMENT + 1 + NUMENT/12
      IF (NUMI+2*NHASH .GT. IMEM(16)-IMEM(5)) THEN
         IER = -3
         GOTO 9900
      ENDIF
      IF (NUMD .GT. IMEM(3)-IMEM(6)) THEN
         IER = -4
         GOTO 9900
      ENDIF

C     Process the type translation table, if needed
      IF (NUMUT .GT. 0) THEN
         DO 120 I=1,191
            IF (LTYPTR(I) .GT. 0) THEN
               IF (LTYPTR(I) .GT. 255) THEN
                  IER = -9
                  GOTO 9900
               ENDIF
               LUTR(I) = LTYPTR(I)
            ELSE IF (LTYPTR(I) .EQ. 0) THEN
               DO 100 J=I,191
                  LUTR(J) = 0
  100          CONTINUE
               GOTO 130
            ELSE IF (LTYPTR(I) .EQ. -1) THEN
               LUTR(I) = 0
            ELSE IF (LTYPTR(I) .EQ. -2) THEN
               DO 110 J=I,191
                  LUTR(J) = -1
  110          CONTINUE
               GOTO 130
            ELSE IF (LTYPTR(I) .EQ. -3) THEN
               LUTR(I) = -1
            ELSE
               IER = -9
               GOTO 9900
            ENDIF
  120    CONTINUE
C        Next line is exit point for short LTYPTR arrays (the usual case)
  130    CONTINUE
      ENDIF

C     Load the pointer list into the hash table and allocate corresponding
C     header blocks (with no data space allocation, initially).  The hash
C     table is located in the free area of IMEM just below the area which
C     might be taken up by new header blocks, and well above the total area
C     needed by all the new entities.  Locations IHASH0 <= i < IHASH0+NHASH
C     hold old pointers and locations IHASH0+NHASH <= i < IHASH0+2*NHASH
C     hold the corresponding new pointers.  About 8% free space is included
C     to cut down the searching time.
      IHASH0 = IMEM(16) - 8*NUMENT - 2*NHASH
      DO 200 I=IHASH0,IHASH0+NHASH-1
         IMEM(I) = 0
  200 CONTINUE
      LENPA = MIN (NUMENL, NUMPA)
      LINUM = 2
      IER = 0
      NUMCRE = 0
      ILIM = NI
      DO 230 I=1,NUMENT,NI
         LINUM = LINUM + 1
         IF (I+NI .GT. NUMENT)  ILIM = NUMENT - I + 1
         READ (LUNIT, FORMI, ERR=9800) (LPTR(J), J=1,ILIM)
         DO 220 J=1,ILIM
            IH = MOD(LPTR(J)/2048, NHASH)
            IF (IMEM(IHASH0+IH) .NE. 0) THEN
C              START LOOP
  210             CONTINUE
                  IH = IH + 1
                  IF (IH .EQ. NHASH)  IH = 0
                  IF (IMEM(IHASH0+IH) .NE. 0) GOTO 210
C                 Since NHASH > NUMENT, this search for an open spot in the
C                 hash table must succeed.
C              END LOOP
            ENDIF
            IMEM(IHASH0+IH) = LPTR(J)
C           Allocate a new header block (for now it looks like a null string)
            CALL D1DEFE (CMEM, IMEM, DMEM, 255, 0, 0, 0, .FALSE.,
     +                   IPTR, IERX)
            IF (IERX .NE. 0) THEN
               IER = -13
               GOTO 9800
            ENDIF
            NUMCRE = NUMCRE + 1
            IMEM(IHASH0+NHASH+IH) = IPTR
C           Record first LENPA new pointers in MEMPA
            IF (I+J-1 .LE. LENPA)  MEMPA(I+J-1) = IPTR
  220    CONTINUE
  230 CONTINUE

C     Now read the user types of the file and finish the translation table
      IF (NUMUT .GT. 0) THEN
         IER = -1
         DO 320 I=1,NUMUT
            LINUM = LINUM + 1
            READ (LUNIT, 50, ERR=9800) ITYP, J, TYPNAM
   50       FORMAT (2I3,A)
            IF (ITYP .LE. 0 .OR. ITYP .GT. 191 .OR. J .LT. 0) THEN
               GOTO 9800
            ELSE IF (LUTR(ITYP) .EQ. 0) THEN
               IF (J .EQ. 0) GOTO 9800
               CALL DTERPT (0)
               IH = 1
               L = 1
C              START LOOP
  300             CONTINUE
                  CALL D2FETN (CMEM, IMEM, DMEM, IH, MEMTYP, K, IERX)
                  IF (IERX .EQ. 0) THEN
                     IF (TYPNAM(1:J) .EQ. MEMTYP(1:K)) THEN
                        LUTR(ITYP) = IH
                        GOTO 310
                     ENDIF
                  ELSE IF (IERX .EQ. -3) THEN
                     IH = 256
                     L = -1
                  ELSE
                     GOTO 310
                  ENDIF
                  IF (L .EQ. -1 .AND. IH .EQ. 192) GOTO 310
                  IF (L .EQ. 1 .AND. IH .EQ. 191) THEN
                     IH = 255
                     L = -1
                  ELSE
                     IH = IH + L
                  ENDIF
                  GOTO 300
C              END LOOP
  310          CONTINUE
               CALL DTERPT (1)
               IF (LUTR(ITYP) .EQ. 0) THEN
                  IER = -12
                  GOTO 9800
               ENDIF
            ELSE IF (LUTR(ITYP) .LT. 0) THEN
               IER = -12
               GOTO 9800
            ENDIF
  320    CONTINUE
      ENDIF

C     Read in the individual entities
      DO 800 I=1,NUMENT
C        Check blank line at begining of entity
         IER = -2
         LINUM = LINUM + 1
         READ (LUNIT, '(A)', ERR=9800) LINE1
         IF (LINE1 .NE. ' ') GOTO 9800

C        Process initial line of entity
         IER = -3
         LINUM = LINUM + 1
         READ (LUNIT, FORMI, ERR=9800) IPTR, NUMC, NUMI, NUMD, ITYP
         IF (IPTR .LE. 0 .OR. NUMC .LT. 0 .OR. NUMI .LT. 0
     +       .OR. NUMD .LT. 0 .OR. ITYP .LE. 0 .OR. ITYP .GT. 255)
     +       GOTO 9800
         IF (ITYP .LE. 191) THEN
            ITYP = LUTR(ITYP)
            IF (ITYP .LE. 0) THEN
               IER = -8
               GOTO 9800
            ENDIF
         ENDIF
C        Find IPTR in hash table and translate to new pointer
         IH = MOD (IPTR/2048, NHASH)
         IF (IMEM(IHASH0+IH) .NE. IPTR) THEN
            IER = -7
C           START LOOP
  400          CONTINUE
               IF (IMEM(IHASH0+IH) .EQ. 0) GOTO 9800
               IH = IH + 1
               IF (IH .EQ. NHASH) IH = 0
               IF (IMEM(IHASH0+IH) .NE. IPTR) GOTO 400
C           END LOOP
         ENDIF
         IPTR = IMEM(IHASH0+NHASH+IH)
         IHDR = IPTR/256
C        Expand new entity to proper size
         CALL D1DEFX (CMEM, IMEM, DMEM, IHDR, NUMC, NUMI, NUMD, .FALSE.,
     +                J, IERX)
         IF (IERX .NE. 0) THEN
            IER = -14
            GOTO 9800
         ENDIF
C        Replace data type info
         IMEM(IHDR+7) = IMEM(IHDR+7) - 255 + ITYP

C        Read entity's character data
         IF (NUMC .GT. 0) THEN
            IER = -4
            L = IMEM(IHDR+1)
            ILIM = 40
            DO 510 J=1,NUMC,40
               LINUM = LINUM + 1
               IF (J+40 .GT. NUMC) ILIM = NUMC - J + 1
               READ (LUNIT, FORMC, ERR=9800) CH, CHX
               DO 500 K=1,ILIM
                  IF (CHX(K:K) .EQ. ' ') THEN
                     CMEM(L:L) = CH(K:K)
                  ELSE IF (CHX(K:K) .EQ. '^') THEN
                     M = INDEX (CCH, CH(K:K)) - 1
                     IF (M .LT. 0) THEN
                        IER = -9
                        GOTO 9800
                     ELSE IF (M .LT. 32) THEN
                        CMEM(L:L) = CHAR (M)
                     ELSE
                        CMEM(L:L) = CHAR (127)
                     ENDIF
                  ELSE IF (CHX(K:K) .EQ. '''') THEN
                     CMEM(L:L) = CHAR (ICHAR (CH(K:K)) + 128)
                  ELSE IF (CHX(K:K) .EQ. '"') THEN
                     M = INDEX (CCH, CH(K:K)) - 1
                     IF (M .LT. 0) THEN
                        IER = -9
                        GOTO 9800
                     ELSE IF (M .LT. 32) THEN
                        CMEM(L:L) = CHAR (M+128)
                     ELSE
                        CMEM(L:L) = CHAR (255)
                     ENDIF
                  ELSE
                     IER = -9
                     GOTO 9800
                  ENDIF
                  L = L + 1
  500          CONTINUE
  510       CONTINUE
         ENDIF

C        Read entity's integer data, translating old pointers to new
         IF (NUMI .GT. 0) THEN
            IER = -5
            L = IMEM(IHDR+2)
            ILIM = NI
            DO 620 J=1,NUMI,NI
               LINUM = LINUM + 1
               IF (J+NI .GT. NUMI) ILIM = NUMI - J + 1
               READ (LUNIT, FORMIA, ERR=9800) 
     +          (LPTR(K), PFLAG(K), K=1,ILIM)
               DO 610 K=1,ILIM
                  IF (PFLAG(K) .EQ. '*' .AND. LPTR(K) .NE. 0) THEN
C                    Find LPTR(K) in hash table and translate to new pointer
                     IH = MOD (LPTR(K)/2048, NHASH)
                     IF (IMEM(IHASH0+IH) .NE. LPTR(K)) THEN
C                       START LOOP
  600                      CONTINUE
                           IF (IMEM(IHASH0+IH) .EQ. 0) THEN
                              IER = -7
                              GOTO 9800
                           ENDIF
                           IH = IH + 1
                           IF (IH .EQ. NHASH) IH = 0
                           IF (IMEM(IHASH0+IH) .NE. LPTR(K)) GOTO 600
C                       END LOOP
                     ENDIF
                     LPTR(K) = IMEM(IHASH0+NHASH+IH)
                  ENDIF
                  IMEM(L) = LPTR(K)
                  L = L + 1
  610          CONTINUE
  620       CONTINUE
         ENDIF

C        Read entity's double precision data
         IF (NUMD .GT. 0) THEN
            IER = -6
            L = IMEM(IHDR+3) - 1
            ILIM = ND
            DO 700 J=1,NUMD,ND
               LINUM = LINUM + 1
               IF (J+ND .GT. NUMD) ILIM = NUMD - J + 1
               READ (LUNIT, FORMD, ERR=9800) (DMEM(L+K), K=1,ILIM)
               L = L + ND
  700       CONTINUE
         ENDIF

  800 CONTINUE
      IER = 0
      RETURN

C     Error Handling

C     Erase all new entities
 9800 CONTINUE
C     Suppress error messages during cleanup
      CALL DTERPT (0)
      DO 9850 I=1,NUMCRE
         IPTR = IMEM(IMEM(17))
         CALL D2ERAS (CMEM, IMEM, DMEM, IPTR, IERX)
 9850 CONTINUE
      CALL DTERPT (1)
      IF (IER .GT. -10)  IER = IER - 10*LINUM

C     General error report
 9900 CONTINUE
      CALL DTERR (1, SUBNAM, IER, 0)
      INQUIRE (LUNIT, OPENED=ISOPEN)
      IF (ISOPEN)  CLOSE (LUNIT)
      LENPA = -1
      RETURN
      END

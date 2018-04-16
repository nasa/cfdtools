      SUBROUTINE D2WRIT (CMEM, IMEM, DMEM, LENPA, MEMPA, SUBROU, FILNAM,
     +                   LUNIT, IER)

C     Write dynamically managed data to a file for subsequent rereading by
C     D2READ.  Specifically, write user-listed data entities together with 
C     all known (recursively) linked entities to a file in D2 File Format.  
C     Alternatively (zero-length list), write all currently active entities
C     to the file.
C
C     In order that D2WRIT may be used to respond to "insufficient-memory"
C     problems, it is designed not to allocate any dynamic memory itself and
C     also to avoid terminating if it can possibly proceed to write a
C     useful file.
C
C     The user must supply a subroutine in the argument SUBROU which
C     identifies the pointers in user-defined data types.  A default
C     subroutine named D2LPUT is provided in the Library.  D2LPUT simply
C     reports unsatisfactory results for any user-defined data type.  The
C     behavior of D2WRIT in such cases is to assume the entity contains
C     no pointers in its integer data space (beyond the unsat. result)
C     The subroutine SUBROU must satisfy a number of severe restrictions,
C     among which is that it cannot itself call any Library subroutines.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     LENPA   Length of list of selected pointers in MEMPA, or, if zero,
C             an indication that all active entities are selected.
C     MEMPA   Array of MEM pointers to the entities to be written to the file.
C     SUBROU  User-supplied subroutine which identifies pointers in user-
C             defined data types.  If there are no user-defined types, or if
C             the user has not created an appropriate subroutine, the default
C             Library routine D2LPUT must be used.  The required calling
C             sequence is:  SUBROUTINE SUBROU (CMEM, IMEM, DMEM, LP),  where
C             LP is declared  INTEGER LP(14)  and initialized as follows:
C                LP(1) = Index of header block of entity to search
C                LP(2) = Data type number of entity
C                LP(3) = 0   (to indicate that this is the initial call
C                             for the entity)
C                LP(4) = One less than the absolute index in IMEM for the
C                        beginning of the entity's integer data space.
C                LP(5) = The absolute index in IMEM for the end of the
C                        entity's integer data space.
C                LP(6 to 14) = Unintialized
C             Upon return from each call, LP(3) is assumed to be the absolute
C             index in IMEM of the beginning of the next block of consecutive 
C             pointer values in the entity, and LP(4) is the absolute index
C             of the end of the block.  If the values in LP(3) and LP(4) do
C             not increase consistently, D2WRIT concludes that SUBROU is not
C             handling that data type satisfactorily and reports a warning
C             error to that effect.  D2WRIT makes no changes to LP between
C             consecutive calls to SUBROU for the same entity.  SUBROU is
C             expected to return an LP(3) value greater than the absolute
C             index of the end of the entity's integer data space to indicate 
C             there are no more pointers in the integer data space of the 
C             entity.  Since D2WRIT temporarily modifies the header blocks of 
C             the selected entities during processing, SUBROU is not permitted 
C             to call any other Library routine because they would conclude 
C             that dynamic memory was corrupted.  
C     FILNAM  File name to use for the file to be written.
C     LUNIT   FORTRAN logical unit number to use during file operations.
C   OUTPUT:
C     IER     Returned error code.  Zero implies no errors.  D2WRIT is
C             designed to proceed to completion without reporting error if
C             at all possible.
C             >0  = Number of user-defined data types for which SUBROU gave
C                   unsatisfactory results.  See warning errors for list.
C             -1  = Given list does not include any valid pointers.
C             -2  = Unable to open given file and/or unit for writing.
C
C   CALLS:
C     D0LPLT
C     D2FETN
C     DTERPT
C     DTERR
C
C   HISTORY:
C     20Sep93  P. Kraushar   Created.
C*****
C
C   Long Name Alias:
C     ENTRY D2_WRITE_DT_FILE (CMEM, IMEM, DMEM, LENPA, MEMPA, SUBROU,
C    +                        FILNAM, LUNIT, IER)
      CHARACTER CMEM*(*), FILNAM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER LENPA, MEMPA(*), LUNIT, IER
      EXTERNAL SUBROU

      INTEGER LINLEN
      PARAMETER (LINLEN=80)

      CHARACTER VERSN*18, FORMI*12, FORMIA*12, FORMD*12, FORMC*5
      CHARACTER CH*40, CHX*40, CCH*33, OCCUR(191), PFLAG(10)
      CHARACTER SUBNAM*6, TYPNAM*74
      INTEGER NUMENL, NUMENT, NUMUT, NUMC, NUMI, NUMD, NI, LIW, ND
      INTEGER LDW, LDD, LDE, IFIRST, ILASTH, IPTR, IHDR, IH, I, ILIM, J
      INTEGER LP(14), LPTR(10), DTJCON, ITYP
      DOUBLE PRECISION DTMCON
      LOGICAL ODDCH
      EXTERNAL DTJCON, DTMCON

      DATA VERSN /'D2 File Format 1.0'/
      DATA FORMC, OCCUR /'(A,A)', 191*' '/
      DATA CCH /'@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_~'/
      DATA SUBNAM /'D2WRIT'/

C****
      IER = 0
      NUMENL = 0
      NUMENT = 0
      NUMUT = 0
      NUMC = 0
      NUMI = 0
      NUMD = 0
      IF (LENPA .GT. 0) THEN
C        Mark listed entities for writing
         IFIRST = 0
         DO 100 I=1,LENPA
            IPTR = MEMPA(I)
            IHDR = IPTR/256
            IF (IHDR .GE. IMEM(16) .AND. IHDR .LT. IMEM(2)) THEN
               IF (IMEM(IHDR) .EQ. IPTR) THEN
C                 Header index in range and not already listed
C                 Add to linked list whose head is IFIRST and uses pointer
C                 check value (first element in header block) for chaining.
                  IF (NUMENL .EQ. 0) THEN
                     IFIRST = IPTR
                  ELSE
                     IMEM(ILASTH) = IPTR
                  ENDIF
                  IMEM(IHDR) = 0
C                 Also mark last element of header block so selected entities
C                 may be easily distinguished from unselected entities while
C                 following links and identifying pointers in the file.
                  IMEM(IHDR+7) = -IMEM(IHDR+7)
                  ILASTH = IHDR
                  NUMENL = NUMENL + 1
               ENDIF
            ENDIF
  100    CONTINUE
         IF (NUMENL .EQ. 0) THEN
C           No valid pointers occurred on list MEMPA
            IER = -1
            GOTO 9900
         ENDIF
C        Now start traversing the selection list and adding entities pointed
C        to by selected entities.  When traversal catches up to the end, it
C        is a complete list.
         NUMENT = NUMENL
         IPTR = IFIRST
C        START LOOP
  110       CONTINUE
            IF (IPTR .EQ. 0) GOTO 150
            IHDR = IPTR/256
            ITYP = MOD (ABS(IMEM(IHDR+7)), 256)
            IF (ITYP .LT. 192)  OCCUR(ITYP) = '*'
C           Accumulate total space used by selected entities
            NUMC = NUMC + ABS(IMEM(IHDR+4))
            NUMI = NUMI + ABS(IMEM(IHDR+5)) + 8
            NUMD = NUMD + ABS(IMEM(IHDR+6))
C           Set up LP for D0LPLT or SUBROU
            LP(1) = IHDR
            LP(2) = ITYP
            LP(3) = 0
            LP(4) = IMEM(IHDR+2) - 1
            LP(5) = LP(4) + ABS(IMEM(IHDR+5))
            ILIM = LP(5)
C           START LOOP
  120          CONTINUE
               IF (ITYP .GE. 192) THEN
                  CALL D0LPLT (CMEM, IMEM, DMEM, LP)
               ELSE
                  J = LP(4)
                  CALL SUBROU (CMEM, IMEM, DMEM, LP)
                  IF (J .GT. LP(3) .OR. LP(3) .GT. LP(4)) THEN
C                    SUBROU did not handle this type satisfactorily
                     LP(3) = ILIM + 1
                  ENDIF
               ENDIF
               IF (LP(3) .GT. ILIM)  GOTO 140
               DO 130 I=LP(3),LP(4)
                  IH = IMEM(I)/256
                  IF (IH .GE. IMEM(16) .AND. IH .LT. IMEM(2)) THEN
C                    IH is in valid header range
                     IF (IMEM(IH) .EQ. IMEM(I)) THEN
C                       Looks like a good pointer not previously seen
C                       Add it to the list
                        IMEM(ILASTH) = IMEM(IH)
                        IMEM(IH) = 0
                        IMEM(IH+7) = -IMEM(IH+7)
                        ILASTH = IH
                        NUMENT = NUMENT + 1
                     ENDIF
                  ENDIF
  130          CONTINUE
               GOTO 120
C           END LOOP
  140       CONTINUE
            IPTR = IMEM(IHDR)
            GOTO 110
C        END LOOP
  150    CONTINUE
      ELSE
C        Select all active entities (i.e. in active list and not deleted)
         IHDR = IMEM(17)
         IFIRST = 0
C        START LOOP
  200       CONTINUE
            IF (IMEM(IHDR) .GT. 0) THEN
C              Active entity, put it on the list, note user type, and mark it.
               IPTR = IMEM(IHDR)
               IMEM(IHDR) = IFIRST
               IFIRST = IPTR
               ITYP = MOD (IMEM(IHDR+7), 256)
               IF (ITYP .LT. 192)  OCCUR(ITYP) = '*'
               IMEM(IHDR+7) = -IMEM(IHDR+7)
C              Accumulate total space used by selected entities
               NUMC = NUMC + ABS(IMEM(IHDR+4))
               NUMI = NUMI + ABS(IMEM(IHDR+5)) + 8
               NUMD = NUMD + ABS(IMEM(IHDR+6))
               NUMENT = NUMENT + 1
            ENDIF
            IHDR = ABS(IMEM(IHDR+7))/256
C           Stop at the entity used to list user type names which was created
C           by D2INIT.  It will always be at the end of the active list.
            IF (IHDR .LT. IMEM(2) - 8) GOTO 200
C        END LOOP
         IF (NUMENT .EQ. 0) THEN
            IER = -1
            GOTO 9900
         ENDIF
      ENDIF

C     Tally number of user types occurring
      DO 300 I=1,191
         IF (OCCUR(I) .EQ. '*')  NUMUT = NUMUT + 1
  300 CONTINUE
C     Calculate and prepare output formats
      LIW = 2 + INT (LOG10 (REAL (DTJCON(2))))
      NI = MIN (LINLEN/(LIW+1), 10)
      WRITE (FORMIA, 10) NI, LIW
   10 FORMAT ('(',I2.2,'(I',I2.2,',A1))')
C     FORMIA ends up looking like '(06(I11,A1))' for 80-char line, 32-bit int
      FORMI = FORMIA
      FORMI(9:10) = '1X'
      LDD = DTJCON(11)
      LDE = 1 + INT (LOG10 (LOG10 (DTMCON(2))))
      LDW = 4 + LDD + LDE
      ND = LINLEN/LDW
      WRITE (FORMD, 20) ND, LDW, LDD, LDE
   20 FORMAT ('(',I1,'G',I2.2,'.',I2.2,'E',I2.2,')')
C     FORMD ends up looking like '(3G23.16E4)' for 80-char line, 64-bit dble
C     Try to open output file
      OPEN (LUNIT, FILE=FILNAM, ERR=9898)
C     Write version and formatting info on first line
      WRITE (LUNIT, 30, ERR=9898) VERSN, LINLEN, LIW, LDW, LDD, LDE
   30 FORMAT(A18,5I4)
C     Write entity counts, user type count and total space requirements
C     on second line
      WRITE (LUNIT, FORMI) NUMENL, NUMENT, NUMUT, NUMC, NUMI, NUMD

C     Write list of entity pointers selected
      IPTR = IFIRST
      I = 0
C     START LOOP
  310    CONTINUE
         I = I + 1
         LPTR(I) = IPTR
         IF (I .EQ. NI) THEN
            WRITE (LUNIT, FORMIA) (LPTR(IH), '*', IH=1,NI)
            I = 0
         ENDIF
         IHDR = IPTR/256
         IPTR = IMEM(IHDR)
         IF (IPTR .NE. 0) GOTO 310
C     END LOOP
      IF (I .GT. 0)  WRITE (LUNIT, FORMIA) (LPTR(IH), '*', IH=1,I)
C     List names of user types
      CALL DTERPT (0)
      DO 320 I=1,191
         IF (OCCUR(I) .EQ. '*') THEN
            CALL D2FETN (CMEM, IMEM, DMEM, I, TYPNAM, ILIM, IER)
            IF (IER .GE. 0) THEN
               WRITE (LUNIT, 40) I, ILIM, TYPNAM(1:ILIM)
   40          FORMAT (2I3,A)
            ELSE
               WRITE (LUNIT, 40) I, 0
            ENDIF
         ENDIF
  320 CONTINUE
      CALL DTERPT (1)

C     Write out the content of each selected entity
      IPTR = IFIRST
C     START LOOP
  400    CONTINUE
         IHDR = IPTR/256
         ITYP = MOD (ABS(IMEM(IHDR+7)), 256)
         NUMC = ABS (IMEM(IHDR+4))
         NUMI = ABS (IMEM(IHDR+5))
         NUMD = ABS (IMEM(IHDR+6))
C        Write blank line, then header line for entity
         WRITE (LUNIT, '(A)')
         WRITE (LUNIT, FORMIA) IPTR, '*', NUMC, ' ', NUMI, ' ', NUMD,
     +          ' ', ITYP
C        Write its character data
         IF (NUMC .GT. 0) THEN
            ILASTH = IMEM(IHDR+1) + NUMC - 1
            ILIM = 40
            DO 420 IH=IMEM(IHDR+1)-1,ILASTH,40
               ODDCH = .FALSE.
               IF (IH+40 .GT. ILASTH)  ILIM = ILASTH - IH
               DO 410 I=1,ILIM
                  J = ICHAR (CMEM(IH+I:IH+I))
                  IF (J .LT. 32) THEN
                     CH(I:I) = CCH(J+1:J+1)
                     CHX(I:I) = '^'
                     ODDCH = .TRUE.
                  ELSE IF (J .LT. 127) THEN
                     CH(I:I) = CMEM(IH+I:IH+I)
                     CHX(I:I) = ' '
                  ELSE IF (J .LT. 128) THEN
                     CH(I:I) = CCH(33:33)
                     CHX(I:I) = '^'
                     ODDCH = .TRUE.
                  ELSE IF (J .LT. 140) THEN
                     CH(I:I) = CCH(J+1:J+1)
                     CHX(I:I) = '"'
                     ODDCH = .TRUE.
                  ELSE IF (J .LT. 255) THEN
                     CH(I:I) = CHAR (J-128)
                     CHX(I:I) = ''''
                     ODDCH = .TRUE.
                  ELSE
                     CH(I:I) = CCH(33:33)
                     CHX(I:I) = '"'
                     ODDCH = .TRUE.
                  ENDIF
  410          CONTINUE
               IF (ODDCH) THEN
                  IF (ILIM .EQ. 40) THEN
                     WRITE (LUNIT, FORMC) CH, CHX
                  ELSE
                     CH(ILIM+1:40) = ' '
                     WRITE (LUNIT, FORMC) CH, CHX(1:ILIM)
                  ENDIF
               ELSE
                  WRITE (LUNIT, FORMC) CH(1:ILIM)
               ENDIF
  420       CONTINUE
         ENDIF
C        Write its integer data
         IF (NUMI .GT. 0) THEN
C           Set up LP for D0LPLT or SUBROU
            LP(1) = IHDR
            LP(2) = ITYP
            LP(3) = 0
            LP(4) = IMEM(IHDR+2) - 1
            LP(5) = LP(4) + NUMI
            ILASTH = LP(5)
            ILIM = NI
            DO 440 IH=IMEM(IHDR+2)-1,ILASTH-1,NI
               IF (IH+NI .GT. ILASTH)  ILIM = ILASTH - IH
               DO 430 I=1,ILIM
                  LPTR(I) = IMEM(IH+I)
                  IF (IH+I .GT. LP(4)) THEN
C                    Identify next block of pointers
                     IF (ITYP .GE. 192) THEN
                        CALL D0LPLT (CMEM, IMEM, DMEM, LP)
                     ELSE
                        CALL SUBROU (CMEM, IMEM, DMEM, LP)
                     ENDIF
                     IF (IH+I .GT. LP(3) .OR. LP(3) .GT. LP(4)) THEN
C                       SUBROU does not handle this type satisfactorily
                        IF (OCCUR(ITYP) .EQ. '*') THEN
C                          Report each problem type only once
                           CALL DTERR (SUBNAM, 0, -1000-ITYP, 0)
                           IER = IER + 1
                           OCCUR(ITYP) = '?'
                        ENDIF
                        LP(3) = ILASTH + 1
                        LP(4) = LP(3)
                     ENDIF
                  ENDIF
                  IF (IH+I .LT. LP(3)) THEN
                     PFLAG(I) = ' '
                  ELSE
                     PFLAG(I) = '*'
                  ENDIF
  430          CONTINUE
               WRITE (LUNIT, FORMIA) (LPTR(I), PFLAG(I), I=1,ILIM)
  440       CONTINUE
         ENDIF
C        Write its double precision data
         IF (NUMD .GT. 0) THEN
            IH = IMEM(IHDR+3)
            ILASTH = IH + NUMD - 1
            WRITE (LUNIT, FORMD) (DMEM(I),I=IH,ILASTH)
         ENDIF
C        Move to next entity
         IPTR = IMEM(IHDR)
         IF (IPTR .NE. 0) GOTO 400
C     END LOOP
C     Close the file
      CLOSE (LUNIT)
C     Restore the entity headers
      IPTR = IFIRST
C     START LOOP
  500    CONTINUE
         IHDR = IPTR/256
         IFIRST = IMEM(IHDR)
         IMEM(IHDR+7) = -IMEM(IHDR+7)
         IMEM(IHDR) = IPTR
         IPTR = IFIRST
         IF (IPTR .NE. 0) GOTO 500
C     END LOOP
      RETURN

C     Error Handlers

C     Open or write error on file
 9898 CONTINUE
      IER = -2
C     General error report
 9900 CONTINUE
      CALL DTERR (SUBNAM, 1, IER, 0)
      RETURN
      END

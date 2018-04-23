cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      PROGRAM PEAKER
c
c     Read a Tecplot file ('wingslices.dat') containing one or more slices of
c     surface data, and output data corresponding to peak temperature across
c     slices as 'wingpeak.dat'.
c
c     The input file is assumed to begin with a title and a list of names
c     for the quantities appearing in each data line.  The data column that
c     corresponds to the name "T" is used for locating peaks.
c 
c     Each slice is preceded by a ZONE line with zone title, and 3 more lines
c     as shown in the sample.
c
c     The number of slices is unknown (read till EOF), and the number of
c     quantities at each point of a slice is unknown (but constant in a
c     given data file).  This number is determined by counting the names.
c
c     For a given slice, the number of points and the number of two-point
c     segments ("elements") is contained in the file.  All quantities for a
c     given point appear on a single line.  Point indices defining elements
c     follow these lines, one pair per line.  (They are just skipped here.)
c
c     Sample input file ('wingslices.dat'):
c
c     TITLE     = ""
c     VARIABLES = "X"
c     "Y"
c     "Z"
c     "p, N/m^2"
c     "T, K"
c     ZONE T="Slc: Y=4"
c      N=167, E=164, ZONETYPE=FELineSeg
c      DATAPACKING=POINT
c      DT=(SINGLE SINGLE SINGLE SINGLE SINGLE )
c      3.06317863E+01 4.00000000E+00 8.26852035E+00 1.03469858E+01 3.6441080E+02
c      3.09121646E+01 4.00000000E+00 8.23717594E+00 9.67996311E+00 3.6260284E+02
c      3.11978912E+01 4.00000000E+00 8.20183944E+00 9.12536144E+00 3.6022982E+02
c      2.94563331E+01 4.00000000E+00 8.36488151E+00 1.27059392E+01 3.67887451+02
c      :::::::::::::::::::::::::::::
c      3.77420311E+01 4.00000000E+00 7.27871274E+00 1.35319226E+03 7.7381402E+02
c      3.80176048E+01 4.00000000E+00 7.34731006E+00 1.29152636E+03 7.6325775E+02
c      1 7
c      2 1
c      3 2
c      4 10
c      ::::
c     ZONE T="Slc: Y=4.2"
c      N=159, E=156, ZONETYPE=FELineSeg
c      DATAPACKING=POINT
c      DT=(SINGLE SINGLE SINGLE SINGLE SINGLE )
c      3.12043075E+01 4.19999980E+00 8.20732402E+00 8.76233387E+00 3.5224945E+02
c      3.03418598E+01 4.19999980E+00 8.29858684E+00 1.02475318E+01 3.5742605E+02
c      :::::::::::::::::::::::::::::
c      1 5
c      8 6
c      ::::
c
c    Corresponding output file ('wingpeak.dat'):
c
c     TITLE     = ""
c     VARIABLES = "X"
c     "Y"
c     "Z"
c     "p, N/m^2"
c     "T, K"
c     ZONE T="LeadingEdgePeakT"
c     I= 40, J=1, K=1,F=POINT
c     DT=(SINGLE SINGLE SINGLE SINGLE SINGLE )
c      3.77420311E+01 4.00000000E+00 7.27871274E+00 1.35319226E+03 7.7381402E+02
c      3.03418598E+01 4.19999980E+00 8.29858684E+00 1.02475318E+01 3.5742605E+02
c      :::::::::::::::::::::::::::::
c
c
c     History:
c
c     ??/??/??  Dinesh Prahbu  Original implementation.
c     03/30/04  David Saunders Generalized to handle variable numbers  of lines
c                              and columns.
c
c     Origin:   ELORET/NASA Ames Research Center, Moffett Field, CA
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IMPLICIT NONE

c     Constants:

      CHARACTER, PARAMETER :: delimiters * 3 = ' =,'

c     Variables:

      INTEGER :: first, i, icol_t, ios, ipeak, islice, last, mark,
     >           nelements, nheader, npoints, nslices, nvars
      INTEGER :: irow_peak(1)

      REAL, ALLOCATABLE, DIMENSION (:,:) :: slice

      CHARACTER :: format_string * 4, line * 80

c     Execution:

      OPEN (1, FILE='wingslices.dat', STATUS='OLD', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (*, '(/, a)') ' Unable to open wingslices.dat for input.'
         GO TO 999
      END IF

      OPEN (2, FILE='wingpeak.dat', STATUS='UNKNOWN', IOSTAT=ios)
      IF (ios /= 0) THEN
         WRITE (*, '(/, a)') ' Unable to open wingpeak.dat for output.'
         GO TO 999
      END IF

c     Transcribe the title and variable names:

      nheader = 0
      DO ! Until we hit "ZONE"
         READ (1, '(A)') line
         IF (line(2:2) == 'T') THEN ! Save "Temperature" variable #
            icol_t = nheader
         ELSE IF (line(1:1) == 'Z') THEN
            EXIT
         END IF

         nheader = nheader + 1
         WRITE (2, '(A)') line(1:LEN_TRIM(line))
      END DO

      nvars = nheader - 1

      WRITE (2, '(A)') 'ZONE T="LeadingEdgePeakT"'

c     Scan the data twice, (a) to determine how many slices, and (b) to
c     avoid storing more than one (variable-length) slice at a time.

      format_string = '(I?)'
      nslices = 0
      DO ! Until EOF

         READ (1, '(A)', IOSTAT=ios) line ! 'N=nnn, E=mmm, ...'
         IF (ios /= 0) THEN
            WRITE (*, '(/, A, I4)')
     >         ' Read error at start of slice', nslices + 1
            GO TO 999
         END IF

c        Decode the 2nd and 4th tokens:

         first = 1;  last = LEN_TRIM (line)

         CALL SCAN2 (line, delimiters, first, last, mark)
         first = mark + 2
         CALL SCAN2 (line, delimiters, first, last, mark)

         WRITE (format_string(3:3), '(I1)') mark - first + 1
         READ (line(first:mark), format_string, IOSTAT=ios) npoints

         first = mark + 2
         CALL SCAN2 (line, delimiters, first, last, mark)
         first = mark + 2
         CALL SCAN2 (line, delimiters, first, last, mark)

         WRITE (format_string(3:3), '(I1)') mark - first + 1
         READ (line(first:mark), format_string, IOSTAT=ios) nelements

         DO i = 1, 2 + npoints + nelements ! Skip to next ZONE
            READ (1, '(A)', IOSTAT=ios)
            IF (ios /= 0) THEN
               WRITE (*, '(/, (A, I6))')
     >            ' Trouble skipping slice', slice,
     >            ' # points:  ', npoints, ' # elements:', nelements
               GO TO 999
            END IF
         END DO

         READ (1, '(A)', IOSTAT=ios) ! ZONE ... or EOF
         IF (ios < 0) EXIT

         nslices = nslices + 1
      END DO ! Next slice

      REWIND (1)
      DO i = 1, nheader ! Skip
         READ (1, '(A)')
      END DO

      WRITE (2, '(A, I3, A)') 'I=', nslices, ', J=1, K=1, F=POINT'

      DO islice = 1, nslices

         READ (1, '(A)') ! ZONE ...
         READ (1, '(A)') line
         first = 1;  last = LEN_TRIM (line)

         CALL SCAN2 (line, delimiters, first, last, mark)
         first = mark + 2
         CALL SCAN2 (line, delimiters, first, last, mark)

         WRITE (format_string(3:3), '(I1)') mark - first + 1
         READ (line(first:mark), format_string, IOSTAT=ios) npoints

         first = mark + 2
         CALL SCAN2 (line, delimiters, first, last, mark)
         first = mark + 2
         CALL SCAN2 (line, delimiters, first, last, mark)

         WRITE (format_string(3:3), '(I1)') mark - first + 1
         READ (line(first:mark), format_string, IOSTAT=ios) nelements

         READ (1, '(A)') ! DATAPACKING=POINT
         READ (1, '(A)') line ! DT=(SINGLE ...)
         IF (islice == 1) WRITE (2, '(A)') line(1:LEN_TRIM(line))

         ALLOCATE (slice(nvars,npoints))

         READ (1, *, IOSTAT=ios) slice
         IF (ios /= 0) THEN
            WRITE (*, '(/, A, I4)')
     >         ' Trouble reading points of slice', islice
            GO TO 999
         END IF

         DO i = 1, nelements
            READ (1, '(A)') ! Skip
         END DO

         irow_peak = MAXLOC (slice(icol_t,:))
         ipeak     = irow_peak(1)

         WRITE (2, '(1P, 12E16.8)') slice(1:nvars,ipeak)

         DEALLOCATE (slice)

      END DO ! Next slice

  999 CONTINUE

      END PROGRAM PEAKER

C+----------------------------------------------------------------------
C
      SUBROUTINE SCAN2 (STRING, SEPS, FIRST, LAST, MARK)
C
C
C     Description and usage:
C
C           Looks for non-blank fields ("tokens") in a string, where the
C        fields are of arbitrary length and separated by any of a set of
C        user-specified separators (e.g., blanks or commas).  The position
C        of the end of the first token is also returned so that this
C        routine may be conveniently used within a loop to process an
C        entire line of text.
C
C           The procedure examines a substring, STRING (FIRST : LAST), which
C        may of course be the entire string (in which case just call SCAN2
C        with FIRST = 1 and LAST = LEN (STRING) ).  The indices returned
C        are relative to STRING itself, not the substring.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING    *         C    I      Text string containing data to be
C                                        scanned.
C        SEPS      *         C    I      String containing the separators.
C                                        Each character in SEPS counts as a
C                                        token delimiter.
C        FIRST               I    I/O    Index of first character of interest
C                                        in STRING.  FIRST >= 1 is assumed.
C                                        Output is index of the beginning of
C                                        the first non-separator, or 0 if no
C                                        token was found.
C        LAST                I    I/O    Index of last character of interest
C                                        in STRING.  LAST <= LEN (STRING) is
C                                        assumed.  Output is index of the end
C                                        of the last non-separator, or 0 if no
C                                        token was found.
C        MARK                I      O    Points to the end of the first token
C                                        found in the specified portion of
C                                        STRING.  Use FIRST = MARK + 2 for
C                                        further searches in the same string.
C                                        MARK is set to 0 if no token was found.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C        (2)  FIRST >= 1 and LAST <= LEN (STRING) are assumed upon entry,
C             because these are the obvious bounds needed by the calling
C             program on the first scan of a given string - no need to
C             evaluate LEN (STRING) here with every call, while use of
C             MIN and MAX tends to disguise programming errors.
C
C        (3)  Entering with FIRST > LAST is an eventual normal consequence
C             of setting FIRST = MARK + 2 for further scans of the same string.
C             The zero DO loop iteration count feature of the language is
C             taken advantage of here - the result is simply to drop
C             out the bottom with FIRST = LAST = MARK = 0.
C
C        (4)  This version adjusts LAST using a BACKward search, which is
C             degenerate for all tokens after the first.  (Doing it in
C             the one forward loop means unnecessary repeated tokenizing
C             to find the end.)
C
C        (5)  Zeroing all of FIRST, LAST, and MARK (rather than just MARK)
C             when no token is found may, in retrospect, be undesirable in
C             that information is lost (esp. LAST), but there are too many
C             applications to change this now.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C         4 Mar 1986    RAK    Variation of SCANNR, which is hard-coded
C                              (for historical reasons and a modest speed
C                              advantage) for separators BLANK, TAB,
C                              COMMA, COLON, and EQUAL.
C
C         5 May  1988   DAS    Reverse search used to find LAST; MAX, MIN,
C                              and LEN also eliminated (see Notes).
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   FIRST, LAST, MARK
      CHARACTER
     >   SEPS * (*), STRING * (*)

C     Local variables.

      INTEGER
     >   HEAD, I, J, NSEPS, TAIL
      LOGICAL
     >   FOUND


C     Execution.
C     ----------

      NSEPS = LEN (SEPS)
      HEAD = FIRST
      TAIL = LAST

      FIRST = 0
      LAST = 0
      MARK = 0
      FOUND = .FALSE.

C     Look at each character in STRING (HEAD : TAIL) from left to right
C     until a token is found.  (Drop through if LAST > FIRST.)

      DO 30, I = HEAD, TAIL

C        Is the character a separator?

         DO 10, J = 1, NSEPS
            IF (STRING (I : I) .EQ. SEPS (J : J)) GO TO 20
   10    CONTINUE

C           Not a separator.  Check for the beginning of the FIRST token.

            IF (.NOT. FOUND) THEN
               FIRST = I
               FOUND = .TRUE.
            END IF
            GO TO 30

   20    CONTINUE

C           We found a separator.

            IF (FOUND) THEN

C              We just passed the "trailing edge" of the first token.

               MARK = I - 1
               GO TO 40
            END IF
   30 CONTINUE

C     We reached the last character while still seeking the first token.
C     Either this is part of the first token, or MARK and LAST are still zero.

      IF (FOUND) THEN
         MARK = TAIL
         LAST = TAIL
      END IF
      GO TO 99


   40 CONTINUE

C     We found the first token but haven't reached the end yet to adjust LAST.

      DO 60, I = TAIL, MARK + 1, -1

C        Keep searching as long as we have a separator.

         DO 50, J = 1, NSEPS
            IF (STRING (I : I) .EQ. SEPS (J : J)) GO TO 60
   50    CONTINUE

C           Not a separator, so we've found LAST.

            LAST = I
            GO TO 99
   60 CONTINUE

C     Dropped through, so there were no more tokens.

      LAST = MARK

C     Termination.
C     ------------

   99 RETURN
      END

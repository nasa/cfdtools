C+----------------------------------------------------------------------
C
      SUBROUTINE RDXYZ (NDIM, LUNDAT, LUNERR, COLUMNS, MAXPTS,
     >                  NPTS, X, Y, Z, FINIS, IER)
C
C  PURPOSE:
C
C        RDXYZ reads one set of (X), (X,Y), or (X,Y,Z) points per call
C     from the indicated columns of a file which is assumed to be open
C     already on unit LUNDAT.
C
C        For historical reasons, the data may take one of two forms (with
C     or without integer counts, basically, with the latter now preferred):
C
C     (1) "SMOOTH" or "PROFILE" format: one or more "curves" with explicit
C         integer count preceding the real values for each "curve."  "NPTS"
C         is identified as a single leading numeric token, with no decimal
C         point.  The value of "NPTS" should match the number of "good"
C         points.
C
C     (2) "Indefinite" format: no explicit counts of the points in each
C         "curve."  A line with "END" for the first signficant characters
C         serves to separate "curves." 
C
C        In both cases, any title line (normally the first line of the
C     file) should be handled externally by the application.  See programs
C     SMOOTH and PROFILE for sample usages.  Other points to note:
C
C        > The file is assumed to be positioned ready to read the next "curve."
C          A value of NPTS = 0 is valid for data form (1).
C
C        > Blank lines and any strings starting with '!' are ignored, meaning
C          data points may be annotated with trailing comments or commented out
C          altogether.
C
C        > EOF may be encountered on the first (significant) read - it is up
C          to the application to know whether this is normal or not.
C
C        > An EOF encountered after the first significant read (causing FINIS
C          to be set .TRUE.) is normal for data form (2) but will be abnormal
C          for data form (1) - this is the only case of a mismatched "NPTS"
C          trapped.
C
C  ARGUMENTS:
C     ARG      DIM   TYPE  I/O/S   DESCRIPTION
C     NDIM      -      I     I     NDIM = 1, 2, or 3 depending on whether
C                                  (X), (X,Y), or (X,Y,Z) are to be read.
C                                  Pass X, X, X or X, Y, Y if NDIM = 1 or 2.
C     LUNDAT    -      I     I     Logical unit number for file being read.
C     LUNERR    -      I     I     Logical unit number for error messages.
C     COLUMNS  NDIM    I     I     Column number(s) for X (,Y (,Z)).
C     MAXPTS    -      I     I     Maximum number of points in any one
C                                  dataset expected by calling program.
C     NPTS      -      I     O     Number of data points found; NPTS >= 0.
C     X,Y,Z   MAXPTS   R     O     Data points found.  See NDIM.
C     FINIS     -      L     O     .TRUE. means EOF was encountered after
C                                  the first significant read - see above.
C     IER       -      I     O     IER=0 means one dataset found normally;
C                                        NPTS may still be 0.
C                                  IER=1 means EOF encountered on the first
C                                        significant read - may be normal.
C                                        FINIS is set too (but redundant).
C                                  IER=2 means an abnormal error, for which
C                                        RDXYZ will have issued a diagnostic:
C                                        out-of-range NPTS; conversion error;
C                                        or missing Y or Z.
C  PROCEDURES:
C     GETLINE    Reads a line as text; handles suppressed pts. & comments
C     TOKENS     Tokenizes a string
C     UPCASE     Converts a string to upper case
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77 with ! comments, IMPLICIT NONE, and
C                some variable names up to 8 characters.
C
C  HISTORY:
C     07/31/90  DAS  Original implementation, when PROFILE's PRREAD was
C                    found inappropriate for reading 2nd derivatives.
C                    Program SMOOTH can use the same capability, since
C                    the <original> data format handled here is that of
C                    <the original> SMOOTH.
C     12/10/91  DAS  Added multi-column and indefinite-number-of-points
C                    capabilities, mainly for SMOOTH purposes.  Upward
C                    compatibility is awkward!  (Handling all cases of
C                    mismatches between "NPTS" and actual number found
C                    for the original data format would be more trouble
C                    than it is worth, so only EOF with too few points
C                    is trapped and warned of properly.)  Maybe the line
C                    buffer should be in/out (as in READCOLS now), but
C                    the present compromise still avoids backspacing
C                    and should be simpler to reuse.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   NDIM, COLUMNS (NDIM), IER, LUNDAT, LUNERR, MAXPTS, NPTS
      REAL
     >   X (MAXPTS), Y (MAXPTS), Z (MAXPTS)
      LOGICAL
     >   FINIS

C     Local constants:

      INTEGER
     >   LENTOK, MAXCOL
      CHARACTER
     >   COMMENT * 1, IFMT * 8, RFMT * 10
      PARAMETER
     >  (COMMENT = '!',
     >   IFMT    = '(BN,I25)',   ! Should match LENTOK
     >   LENTOK  = 25,           ! Largest # characters in any numeric value
     >   MAXCOL  = 20,           ! Largest column number handled
     >   RFMT    = '(BN,F25.0)') ! Should match LENTOK

C  *  Local variables:

      INTEGER
     >   I, IOS, J, LAST, LASTCOL, NGIVEN, NUMBER
      LOGICAL
     >   HAVELINE
      CHARACTER
     >   LINE * 132, LIST (MAXCOL) * (LENTOK)

C  *  Procedures:

      EXTERNAL
     >   GETLINE, TOKENS, UPCASE

C     Execution:

      IER = 0
      NPTS = 0
      NGIVEN = 0
      FINIS = .FALSE.

C     Look for the start of a "curve."  The original "NPTS" dataset form
C     forces this special treatment at the beginning.

  100 CALL GETLINE (LUNDAT, COMMENT, LINE, LAST, IOS)
      IF (IOS .LT. 0) THEN                      ! Probably a normal EOF
         IER = 1
         FINIS = .TRUE.                         ! Probably redundant
         GO TO 999
      END IF
      IF (IOS .NE. 0) GO TO 900                 ! Read error
      IF (LAST .EQ. 0) GO TO 100                ! Insignificant line

      LASTCOL = COLUMNS (1)
      IF (NDIM .GE. 2) LASTCOL = MAX (LASTCOL, COLUMNS (2))
      IF (NDIM .EQ. 3) LASTCOL = MAX (LASTCOL, COLUMNS (3))

      NUMBER = LASTCOL
      CALL TOKENS (LINE (1 : LAST), NUMBER, LIST)

      HAVELINE = .TRUE.
      IF (NUMBER .EQ. 1) THEN                   ! Must be NPTS unless NDIM = 1

         IF (INDEX (LIST (1), '.') .EQ. 0) THEN    ! No decimal point
            READ (LIST (1), IFMT, ERR=905) NGIVEN     ! Try to read NPTS

            IF (NGIVEN .LT. 0 .OR. NGIVEN .GT. MAXPTS) GO TO 915  ! Bad NPTS
            IF (NGIVEN .EQ. 0) GO TO 999                          ! NPTS = 0 OK
            HAVELINE = .FALSE.            ! Need to get the first data proper
         END IF

      ELSE              ! More than one token found - assumed to be data proper
      END IF

C     Start of loop over the data proper.

      I = 0

  200 CONTINUE
         IF (HAVELINE) THEN    ! Annoyance caused by original "NPTS" case
            HAVELINE = .FALSE.
         ELSE
            CALL GETLINE (LUNDAT, COMMENT, LINE, LAST, IOS)

            IF (IOS .LT. 0) THEN                   ! EOF - probably OK
               FINIS = .TRUE.                      ! Leave IER = 0
               IF (I .LT. NGIVEN) GO TO 925        ! "NPTS" was too big
               GO TO 999
            END IF

            IF (IOS .NE. 0) GO TO 900              ! Read error
            IF (LAST .EQ. 0) GO TO 200             ! Insignificant line

            NUMBER = LASTCOL
            CALL TOKENS (LINE (1 : LAST), NUMBER, LIST)
         END IF

C        Look for end of "curve" (indefinite data format case):

         CALL UPCASE (LIST (1) (1 : 3))
         IF (LIST (1) (1 : 3) .EQ. 'END') GO TO 999

         I = I + 1
         IF (I .GT. MAXPTS) GO TO 915              ! No more room in X/Y/Z
         IF (NUMBER .LT. LASTCOL) GO TO 920        ! Missing value

C        At last we can do the actual (internal) read:

         J = COLUMNS (1)
         READ (LIST (J), RFMT, ERR=910) X (I)
         IF (NDIM .GE. 2) THEN
            J = COLUMNS (2)
            READ (LIST (J), RFMT, ERR=910) Y (I)
            IF (NDIM .EQ. 3) THEN
               J = COLUMNS (3)
               READ (LIST (J), RFMT, ERR=910) Z (I)
            END IF
         END IF
         NPTS = I
         IF (I .NE. NGIVEN)         ! This works for both data formats
     >GO TO 200      

      GO TO 999


C     Error handling:

  900 WRITE (LUNERR, 1001) 'Bad return from GETLINE.'
      GO TO 990

  905 WRITE (LUNERR, 1002) 'Error reading an integer from ', LIST (1)
      GO TO 990

  910 WRITE (LUNERR, 1002) 'Error reading a real value from ', LIST (J)
      GO TO 990

  915 WRITE (LUNERR, 1001) 'Bad number of points found:', NPTS,
     >                     'Maximum provided for:', MAXPTS
      GO TO 990

  920 WRITE (LUNERR, 1001)
     >   'Real value missing from current dataset.  Line #:', I
      GO TO 990

  925 WRITE (LUNERR, 1001)
     >   'Warning: EOF before indicated # pts. found:', NGIVEN,
     >   '         Proceeding with # pts. read:      ', I
      GO TO 999

  990 IER = 2
C*****GO TO 999

  999 RETURN

C     Formats:

 1001 FORMAT (/, ' *** RDXYZ: ', A, I8)
 1002 FORMAT (/, ' *** RDXYZ: ', A, A)

      END

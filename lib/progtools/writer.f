C+----------------------------------------------------------------------
C
      SUBROUTINE WRITER
C
C  Description and usage:
C
C        WRITER is an error message utility modeled after prompting utility
C     READER.  This implementation sends messages to both the logical unit
C     passed by argument (normally a log file) and to Standard Output
C     (logical unit 6, normally the user's console).  If LUNERR = 6, then
C     the message is written only once.  Message text is limited to one
C     line (plus optional numerical values).  Use WRITER repeatedly for
C     longer messages.
C
C        WRITER was first developed for the CALibration subsystem of FMDAS
C     (data acquisition and the display software for the Fluid Mechanics
C     Laboratory at NASA Ames).
C
C        WRITER uses different entry points to handle I*2 and I*4 integers,
C     R*4 and R*8 reals, and C*(*) character values along with the text
C     of the message.
C
C        Typical calls (unrelated) follow:
C
C              CALL WRITI4 (LUNERR, SUBNAM, 1, IER,
C             >   'Premature return from CBREAD.  IER = ')
C
C              CALL WRITR4 (LUNERR, SUBNAM, NCOEFS, COEFS,
C             >   'Coefficients found:')
C
C        Corresponding messages:
C
C              CAL: Premature return from CBREAD.  IER =  6
C              CBXYZ: Coefficients found:
C                1.2345E-01  6.7890E-20 -2.4680E-29
C
C  Arguments:
C
C     Name  Dimension  Type  I/O/S  Description
C     LUNLOG    -        I     I    Logical unit number for output to a
C                                   log file.  Use LUNLOG = 6 to suppress
C                                   log output.  (Messages always go to
C                                   unit 6, and to LUNLOG if it is not 6.)
C     SUBNAM    -      C*(*)   I    Name of routine calling WRITER
C                                   (variable number of characters).
C     NVALS     -        I     I    Number of values to be displayed along
C                                   with the message.  NVALS = 0 is valid.
C                                   NVALS < 0 or NVALS > MXVALS produce a
C                                   warning, with the value 1 used.
C
C     One of the following applies to each entry point (ignored if NVALS=0):
C
C     I2        *       I*2    I    Integer*2 value(s)
C     I4        *       I*4    I    Integer*4 value(s)
C     R4        *       R*4    I    Real*4 value(s)
C     R8        *       R*8    I    Real*8 value(s)
C     L1        *       L*1    I    Logical*1 value(s)
C     L4        *       L*4    I    Logical*4 value(s)
C     CH        *       C*(*)  I    Character value(s)
C
C     TEXT      -       C*(*)  I    Text of message.  Since the
C                                   message should fit on one line,
C                                   and WRITER inserts the calling
C                                   module's name, the use should ensure
C                                   LEN(TEXT) + LEN(SUBNAM) + 2 <= 80
C
C  External files used:
C
C     Unit    I/O/S   Description
C     LUNCRT=6  O     Normally the operator's console.
C     LUNLOG    O     Normally a log file.  See argument description.
C
C  Method:
C
C        Separate ENTRY points are employed for writing INTEGER*2s,
C     INTEGER*4s,  REAL*4s,  REAL*8s, or  CHARACTER strings.  Other
C     entry points can be added if the need arises. As with READER,
C     the top of WRITER is NOT a valid entry point.
C
C        If a single diagnostic value will fit on the same line as the
C     text, this is arranged for, with 80 characters chosen as the
C     limit (not counting carriage control).
C
C  Error Handling:
C
C        The given NVALS is range-checked to avoid voluminous output,
C     but cannot be overwritten because NVALS is not necessarily a true
C     variable in the calling program.  The TEXT argument is truncated
C     if it would exceed one line.
C
C  Notes:
C
C     1.  I2, I4, R4, R8, L1, L4, and CH are not dimensioned with NVALS
C         below, because NVALS may be zero.
C
C  Environment:  DEC VAX/VMS, FORTRAN 77.
C
C  System dependencies:
C     1.  IMPLICIT NONE is nonstandard.
C     2.  END DO avoids a bunch of labels.
C     3.  Variable names longer than six characters (maximum of
C         eight here) are nonstandard.
C     4.  INTEGER*2, LOGICAL*1, etc. may not be supported by some systems.
C
C  Author:  David Saunders, Informatics, Palo Alto, CA.
C
C  History:
C
C     05/16/84  DAS  Initial implementation (READER style; RSX/PDP-11;
C                    output to LUNERR argument only). (No true CHARACTER
C                    variables, so SUBNAM was fixed at 6 characters.)
C     10/08/85  CLH  Adapted for VAX/VMS (no need for end-of-string
C                    delimiter now; changed from INTEGER or REAL to
C                    I*2, I*4, R*4, or R*8).
C     10/13/88  CWB  Provided for message output to Standard Output
C                    (unit 6) and to optional log file (LUNLOG).  (No
C                    way of avoiding use of 6 here, since changing the
C                    calling sequence is out of the question.)  SUBNAM
C                    is now variable length; added WRITCH entry point.
C     08/04/89  CLH  Excess trailing blanks in text and module name
C                    section of output (not value section) suppressed.
C                    If one I*2 or I*4 integer on same line with
C                    text string, format is tailored to integer's
C                    magnitude.  Added LENSUB, LENTXT, SUBLEN, and
C                    TXTLEN.
C     12/31/90  DAS  Handled Logical*1 and Logical*4 values; avoided
C                    nonstandard I<IFMTLEN> in FORMAT 1110 via VARFMT;
C                    commented-out the code that suppresses trailing
C                    blanks from SUBNAM and TEXT, since this prevents
C                    alignment of outputs from separate calls.
C                       
C-----------------------------------------------------------------------

      IMPLICIT  NONE

C ... Arguments:

      INTEGER   LUNLOG, NVALS
      INTEGER*2 I2(*)
      INTEGER*4 I4(*)
      REAL*4    R4(*)
      REAL*8    R8(*)
      LOGICAL*1 L1(*)
      LOGICAL*4 L4(*)
      CHARACTER CH(*)*(*), SUBNAM*(*), TEXT*(*)

C ... Local constants:

      INTEGER   LUNCRT, MXLINE, MXVALS
      CHARACTER ISTAR2*2, ISTAR4*2, RSTAR4*2, RSTAR8*2, LSTAR1*2,
     >          LSTAR4*2, CSTARS*2

      PARAMETER (LUNCRT=6,   MXLINE=80,   MXVALS=100,
     >          ISTAR2='I2', ISTAR4='I4', RSTAR4='R4', RSTAR8='R8',
     >          LSTAR1='L1', LSTAR4='L4', CSTARS='A ')

C ... Local variables:

      INTEGER   I, INC, IFMTLEN, LENSUB, LENTXT, LUN, MSGLEN, NCHARS,
     >          NPRINT, SUBLEN, TXTLEN
      CHARACTER CASE*2, VARFMT*17

C ... Storage:

      DATA      VARFMT /'(1X,A,'': '',A,Inn)'/
C                        1234567 890 1234567

C ... Execution:

      STOP 'Illegal entry to subroutine WRITER.'

C-----------------------------------------------------------------------

      ENTRY WRITI2 (LUNLOG, SUBNAM, NVALS, I2, TEXT)

C-----------------------------------------------------------------------

      CASE = ISTAR2
      GO TO 200

C-----------------------------------------------------------------------

      ENTRY WRITI4 (LUNLOG, SUBNAM, NVALS, I4, TEXT)

C-----------------------------------------------------------------------

      CASE = ISTAR4
      GO TO 200

C-----------------------------------------------------------------------

      ENTRY WRITR4 (LUNLOG, SUBNAM, NVALS, R4, TEXT)

C-----------------------------------------------------------------------

      CASE = RSTAR4
      GO TO 200

C-----------------------------------------------------------------------
         
      ENTRY WRITR8 (LUNLOG, SUBNAM, NVALS, R8, TEXT)

C-----------------------------------------------------------------------

      CASE = RSTAR8
      GO TO 200

C-----------------------------------------------------------------------

      ENTRY WRITL1 (LUNLOG, SUBNAM, NVALS, L1, TEXT)

C-----------------------------------------------------------------------

      CASE = LSTAR1
      GO TO 200

C-----------------------------------------------------------------------
         
      ENTRY WRITL4 (LUNLOG, SUBNAM, NVALS, L4, TEXT)

C-----------------------------------------------------------------------

      CASE = LSTAR4
      GO TO 200

C-----------------------------------------------------------------------

      ENTRY WRITCH (LUNLOG, SUBNAM, NVALS, CH, TEXT)
                                                     
C-----------------------------------------------------------------------

      CASE = CSTARS
C*****GO TO 200

  200 CONTINUE

      INC = LUNLOG - LUNCRT
      IF (INC .EQ. 0) INC = 1

C ... Find the "end" of the text string:

      TXTLEN = LEN (TEXT)
      LENTXT = TXTLEN

CCC      DO I=TXTLEN, 1, -1                   ! NO:  Trailing blanks should
CCC         IF (TEXT(I:I).NE.' ') GO TO 240   !      be allowed for alignment.
CCC         LENTXT = LENTXT - 1
CCC      END DO
CCC
CCC  240 CONTINUE

C ... Find the "end" of the subroutine string:

      SUBLEN = LEN (SUBNAM)
      LENSUB = SUBLEN

CCC      DO I=SUBLEN, 1, -1                     ! See above
CCC         IF (SUBNAM(I:I).NE.' ') GO TO 250
CCC         LENSUB = LENSUB - 1
CCC      END DO
CCC
CCC  250 CONTINUE

C ... Constrain length of message (truncate if needed to fit on one
C     line with SUBNAM and ': '):

      MSGLEN = LENSUB + 2
      NCHARS = MIN (LENTXT, MXLINE - MSGLEN)
      MSGLEN = MSGLEN + NCHARS

      NPRINT = NVALS
      IF (NPRINT .GT. MXVALS  .OR. NPRINT .LT. 0) THEN
         IF (ABS (NVALS) .LT. 10) THEN
            IFMTLEN = 2
         ELSE IF (ABS (NVALS) .LT. 100) THEN
            IFMTLEN = 3
         ELSE IF (ABS (NVALS) .LT. 10000) THEN
            IFMTLEN = 5
         ELSE
            IFMTLEN = 12
         END IF
         WRITE (VARFMT(15:16), '(I2)') IFMTLEN
         DO LUN=LUNCRT, LUNLOG, INC
            WRITE (LUN, VARFMT) SUBNAM(1:LENSUB),
     >         'Warning! Bad NVALS = ', NVALS
         END DO
         NPRINT = 1
      END IF


      IF (NPRINT .EQ. 0) THEN

C ...    No need to worry about appending diagnostic values:

         DO LUN=LUNCRT, LUNLOG, INC
            WRITE (LUN, 1150, ERR=999) SUBNAM(1:LENSUB),
     >         TEXT(1:NCHARS)
         END DO


      ELSE IF (CASE .EQ. ISTAR2) THEN

         IF (ABS (I2(1)) .LT. 10) THEN
            IFMTLEN = 2
         ELSE IF (ABS (I2(1)) .LT. 100) THEN
            IFMTLEN = 3
         ELSE
            IFMTLEN = 6
         END IF
         IF (NPRINT .EQ. 1  .AND.  MSGLEN + IFMTLEN .LE. MXLINE) THEN

C ...       The one integer fits on the same line as the message.

            WRITE (VARFMT(15:16), '(I2)') IFMTLEN
            DO LUN=LUNCRT, LUNLOG, INC
               WRITE (LUN, VARFMT, ERR=999) SUBNAM(1:LENSUB),
     >            TEXT(1:NCHARS), I2(1)
            END DO

         ELSE

C ...       Put all values on the next line:

            DO LUN=LUNCRT, LUNLOG, INC
               WRITE (LUN, 1120, ERR=999) SUBNAM(1:LENSUB),
     >            TEXT(1:NCHARS), (I2(I), I = 1, NPRINT)
            END DO
         END IF


      ELSE IF (CASE .EQ. ISTAR4) THEN

         IF (ABS (I4(1)) .LT. 10) THEN
            IFMTLEN = 2
         ELSE IF (ABS (I4(1)) .LT. 100) THEN
            IFMTLEN = 3
         ELSE IF (ABS (I4(1)) .LT. 1000) THEN
            IFMTLEN = 4
         ELSE
            IFMTLEN = 11
         END IF
         IF (NPRINT .EQ. 1  .AND.  MSGLEN + IFMTLEN .LE. MXLINE) THEN
            WRITE (VARFMT(15:16), '(I2)') IFMTLEN
            DO LUN=LUNCRT, LUNLOG, INC
               WRITE (LUN, VARFMT, ERR=999) SUBNAM(1:LENSUB),
     >            TEXT(1:NCHARS), I4(1)
            END DO
         ELSE
            DO LUN=LUNCRT, LUNLOG, INC
               WRITE (LUN, 1140, ERR=999) SUBNAM(1:LENSUB),
     >            TEXT(1:NCHARS), (I4(I), I = 1, NPRINT)
            END DO
         END IF


      ELSE IF (CASE .EQ. RSTAR4) THEN

         IF (NPRINT .EQ. 1  .AND.  MSGLEN + 12 .LE. MXLINE) THEN
            DO LUN=LUNCRT, LUNLOG, INC
               WRITE (LUN, 1150, ERR=999) SUBNAM(1:LENSUB),
     >            TEXT(1:NCHARS), R4(1)
            END DO
         ELSE
            DO LUN=LUNCRT, LUNLOG, INC
               WRITE (LUN, 1160, ERR=999) SUBNAM(1:LENSUB),
     >            TEXT(1:NCHARS), (R4(I), I = 1, NPRINT)
            END DO
         END IF

                          
      ELSE IF (CASE .EQ. RSTAR8) THEN

         IF (NPRINT .EQ. 1  .AND.  MSGLEN + 23 .LE. MXLINE) THEN
            DO LUN=LUNCRT, LUNLOG, INC
               WRITE (LUN, 1170, ERR=999) SUBNAM(1:LENSUB),
     >            TEXT(1:NCHARS), R8(1)
            END DO
         ELSE
            DO LUN=LUNCRT, LUNLOG, INC
               WRITE (LUN, 1180, ERR=999) SUBNAM(1:LENSUB),
     >            TEXT(1:NCHARS), (R8(I), I = 1, NPRINT)
            END DO
         END IF

      ELSE IF (CASE .EQ. LSTAR1) THEN

         IF (NPRINT .EQ. 1  .AND.  MSGLEN + 1 .LE. MXLINE) THEN
            DO LUN=LUNCRT, LUNLOG, INC
               WRITE (LUN, 1190, ERR=999) SUBNAM(1:LENSUB),
     >            TEXT(1:NCHARS), L1(1)
            END DO
         ELSE
            DO LUN=LUNCRT, LUNLOG, INC
               WRITE (LUN, 1200, ERR=999) SUBNAM(1:LENSUB),
     >            TEXT(1:NCHARS), (L1(I), I = 1, NPRINT)
            END DO
         END IF

      ELSE IF (CASE .EQ. LSTAR4) THEN

         IF (NPRINT .EQ. 1  .AND.  MSGLEN + 1 .LE. MXLINE) THEN
            DO LUN=LUNCRT, LUNLOG, INC
               WRITE (LUN, 1190, ERR=999) SUBNAM(1:LENSUB),
     >            TEXT(1:NCHARS), L4(1)
            END DO
         ELSE
            DO LUN=LUNCRT, LUNLOG, INC
               WRITE (LUN, 1200, ERR=999) SUBNAM(1:LENSUB),
     >            TEXT(1:NCHARS), (L4(I), I = 1, NPRINT)
            END DO
         END IF

      ELSE IF (CASE .EQ. CSTARS) THEN

         IF (NPRINT .EQ. 1  .AND.
     >        MSGLEN + LEN (CH(1)) .LT. MXLINE) THEN     !.LT. allows for 1X.
            DO LUN=LUNCRT, LUNLOG, INC
               WRITE (LUN, 1210, ERR=999) SUBNAM(1:LENSUB),
     >            TEXT(1:NCHARS), CH(1)
            END DO
         ELSE
            DO LUN=LUNCRT, LUNLOG, INC
               WRITE (LUN, 1220, ERR=999) SUBNAM(1:LENSUB),
     >            TEXT(1:NCHARS), (CH(I), I = 1, NPRINT)
            END DO
         END IF

      END IF

  999 RETURN

C ... Formats:

CCCCC 1110 FORMAT (1X, A, ': ', A, I<IFMTLEN>)
 1120 FORMAT (1X, A, ': ', A, /, (1X, 10I8))
 1140 FORMAT (1X, A, ': ', A, /, (1X, 6I12))
 1150 FORMAT (1X, A, ': ', A, 1P, E12.4)
 1160 FORMAT (1X, A, ': ', A, /, (1X, 1P, 6E12.4))
 1170 FORMAT (1X, A, ': ', A, 1P, E23.15)
 1180 FORMAT (1X, A, ': ', A, /, (1X, 1P, 3E23.15))
 1190 FORMAT (1X, A, ': ', A, L1)
 1200 FORMAT (1X, A, ': ', A, /, (40L2))
 1210 FORMAT (1X, A, ': ', A, A)
 1220 FORMAT (1X, A, ': ', A, /, (1X, A))

      END

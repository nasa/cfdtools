C+----------------------------------------------------------------------
C
      SUBROUTINE RDQPL ( MXPTS, TITLE, NPTS, X, Y, LUNRD )
C
C PURPOSE: RDQPL is the complement of SAVQPL (both of which are really
C          obsolete now).  RDQPL reads one curve (with a title) from a
C          file in QPLOT format.  No attempt is made to handle QPLOT's
C          optional namelist. In fact, any text records other than the
C          first record are ignored.  Use RDQPL with discretion, since
C          any curve other than the first in a file is not  likely  to
C          have a useful title associated with  it.
C
C METHOD:  *  Read a title as the first error.  Any error is fatal.
C
C          *  Read (x,y) pairs in a counting loop (list-directed I/O).
C             Use ERR=... to skip unwanted labels and legend info etc.
C             until no error is encountered,  then  keep reading until
C             either an EOF or another error is encountered.  No BACK-
C             SPACE is done at such an error in case it was due to  an
C             END FRAME record (meaning valid title probably follows).
C
C ARGUMENTS:
C
C    ARG     DIM     TYPE I/O/S DESCRIPTION
C  MXPTS      -        I    I   Max. no. of pts. allowed in one curve
C  TITLE      -      C*(*)  O   Title found for curve
C  NPTS       -        I    O   No. of data pts. found on this call
C  X, Y      NPTS      R    O   Coordinates found as one curve.
C  LUNRD      -        I    I   Logical unit number for QPLOT file.
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C AUTHOR: David Saunders, Sterling Software, Palo Alto, CA.  12/05/83
C         <Resurrected 08/25/88 - missed being moved to OBSOLETE lib.>
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

C     Arguments:

      INTEGER   LUNRD, MXPTS, NPTS
      REAL      X(MXPTS), Y(MXPTS)
      CHARACTER TITLE*(*)

C     Execution.

C     Look for a title record:

      READ ( LUNRD, '(A)', ERR=70 ) TITLE

   30 READ ( LUNRD, *, ERR=30, END=80 ) X(1), Y(1)

      NPTS = 1

   40 CONTINUE
         NPTS = NPTS + 1
         IF ( NPTS .GT. MXPTS ) GO TO 90

         READ ( LUNRD, *, ERR=50, END=50 ) X(NPTS), Y(NPTS)
      GO TO 40

   50 NPTS = NPTS - 1
      IF ( NPTS .GE. 1 ) RETURN

      STOP 'RDQPL: No data points found.'
   70 STOP 'RDQPL: Unable to read first record as text.'
   80 STOP 'RDQPL: Unexpected EOF.'
   90 STOP 'RDQPL: Too many data points found.'

      END

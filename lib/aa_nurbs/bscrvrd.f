C+------------------------------------------------------------------------------
C
      SUBROUTINE BSCRVRD (LUNIT, TITLE, MAXC, C, NC, IER)
C
C     Purpose:
C
C           BSCRVRD reads one B-spline curve from a file containing one
C        or more such curves in the "C-array" format written by DTRDWR.
C        The curve may be nonrational or rational, with any number of
C        dimensions.
C
C           The 5-word header of each curve following its title must be
C        readable list-directed as now ensured by DTRDWR's write option.
C        (DTRDWR's read option cannot handle more than one curve or
C        surface per file because it reads to EOF in a single read.)
C
C           The file is assumed to be opened ready for reading the title
C        of the next curve.
C
C     Arguments:
C
C        Name   Dim   Type   Description
C        LUNIT         I     Logical unit number for the file.
C        TITLE       C*(*)   Descriptive title.  The IGES format has room
C                            for up to 72 characters.
C        MAXC          I     Dimension of C-array >= NC.
C        C     MAXC    R     C-array containing DTRC format B-spline data
C        NC            O     Number of meaningful elements in the C-array:
C                            NC = 5 + # knots + # control pts. * NDIM, where
C                            NDIM = 2 for non-rational (X, Y) curve;
C                                 = 3 for rational (X, Y) curve, etc.
C        IER           I     Error flag:
C                             0  no error detected;
C                            -1  EOF encountered on first read (may be OK);
C                            -2  not one independent variable;
C                            -3  NC is greater than the specified
C                                size of the C-array;
C                            -4  bad number or unexpected EOF
C
C     History:
C
C     02/19/93   D.Saunders   Adapted BSCRVRD from DTRDWR after the latter's
C                             write option was revised to make BSCRVRD viable.
C     02/01/94       "        Dispensed with DOUBLE PRECISION.  Compile with
C                             /REAL_LENGTH=64 or -r8 switches as necessary.
C
C     Author:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      CHARACTER
     >   TITLE * (*)
      INTEGER
     >   LUNIT, MAXC, NC, IER
      REAL
     >   C (MAXC)

C     Local variables:

      INTEGER
     >   I, KORDER, NCPTS, NDIM, NKNOTS

C     Execution:


      READ (LUNIT, '(A)', END=800) TITLE

C     Each header must be on separate lines from the rest of the curve:

      READ (LUNIT, *, ERR=888) (C (I), I = 1, 5)

C     Check the header:

      IF (INT (C (1)) .NE. 1) THEN  ! Should be 1 independent variable
         IER = -2
         GO TO 999
      END IF

      NDIM = INT (C (2))            ! No. of dependent variables
      NDIM = ABS (NDIM)             ! Rational or nonrational
      KORDER = INT (C (3))          ! Degree + 1
      NCPTS  = INT (C (4))          ! No. of control points
      NKNOTS = KORDER + NCPTS       ! No. of knots

C     Check size of C-array:

      NC = 5 + NKNOTS + NCPTS * NDIM

      IF (NC .GT. MAXC) THEN
         IER = -3
         GO TO 999
      END IF

C     Read the rest of the curve:

      READ (LUNIT, *, ERR=888) (C (I), I = 6, NC)
      IER = 0
      GO TO 999


  800 IER = -1   ! No more curves
      GO TO 999

  888 IER = -4   ! Bad number encountered or unexpected EOF

  999 RETURN
      END

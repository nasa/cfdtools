C+------------------------------------------------------------------------------
C
      SUBROUTINE DTRDWR (MODE, LUNIT, TITLE, MAXC, C, NC, IER)
C
C     Purpose:
C
C           Read or write a DT-format C-array B-spline curve or surface.
C        There may be only one curve or one surface in the file, which
C        is assumed to be opened ready for (list-directed) I/O.
C
C           This is a straightforward merge of DTRDCA and DTWRCA from the
C        DT_NURBS library.
C
C     Arguments:
C
C        Name   Dim   Type   Description
C        MODE          I     MODE = 'R' means read from the specified file;
C                                 = 'W' means write to it.
C        LUNIT         I     Logical unit number for the file.
C        TITLE       C*(*)   Descriptive title.  The IGES format has room
C                            for up to 72 characters.
C        MAXC          I     Dimension of C-array >= NC.
C        C     MAXC          C-array containing DTRC format B-spline data
C        NC            I     Number of meaningful elements in the C-array:
C                            Curve:   # control pts. * NDIM + # knots + 5 where
C                                     NDIM = 2 for non-rational (X, Y) curve;
C                                          = 3 for rational (X, Y) curve, etc.
C                            Surface: # U ctrl. pts. * # V ctrl. pts. * NDIM +
C                                     # U knots + # V knots + 8
C        IER           I     Error flag:
C                             0  no error detected;
C                            -1  more than 2 independent variables;
C                            -2  incorrect number of dimensions:
C                                second element of C-array > 3 or < -4;
C                            -3  NC is greater than the specified
C                                size of the C-array;
C                            -4  error reading the file.
C
C     History:
C
C        Sep. 1989  C.P.Chi      Original code (DT_NURBS library).
C
C        04/27/92   D.Saunders,  Merged DTRDCA, DTWRCA as DTRDWR.
C                   Sterling/    Clarified NC.  Cleaned up the code
C                   NASA Ames    for use in a design application.
C
C        04/28/92   D.A.S.       Eliminated the file name argument and the
C                                open/close, which were incompatible with
C                                conversion program DT2IGES.
C
C        05/08/92   "  "         Introduced the TITLE argument.
C
C        02/19/93   "  "         Made written curves more readable, at the
C                                expense of more than one WRITE.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      CHARACTER
     >   MODE * 1, TITLE * (*)
      INTEGER
     >   LUNIT, MAXC, NC, IER
      DOUBLE PRECISION
     >   C (MAXC)

C     Local variables:

      INTEGER
     >   I, NDEGL, NDEGP, NDIM, NDV, NIV, NKU, NKV, NLNS, NPTS
      LOGICAL
     >   READ

C     Execution:

      IER = 0
      READ = MODE .EQ. 'R'

      IF (READ) THEN
         READ (LUNIT, '(A)', ERR=888) TITLE
         READ (LUNIT, *, ERR=888, END=100) C
      END IF

C     Analyze C-array.

  100 NIV = INT (C (1))             ! No. of independent variables

      IF (NIV .LT. 1 .OR. NIV .GT. 2) THEN
         IER = -1
         GO TO 999
      END IF

      NDV = INT (C (2))             ! No. of dependent variables

      IF (NDV .GT. 3 .OR. NDV .LT. -4) THEN
         IER = -2
         GO TO 999
      END IF

      NDIM = ABS (NDV)              ! Rational or nonrational

      IF (NIV .EQ. 1) THEN          ! Curve

         NDEGP = INT (C (3)) - 1    ! Degree
         NPTS  = INT (C (4))        ! No. of control points
         NKU   = NDEGP + NPTS + 1   ! No. of knots

C ...    Check size of C-array

         NC = NPTS * NDIM + NKU + 5

         IF (NC .GT. MAXC) THEN
            IER = -3
            GO TO 999
         END IF

      ELSE                          ! Surface

         NDEGL = INT (C (3)) - 1    ! Degrees
         NDEGP = INT (C (4)) - 1
         NLNS  = INT (C (5))        ! Control points
         NPTS  = INT (C (6))
         NKU   = NDEGL + NLNS + 1   ! Knots
         NKV   = NDEGP + NPTS + 1

         NC = (NLNS * NPTS) * NDIM + NKU + NKV + 8

         IF (NC .GT. MAXC) THEN
            IER = -3
            GO TO 999
         END IF

      END IF

      IF (.NOT. READ) THEN
         WRITE (LUNIT, '(A)') TITLE

         IF (NIV .EQ. 1) THEN          ! Make a curve more readable
            WRITE (LUNIT, *) (C (I), I = 1, 5)
            WRITE (LUNIT, *)
            WRITE (LUNIT, *) (C (I), I = 6, 5 + NKU)
            WRITE (LUNIT, *)
            WRITE (LUNIT, *) (C (I), I = 6 + NKU, 5 + NKU + NPTS)
            WRITE (LUNIT, *)
            WRITE (LUNIT, *) (C (I), I = 6 + NKU + NPTS, NC)
         ELSE
            WRITE (LUNIT, *) (C (I), I = 1, NC)
         END IF
      END IF

      GO TO 999

  888 IER = -4

  999 RETURN
      END

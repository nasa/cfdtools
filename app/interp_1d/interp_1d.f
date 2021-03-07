C+----------------------------------------------------------------------
C
      PROGRAM INTERP_1D
C
C
C     PURPOSE:
C
C        INTERP_1D provides for interactive interpolating/extrapolating
C     within a monotonic set of data points (at least two points).   In
C     particular,  it  answers the common case of interpolating a third
C     point from a given pair of points,  since  using  a calculator is
C     error-prone.
C
C        Both linear and local cubic spline interpolation are  provided
C     here, possibly for comparison purposes.  64-bit arithmetic may be
C     advisable (as for careful capture of the end of a trajectory).
C
C        Also calculated are the straight line slope and intercept  for
C     the appropriate interval in the table,  along with the derivative
C     from the spline fit.
C
C
C     METHOD:
C
C        After the data file  ("table")  and spline method are entered,
C     the inner loop over abscissas of interest prompts for each target
C     abscissa and displays the interpolated (or extrapolated) results.
C
C
C     DATA FORMAT:
C
C        (X, Y) data in two columns, read until EOF.
C
C     EXTERNAL REFERENCES:
C
C        TABLE1    Linear interpolation table look-up utility which
C                  returns the relevant interval - otherwise, we could
C                  use the linear option in LCSFIT.
C        LCSFIT    Local cubic spline utility (4-point method).
C
C
C     HISTORY:
C
C     DAS   10/21/83   Test program for TABLE1.
C     DAS   08/08/86   Slight clean up as INTERP1D.
C     DAS   04/07/87   More tidy up for publication.
C   DAS/RAK 03/04/88   Slope and intercept displayed now; table now
C                      entered as an indefinite list of pairs with
C                      redirection to a file permitted via RDTUPLES.
C     DAS   05/23/02   INTERP_1D version intended for higher precision
C                      data, prompted by Traj_opt gradient difficulties.
C                      Read the data from a file instead of prompting,
C                      but prompt for target abscissas.  See SMOOTH
C                      for much more complete handling of 32-bit data.
C
C     AUTHOR:  David Saunders, ELORET/NASA Ames Research Center.
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

C     Local constants:

      INTEGER, PARAMETER ::
     >   LUNCRT = 6, LUNDAT = 1, LUNKBD = 5, MXTABLE = 500

C     Local variables:

      INTEGER
     >   INDEX, IOS, NTABLE

      REAL
     >   SLOPE, X, XTABLE (MXTABLE), Y, YINTERCEPT, YTABLE (MXTABLE),
     >   YSPLINE, YPSLINE

      LOGICAL
     >   CR, EOF

      CHARACTER
     >   FILENAME * 64, METHOD * 1

C     Procedures:

      REAL
     >   TABLE1

      EXTERNAL
     >   LCSFIT, TABLE1

C     Execution:

      WRITE (LUNCRT, '(/, 2A)')
     >   ' Interactive 1-D interpolation or extrapolation. ',
     >   ' ^D quits loop over target Xs.'

      WRITE (LUNCRT, '(/, A)', ADVANCE='NO')
     >   ' LCSFIT method (M=mono, B=loose, L=linear): '
      READ  (LUNKBD, *) METHOD

      IF (METHOD == 'm') METHOD = 'M'
      IF (METHOD == 'b') METHOD = 'B'
      IF (METHOD == 'l') METHOD = 'L'

      WRITE (LUNCRT, '(/, A)', ADVANCE='NO') ' Data file name: '
      READ  (LUNKBD, '(A)') FILENAME

      OPEN  (LUNDAT, FILE=FILENAME, STATUS='OLD')

      NTABLE = 0
      DO INDEX = 1, MXTABLE
         READ (LUNDAT, *, IOSTAT=IOS) XTABLE(INDEX), YTABLE(INDEX)
         IF (IOS /= 0) THEN
            CLOSE (LUNDAT)
            WRITE (LUNCRT, '(A, I4)') ' # data points found: ', NTABLE
            EXIT
         END IF
         NTABLE = NTABLE + 1
      END DO

      WRITE (LUNCRT, '(/)')
      INDEX = 1

      DO ! Until EOF
         WRITE (LUNCRT, '(A)', ADVANCE='NO') ' Enter X: '
         READ  (LUNKBD, *, IOSTAT=IOS) X
         IF (IOS /= 0) EXIT

         Y = TABLE1 (NTABLE, XTABLE, YTABLE, INDEX, X, IOS)

         SLOPE = (YTABLE (INDEX + 1) - YTABLE (INDEX)) /
     >           (XTABLE (INDEX + 1) - XTABLE (INDEX))
         YINTERCEPT = YTABLE (INDEX) - SLOPE * XTABLE (INDEX)

         CALL LCSFIT (NTABLE, XTABLE, YTABLE, .TRUE., METHOD,
     >                1, X, YSPLINE, YPSLINE)

         WRITE (LUNCRT, '(/, (A, 1P, 2E25.15, A, E25.15))')
     >      ' Linear Y, Y'': ', Y, SLOPE,
     >      '     Intercept: ', YINTERCEPT,
     >      ' Spline Y, Y'': ', YSPLINE, YPSLINE

      END DO ! Another target abscissa?

      END PROGRAM INTERP_1D

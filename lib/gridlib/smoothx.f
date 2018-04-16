C+------------------------------------------------------------------------------
C
      subroutine SMOOTHX (xinout, NumX, XFixed, NumFixed, Interval,
     &                    work, status)
C
C Purpose:
C     Given an initial guess at an abscissa distribution in xinout(), and
C     another list with the abscissas that are fixed, XFixed(), {both must
C     be monotonically increasing} generate a new list of abscissas which 
C     have a smoothly varying distribution, with the fixed points unchanged and 
C     the points within Interval points either side of the fixed points moved
C     as necessary to get smoothness.  The first and last original points are
C     not moved, and do not need to appear in the list of fixed points.
C
C     Program ADJUST is available to experiment with SMOOTHX's behavior.
C
C Method:
C     Merge the fixed points into the input list, and blank out the input
C     points on either side of each fixed point, using Interval to 
C     determine how many points to blank out.  Fit a monotone spline through
C     this blanked-out list (X vs. I), and evaluate it at the blanked out
C     points.  The blanking-out allows the curve to change shape around the
C     fixed points when it gets splined.  If the curve comes out jagged, or
C     irregular, then Interval should be made larger.
C
C Arguments:
C     NumX        Integer   in:  length of xinout vector
C     NumFixed       "      in:  length of XFixed vector
C     Interval       "      in:  # of elements to blank either side of fp's
C     status         "      out: completion status (=0 for success)
C
C     xinout(*)    Real     in/out:  abscissa values
C     XFixed(*)      "      in:      abscissae that can't be changed
C     work(NumX*2)   "      scratch: used for the curve to be fit
C
C Error handling:
C     If the end points of XFixed() are outside the range of xinout(), then 
C     abort with status = 1 for lower endpoint, = 2 for upper endpoint.
C
C Notes:
C     Indices are real because they're used as abscissas to the spline 
C     routines.
C     If two fixed abscissas are close to the same point in the input array,
C     the index of the second will correspond to the next higher input point,
C     even if it is further away from that point.
C
C Procedures:
C     LCSFIT  Local cubic spline with monotone option
C   
C System dependencies:
C     Lower case characters in code is non-ANSI 77 standard
C     "!" comments and names longer than 6 characters are also non-standard
C
C Coding style:
C     All upper case symbols are EXTERNAL references
C     Mixed case symbols are constants (either parameter or input arguments)
C     All lower case symbols are read/write variables
C 
C Author:  David B. Serafini, Sterling Software/NASA Ames Research Center
C
C History:
C   3/86   <dbs>    Initial design and coding for GRID3D program
C   3/86   <dbs>    Added endpoint slope calculation because CSFIT was
C                   generating non-monotonically increasing curves for
C                   the x-values.
C   2/87   <dbs>    Changed to MSFIT for spline fit, and eliminated
C                   monotonicity correction.
C   2/87   <dbs>    Changed from fitting only the fixed points to
C                   blanking-out the input points in the neighborhood
C                   of the fp's and fitting all the rest.
C   9/87 D.Saunders Forced inclusion of last point (NumX) no matter what.
C                   Otherwise, a fixed point too near the high end means
C                   the splined points don't cover the data range.
C   2/96     "      LCSFIT replaces MSFIT/CSEVAL. Regularized some of the style.
C  11/98     "      A new compiler accessed XFixed(NumFixed+1) in the following:
C                      if (fptr <= NumFixed .AND. xinout(iptr) >= XFixed(fptr))
C                   so the test is safer now.
C-------------------------------------------------------------------------------

      implicit none

C     >>arguments
      integer NumX,                 !in:  length of xinout vector
     &        NumFixed,             !in:  length of XFixed vector
     &        Interval,             !in:  # of elements to blank around fp's
     &        status                !out: completion status (=0 for success)

      real    xinout(*),            !in/out:  abscissa values
     &        XFixed(*),            !in:      abscissae that can't be changed
     &        work(NumX,2)          !scratch: used for curve to be fit

C     >>local constants
      integer XFit, IFit
      parameter (XFit=1, IFit=2)

C     >>local variables
      integer iptr,                 !points to current value in xinout()
     &        optr,                 !points to current space in work(n,XFIT)
     &        fptr,                 !points to current value in XFixed()
     &        lastfp,               !value of iptr of last XFixed() merged
     &        nfit                  !number of points on curve to spline

      real    cutoff                !index of point before a fixed point
                                    !that marks beginning of blank interval
      real    deriv                 !play safe with LCSFIT
      logical check                 !avoids a compound if test

C     >>procedures
      external  LCSFIT
      intrinsic REAL


C-----------------------------------------------------------------------


C     >>execution

C     >>check endpoints
      if (XFixed(1) .LT. xinout(1)) then
          status = 1
          go to 999
      else if (XFixed(NumFixed) .GT. xinout(NumX)) then
          status = 2
          go to 999
      end if

C     >>initialize
      iptr = 2
      optr = 2
      fptr = 1
      lastfp = -Interval
      work(1,XFit) = xinout(1) ! Always include first point in output
      work(1,IFit) = 1.

C     >>LOOP over input points (iptr) starting at 2
   10 continue

C         >>check to see if there is a fixed point between this input point
C           and the previous one

          check = fptr .LE. NumFixed
          if (check) check = xinout(iptr) .GE. XFixed(fptr)

          if (check) then

C             >>if the fixed point is closer to the previous input point,
C               set iptr back one, except when the previous point was fixed

              if (XFixed(fptr) - xinout(iptr-1) .LE.
     &            xinout(iptr) - XFixed(fptr)
     &            .AND. lastfp .LT. iptr - 1)   iptr = iptr - 1

C             >>blank-out the output arrays back to the beginning of an
C               interval around the fp, or back to the previous fp, or the
C               lower endpoint, whichever is closer

              cutoff = REAL (MAX (lastfp, iptr-Interval-1, 0))

C             >>LOOP backward from optr until the cutoff point is reached

   20         if (optr .GT. 1) then
                  if (work(optr-1,IFit) .GT. cutoff) then
                      optr = optr - 1
                      go to 20
                  end if
              end if

C             { work(optr-1,IFit) .EQ. cutoff .OR. optr .EQ. 1 }

C             >>put the fixed point in output, and reset pointers

              work(optr,XFit) = XFixed(fptr)
              work(optr,IFit) = REAL (iptr)

              optr = optr + 1
              fptr = fptr + 1
              lastfp = iptr

          else !{ xinout(iptr) .LT. XFixed(fptr) .OR. fptr .GT. NumFixed }
C             >>check lastfp to see if we're in an interval after a fp,
C               and if not, store the input point in output arrays.
C               Be sure to include very last point so spline covers full range.

              if (iptr .GT. lastfp+Interval .OR. iptr .EQ. NumX) then
                  work(optr,XFit) = xinout(iptr)
                  work(optr,IFit) = REAL (iptr)
                  optr = optr + 1

C             { else } !do nothing, skip this input point

              end if

          end if

C     >>ENDLOOP:
      iptr = iptr + 1
      if (iptr .LE. NumX) go to 10

C     >>decrement optr so it equals the number of points stored
      nfit = optr - 1

C     >>monotonic local spline interpolation at non-fit points
      optr = 1
      do 30 iptr = 1, NumX

          if (INT (work(optr,IFit)) .EQ. iptr) then
C             >>this index was fit, so don't evaluate

              xinout(iptr) = work(optr,XFit)
              optr = optr + 1

          else
C             >>not fit, so interpolate

              call LCSFIT (nfit, work(1,IFit), work(1,XFit), .TRUE.,
     &                     'M', 1, REAL (iptr), xinout(iptr), deriv)

          end if

   30 continue

C     >>done
      status = 0

 999  return

      end subroutine SMOOTHX

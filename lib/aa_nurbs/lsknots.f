C+------------------------------------------------------------------------------
C
      SUBROUTINE LSKNOTS (NDEG, NDATA, NCPTS, UBAR, NU, U, IER)
C
C  ONE-LINER: Generates interior knots for a least squares B-spline curve fit
C
C  DESCRIPTION:
C
C        LSKNOTS derives a set of interior knots from the coordinates of NDATA
C     discrete points defining a curve, as needed for least squares (or
C     interpolatory) B-spline curve fitting with the indicated number of
C     control points and the indicated degree.  (The first and last NDEG + 1
C     knots are typically established at a higher level to ensure matching
C     at the end points.  If an interior point among the data points is also
C     to be constrained, this routine could be called once for each subset
C     of data points; the interior knot multiplicity and the distribution of
C     control points among the subcurves would be handled at the higher level.)
C
C        The algorithm for the interpolation case (NCPTS = NDATA) appears on
C     p.160 of lecture notes by Dr. Wayne Tiller, "Theory and Implementation
C     of NURBS Curves and Surfaces" (c. 1992).  It is outlined as follows:
C
C        If UBAR(J), J = 1, NDATA are the cumulative chord lengths associated
C     with the data points, D = UBAR(NDATA), and K = NDEG + 1 is the order of
C     the B-spline of interest, the knot vector for the basic, interpolatory
C     case is:
C
C        U = {0, ..., 0, U(1), ..., U(NDATA-K), D, ..., D}
C        ~        K                                K
C
C     where U(J) is the average of the NDEG values UBAR(J+1:J+NDEG).
C
C        For the least squares case (NCPTS < NDATA), the knots are chosen so
C     that each knot span contains roughly the same number of UBAR(*) values.
C     The number of interior knots is NU = NCPTS - K in either case.  The
C     heuristic algorithm here involves averaging NU groups of UBAR(*) values
C     of size 2 * (NDATA / (NU + 1)) which overlap by half.  The excess, M,
C     from this truncated divide is distributed among the first M interior
C     knots, where M <= NU, rather than having all of them absorbed by
C     the last group averaged (as was done originally).
C     
C
C  ARGUMENTS:
C
C     Name    Type/Dimension  I/O/S  Description
C
C     NDEG      I             I      Degree of intended B-spline. NDEG >= 1.
C
C     NDATA     I             I      Number of data points. NDATA >= NDEG + 1.
C
C     NCPTS     I             I      Number of control points (including
C                                    end-points).  NDEG + 1 <= NCPTS <= NDATA.
C
C     UBAR      R (NDATA)     I      Cumulative chord lengths associated with
C                                    the data points (possibly normalized).
C
C     NU        I               O    Output as the number of interior knots.
C                                    NU = NCPTS - (NDEG + 1) >= 0.
C
C     U         R (NU)          O    Desired interior knots.  (May be empty.)
C
C     IER       I               O    0  means normal return;
C                                    -1 means NCPTS > NDATA;
C                                    -2 means NCPTS < NDEG + 1.
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77 with minor extensions:
C                IMPLICIT NONE
C                Variables up to 8 characters
C  HISTORY:
C
C     13 Mar. 1992   DAS   Initial implementation, to facilitate use of
C                          the DTNURBS library's curve/surface fitting.
C     18 Mar. 1992   DAS   Averaging NDEG + (NDATA - NCPTS) of the UBAR(*)
C                          values (by analogy with the interpolation case)
C                          is not good: the knots bunch nearer the middle
C                          as NCPTS is reduced.
C     14 Sep. 1992   DAS   The last knot originally absorbed all the slack
C                          from the truncated divide.  Now, one extra point
C                          is added to each group averaged until all the
C                          remainders have been absorbed.  This distributes
C                          the knots a little better (but is still not as
C                          symmetric as one would like, and the case of
C                          group size 2 requires special treatment).
C                          Unusually irregular spacing of data points can
C                          still mean loss of interlacing (knot spans not
C                          containing at least one data point), especially
C                          if the first or last data interval is too large.
C     01 Feb. 1994   DAS   Dispensed with DOUBLE PRECISION.  Compile with
C                          /REAL_LENGTH=64 or -r8 switches as necessary.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   NCPTS, NDATA, NDEG, NU, IER
      REAL
     &   U (*), UBAR (NDATA)

C     Local constants.

      REAL
     &   ZERO, ONE
      PARAMETER
     &  (ZERO = 0.E+0,
     &   ONE  = 1.E+0)

C     Local variables.

      INTEGER
     &   I, I1, I2, J, K, NAVG, NHALF, NREMAIN
      REAL
     &   RN, RNAVG, RNAVGP1, SUM

C     Execution.

      K = NDEG + 1 

      IF (NCPTS .LT. K) THEN
         IER = -2
         GO TO 99
      ELSE IF (NCPTS .GT. NDATA) THEN
         IER = -1
         GO TO 99
      ELSE
         IER = 0
      END IF

      NU = NCPTS - K

      IF (NCPTS .EQ. NDATA) THEN   ! Interpolation case

         NAVG = NDEG
         RNAVG = ONE / DBLE (NAVG)

         DO 30, J = 1, NU
            SUM = ZERO
            DO 20, I = J + 1, J + NAVG
               SUM = UBAR (I) + SUM
   20       CONTINUE
            U (J) = SUM * RNAVG
   30    CONTINUE

      ELSE                         ! Least squares case

         NHALF = NDATA / (NU + 1)  ! Half the number in each group averaged,
         NREMAIN = MOD (NDATA, NU + 1)  ! if the divide is exact;
         NAVG = 2 * NHALF               ! NAVG may be 1 greater if not.
         RNAVG = ONE / DBLE (NAVG)
         RNAVGP1 = ONE / DBLE (NAVG + 1)

         I1 = 0
         DO 50, J = 1, NU
            SUM = ZERO
            I2 = I1 + NAVG

            IF (J .LE. NREMAIN) THEN ! Fudge the size of the group averaged
               I2 = I2 + 1           ! to spread any excess among the first
               RN = RNAVGP1          ! so many knots.
            ELSE
               RN = RNAVG
            END IF

            DO 40, I = I1 + 1, I2
               SUM = UBAR (I) + SUM
   40       CONTINUE

            U (J) = SUM * RN
            I1 = I2 - NHALF
   50    CONTINUE

C        The following avoids averaging only two at the right end, which
C        would give no data point between the last interior knot and the
C        end-point knot.  Other losses of interlacing are still possible.

         IF (NHALF .EQ. 1) THEN   ! Redo the last few (backwards)
            I1 = NDATA
            DO 70, J = NU, NREMAIN + 1, -1
               SUM = ZERO
               DO 60, I = I1, I1 - 2
                  SUM = UBAR (I) + SUM
   60          CONTINUE
               U (J) = SUM * RNAVGP1
               I1 = I1 - 1
   70       CONTINUE
         END IF
      END IF

   99 RETURN
      END

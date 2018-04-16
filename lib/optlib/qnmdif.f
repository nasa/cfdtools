*DECK QNMDIF
C+----------------------------------------------------------------------
C
      SUBROUTINE QNMDIF (N, NLDIM, NFTOTL, NITER, NTYPE, LWRITE,
     >   X, F, FM, G, H, EL, D, ETA, TOL, STEPMX, EPSMCH, EPSOBJ,
     >   UNITL, LOCAL, CONV, CONTIN, RESCUE, PRINT, FUN, USER)
C
C
C     Description and usage
C     ---------------------
C
C        QNMDIF uses a quasi-Newton method with difference
C     approximations to derivatives to minimize a function F(X)
C     of the N independent variables X1,X2,...,XN.  The method is
C     iterative, and so requires an initial estimate of the position
C     of the minimum.  It is intended for functions which are 
C     continuous and have continuous first and second derivatives in
C     the region between the initial estimate and the minimum, although
C     it may work if the derivatives have occasional discontinuities.
C
C        The user is required to provide a subroutine to evaluate
C     the function to be minimized.  The first derivatives are
C     calculated by differencing function values.
C
C        QNMDIF is a revised version of the program given in the NPL
C     report by Gill, Murray, and Pitfield, "The Implementation of
C     Two Revised Quasi-Newton Algorithms for Unconstrained
C     Optimization", National Physical Laboratory, Div. of Numerical
C     Analysis and Computing, Report No. DNAC 11, 1972.
C
C        From an initial point it generates a sequence of points which
C     is intended to converge to a local minimum of F(X).  A typical
C     iteration starts at the point X, with an estimate G of the
C     gradient vector, and a lower triangular matrix EL and a diagonal
C     matrix D such that EL * D * ELT is a positive definite approximation
C     to the current matrix of second derivatives (Hessian matrix).
C     (note that ELT denotes the transpose of EL).  The linear equations
C
C              EL * D * ELT  *  P = -G
C
C     are solved to give an estimate P of the direction of search.
C     A scalar ALPHA is found so that X + ALPHA * P is approximately a
C     minimum with respect to ALPHA, and then X is replaced by
C     X + ALPHA * P.  The matrices EL and D are then updated according
C     to the complementary Davidon-Fletcher-Powell update formula.
C     Most iterations use forward difference approximations to compute
C     G, but central differences are used when they seem necessary.
C
C        QNMDIF is the main subroutine in a set of routines that
C     includes the following six subroutines:  QNMDIF, LINSCH, OUTPUT,
C     APROXG, UPCHOL, and C1MDCH.  In addition, the user must supply
C     a subroutine called USER which may be a dummy - see below.
CRAK  A second group of subroutines comprised of CENDIF, FDSTEP, and
CRAK  DIFFER may be used to estimate good finite difference step sizes.
C
C        The subroutines QNMDIF, LINSCH, and UPCHOL contain locally
C     dimensioned arrays, which must all be of length at least N.  In
C     the current version, these arrays are of length MAXN, a PARAMETER.
C     QNMDIF and companion routine CENDIF both check that N <= MAXN.
C
CRAK     NOTE: This particular version of QNMDIF has been revised for
CRAK  expensive objective functions which cannot be evaluated to full
CRAK  machine precision.
C
C
C
C     Parameters
C     ----------
C
C     N - The number of independent variables, an integer constant.
C        Present maximum value is 40, set by parameter MAXN, below.
C
C     NLDIM - The value N * (N-1)/ 2, which is the number of elements
C        in the matrix EL, stored as a vector.  Note, however, that
C        NLDIM must always be at least 1, so that if QNMDIF is called
C        with N = 1, NLDIM should be set to 1.
C
C     NFTOTL - An integer, which is set by QNMDIF on exit to the total
C        number of calls to the objective function executed during the
C        location of the minimum.  It need not be initialized.
CRAK
CRAK     MODIFICATION NOTE:  NFTOTL must now be initialized by the
CRAK     user, even if only to zero.  When F, G, EL, D are supplied
CRAK     (i.e. UNITL = .FALSE.), NFTOTL may be initialized to reflect
CRAK     the work done prior to calling QNMDIF.
CRAK
C
C     NITER - An integer, which is set by QNMDIF to the number of
C        iterations executed during the location of the minimum,
C        where an iteration is the procedure described in the comments
C        at the beginning of QNMDIF.  It need not be initialized.
CRAK
CRAK     MODIFICATION NOTE:  NITER may now be used to limit the number
CRAK     of opimization iterations to be performed by QNMDIF.  Set
CRAK     NITER to a large positive value in calling routine to disable
CRAK     the test (and hence revert to original mode).
CRAK
CRAK
CRAK  NTYPE - An integer flag set by QNMDIF to indicate whether the
CRAK     gradient evaluated before returning used forward (NTYPE = 0)
CRAK     or central differences (NTYPE = 1).  On input, must be set to
CRAK     correspond to the user-supplied gradient if UNITL = .FALSE.
CRAK     or to the desired initial gradient type if starting from
CRAK     scratch (UNITL = .TRUE.).
CRAK
CRAK
CRAK  LWRITE - Required integer input denoting the unit number of
CRAK     the primary output device (often device 6).
CRAK
C
C     X - The array of independent variables, a vector of length N.
C        On entry, it should be set to the initial estimate of the
C        location of the minimum.  On output it will contain the best
C        estimate found of the minimum.
C
CRAK
CRAK  F - A real scalar.  On entry, it must be set to initial function
CRAK     value.  On exit it will be the lowest value found of the
CRAK     objective function.
CRAK
CRAK
CRAK  FM - A real scalar.  On entry, it should be set to a lower
CRAK     bound on the value of the function at the minimum, or to
CRAK     a large negative value if no estimate is available.
CRAK
CRAK
CRAK  G - A real array of length N, containing the gradient, which
CRAK     must be initialized if UNITL = .FALSE., and otherwise not.
CRAK     On exit, equals the last gradient evaluated, which may be
CRAK     useful.
CRAK
C
C     H - A difference array, of length N, where H (I) is the interval
C        used to compute the partial derivative of F with respect to
C        X (I).  Great care should be used in choosing the H array, and
C        the user should consider that the smaller the value of H, the
C        less error is incurred in computing the derivative when exact
C        arithmetic is used, since the error in the forward difference
C        approximation is of order H (using the Taylor series expansion);
C        on the other hand, the smaller H is, the more cancellation
C        error will occur, since values of the function computed at close
C        points will have many matching figures, which results in a loss
C        of accuracy in the computed derivative when these close values
C        are subtracted.  Each H (I) should be chosen small enough for the
C        difference estimates to be close to the true derivatives, but
C        not so small that cancellation error destroys the accuracy.  It
C        is suggested that the user experiment with different values of
C        H (I) before using this subroutine, to see how the function
C        values behave for various sizes of H (I).
C
C
C     EL - A vector of length NLDIM (at least 1), used to hold the
C        lower triangular factor of the approximate Hessian.  If the
C        user wishes to provide an initial estimate of the Hessian, EL
C        should be filled with the elements of the lower triangle of the
C        Cholesky factorization EL * D * ELT, where ELT is EL tranpose,
C        stored by rows.  Otherwise, set UNITL to .TRUE., and QNMDIF will
C        use the identity matrix as the initial approximation (i.e., the
C        first iteration will be steepest descent).  On exit, EL contains
C        the lower triangle of the factors of the final approximate Hessian.
C
C
C     D - A vector of length N, used to hold the diagonal matrix of the
C        factors of the approximate Hessian.  If the user wishes to
C        provide an initial estimate of the Hessian, D should be set in
C        the same way as EL (described above).  Otherwise, when UNITL is
C        set to .TRUE., the identity is used for the first iteration, and
C        the user need not set EL and D. On exit, D contains the diagonal
C        factor of the final approximate Hessian.
C
C
C     ETA - A scalar constant, set by the user to a number strictly
C        between 0 and 1.  It controls the accuracy of the search for the
C        minimum of the function along the current direction of search,
C        in the following way:  if the exact minimum along the line were
C        found at each step, the projected gradient (gradient of the
C        function along the line) would be zero.  However, it has been
C        found that finding the exact minimum at each step is not
C        computationally efficient, in that it requires many extra
C        function evaluations, and does not in general assist the speed
C        of convergence dramatically.  The parameter ETA is multipled by
C        the projected gradient at the start of each minimization along a
C        line.  When the projected gradient at a point is less than ETA
C        times the initial value, the linear search is terminated.  Thus,
C        a small value of ETA means that the projected gradient must be
C        very small, i.e., close to the minimum along the line.  A value
C        of ETA close to 1 means that the projected gradient need
C        decrease only slightly from the starting point, so that the
C        resulting point is only a rough approximation to the minimum
C        along the line.  A guide to its value can be found in the NPL
C        report mentioned above, but the user may wish to try different
C        values on experimental runs to find the best value for his problem.
C        A value of 0.2E+0 is generally acceptable for QNMDIF.
C
C
C     TOL - A real positive number, which specifies the accuracy to
C        which the minimum is to be found.  QNMDIF terminates when the
C        norm of the approximate gradient is less than TOL ** (1/ 3) and
C        the distance moved during the last iteration is less than a
C        number close to TOL.  TOL should be chosen as approximately
C        the Euclidean distance which can be allowed between the point
C        returned by QNMDIF and the actual minimum.  Its value should
C        be chosen with care, keeping in mind that the error in the
C        approximate gradient is of order MAXH (I) ** 2 at best (when
C        central differences are used), where MAXH is the maximum
C        component of H in absolute value, that the function may be
C        known only to a given number of figures, and that the
C        parameters may be known only to a given number of figures
C        (for example, when the problems arises from a physical
C        situation, or a problem involving measured data).  A lower
C        bound for an acceptable value of TOL for most problems
C        is SQRT (EPSOBJ), but a larger value may be required.
C
C
C     STEPMX - The maximum allowed step from any point to the minimum
C        along a line from that point.  this parameter should be used
C        when the user is confident that the parameter vector should
C        not move more than the distance STEPMX from any particular
C        location.  If the linear search fails to find a minimum closer
C        to the current point than STEPMX, then the distance moved is
C        STEPMX.  It should be set to an upper bound on the distance
C        between the starting point and the minimum.  If the user does
c        not have a reasonable guess for this value, it can simply be
C        set to a very large number (say 1.0E+10), and then will have
C        essentially no effect on the linear search.
C
C
C     EPSMCH -  The machine precision of the computer to be used, i.e.,
C        the smallest value such that, in floating point arithmetic,
C        1 + EPSMCH > 1.  on the IBM 360, EPSMCH is 16.0 ** (-13).
C
CRAK
CRAK  EPSOBJ - Must be set by calling routine to value of smallest
CRAK     significant (absolute) change in objective function.  QNMDIF
CRAK     resets EPSOBJ to EPSMCH if input value is less than EPSMCH.
CRAK
C
C     UNITL - A logical parameter, set by the user.  In general, the
C        user will not know an approximation to the Hessian matrix
C        at the starting point, and in this event UNITL should be set
C        to .TRUE., so that the identity matrix is used as the starting
C        approximation to the Hessian.  If, however, the user does have
C        an idea of the approximate Hessian, UNITL should be set to
C        .FALSE., and the EL and D arrays should be set to the Cholesky
C        factors of the initial approximation (see comments under EL
C        and D, above).
CRAK
CRAK     MODIFICATION NOTE:  UNITL = .FALSE now implies that all of the
CRAK     information needed to begin an iteration (F, X, G, EL, and D)
CRAK     has been supplied (not merely EL and D).
CRAK
C
C     LOCAL - A logical parameter, to be set to .TRUE. or .FALSE. by the
C        user.  Gill and Murray have developed a sophisticated (?) "local
C        search" which follows the normal minimization, and which checks
C        through a random search whether a point can be found lower than
C        the estimated minimum.  The purpose of this search is to avoid
C        convergence to a saddle point, where the gradient is zero, but
C        which is not a minimum.  This search requires a fair number of
C        extra function evaluations, but will essentially guarantee
C        convergence to a true minimum.  If the user feels that there is
C        no danger of convergence to a saddle point, and in particular,
C        if the function is extremely expensive to evaluate, then LOCAL
C        should be set to .FALSE., and the local search will not be
C        executed.
C
C
C     CONV - A logical flag set by QNMDIF to .TRUE. if there has been
C        satisfactory convergence to a minimum, and to .FALSE. if not.
C        If LOCAL is .TRUE., then CONV = .TRUE. means that the local
C        search failed to find a significantly lower point.  if LOCAL
C        was set to .FALSE., then CONV = .TRUE. means that the gradient
C        at the final point was essentially zero, and that the steps to
C        the final point had converged.  Note that if LOCAL is set to
C        .FALSE., CONV = .TRUE. can happen if there is convergence to
C        a saddle point, since there is no way to check the second order
C        conditions to distinguish a minimum from a saddle point.
C
CRAK
CRAK  CONTIN - A logical input flag which must be .TRUE. if QNNDIF is
CRAK     to update the gradient and Hessian before returning due to
CRAK     iteration count.  Efficiency is not compromised if the final
CRAK     values are saved and later used by the calling routine for
CRAK     restarting.
CRAK
CRAK
CRAK  RESCUE - A logical input flag which may be set .TRUE. if the user
CRAK     wishes to begin the run with the local search procedure, which
CRAK     counts as two iterations.  Normally, use RESCUE = .FALSE. and set
CRAK     LOCAL = .TRUE. if some assurance against premature convergence
CRAK     is desired.  The RESCUE option may be useful if the objective
CRAK     function routine has blundered into a bad region.
CRAK
C
C     PRINT - A logical flag, set by the user to control printout.
C        for PRINT = .FALSE., most printing is suppressed.
C
C
C     FUN - A subroutine provided by the user, which must be declared
C        as EXTERNAL in the calling subprogram.  Its call must be of
C        the form
C
C        CALL FUN (N, X, FVAL)
C
C        where N is the number of variables, X is the vector of
C        variables and FVAL is the computed value of the function.
C
CRAK
CRAK  USER - A subroutine provided by the user, which must be declared
CRAK     as EXTERNAL in the calling subprogram.  Its call must be of
CRAK     the form
CRAK
CRAK     CALL USER (IENTRY, N, NLDIM, X, F, G, H, EL, D, NFTOTL,
CRAK    >   NITER, NTYPE, CONTIN)
CRAK
CRAK     where the variables are the current values of those described
CRAK     above.  QNMDIF calls USER at the end of each successful line
CRAK     search with IENTRY = 1, at which point X, F, H, NFTOTL, and
CRAK     NITER are current, but G, EL, and D have not yet been updated.
CRAK     This permits certain housekeeping chores to be performed before
CRAK     the calculations for the next iteration are begun.  QNMDIF sets
CRAK     IENTRY = 2 and calls USER again following calculation of gradient
CRAK     and Hessian updates, as well as after the local search if
CRAK     requested.  Disk files for saving restart information should be
CRAK     written when IENTRY = 2, with an initial REWIND command.
CRAK     Different applications may dictate other choices for parameter
CRAK     list and calling points.  Note that either or both of these
CRAK     options may be simply turned off by providing subroutine USER
CRAK     with RETURN statements for the unneeded cases.
CRAK
CRAK
C
C     Development history
C     -------------------
C
C        Date           Initials   Description
C            ?   1972              Initial coding in ALGOL at the NPL
C                                     (Gill, Murray, and Pitfield)
C            ?                     FORTRAN translation, Stanford Univ.
C                                     (Margaret Wright)
C        20 Sep. 1982   RAK        Version 1.1, modified as noted for
C                                     NASA release
C        13 Aug. 1986   RAK        Version 1.2 - mostly cosmetic mods.
C                                     Included CENDIF/FDSTEP/DIFFER for
C                                     stepsize selection (though they must
C                                     be called separately).  Added PARAMETER
C                                     statements for local array dimensions,
C                                     with test that N <= MAXN on entry.
C                                     Increase MAXN from 30 to 40.  Some
C                                     protection added in LINSCH.
C
C
C     Author
C     ------
C
C        Robert Kennelly
C        NASA - Ames Research Center
C        Mail Stop 227-6
C        Moffett Field, CA  94035
C
C        Telephone:  (415) 694-5944
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT REAL (A - H, O - Z)

C     Constants.

      INTEGER
     >   MAXN
      PARAMETER
     >   (MAXN = 40)

C     Arguments.

      LOGICAL
     >   CONTIN, CONV, COUNT, FINAL, LOCAL, PRINT, RESCUE, SUCCES, UNITL
      REAL
     >   D(N), EL(NLDIM), G(N), H(N), X(N)

C     Local variables.

      LOGICAL
     >   DONE
      REAL
     >   OLDG(MAXN), P(MAXN), PP(MAXN), Y(MAXN)

C     Procedures.

      EXTERNAL
     >   FUN, USER


C     Execution.
C     ----------

      IF (N .GT. MAXN) STOP 'QNMDIF:  ERROR - too many variables!'


C     Initialization.
C     ---------------

CRAK  Save input value of NITER for use as limit.

      NITMAX = NITER
      NITER = 0

C     If UNITL is .TRUE., use identity matrix as approximation to Hessian
C     for the first iteration (steepest descent).

      FNEW = F
      OLDF = F
      IF (.NOT.UNITL) GO TO 30
         DO 10 I = 1, N
            G(I) = 0.0E+0
            D(I) = 1.0E+0
   10    CONTINUE
         DO 20 I = 1, NLDIM
            EL(I) = 0.0E+0
   20    CONTINUE
         IF (.NOT.RESCUE) CALL APROXG (N, NTYPE, X, G, H, PP, FUN,
     >      FNEW, NFTOTL, PRINT, LWRITE)
   30 CONTINUE
      IF (EPSOBJ .LT. EPSMCH) EPSOBJ = EPSMCH
      RTEPS = SQRT (EPSOBJ)
      TINY = EPSMCH ** 2
      CONV = .TRUE.
      FINAL = .FALSE.
      COUNT = .FALSE.
      BOUNDK = 0.01E+0/ (SQRT (REAL (N)) * EPSMCH)

C     Calculate maximum component of vector H and ALPMIN, which
C     is used as a measure of the smallest step that should be
C     taken along any direction.

      HMAX = 0.0E+0
      DO 40 I = 1, N
         U = ABS (H(I))
         IF (U .GT. HMAX) HMAX = U
         P(I) = 0.0E+0
   40 CONTINUE
      ALPMIN = HMAX ** (2.0E+0/ 3.0E+0)

CRAK  Start with local search if rescue option selected.

      IF (RESCUE) GO TO 290


C     Main iteration loop.
C     --------------------

   50 CONTINUE

CRAK  Save current variable and gradient vectors.

      OLDF = FNEW
      DO 60 I = 1, N
         OLDG(I) = G(I)
         Y(I) = X(I)
   60 CONTINUE


C     Search direction.
C     -----------------

C     Calculate the direction of search, P, by solving the system of
C     linear equations EL * D * ELT * P = -G, where EL is a unit lower
C     triangular matrix , ELT is its transpose, and D is a diagonal
C     matrix represented in the program by the vector D.

C     Forward solution.

      P(1) = - OLDG(1)
      IF (N .EQ. 1) GO TO 90
         IR = 1
         DO 80 I = 2, N
            SUM = - OLDG(I)
            IT  = I - 1
            DO 70 K = 1, IT
               SUM = SUM - P(K) * EL(IR)
               IR  = IR + 1
   70       CONTINUE
            P(I) = SUM
   80    CONTINUE
   90 CONTINUE

C     Back substitution.

      P(N) = P(N)/ D(N)
      IF (N .EQ. 1) GO TO 120
         DO 110 I = 2, N
            IN = N + 1 - I
            IR = IR - 1
            IS = IR
            IT = IN + 1
            SUM = P(IN)/ D(IN)
            DO 100 K = IT, N
               KN = N + IT - K
               SUM = SUM - P(KN) * EL(IS)
               IS = IS + 2 - KN
  100       CONTINUE
            P(IN) = SUM
  110    CONTINUE
  120 CONTINUE

C     Compute the norm of P (the direction of search), the norm of
C     G (the approximate gradient vector), and the dot product
C     of P and G (which should be negative if P is a descent
C     direction).

      SUM = 0.0E+0
      GTP = 0.0E+0
      GNM = 0.0E+0
      DO 130 I = 1, N
         SUM = SUM + P(I) ** 2
         GTP = GTP + G(I) * P(I)
         GNM = GNM + G(I) ** 2
  130 CONTINUE
      PNORM = SQRT (SUM) + TINY
      GNM = SQRT (GNM) + TINY

CRAK  Protect GTP from being unreasonably small.

      IF (ABS (GTP) .LT. TINY) GTP = SIGN (TINY, GTP)


C     Line search.
C     ------------

C     Initialize ALPHA, the step to be taken along the direction of
C     search P.  The initial value is input to the linear search
C     subroutine, which estimates the approximate minimum of the
C     objective function along P.

CRAK  EPSOBJ has been substituted for SQRT (EPSMCH) as a tolerance
CRAK  here.  In all other appearances of EPSOBJ, it was substituted
CRAK  directly for EPSMCH.

      DELTAF = MAX (ABS (FNEW - FM), EPSOBJ)
      ALPHA = MIN (- 2.0E+0  *  DELTAF/  GTP, 1.E+0)

CRAK  Line search tolerance (T) is scale dependent.  For problems with
CRAK  second derivatives very different from order one (especially
CRAK  where desired accuracy is limited by available objective function
CRAK  precision), it may be desirable to try other values.  The original
CRAK  version of QNMDIF used T = RTEPS/ PNORM; we use EPSOBJ/ PNORM to
CRAK  to force LINSCH to work a bit harder.

      T = EPSOBJ/ PNORM
      SFTBND = 0.E+0
      IF (NTYPE .EQ. 0) SFTBND = ALPMIN/ PNORM

C     Find suitably lower point along direction of search.  If the linear
C     search is successful, the vector X and the function value FNEW are
C     modified within the subroutine LINSCH.  If the search fails, they
C     remain unchanged.

      CALL LINSCH (N, FUN, RTEPS, T, ETA, SFTBND, STEPMX/ PNORM, P,
     >    X, FNEW, ALPHA, GTP, NFTOTL, SUCCES, PRINT, LWRITE)
      IF (SUCCES .OR. NTYPE.EQ.1) NITER = NITER + 1
      COUNT = NITER .GE. NITMAX
      IF (SUCCES) GO TO 150
         IF (PRINT) WRITE (LWRITE, 1000)
 1000    FORMAT (34H0QNMDIF:  UNSUCCESSFUL LINE SEARCH)

C        If central differences were used to approximate the gradient,
C        branch to failure exit.

         IF (NTYPE .EQ. 1) GO TO 290

CRAK     Start over using central differences if first line search based
CRAK     on forward difference input derivative data failed (since array
CRAK     PP not initialized).

         IF (NITER.GT.0 .OR. UNITL) GO TO 140
            NTYPE = 1
            CALL APROXG (N, NTYPE, X, G, H, PP, FUN, FNEW, NFTOTL,
     >         PRINT, LWRITE)
            GO TO 50
  140    CONTINUE

C        Otherwise, switch to central difference approximation (using
C        intermediate values stored in array PP) and try again.

         NTYPE = 1
         CALL APROXG (N, 2, X, G, H, PP, FUN, FNEW, NFTOTL,
     >      PRINT, LWRITE)
         GO TO 50
  150 CONTINUE

C     The linear search was successful - print if required.

      IF (.NOT.PRINT) GO TO 160
         WRITE (LWRITE, 1010)
 1010    FORMAT (1H1)
         CALL OUTPUT (N, NLDIM, Y, G, H, P, OLDF, ALPHA, D,
     >      EL, NITER, NFTOTL, FINAL, LWRITE)
         WRITE (LWRITE, 1010)
  160 CONTINUE

CRAK  Invoke user-supplied routine for housekeeping, then jump out if
CRAK  maximum number of iterations have been performed and final
CRAK  derivative information was not requested.  Note that the values
CRAK  printed out in "final information" will not be up-to-date.

      IENTRY = 1
      CALL USER (IENTRY, N, NLDIM, X, FNEW, G, H, EL, D, NFTOTL,
     >    NITER, NTYPE, CONTIN)
      IF (COUNT .AND. .NOT.CONTIN) GO TO 290


C     Calculate gradient.
C     -------------------

C     Check whether central differences were used for the gradient.
C     If so, and if the step taken was sufficiently large, switch
C     back to forward differences.

      IF (NTYPE .EQ. 1 .AND. ALPHA .GT. ALPMIN/ PNORM)  NTYPE = 0

C     Calculate the new approximate gradient.

      CALL APROXG (N, NTYPE, X, G, H, PP, FUN, FNEW, NFTOTL,
     >   PRINT, LWRITE)
      GNMSQ = 0.0E+0
      DO 170 I = 1, N
         GNMSQ = GNMSQ + G(I) ** 2
  170 CONTINUE

C     If we are using forward differences, and the norm of the gradient
C     is small, recalculate the gradient using central differences.

      IF (GNMSQ .GT. HMAX .OR. NTYPE .NE. 0) GO TO 190
         NTYPE = 1
         CALL APROXG (N, 2, X, G, H, PP, FUN, FNEW, NFTOTL,
     >      PRINT, LWRITE)
         GNMSQ = 0.0E+0
         DO 180 I = 1, N
            GNMSQ = GNMSQ + G(I) ** 2
  180    CONTINUE
  190 CONTINUE


C     Update Hessian.
C     ---------------

C     Update the Cholesky factors of the Hessian approximation using
C     the rank-two COMDFP (Complementary Davidon-Fletcher-Powell)
C     update formula, providing that new GTP > old GTP.

      GTPFIN = 0.E+0
      DO 200 I = 1, N
         GTPFIN = GTPFIN + G(I) * P(I)
  200 CONTINUE
      V = ALPHA * (GTPFIN - GTP)
      IF (V .GT. 0.E+0) GO TO 210
         IF (PRINT) WRITE (LWRITE, 1020) NITER, GTP, GTPFIN
 1020    FORMAT (37H0QNMDIF:  WARNING - ON ITERATION NO. , I6,
     >      29H, ABS (GTP) DID NOT DECREASE./
     >      10X, 15HGTP(INITIAL) = , E25.15/
     >      10X, 15HGTP(FINAL)   = , E25.15/
     >      10X, 37HUPDATES TO THE HESSIAN WERE BYPASSED.)
         GO TO 270
  210 CONTINUE
         V = SQRT (V)
         DO 220 I = 1, N
            P(I) = (G(I) - OLDG(I))/ V
  220    CONTINUE

C        Update the Cholesky factors after a positive rank-one change.

         CALL C1MDCH (N, NLDIM, D, EL, 1.0E+0, P, IFAIL)
         IF (IFAIL .NE. 0) GO TO 290
         V = SQRT (ABS (GTP))
         DO 230 I = 1, N
            OLDG(I) = - OLDG(I)/ V
  230    CONTINUE

C        Update Cholesky factors following negative rank-one correction.

         CALL UPCHOL (N, NLDIM, EPSMCH, D, EL, OLDG)

C        Check whether the diagonal elements of the updated factors are
C        sufficiently large.  If not, modify them appropriately to ensure
C        that the approximate Hessian is sufficiently positive definite.

         DMAX = D(1)
         IF (N .EQ. 1) GO TO 250
            DO 240 I = 2, N
               IF (D(I) .GT. DMAX)  DMAX = D(I)
  240       CONTINUE
  250    CONTINUE
         BL = DMAX/ BOUNDK
         DO 260 I = 1, N
            IF (D(I) .LT. BL)  D(I) = BL
  260    CONTINUE
  270 CONTINUE

CRAK  Invoke user-supplied routine for housekeeping.

      IENTRY = 2
      CALL USER (IENTRY, N, NLDIM, X, FNEW, G, H, EL, D, NFTOTL,
     >   NITER, NTYPE, CONTIN)


C     Check convergence.
C     ------------------

      SUM = 0.E+0
      DO 280 I = 1, N
         SUM = SUM + X(I) ** 2
  280 CONTINUE

C     Check whether norm of gradient and last step taken were sufficiently
C     small.  If not, verify that objective decreased and check iteration
C     count before continuing.

      IF (GNMSQ.LT.TOL ** .667E+0 .AND.
     >   ALPHA * PNORM.LT.SQRT (EPSMCH * SUM) + TOL) GO TO 300
      IF (OLDF.GT.FNEW .AND. .NOT.COUNT) GO TO 50

CRAK  Exits from main loop:
CRAK
CRAK     290 - Failure, reset CONV
CRAK     300 - Success, leave CONV set to .TRUE.

  290 CONTINUE
      CONV = .FALSE.
  300 CONTINUE

CRAK  Perform local search if requested, but check iteration count.

      IF (.NOT.(LOCAL.OR.RESCUE) .OR. COUNT) GO TO 500


C     Local search.
C     -------------

C     Carry out local search, which is a procedure of searching along
C     random orthogonal directions to see if a lower point can be found.

      IF (PRINT) WRITE (LWRITE, 1030)
 1030 FORMAT (38H0QNMDIF:  BEGIN LOCAL SEARCH PROCEDURE)
      U = ALPMIN
      EPS10 = 10.E+0 * EPSOBJ * ABS (FNEW)

C     Modify independent variable vector with a small change.

  310 CONTINUE
      DO 320 I = 1, N
         Y(I) = X(I) + U
  320 CONTINUE
      CALL FUN (N, Y, V)
      NFTOTL = NFTOTL + 1

C     Check whether a sufficiently large change in the function occurred.
C     If not, increase step size and try again.

      IF (ABS (V - OLDF) .GE. EPS10) GO TO 330
         U = 5.0E+0 * U
         GO TO 310
  330 CONTINUE

C     Calculate orthogonal direction at the modified point, W.

      P(1) = U
      IF (N .EQ. 1) GO TO 350
         DO 340 I = 2, N
            P(I) = - P(I - 1)
  340    CONTINUE
         IF (MOD(N, 2) .EQ. 1)  P(N) = 0.0E+0
  350 CONTINUE
      U = SQRT (ALPMIN)
      EPS10 = EPS10 * ABS (V)

CRAK  Begin first "iteration" (rewritten without loop - RAK, 16 June 81)

      DO 360 I = 1, N
         OLDG(I) = Y(I) + U * P(I)
  360 CONTINUE
      CALL FUN (N, OLDG, OLDF)
      NFTOTL = NFTOTL + 1
      GTP = OLDF - V

C     Ensure that the step has made a significant change in the function
C     value - if not, increase step size.

      IF (ABS (GTP) .GE. EPS10) GO TO 370
         U = 5.0E+0 * U
         GO TO 310
  370 CONTINUE
      GTP = - ABS (GTP/ U) - 1.0E+0
      ALPHA = U
      IF (OLDF .LE. V) GO TO 390
         V = OLDF
         DO 380 I = 1, N
            Y(I) = OLDG(I)
            P(I) = - P(I)
  380    CONTINUE
  390 CONTINUE
      SUM = 0.0E+0
      DO 400 I = 1, N
         SUM = SUM + P(I) ** 2
  400 CONTINUE
      PNORM = SQRT (SUM) + TINY
      SFTBND = 0.0E+0
      T = RTEPS/ PNORM
      CALL LINSCH (N, FUN, RTEPS, T, EPSMCH, SFTBND, STEPMX/ PNORM,
     >   P, Y, V, ALPHA, GTP, NFTOTL, SUCCES, PRINT, LWRITE)

CRAK  Begin second iteration.

      DO 410 I = 1, N
         P(I) = X(I) - Y(I)
         Y(I) = X(I)
  410 CONTINUE
      IF (V .LT. FNEW) GO TO 440
         SUM = 0.0E+0
         DO 420 I = 1, N
            SUM = SUM + P(I) ** 2
  420    CONTINUE
         PNORM = SQRT (SUM) + TINY
         U = HMAX/ PNORM
         DO 430 I = 1, N
            OLDG(I) = X(I) + U * P(I)
  430    CONTINUE
         CALL FUN (N, OLDG, OLDF)
         NFTOTL = NFTOTL + 1
         GTP = - ABS ((OLDF - FNEW)/ U)
         GO TO 450
  440 CONTINUE
         GTP = (V - FNEW) * 4.0E+0/ 3.0E+0
  450 CONTINUE
      ALPHA = 1.0E+0
      IF (V.GE.FNEW .AND. OLDF.LE.FNEW) GO TO 470
         DO 460 I = 1, N
            P(I) = -P(I)
  460    CONTINUE
  470 CONTINUE
      V = FNEW
      SFTBND = 0.0E+0
      T = RTEPS/ PNORM
      CALL LINSCH (N, FUN, RTEPS, T, EPSMCH, SFTBND, STEPMX/ PNORM,
     >   P, Y, V, ALPHA, GTP, NFTOTL, SUCCES, PRINT, LWRITE)
      NITER = NITER + 2
      COUNT = NITER .GE. NITMAX

CRAK  Invoke user-supplied routine for housekeeping.  The second call
CRAK  here (with IENTRY = 2) will usually be repeated below with more
CRAK  complete information, if requested.

      IENTRY = 1
      CALL USER (IENTRY, N, NLDIM, Y, V, G, H, EL, D, NFTOTL,
     >    NITER, NTYPE, CONTIN)
      IENTRY = 2
      CALL USER (IENTRY, N, NLDIM, Y, V, G, H, EL, D, NFTOTL,
     >    NITER, NTYPE, CONTIN)
      IF (PRINT) WRITE (LWRITE, 1040)
 1040 FORMAT (26H0QNMDIF:  END LOCAL SEARCH)

C     Terminate if a lower point was not found.  No new gradient will
C     be calculated, even if CONTIN = .TRUE., but previous values are
C     still current.  In special RESCUE = .TRUE. case, the gradient
C     will be either the input value or set to zero, depending on UNITL.

      IF (V .GE. FNEW) GO TO 500

C     If a lower point was found, the norm of the step to the lower
C     point is small, and the function decreased during the last
C     iteration, or iteration limit has been reached, and no final
C     gradient is requested, terminate.

      FNEW = V
      SUM = 0.0E+0
      DO 490 I = 1, N
         X(I) = Y(I)
         SUM = SUM + Y(I) ** 2
  490 CONTINUE
      DONE = (ALPHA * PNORM .LT. SQRT (EPSMCH * SUM) + TOL) .AND. CONV
      IF ((DONE .OR. COUNT) .AND. .NOT.CONTIN) GO TO 500

CRAK  Calculate gradient, call user interface to permit saving new
CRAK  information, and quit if through.  Note that stopping here leaves
CRAK  the Hessian as it was, i.e. no update is used.  This is a
CRAK  temporary expedient for this special case - the matrix update
CRAK  procedure will be modularized and called following APROXG in
CRAK  a future version.

      NTYPE = 1
      IF (RESCUE) NTYPE = 0
      CALL APROXG (N, NTYPE, X, G, H, PP, FUN, FNEW, NFTOTL,
     >   PRINT, LWRITE)
      IENTRY = 2
      CALL USER (IENTRY, N, NLDIM, X, V, G, H, EL, D, NFTOTL,
     >    NITER, NTYPE, CONTIN)
      IF (DONE .OR. COUNT) GO TO 500
      RESCUE = .FALSE.
      CONV = .TRUE.
      GO TO 50


C     Termination.
C     ------------

  500 CONTINUE
      F = FNEW

CRAK  This version always prints final results - restore the test if desired.
CRAK  IF (.NOT.PRINT) GO TO 510

C        Print final optimization results.

         WRITE (LWRITE, 1010)
         IF (COUNT) WRITE (LWRITE, 1050)
 1050    FORMAT (37H0MAXIMUM NUMBER OF ITERATIONS REACHED)
         IF (CONV) WRITE (LWRITE, 1060)
 1060    FORMAT (44H0THE TEST FOR CONVERGENCE HAS BEEN SATISFIED)
         IF (.NOT. CONV) WRITE (LWRITE, 1070)
 1070    FORMAT (50H0QNMDIF HAS FAILED TO SATISFY THE CONVERGENCE TEST)
         FINAL = .TRUE.
         CALL OUTPUT (N, NLDIM, X, G, H, P, FNEW, ALPHA, D, EL, NITER,
     >      NFTOTL, FINAL, LWRITE)
  510 CONTINUE

      RETURN
      END
*DECK APROXG
C+----------------------------------------------------------------------
C
      SUBROUTINE APROXG (N, NTYPE, X, G, H, WORK, FUN, FNEW, NFTOTL,
     >   PRINT, LWRITE)
C
C
C        This routine is called by subroutine QNMDIF.  It assigns a
C     finite-difference approximation to the gradient of the function
C     F (X) to the 1 * N array G (I). The intervals for differencing F
C     along each of the coordinate directions are given in array H (I).
C     The integer variable NTYPE determines which type of approximation
C     is given.
C
C        NTYPE=0:  Forward differences are obtained and the function
C           values are calculated at N forward points and stored in
C           the 1 * N array WORK.
C
C        NTYPE=1:  Central differences requiring the full 2N function
C           evaluations are computed.
C
C        NTYPE=2:  Central differences are evaluated using the forward
C           points previously stored in array WORK.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT REAL (A - H, O - Z)

C     Arguments.

      LOGICAL
     >   PRINT
      REAL
     >   G(N), H(N), WORK(N), X(N)
      EXTERNAL
     >   FUN

C     Execution.
C     ----------

      IF (PRINT .AND. NTYPE .EQ. 0) WRITE (LWRITE, 1000) NTYPE
 1000 FORMAT (51H0APROXG:  CALCULATE GRADIENT BY FORWARD DIFFERENCES,
     >   10H (NTYPE = , I1, 1H))
      IF (PRINT .AND. NTYPE .GT. 0) WRITE (LWRITE, 1010) NTYPE
 1010 FORMAT (51H0APROXG:  CALCULATE GRADIENT BY CENTRAL DIFFERENCES,
     >   10H (NTYPE = , I1, 1H))

      IF (NTYPE .EQ. 2) GO TO 20

C        Compute forward points.

         DO 10 I = 1, N
            XI = X(I)
            X(I) = XI + H(I)
            CALL FUN (N, X, FUNX)
            WORK(I) = FUNX
            X(I) = XI
   10    CONTINUE
         NFTOTL = NFTOTL + N
   20 CONTINUE
      IF (NTYPE .GT. 0) GO TO 40

C        Compute forward difference gradient.

         DO 30 I = 1, N
            G(I) = (WORK(I) - FNEW)/ H(I)
   30    CONTINUE
         GO TO 60
   40 CONTINUE

C        Compute backward points and central difference gradient.

         DO 50 I = 1, N
            XI = X(I)
            HI = H(I)
            X(I) = XI - HI
            CALL FUN (N, X, FUNX)
            G(I) = 0.5E+0 * (WORK(I) - FUNX)/ HI
            X(I) = XI
   50    CONTINUE
         NFTOTL = NFTOTL + N
   60 CONTINUE


C     Termination.
C     ------------

      RETURN
      END
*DECK C1MDCH
C+----------------------------------------------------------------------
C
      SUBROUTINE C1MDCH (N, NLDIM, D, EL, ALPHA, Z, IFAIL)
C
C
C        The subroutine C1MDCH updates the Cholesky factorization of the
C     matrix EL * D * ELT  + ALPHA * Z * ZT, where EL is a unit lower
C     triangular matrix, D is a diagonal matrix with positive elements,
C     ELT is the transpose of EL, ALPHA is a positive scalar, Z is a
C     vector, and ZT is the transpose of Z, so that ALPHA * Z * ZT is
C     a positive rank-one change.
C
C        The algorithm used is described in Gill, Golub, Murray, and
C     Saunders, "Methods for Modifying Matrix factorizations", Stanford
C     Computer Science Report CS-72-322, Stanford University, Dec. 1972.
C
C        N is the dimension of the matrices and vectors, NLDIM is
C     N * (N-1)/ 2, the length of the vector EL in which the strict
C     lower triangle of the matrix EL is stored by rows, D is a
C     vector of length N containing the elements of the matrix D,
C     EL is a vector of length NLDIM containing the elements of EL,
C     ALPHA is the positive scalar multiple of the rank-one change,
C     Z is the vector involved, and IFAIL is an integer flag set
C     to zero if there are no errors, and to 1 if overflow occurred.
C
C        Both EL and D are overwritten with the modified factors of the
C     altered matrix, and the values of ALPHA and Z are not retained.
C
C        The subroutine sets IFAIL to zero if there are no problems.
C     IFAIL is set to 1 if any element of the diagonal formed is zero
C     or negative, or when the ratio of the new diagonal element to
C     the old is too large and would cause overflow.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT REAL (A - H, O - Z)

C     Constants.

C     RMAX should be set to the largest positive floating point value such
C     that +RMAX and -RMAX can both be represented in the computer.  For
C     convenience, the smallest useful value (VAX) has been used - it is
C     unlikely to matter, but could be increased on another system.

      REAL
     >   RMAX
      PARAMETER
     >   (RMAX = 1.E+38)

C     Arguments.

      REAL
     >   D(N), EL(NLDIM), Z(N)


C     Execution.
C     ----------

      A = ALPHA
      K = 0
      IFAIL = 1
      DO 70 I = 1, N
         P1 = Z(I)
         DI = D(I)
         T = A * P1
         D(I) = DI + T * P1
         DB = D(I)

C        Exit with IFAIL = 1 if any of the new diagonal elements is
C        non-positive, or if the ratio of diagonals will overflow.

         IF (DB .GE. 1.0E+0) GO TO 10
            IF (DB .LT. 0.0E+0 .OR. DI .GT. (RMAX * DB))  RETURN
   10    CONTINUE
         GAMMA = DI/ DB
         BETA = T/ DB
         K = K + I
         J = K
         A = A * GAMMA
         IF (I .EQ. N) GO TO 60
            IP1 = I + 1
            IF (GAMMA .GT. 0.25E+0) GO TO 30

C              If ratio of diagonals is less than 4.0, proceed with
C              normal update.

               DO 20 IB = IP1, N
                  T = EL(J)
                  EL(J) = T * GAMMA + BETA * Z(IB)
                  Z(IB) = Z(IB) - P1 * T
                  J = J + IB - 1
   20          CONTINUE
               GO TO 50
   30       CONTINUE

C              Use alternative update if ratio of diagonals exceeds 4.0.

               DO 40 IB = IP1, N
                  T = EL(J)
                  Z(IB) = Z(IB) - P1 * T
                  W = Z(IB)
                  EL(J) = BETA * W + T
                  J = J + IB - 1
   40          CONTINUE
   50       CONTINUE
   60    CONTINUE
   70 CONTINUE
      IFAIL = 0


C     Termination.
C     ------------

      RETURN
      END
*DECK LINSCH
C+----------------------------------------------------------------------
C
      SUBROUTINE LINSCH (N, FUN, EPS, T, ETA, SFTBND, STEPMX,
     >   P, X, F, ALPHA, GTP, NFTOTL, SUCCES, PRINT, LWRITE)
C
C
C        Called by QNMDIF to execute a linear search along a given
C     direction P, starting from a given point X, to locate an
C     approximation ALPHA to the point at which the objective function
C     attains its minimum value along the direction P.  The method used
C     is that of successive quadratic interpolation with safeguards.
C
C
C     History:
C
C        12 Aug. 1986    RAK    Print and test T on entry to monitor
C                               conflicts with SFTBND and STEPMX.  Check
C                               D and GTP for before main loop.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT REAL (A - H, O - Z)

C     Constants.

      INTEGER
     >   MAXN
      PARAMETER
     >   (MAXN = 40)

C     Arguments.

      LOGICAL
     >   PRINT, SUCCES
      REAL
     >   P(N), X(N)

C     Local variables.

      REAL
     >   Z(MAXN)

C     Procedures.

      EXTERNAL
     >   FUN

C     Execution.
C     ----------

      IF (PRINT) WRITE (LWRITE, 1000) SFTBND, STEPMX, T
 1000 FORMAT (36H0LINSCH:  BEGIN SEARCH, MIN. STEP = , E15.8/
     >   24X, 12HMAX. STEP = , E15.8/
     >   24X, 12HTOLERANCE = , E15.8)

      IF (T .GT. STEPMX) THEN
         WRITE (LWRITE, 1010)
 1010    FORMAT (51H0LINSCH:  MAX. STEP MUST BE GREATER THAN TOLERANCE!)
         SUCCES = .FALSE.
         GO TO 330
      END IF
         

C     Initialization.
C     ---------------

C     NLIN is the number of function evaluations during the linear
C     search.  It is used locally as a flag to compute the parabolic
C     step the first time using GTP, the directional derivative along
C     the search vector.

CRAK  MAXLIN is a (fairly loose) upper limit on the number of function
CRAK  evaluations permitted before declaring the search a failure.

      NLIN = 0
      MAXLIN = 20
      V = 0.0E+0
      W = 0.0E+0
      XMIN = 0.0E+0
      FA = F
      FV = F
      FW = F
      FMIN = F
      FOLD = F
      D = ALPHA
      TOL = T
      T2 = 2.0E+0 * TOL

C     All points in the linear search are scaled (shifted) so that
C     the currently lowest point (XMIN) is the origin.  A and B define
C     the interval of uncertainty, W is the last value of XMIN, and V
C     is the highest of the three points through which a parabola may
C     be fitted.

      A = 0.0E+0
      B = STEPMX + EPS * ABS (STEPMX) + T
      E = B
      SCXBD = STEPMX

C     NPCNT is a local count of the number of non-parabolic interpolation
C     steps used in LINSCH.  It may be useful for diagnostic purposes if
C     there are problems locating the minimum of a particular function.
C     GTEST1 and GTEST2 are used to test convergence.

      NPCNT = 0
      B1 = B
      A1 = A
      GTEST1 = - 1.0E-4 * GTP
      GTEST2 = - ETA * GTP
      ALPHA = 0.0E+0

C     D is the estimate of the next point at which the function is to be
C     evaluated, with respect to the current origin.  Make sure that it,
C     and GTP, are sensible before proceding.

      IF (D .LE. 0.0E+0 .OR. GTP .GE. 0.0E+0) THEN
         WRITE (LWRITE, 1020) D, GTP
 1020    FORMAT (32H0LINSCH:  BAD INITIAL DATA, D = , E25.15,
     >      8H, GTP = , E25.15)
         SUCCES = .FALSE.
         GO TO 330
      END IF


C     Begin main loop.
C     ----------------

   10 CONTINUE
      IF (D .LT. SCXBD) GO TO 20

C        D exceeds the shifted bound on the minimum, so adjust D and
C        the bound.

         D = SCXBD
         SCXBD = (SCXBD - TOL)/ (1.0E+0 + EPS)
   20 CONTINUE

C     U is the point at which the function will actually be evaluated,
C     scaled with respect to the shifted origin.  Thus XMIN, the actual
C     step from the original starting point, must be added to obtain
C     the true step.  The estimate D will be used as the value for U
C     only if it is sufficiently far from zero (the current minimum),
C     since the function is not allowed to be evaluated at points closer
C     together than TOL.

      U = D
      IF (ABS (D) .LT. TOL) U = SIGN (TOL, D)
      R = XMIN + U
      DO 30 I = 1, N
         Z(I) = X(I) + R * P(I)
   30 CONTINUE


C     Evaluate function at new point.
C     -------------------------------

      IF (PRINT) WRITE (LWRITE, 1030) R
 1030 FORMAT (17H0LINSCH:  STEP = , E25.15)
      CALL FUN (N, Z, FU)
      NLIN = NLIN + 1

C     Update A, B, V, W, and XMIN.  Check whether new function value is
C     lower than previous minimum.

      IF (FU .GT. FMIN) GO TO 60

C     The following code treats the case where FU is the new lowest
C     point, so that the new point, U, becomes the next origin, and the
C     other points are shifted accordingly.

C     Shift left or right endpoint depending on which side of origin
C     the new point is.

      IF (U .LT. 0.0E+0) GO TO 40
         A = 0.0E+0
         FA = FMIN
         GO TO 50
   40 CONTINUE
         B = 0.0E+0
   50 CONTINUE


C     Shift all points with respect to new minimum.
C     ---------------------------------------------

      V = W
      FV = FW
      W = 0.0E+0
      FW = FMIN
      FMIN = FU
      XMIN = XMIN + U
      A = A - U
      B = B - U
      V = V - U
      W = W - U
      SCXBD = SCXBD - U
      TOL = EPS * ABS (XMIN) + T
      T2 = 2.0E+0 * TOL
      GO TO 110

C     In the following code, the new function value was greater than
C     the previously lowest.  Check the relationship of the new point
C     to the other values (next lowest, etc.).  The origin remains
C     unchanged, but other points may be interchanged.

C     Shift either the left or right endpoint of the interval of
C     uncertainty, depending on whether the new point, U, was to the
C     left or right of the origin.

   60 IF (U .GT. 0.0E+0) GO TO 70
         A = U
         FA = FU
         GO TO 80
   70 CONTINUE
         B = U
   80 CONTINUE
      IF (FU .GT. FW .AND. W .NE. 0.0E+0) GO TO 90

C     If FU is less than or equal previous second best point, or W = 0,
C     i.e., this is the first time through this section of code.

      V = W
      FV = FW
      W = U
      FW = FU
      GO TO 100
   90 IF (FU .GT. FV .AND. V .NE. 0.0E+0 .AND. V .NE. W) GO TO 100

C     FU .LE. FW, or V = 0, or V = W (first or second time through).

      V = U
      FV = FU
  100 U = 0.0E+0

C     Compute midpoint of interval of uncertainty.

  110 XM = 0.5E+0 * (A + B)


C     Check stopping criteria.
C     ------------------------

C     Stop if interval of uncertainty is sufficiently small, or
C     if the best point is less than the required bound, or
C     if function value has decreased sufficiently, or
C     if MAXLIN function evaluations have already been performed.

      IF ((ABS (XM) .LE. (T2 - 0.5E+0 * (B - A))) .OR.
     >   ((XMIN + B) .LE. SFTBND) .OR.
     >   ((FA - FMIN).LE.(ABS (A) * GTEST2) .AND. FMIN.LT.FOLD) .OR.
     >   (NLIN .GE. MAXLIN)) GO TO 280

C     If stopping criteria are not met, continue.

      S = 0.0E+0
      Q = 0.0E+0
      R = 0.0E+0
      IF (ABS (E) .LE. TOL) GO TO 170

C     Otherwise, fit parabola through XMIN, V, and W.

      IF (NLIN .GT. 1) GO TO 130

C     If NLIN = 1, this is the first parabolic fit and the (known)
C     approximate gradient at the initial point can be used.

      Q = 2.0E+0 * (FW - FMIN - W * GTP)
      IF (XMIN .EQ. 0.0E+0) GO TO 120
      S = (2.0E+0 * (FMIN - FW) + W * GTP) * W
      GO TO 140
  120 S = GTP * W * W
      GO TO 140

C     NLIN is greater than 1, so the fit uses function values only.

  130 R = W * (FV - FMIN)
      Q = V * (FW - FMIN)
      S = R * W - Q * V
      Q = 2.0E+0 * (Q - R)
  140 IF (Q .LE. 0.0E+0) GO TO 150
      S = - S
      GO TO 160
  150 Q = - Q
  160 R = E
      IF (D .NE. B1 .OR. B .LE. SCXBD)  E = D

C     Construct an artificial bound on the estimated step length.

  170 A1 = A
      B1 = B
      IF (XMIN .NE. 0.0E+0) GO TO 180
      D = XM
      GO TO 230
  180 IF (B .LE. SCXBD) GO TO 190

C     Expand interval by 4 if minimum is still not bracketed.

      D = - 4.0E+0 * A
      GO TO 230

C     B is less than or equal to SCXBD.

  190 D1 = A

C     Determine interval of length D2 in which to set artificial bound.

      D2 = B
      IF (ABS (D2) .GT. TOL .AND.
     >  ((W .LE. 0.0E+0) .OR. (ABS (D1) .LE. TOL))) GO TO 200

C     ABS (D2) is less than or equal to TOL, or W is positive and
C     ABS (D1) is greater than TOL.  In either case, interchange D1, D2.

      U = D1
      D1 = D2
      D2 = U

C     Use parabolic interpolation only if new point lies in (A1,B1).

  200 U = - D1/ D2
      IF (U .LT. 1.0E+0) GO TO 210

C     Extrapolation - U exceeds 1.

      FACT = 5.0E+0 * (0.1E+0 + 1.0E+0/ U)/ 11.0E+0
      GO TO 220

C     Interpolation step - U less than 1.

  210 FACT = 0.5E+0 * SQRT (U)
  220 D = D2 * FACT

C     If D > 0 then B1 = D else A1 = D.

  230 IF (D .LE. 0.0E+0) GO TO 240
      B1 = D
      GO TO 250
  240 A1 = D
  250 IF (ABS (S) .GE. ABS (0.5E+0 * Q * R) .OR. (S .LE. (Q * A1))
     >  .OR. (S .GE. (Q * B1))) GO TO 260

C     A parabolic interpolation step.

      D = S/ Q

C     F must not be evaluated too close to A or B.

      IF ((D - A) .GE. T2 .AND. (B - D) .GE. T2) GO TO 10

C     Otherwise, set D to plus or minus TOL.

      D = SIGN (TOL, XM)
      GO TO 10

C     A non-interpolation step.

  260 NPCNT = NPCNT + 1
      IF (XM .LE. 0.0E+0) GO TO 270
      E = B
      GO TO 10
  270 E = A
      GO TO 10


C     Check safeguards.
C     -----------------

  280 CONTINUE
      SUCCES = .FALSE.

C     Check that new point satisfies safeguard conditions.  If the
C     function did not decrease, or step to the minimum was less than
C     SFTBND, LINSCH has failed to locate an acceptable minimum.

  290 CONTINUE
         IF (XMIN + B .LE. SFTBND) GO TO 330
         IF (FOLD - FMIN .GT. GTEST1 * XMIN) GO TO 310
         IF (NLIN .GE. MAXLIN) GO TO 330

C        A sufficiently lower point was not found - try halving step.

         XMIN = XMIN * 0.5E+0
         IF (XMIN .LE. T) GO TO 330
         DO 300 I = 1, N
            Z(I) = X(I) + XMIN * P(I)
  300    CONTINUE
         IF (PRINT) WRITE (LWRITE, 1030) XMIN
         CALL FUN (N, Z, FMIN)
         NLIN = NLIN + 1
         GO TO 290

C     A sufficiently lower point was found - set output values.

  310 CONTINUE
      SUCCES = .TRUE.
      ALPHA = XMIN
      IF (SCXBD .LE. 0.0E+0) ALPHA = STEPMX
      DO 320 I = 1, N
         X(I) = X(I) + ALPHA * P(I)
  320 CONTINUE
      F = FMIN


C     Termination.
C     ------------

  330 CONTINUE
      NFTOTL = NFTOTL + NLIN

      RETURN
      END
*DECK OUTPUT
C+----------------------------------------------------------------------
C
      SUBROUTINE OUTPUT (N, NLDIM, X, G, H, P, FNEW, ALPHA, D,
     >   EL, NITER, NFTOTL, FINAL, LWRITE)
C
C
C        This routine prints details of an iteration by QNMDIF.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT REAL (A - H, O - Z)

C     Arguments.

      LOGICAL
     >   FINAL
      REAL
     >   D(N), EL(NLDIM), G(N), H(N), P(N), X(N)


C     Execution.
C     ----------

      IF (FINAL) GO TO 30
         WRITE (LWRITE, 1000) NITER, FNEW
 1000    FORMAT (13H0AT ITERATION, I5, 22H THE FUNCTION VALUE IS,
     >      E25.15)
         WRITE (LWRITE, 1010)
 1010    FORMAT (44H0          CURRENT SOLUTION         GRADIENT,
     >      50H                 SEARCH DIRECTION         STEPSIZE)
         DO 10 I = 1, N
            WRITE (LWRITE, 1020) I, X(I), G(I), P(I), H(I)
 1020       FORMAT (1H , I5, 4E25.15)
   10    CONTINUE
         WRITE (LWRITE, 1030) ALPHA, NFTOTL
 1030    FORMAT (19H0LINEAR SEARCH STEP, E25.15, 5X,
     >      50HNUMBER OF FUNCTION EVALUATIONS AT END OF ITERATION, I9)
         WRITE (LWRITE, 1040)
 1040    FORMAT (40H0CHOLESKY FACTORS OF APPROXIMATE HESSIAN/
     >      31H0   ELEMENTS OF DIAGONAL MATRIX)
         WRITE (LWRITE, 1050) (D(I), I = 1, N)
 1050    FORMAT (1H , 5E25.15)
         IF (N .EQ. 1) GO TO 60
         WRITE (LWRITE, 1060)
 1060    FORMAT (27H0   LOWER TRIANGULAR FACTOR)
         K2 = 0
         NM1 = N - 1
         DO 20 I = 1, NM1
            K1 = K2 + 1
            K2 = K2 + I
            WRITE (LWRITE, 1050) (EL(J), J = K1, K2)
   20    CONTINUE
         GO TO 60
   30 CONTINUE

CRAK     Final iteration gets special treatment.

         WRITE (LWRITE, 1070) NITER, FNEW
 1070    FORMAT (/16H0FINAL ITERATION, I6, 5X, 14HFUNCTION VALUE,
     >      E25.15)
         WRITE (LWRITE, 1080)
 1080    FORMAT (44H0          FINAL SOLUTION           GRADIENT)
         DO 40 I = 1, N
            WRITE (LWRITE, 1020) I, X(I), G(I)
   40    CONTINUE
         WRITE (LWRITE, 1090) NFTOTL
 1090    FORMAT (31H0NUMBER OF FUNCTION EVALUATIONS, I9)
         WRITE (LWRITE, 1040)
         WRITE (LWRITE, 1050) (D(I), I = 1, N)
         IF (N .EQ. 1) GO TO 60
         WRITE (LWRITE, 1060)
         K2 = 0
         NM1 = N - 1
         DO 50 I = 1, NM1
            K1 = K2 + 1
            K2 = K2 + I
            WRITE (LWRITE, 1050) (EL(J), J = K1, K2)
   50    CONTINUE
   60 CONTINUE

CRAK  Compute and print estimate of condition number of Hessian and
CRAK  norm of gradient (all iterations).

      DMAX = D(1)
      DMIN = D(1)
      SUM = 0.E+0
      DO 70 I = 1, N
         IF (D(I) .GT. DMAX) DMAX = D(I)
         IF (D(I) .LT. DMIN) DMIN = D(I)
         SUM = SUM + G(I) ** 2
   70 CONTINUE
      BOUNDK = DMAX/ DMIN
      GNORM = SQRT (SUM)
      WRITE (LWRITE, 1100) BOUNDK, GNORM
 1100 FORMAT (43H0LOWER BOUND ON CONDITION NUMBER OF HESSIAN, E25.15/
     >   17H0NORM OF GRADIENT, E25.15)


C     Termination.
C     ------------

      RETURN
      END
*DECK UPCHOL
C+----------------------------------------------------------------------
C
      SUBROUTINE UPCHOL (N, NLDIM, EPSMCH, D, EL, Z)
C
C
C        The subroutine UPCHOL forms the updated Cholesky factorization
C     of the matrix EL * D * ELT - Z * ZT, where EL is a unit lower
C     triangular matrix, D is diagonal matrix with positive elements,
C     ELT is the transpose of EL, Z is a vector, and ZT is its transpose
C     (so that -Z * ZT is a negative rank-one correction).
C
C        The algorithm used is described in Gill, Golub, Murray, and
C     Saunders, "Methods for Modifying Matrix Factorizations", Stanford
C     Computer Science Report CS-72-322, Stanford University, Dec. 1972.
C
C        N is the dimension of the matrices and vectors, NLDIM is
C     N * (N-1)/ 2, the length of the vector EL in which the strict
C     lower triangle of the matrix EL is stored by rows, EPSMCH is
C     "machine epsilon", used to ensure that the resulting matrix
C     is sufficiently positive definite, D is a vector of length N
C     containing the elements of the matrix D, EL is a vector of length
C     NLDIM containing the elements of EL, and Z is the vector to be
C     used in the update.
C
C       The vectors EL and D are overwritten with the updated Cholesky
C    factors, and the elements of Z are overwritten.  UPCHOL modifies
C    the Cholesky factors of a matrix that could be indefinite, but it
C    ensures that the new matrix is positive definite.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT REAL (A - H, O - Z)

C     Constants.

      INTEGER
     >   MAXN
      PARAMETER
     >   (MAXN = 40)

C     Arguments.

      REAL
     >   D(N), EL(NLDIM), Z(N)

C     Local variables.

      REAL
     >   P(MAXN)


C     Execution.
C     ----------

C     Solve EL * P = Z.

      J = 1
      P(1) = Z(1)
      GAMMA = P(1) ** 2/ D(1)
      IF (N .EQ. 1) GO TO 30
         DO 20 I = 2, N
            K = I - 1
            T = Z(I)
            DO 10 IB = 1, K
               T = T - P(IB) * EL(J)
               J = J + 1
   10       CONTINUE
            P(I) = T
            GAMMA = GAMMA + T * T/ D(I)
   20    CONTINUE
   30 CONTINUE

C     If 1.0 - GAMMA < EPSMCH, then the modified matrix is not sufficiently
C     positive definite.  GAMMA is replaced by a quantity which ensures
C     that the modified matrix is positive definite regardless of subsequent 
C     rounding error.

      GAMMA = 1.0E+0 - GAMMA
      IF (GAMMA .GT. EPSMCH) GO TO 60
         IF (- GAMMA .GT. EPSMCH) GO TO 40
            GAMMA = EPSMCH
            GO TO 50
   40    CONTINUE
            GAMMA = - GAMMA
   50    CONTINUE
   60 CONTINUE
      K = J + N + N

C     Solve D * ELTRANSPOSE * Z = P.

      DO 90 JJ = 1, N
         J = N + 1 - JJ
         PJ = P(J)
         DJ = D(J)
         T = PJ/ DJ
         Z(J) = PJ
         BETA = - T/ GAMMA
         G = GAMMA + PJ * T
         D(J) = DJ * GAMMA/ G
         GAMMA = G
         K = K - J - 1
         IQ = K
         IF (J .EQ. N) GO TO 80
            JP1 = J + 1
            DO 70 IB = JP1, N
               T = EL(IQ)
               EL(IQ) = T + BETA * Z(IB)
               Z(IB) = Z(IB) + PJ * T
               IQ = IQ + IB - 1
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE


C     Termination.
C     ------------

      RETURN
      END

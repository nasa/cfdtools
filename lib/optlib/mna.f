C+----------------------------------------------------------------------
C
      SUBROUTINE MNA (N, NLDIM, F, X, L, D, CONV, FGH, APROXH, MXITER)
C
C     PURPOSE:
C
C     MNA implements a Modified Newton Algorithm for the unconstrained
C     optimization of a function of N variables.  It requires analytic
C     first derivatives of the function via a user-supplied subroutine
C     and it is also expected that the Hessian matrix of second deriv-
C     atives is known. However, subroutine HESS is supplied here which
C     provides a good approximation to the Hessian if the exact one is
c     unknown, or difficult to compute.
C
C     This is an adaptation of ALGOL procedure MNA,  documented in the
C     report below, to which the user is referred:
C
C     REFERENCE:  "The Implementation of Two Modified Newton
C                 Algorithms for Unconstrained Optimization"
C                 by GILL, Murray, and Picken (August 1972).
C                 (Report DNAC 24 of the  National  Physical
C                 Laboratory, England.)
C
C     The routine seeks the point X at which the twice continuously
C     differentiable function F(X) of N independent variables takes
C     its minimum value. Ideally, the variables should be scaled so
C     that the Hessian matrix at the solution is approximately row-
C     equilibrated,  with the minimum of F(X) having a value in the
C     range (-1, +1), and the maximum of F(X) in a unit sphere sur-
C     rounding the minimum in the range (-2, +2).  It  may  not  be
C     possible to fulfill either of these requirements.
C
C     Given an initial approximation to the minimum, MNA computes a
C     lower function value at each iteration.  When the convergence
C     criteria are satisfied, the routine gives the estimated posi-
C     tion of the minimum, the final function value, and  the final
C     Cholesky factorization of the approximate Hessian matrix.
C
C     The arguments and COMMON variables  of  subroutine  MNA   are
C     initialized in the user's driving program.  A description  of
C     these is given below.
C
C     INPUT TO MNA:
C
C     N     : The number of variables X(*) of the function F(X).
C     NLDIM : Dimension of array  L  used to store (row-by-row) the
C             lower triangle of the Hessian matrix. NLDIM=N*(N-1)/2
C     FGH   : A subroutine which must compute function value  F and
C             gradient vector G.  It may also compute exact  second
C             derivatives, but see APROXH below.
C             The call  FGH ( N, X, F, G, L, D, NTYPE ) will return
C             as follows:
C             If  NTYPE = 2,  only function and gradient need to be
C                             computed.
C             If  NTYPE = 3,  either the Hessian will  be  returned
C                             in arrays L and D, or a call will  be
C                             made to subroutine HESS (q.v.).
C     APROXH: This logical variable should be set to .TRUE.  if  no
C             2nd derivatives are available.   The subroutine  HESS
C             will then be invoked to provide approximations.  Note
C             that experience has shown that results will be almost
C             the same either way.
C             If APROXH = .FALSE.,  then subroutine FGH must do the
C             same as calling HESS(N,NLDIM,X,L,D) which assigns the
C             lower triangle of the matrix of  second  derivatives,
C             stored by rows, to the array L(I), I=1,2,..,N*(N-1)/2
C             and its diagonal is assigned to the array D(*).
C     X     : An initial estimate of the solution.
C     L, D  : Arrays not initialized by the user.
C             Note that all other local arrays, here and elsewhere,
C             should be dimensioned (at least N) by the user.
C     MXITER: Maximum number of iterations acceptable by user.
C     DEPS  : The relative machine precision.
C     ETA   : The termination criterion for the linear search.   It
C             should have a value in the range 0. to 1.  The closer
C             to zero it is, the greater the number of  evaluations
C             of the function that will  be  performed,  while  the
C             closer to 1 it is, the greater the number  of  itera-
C             tions likely.  A value of 0.9 is suggested  to  start
C             with, though if  approximate  second  derivatives are
C             used, the optimum choice decreases with increasing N.
C             Recommended values of ETA are 0.5 for N < 10, 0.1 for
C             10 <= N < 20, and 0.01 for N >= 20.
C     TOL   : The overall termination criterion (norm of the gradi-
C             ent).  A typical value is 1.E-6.  (A good estimate is
C             approximately SQRT (DEPS),  but this  can be  relaxed
C             (increased) if fewer significant figures are  accept-
C             able.)
C     STEPMX: An upper bound on the step allowed along a  direction
C             of search.  This can be used to prevent  overflow  in
C             in the computation.  If an approximate solution isn't
C             known,  and overflow in computing the function is un-
C             likely, then STEPMX can be set very large (1.E+11) so
C             that it will not influence the algorithm  at all.  If
C             the solution is known to be within a certain range of
C             of the initial estimate,  then STEPMX can be lowered.
C     NFMAX : Max. no. of fn. evals. permitted on any 1 line search
C     PRINT : Set this logical variable to true if complete  output
C             is desired after every iteration.
C     BRIEF : Set this .TRUE. if more concise output is preferred.
C
C     OUTPUT FROM MNA:
C
C     NFEVAL: The total number of function evaluations used.
C     NGEVAL: The total number of gradient evaluations used.
C     ICOUNT: The total number of linear searches performed.
C     CONV  : A logical variable set to .TRUE. if termination  occurs
C             with the convergence criteria satisfied, and .FALSE. if
C             a lower point cannot be found along a particular direc-
C             tion of search.
C     Also output are final values for F, X(*), D(*), and L(*).
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL   (A-H, L, O-Z )
C
      DIMENSION       X(N), L(NLDIM), D(N)
C
      LOGICAL         PRINT, BRIEF, CONV, SUCCES, APROXH
      COMMON /CRITER/ ETA, TOL, STEPMX, DEPS, NFMAX
      COMMON /NUMBRS/ NUMF, NUMG
      COMMON /COUNTS/ NFEVAL, NGEVAL, ICOUNT
      COMMON /LOGICL/ BRIEF, PRINT
C
      INTEGER         LUNPR
      PARAMETER      (LUNPR=6, MAXN=8, MXASCH=5)
      DIMENSION       G(MAXN), P(MAXN)
      EXTERNAL        FGH
C
C
      RTEPS = SQRT (DEPS)
      DEPSQ = DEPS*DEPS
      CONV  = .TRUE.
      NFEVAL = 0
      NGEVAL = 0
      ICOUNT = 0
      IASCH  = 0
      NUMF   = 0
      ALPHA  = 1.
      NP1    = N+1
C
C                 START OF MAIN ITERATION:
    1 CONTINUE
      IF (NUMF.EQ.1)  GOTO 5
C
      CALL FGH (N, X, FNEW, G, L, D, 3)
C
      NFEVAL = NFEVAL+1
      NGEVAL = NGEVAL+1
    5 IF (.NOT. APROXH)  GOTO 12
C
      DO 10 I=1,N
         D(I) = G(I)
   10 CONTINUE
C
      CALL HESS (N, NLDIM, X, L, D, FGH)
C
   12 CALL FORMCH (N, NLDIM, L, D, DEPS, IQ)
C
      GNORM = 0.
      DO 15 I=1,N
         GNORM = GNORM + G(I)**2
   15 CONTINUE
C
      GNORM = SQRT (GNORM + DEPSQ)
C
C             OVERALL CONVERGENCE CRITERION:
C
      IF (ICOUNT.GE.MXITER .AND. GNORM.GT.TOL*10.) GO TO 998
      IF (GNORM.GT.TOL) GO TO 20
      IF (IQ.EQ.0)      GO TO 999
      GO TO 100
C
C  *  SOLVE  LDL@P = -G  FOR DIRECTION OF SEARCH, P:
   20 P(1) = -G(1)
      IR = 1
      DO 50 I=2,N
         SUM = -G(I)
         IT  = I-1
         DO 40 K=1,IT
            SUM = SUM - P(K)*L(IR)
            IR  = IR+1
   40    CONTINUE
         P(I) = SUM
   50 CONTINUE
C
      P(N) = P(N)/D(N)
      DO 70 I=2,N
         IN = NP1-I
         IR = IR-1
         IS = IR
         IT = IN+1
         ITN= IT+N
         SUM= P(IN)/D(IN)
         DO 60 K=IT,N
            KN = ITN-K
            SUM= SUM - P(KN)*L(IS)
            IS = IS+2-KN
   60    CONTINUE
         P(IN) = SUM
   70 CONTINUE
C
C           COMPUTE THE NORM OF P AND DO LINEAR SEARCH:
   75 PNORM = 0.0
      DO 80 I=1,N
         PNORM = PNORM + P(I)**2
   80 CONTINUE
      PNORM = SQRT (PNORM + DEPSQ)
C
      IF (BRIEF) WRITE (LUNPR,91)
     +   ICOUNT,FNEW,X(1),X(N),GNORM,PNORM,ALPHA,NFEVAL
   91 FORMAT (1X, I5, 1P, 3E24.15, 3E11.2, I5)
C
      IF (PRINT) CALL DETAILS (N,NLDIM,X,G,P,D,L,FNEW,ALPHA)
C
      ALPHA = 1.
C
      CALL LNSRCH (N, NLDIM, FNEW, ALPHA, ETA, RTEPS/PNORM,
     +             STEPMX/PNORM, RTEPS, NFMAX, SUCCES,
     +             P, X, G, FGH, L, D)
C
      NFEVAL = NFEVAL + NUMF
      NGEVAL = NGEVAL + NUMG
      ICOUNT = ICOUNT + 1
      IF (.NOT.SUCCES)  GOTO 998
      GOTO 1
C
  100 IASCH = IASCH+1
      IF (IASCH .LE. MXASCH)  GOTO 101
      WRITE (LUNPR,103)
      GOTO 998
  101 WRITE (LUNPR,102)
  102 FORMAT ('0ALTERNATIVE SEARCH STARTED')
  103 FORMAT ('0BALE OUT - - TOO MANY ALTERNATIVE SEARCHES')
C
      IF (BRIEF.OR.PRINT) CALL DETAILS (N,NLDIM,X,G,P,D,L,FNEW,ALPHA)
      IT = IQ*(IQ-1)/2 + 1
      DO 110 IN=1,N
         I = NP1-IN
         IF ( I.GT.IQ )  GOTO 109
         IS = IT
         IT = IT-1
         SUM = 0.0
         IF ( I.EQ.IQ )  SUM = 1.
         I1 = I+1
         IF ( I1.GT.IQ )  GOTO 108
         DO 105 IK=I1,IQ
            K = IQ+I1-IK
            SUM = SUM - P(K)*L(IS)
            IS = IS+2-K
  105    CONTINUE
  108    P(I) = SUM
         GOTO 110
  109    P(I) = 0.0
  110 CONTINUE
C
      GTP = 0.0
      DO 120 I=1,N
         GTP = GTP + G(I)*P(I)
  120 CONTINUE
C
C        ASCERTAIN DOWNHILL DIRECTION FOR LINEAR SEARCH:
      IF ( GTP.LE.0.0 )  GOTO 75
      DO 125 I=1,N
         P(I) = -P(I)
  125 CONTINUE
      GOTO 75
C
C                 FAILED TO CONVERGE:
  998 CONV = .FALSE.
C
  999 F = FNEW
      IF (BRIEF.OR.PRINT) CALL DETAILS (N,NLDIM,X,G,P,D,L,FNEW,ALPHA)
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE CKGRAD (N, X, FGH)
C
C     PURPOSE: CKGRAD checks possible errors in the FGH-type routines
C              (function, gradient, Hessian) used by routine MNA.  It
C              does the necessary finite differencing to compare with
C              the analytic derivatives.
C     INPUTS:
C        N     Number of variables
C        X     Current values of the variables
C        FGH   User-supplied routine for calculating function, etc.
C
C     AUTHOR:  David Saunders, Informatics, 07/08/75.
C
C-----------------------------------------------------------------------
C
      REAL          X(N)
      PARAMETER    (LUNPR=6, MAXN = 8, MAXL = MAXN*(MAXN-1)/2)
      REAL          D(MAXN),G(MAXN),GXH(MAXN),L(MAXL)
C
      LOGICAL       LZERO
C
      EXTERNAL      FGH
C
      DATA          H/1.E-7/
C                   H is machine dependent.  Try EPS**0.5 or EPS**0.33,
C                   where EPS is the machine precision.
C
      CALL FGH (N, X, F, G, L, D, 2)
C
      WRITE (LUNPR,20) F
 20   FORMAT('0 F(X) =',1P,E23.15,/'0 I',11X,'X(I)',20X,'X(I)+H',19X,
     +       'F(X+H)',18X,'G-TAYLOR',16X,'G-ROUTINE')
C
      DO 50 I=1,N
         XI = X(I)
         LZERO = ABS(XI).LT.H
         X(I) = XI*(1.+H)
         IF (LZERO) X(I) = XI+H
C
         CALL FGH (N, X, FXH, GXH, L, D, 2)
C
         DEL = XI*H
         IF (LZERO) DEL = H
         GTAYLR = (FXH - F) / DEL
         WRITE (LUNPR,40) I,XI,X(I),FXH,GTAYLR,G(I)
 40      FORMAT(1X,I2,1P,5E25.15)
         X(I) = XI
 50   CONTINUE
C
C  *  CHECK 2ND DERIVATIVES ALSO:

      NLDIM = N*(N-1)/2
      CALL FGH (N, X, F, G, L, D, 3)
C
      WRITE (LUNPR,60)
      WRITE (LUNPR,70) (D(I),I=1,N)
      WRITE (LUNPR,80) (L(I),I=1,NLDIM)
C
      DO 55 I=1,N
         D(I)=G(I)
 55   CONTINUE
C
      CALL HESS (N, NLDIM, X, L, D, FGH)
C
      WRITE (LUNPR,70) (D(I),I=1,N)
      WRITE (LUNPR,80) (L(I),I=1,NLDIM)
C
 60   FORMAT('02ND DERIVS BY FUNCTION AND BY HESS ROUTINES....')
 70   FORMAT('0D =',/(1P,5E25.15))
 80   FORMAT('0L =',/(1P,5E25.15))
C
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE FORMCH (N, NLDIM, L, D, DEPS, NSPHI)
C
C     FORMCH (FORM CHolesky factors) is called by subroutine MNA.
C     Given an N*N symmetric matrix G, it forms the Cholesky  factor-
C     ization  LDL'  of a positive definite matrix  G + E, where L is
C     unit lower triangular,  stored with the unit diagonal  omitted,
C     and D is a diagonal matrix.  Matrix E is diagonal and is ident-
C     ically zero when G is sufficiently positive-definite.  Matrix G
C     is stored with its diagonal in the 1*N array D(I) and its lower
C     triangle is stored row-wise in the 1*N(N-1)/2 array L(I).   The
C     factorization is over-written on L and D.
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL   ( A-H, L, O-Z )
      DIMENSION       L(NLDIM), D(N)
C
C  *  Find "BOUND", which is the max. of ABS (L(I))/N and ABS (D(I)):
      J = N*(N-1)/2
      BOUND = -1.
      DO 20 I=1,J
         BOUND = AMAX1 (ABS (L(I)), BOUND)
 20   CONTINUE
      BOUND = BOUND/N
      DO 30 I=1,N
         BOUND = AMAX1 (ABS (D(I)), BOUND)
 30   CONTINUE
C
      Z    = 1./BOUND
      M    = 1
      DMAX = 0.
      NSPHI= 0
C
      DO 100 J=1,N
         U = D(J)
         IQ = M
         IF ( J.EQ.1 )  GOTO 50
         IT = J-1
         DO 40 K=1,IT
            V = L(M)
            W = V/D(K)
            L(M) = W
            M = M+1
            U = U - W*V
   40    CONTINUE
C
C  *     Intermediate D(J) are modified such that D(J)>=DEPS:
   50    W = U
         IF ( U.GE.DEPS )  GOTO 60
         W = ABS(U)
         IF ( W.LT.DEPS )  W = DEPS
         V = ABS( W - U )
         IF ( V.LE.DMAX )  GOTO 60
         DMAX = V
         NSPHI = J
   60    V = 0.0
         IS = M
         IF ( J.EQ.N )  GOTO 99
         JP1 = J+1
         DO 80 I=JP1,N
            U = 0.0
            IR = IQ
            IF ( J.EQ.1 )  GOTO 71
            IT = J-1
            DO 70 K=1,IT
               U = U - L(IR)*L(IS)
               IR = IR+1
               IS = IS+1
   70       CONTINUE
   71       U = U + L(IS)
            L(IS) = U
            IS = IS + I - J
            U = ABS(U)
            IF ( U.GT.V )  V = U
   80    CONTINUE
C
C  *     D(J) are modified such that ABS (SQRT (D(J)) * any element
C  *     of the Jth column of L ) <= SQRT (BOUND):
   99    V = V*V*Z
         D(J) = AMAX1(W,V)
  100 CONTINUE
C
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE HESS (N, NLDIM, X, L, D, FGH)
C
C     Subroutine HESS  may be used to approximate the  N * N  Hessian
C     matrix at the point X, required by subroutine MNA. The approxi-
C     mation is obtained by using finite differences.    Notice  that
C     the gradient at the point X is assumed to be  stored  in  array
C     D(I).    After the call, the approximate Hessian is stored with
C     its diagonal in D(I),  and its lower triangle in array L(I), I=
C     1,2,...,N*(N-1)/2, (stored row-by-row). The array X(I) is  left
C     as found on entry.  H = SQRT(EPS) approx. is machine dependent.
C     The quantity NGEVAL output by MNA does not include the gradient
C     evaluations performed in HESS.  The total number of such evalu-
C     ations used is actually  NGEVAL + N*NITERS.
C
C-----------------------------------------------------------------------
C
      REAL          X(N),L(NLDIM),D(N)
      PARAMETER     (MAXN=8)
      REAL          GXH(MAXN)
C
      DATA          H/1.E-7/
C                   H is machine-dependent.  Try SQRT(machine epsilon).
C
      IP = 0
      HH = 1./(H+H)
      DO 50 I=1,N
         IP = IP+I
         IQ = IP
         XI = X(I)
         X(I) = XI+H
         CALL FGH (N, X, F, GXH, L, D, 2)
         X(I) = XI
         GI = D(I)
         D(I) = ( GXH(I) - GI )/H
         IP1 = I+1
         IF ( IP1.GT.N )  GOTO 31
         DO 30 J=IP1,N
            L(IQ) = GXH(J) - GI
            IQ = IQ-1+J
   30    CONTINUE
         IF ( I.EQ.1 )  GOTO 50
   31    J1 = I-1
         J2 = IP-J1-I
         DO 40 J=1,J1
            J2J = J2 + J
            L(J2J) = ( GXH(J) - GI + L(J2J) )*HH
   40    CONTINUE
   50 CONTINUE
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE LNSRCH (N, NLDIM, F, ALPHA, ETA, TAU, STEPMX,
     +                   RTEPS, NFMAX, SUCCES, P, X, G, FGH, L, DIAG)
C
C     Line search from Gill and Murray, translated by Margaret Wright.
C     Finds an approximation ALPHA to the point at which the function
C     F (X + ALPHA*P)  attains its minimum value along direction  P.
C     Uses successive cubic interpolation with safeguards.
C     See NPL Report NAC 37 (1974).
C
C     The values A and B define the interval of uncertainty within
C     which the minimum is located. XMIN and W are the points with
C     the lowest and second lowest function values obtained so far.
C     A, XMIN and W are initialized at the origin.
C
C-----------------------------------------------------------------------
C
      REAL          P(N), X(N), G(N), L(NLDIM), DIAG(N)
C
      COMMON/NUMBRS/NUMF,NUMG
C
      PARAMETER     (LUNPR=6, MAXN=8)
      REAL          X1(MAXN), G1(MAXN)
      LOGICAL       SUCCES
C
C  *  Initialize various values:
C
      NUMF=0
      OLDF=F
      FMIN=F
      FW=F
      GW=0.
      T =TAU*0.1
      DO 50 I=1,N
         GW = GW + G(I)*P(I)
 50   CONTINUE
      IF (GW .GT. 0.)  WRITE (LUNPR,51) GW
 51   FORMAT('0LNSRCH,51: Gradient does not correspond to a',
     +       ' descent direction.  G^P =', 1P, E12.2)
C
      GMIN=GW
      GTEST1=-1.E-4*GMIN
      GTEST2=-(ETA+RTEPS)*GMIN
      XMIN=0.
      W =0.
      RR=0.
      A =0.
      D =ALPHA
C
C  *  The initial step must be in a positive direction along p:
      IF ( D.LE.0. )  GOTO 500
      TOL=T
      B  =STEPMX + RTEPS*ABS(STEPMX) + T
      B1 =B
      E  =B
      SCXBD=STEPMX
C
C  *  Begin iteration loop.   D is the tentative step
C  *  to the point at which the function is to be evaluated.
 100  CONTINUE
      D = AMIN1 (D, SCXBD)
      IF (ABS (D) .LT. TOL)  GOTO 110
      U = D
      GOTO 120
 110  U = SIGN (TOL, D)
C
C  *  Evaluate function at new point.
C  *  Check to see if it has been evaluated twice at same point.
 120  OLDR = RR
      RR   = XMIN + U
      IF ( OLDR.EQ.RR )  GOTO 500
      DO 130 I=1,N
         X1(I) = X(I) + RR*P(I)
 130  CONTINUE
C
      IF ( NUMF.GE.NFMAX ) GOTO 560
C
      NTYPE = 2
      IF ( NUMF.EQ.0 )  NTYPE = 3
C
      CALL FGH (N, X1, FU, G1, L, DIAG, NTYPE)
      NUMF = NUMF+1
C
C  *  Compute projected gradient at new point:
      GU = 0.
      DO 140 I=1,N
         GU = GU + G1(I)*P(I)
 140  CONTINUE
C
C  *  Update A,B,W,XMIN.
C  *  Check if new value is lower than previous best:
      IF ( FU.GT.FMIN )  GOTO 190
C
C  *  The new value is less than or equal to FMIN.
C  *  The new point becomes the next origin, and all
C  *  other points are shifted accordingly.
      FW  =FMIN
      FMIN=FU
      GW  =GMIN
      GMIN=GU
      XMIN=XMIN+U
      A=A-U
      B=B-U
      W=-U
      SCXBD=SCXBD-U
      IF ( GU.GE.0. )  GOTO 160
      A=0.
      GOTO 170
 160  B=0.
 170  TOL=ABS(XMIN)*RTEPS + T
      DO 180 I=1,N
         G(I) = G1(I)
 180  CONTINUE
      GOTO 230
C
C  *  In the following, the new point exceeds the previous best.
C  *  The origin remains unchanged, but the new point may qualify
C  *  as W, the second best point.
 190  CONTINUE
      IF ( U.GE.0. )  GOTO 200
      A=U
      GOTO 210
 200  B=U
 210  IF ( FU.GT.FW .AND. W.NE.0. )  GOTO 220
      W =U
      FW=FU
      GW=GU
 220  U =0.
C
C  *  Check termination criteria:
 230  CONTINUE
      XM=(A+B)*0.5
      IF ( (B-A).LE.TOL+TOL )  GOTO 500
      IF ( ABS(GMIN).LE.GTEST2 .AND. FMIN.LT.OLDF )  GOTO 500
      R=0.
      Q=0.
      S=0.
      IF ( ABS(E).LE.TOL )  GOTO 300
C
C  *  Fit cubic through U and W:
      R=3.*(FMIN-FW)/W + GMIN + GW
      Q=R*R - GW*GMIN
      IF ( Q.GE.0. )  GOTO 250
      R=0.
      Q=0.
      GOTO 300
C
C  *  Compute minimum of fitted cubic:
 250  Q = SIGN (SQRT (Q), W)
      S = (GMIN-R-Q)*W
      Q = GW-GMIN+(Q+Q)
      IF ( Q.LE.0. )  GOTO 280
      S =-S
      GOTO 290
 280  Q =-Q
 290  R = E
      IF ( B1.NE.D .OR. B.LE.SCXBD )  E = D
C
C  *  Construct an artificial bound on the estimated step length:
 300  A1=A
      B1=B
      IF ( B.LE.SCXBD )  GOTO 310
C
C  *  B EXCEEDS SCXBD -- EXPAND INTERVAL:
      D =-4.*W
      IF ( D.GE.B )  D = SCXBD
      GOTO 370
C
C  *  Here, B is less than or equal to SCXBD.
 310  IF ( ( A.NE.0. .OR. W.GE.0. ) .AND.
     +     ( B.NE.0. .OR. W.LE.0. ) )  GOTO 360
C
C  *  The minimum is bracketed by the points W and 0:
      D1=W
      IF ( A.EQ.0. )  GOTO 320
      D2=A
      GOTO 330
 320  D2=B
 330  U =-D1/D2
      IF ( U.GE.1. )  GOTO 340
      UINC = SQRT(U)*0.5
      GOTO 350
 340  UINC = ( 0.5 + 5./U )/11.
 350  D = D2*UINC
      GOTO 370
C
 360  D =XM
 370  IF ( D.LE.0. )  GOTO 380
      B1=D
      GOTO 390
 380  A1=D
 390  CONTINUE
C
C  *  Use cubic interpolation only if the new point is in (A1,B1):
      IF (ABS (S) .GE. ABS (0.5*Q*R))   GOTO 410
      IF (S .LE. Q*A1 .OR. S .GE. Q*B1) GOTO 410
C
C  *  A cubic interpolation step:
      D = S/Q
C
C  *  The function must not be evaluated closer than TOL to A or B:
      IF ( D-A.GE.TOL .AND. B-D.GE.TOL )  GOTO 100
      D = SIGN(TOL,XM)
      GOTO 100
C
 410  E = B-A
      GOTO 100
C
C  *  Termination of search.
C  *  Check decrease in function value:
 500   IF ( OLDF-FMIN.LE.GTEST1*XMIN )  GOTO 540
C
C  *  A sufficiently lower point has been found:
      IF ( SCXBD.GT.0. )  GOTO 510
      ALPHA = STEPMX
      XMIN  = STEPMX
      GOTO 520
C
 510  ALPHA = XMIN
 520  DO 530 I=1,N
         X(I) = X(I) + XMIN*P(I)
 530  CONTINUE
      F = FMIN
      SUCCES = .TRUE.
      RETURN
C
C
C  *  Function has not decreased enough.
C  *  Try to find a suitable point by halving the step size:
 540  XMIN = XMIN*0.5
      IF ( XMIN.LE.T )  GOTO 560
      IF ( NUMF.GE.NFMAX ) GOTO 560
C
      DO 550 I=1,N
         X1(I) = X(I) + XMIN*P(I)
 550  CONTINUE
C
      CALL FGH (N, X1, FMIN, G, L, DIAG, 2)
      NUMF = NUMF+1
      GOTO 500
C
C  *  Too many function evaluations, or
C  *  XMIN is less than tolerance -- failure:
 560  SUCCES = .FALSE.
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE DETAILS (N, NLDIM, X, G, P, D, L, FNEW, ALPHA)
C
C     This routine outputs details of an iteration by MNA.
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL   ( A-H, L, O-Z )
      DIMENSION       X(N), G(N), P(N), D(N), L(NLDIM)
      COMMON /COUNTS/ NFEVAL, NGEVAL, ICOUNT
      COMMON /LOGICL/ BRIEF,PRINT
      LOGICAL         BRIEF,PRINT
      INTEGER         LUNPR
      PARAMETER       ( LUNPR=6 )
C
C
      WRITE (LUNPR,11) ICOUNT
      WRITE (LUNPR,12)
      WRITE (LUNPR,13) ( X(I),G(I),P(I), I=1,N )
      WRITE (LUNPR,14) FNEW,ALPHA,NFEVAL
      WRITE (LUNPR,16)
      WRITE (LUNPR,18) ( D(I), I=1,N )
      IF ( .NOT.PRINT )  RETURN
      WRITE (LUNPR,17)
      K2 = 0
      NM1= N-1
      DO 10 I=1,NM1
         K1 = K2+1
         K2 = K2+I
         WRITE (LUNPR,18) ( L(J), J=K1,K2 )
   10 CONTINUE
C
   11 FORMAT ('0STATUS AT ITERATION', I5)
   12 FORMAT ('0     CURRENT SOLUTION            GRADIENT',
     +        '             DIRECTION OF SEARCH')
   13 FORMAT (1P, 3E25.15)
   14 FORMAT ('0FUNCTION VALUE=',1P,E23.15,'   STEP=',E15.7,I8,
     +        ' FUNCTION EVALS')
   16 FORMAT ('0DIAGONAL ELEMENTS= ')
   17 FORMAT ('0SUBDIAGONALS OF L= ')
   18 FORMAT (1X,1P,10E12.3)
C
      RETURN
      END

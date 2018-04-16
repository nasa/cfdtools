C-----------------------------------------------------------------------
C
      SUBROUTINE CALERF (ARG, JINT, RESULT)
C
C     This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
C     for a real argument  x.  It contains three FUNCTION type
C     subprograms: ERF, ERFC, and ERFCX, and one SUBROUTINE type
C     subprogram, CALERF.  The calling statements for the primary
C     entries are:
C
C                   Y = ERF (X)
C
C                   Y = ERFC (X)
C     and
C                   Y = ERFCX (X)
C
C     The routine  CALERF  is intended for internal packet use only,
C     all computations within the packet being concentrated in this
C     routine.  The function subprograms invoke  CALERF  with the
C     statement
C
C          CALL CALERF (ARG, JINT, RESULT)
C
C     where the parameter usage is as follows:
C
C      Function                 Parameters for CALERF
C       call               ARG            JINT      Result
C
C     ERF (ARG)      ANY REAL ARGUMENT      0      ERF (ARG)
C     ERFC (ARG)     |ARG| < XBIG           1      ERFC (ARG)
C     ERFCX (ARG)    XNEG < ARG < XMAX      2      ERFCX (ARG)
C
C     The main computation evaluates near-minimax approximations
C     from "Rational Chebyshev approximations for the error function"
C     by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
C     transportable program uses rational functions that theoretically
C     approximate  erf(x)  and  erfc(x)  to at least 18 significant
C     decimal digits.  The accuracy achieved depends on the arithmetic
C     system, the compiler, the intrinsic functions, and proper
C     selection of the machine-dependent constants.
C
C-----------------------------------------------------------------------
C
C     Explanation of machine-dependent constants:
C
C     XMIN   = the smallest positive floating-point number.
C     XINF   = the largest positive finite floating-point number.
C     XNEG   = the largest negative argument acceptable to ERFCX;
C              the negative of the solution to the equation
C              2*exp(x*x) = XINF.
C     XSMALL = argument below which erf(x) may be represented by
C              2*x/sqrt(pi)  and above which  x*x  will not underflow.
C              A conservative value is the largest machine number X
C              such that   1.0 + X = 1.0   to machine precision.
C     XBIG   = largest argument acceptable to ERFC;  solution to
C              the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
C              W(x) = exp(-x*x)/[x*sqrt(pi)].
C     XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
C              machine precision.  A conservative value is
C              1/[2*sqrt(XSMALL)]
C     XMAX   = largest acceptable argument to ERFCX; the minimum
C              of XINF and 1/[sqrt(pi)*XMIN].
C
C     Approximate values for some important machines are:
C
C                             XMIN       XINF        XNEG     XSMALL
C
C     CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
C     CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
C     IEEE (IBM/XT,
C       SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
C     IEEE (IBM/XT,
C       SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
C     IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
C     UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
C     VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
C     VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
C
C
C                             XBIG       XHUGE       XMAX
C
C     CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
C     CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
C     IEEE (IBM/XT,
C       SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
C     IEEE (IBM/XT,
C       SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
C     IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
C     UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
C     VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
C     VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
C
C-----------------------------------------------------------------------
C
C     Error returns:
C
C     The program returns  ERFC = 0      for  ARG >= XBIG;
C
C                          ERFCX = XINF  for  ARG <  XNEG;
C        and
C                          ERFCX = 0     for  ARG >= XMAX.
C
C
C     Author: W. J. Cody
C             Mathematics and Computer Science Division
C             Argonne National Laboratory
C             Argonne, IL 60439
C
C     Latest modification by the author: March 19, 1990
C
C     Adjustments by D.A.Saunders, ELORET/NASA Ames, 02/02/01:
C
C        Source code was obtained from www.netlib.org/specfun/erf and
C        stripped of the CS/CD comments.  The REAL form here can be
C        compiled with the appropriate switch to give 32- or 64-bit
C        precision.  Machine-dependent constants are now supplied by
C        Fortran 90 intrinsics.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C-----------------------------------------------------------------------
C     Arguments:
C-----------------------------------------------------------------------

      REAL,    INTENT (IN)  ::
     >   ARG

      INTEGER, INTENT (IN)  ::
     >   JINT

      REAL,    INTENT (OUT) ::
     >   RESULT

C-----------------------------------------------------------------------
C     Local variables:
C-----------------------------------------------------------------------

      INTEGER
     >   I

      REAL
     >   DEL, FOUR, HALF, ONE, SIXTEN, SQRPI, TWO, THRESH,
     >   X, XDEN, XNUM, Y, YSQ, ZERO

      REAL
     >   A(5), B(4), C(9), D(8), P(6), Q(5)

      REAL, SAVE ::
     >   XBIG, XHUGE, XINF, XMAX, XMIN, XNEG, XSMALL

      LOGICAL, SAVE ::
     >   FIRST = .TRUE.

C     Extra variables needed by the Newton iteration for XBIG:

      INTEGER
     >   ITER, ITMAX
      REAL
     >   ALPHA, DV, F, FNORM, FNORM0, ROOTPI, T1, T2, TOL, V, VLAST, VSQ
      LOGICAL
     >   CONVGD, FAIL

C-----------------------------------------------------------------------
C     Mathematical constants:
C-----------------------------------------------------------------------

      DATA
     >   FOUR, ONE, HALF, TWO, ZERO /4.0E0, 1.0E0, 0.5E0, 2.0E0, 0.0E0/,
     >   SQRPI /5.6418958354775628695E-1/, THRESH /0.46875E0/,
     >   SIXTEN /16.0E0/

C-----------------------------------------------------------------------
C     Coefficients for approximation to  erf  in first interval:
C-----------------------------------------------------------------------

      DATA
     >   A /3.16112374387056560E00, 1.13864154151050156E02,
     >      3.77485237685302021E02, 3.20937758913846947E03,
     >      1.85777706184603153E-1/,
     >   B /2.36012909523441209E01, 2.44024637934444173E02,
     >      1.28261652607737228E03, 2.84423683343917062E03/

C-----------------------------------------------------------------------
C     Coefficients for approximation to  erfc  in second interval:
C-----------------------------------------------------------------------

      DATA
     >   C /5.64188496988670089E-1, 8.88314979438837594E0,
     >      6.61191906371416295E01, 2.98635138197400131E02,
     >      8.81952221241769090E02, 1.71204761263407058E03,
     >      2.05107837782607147E03, 1.23033935479799725E03,
     >      2.15311535474403846E-8/,
     >   D /1.57449261107098347E01, 1.17693950891312499E02,
     >      5.37181101862009858E02, 1.62138957456669019E03,
     >      3.29079923573345963E03, 4.36261909014324716E03,
     >      3.43936767414372164E03, 1.23033935480374942E03/

C-----------------------------------------------------------------------
C     Coefficients for approximation to  erfc  in third interval:
C-----------------------------------------------------------------------

      DATA
     >   P /3.05326634961232344E-1, 3.60344899949804439E-1,
     >      1.25781726111229246E-1, 1.60837851487422766E-2,
     >      6.58749161529837803E-4, 1.63153871373020978E-2/,
     >   Q /2.56852019228982242E00, 1.87295284992346047E00,
     >      5.27905102951428412E-1, 6.05183413124413191E-2,
     >      2.33520497626869185E-3/

C-----------------------------------------------------------------------
C     Execution:
C-----------------------------------------------------------------------

      IF (FIRST) THEN

C        Calculate precision-dependent quantities
C        (a minor price to pay for portability):

         FIRST  = .FALSE.
         XMIN   =  TINY (XMIN)
         XINF   =  HUGE (XINF)
         XNEG   = -SQRT (LOG (HALF * XINF))
         XSMALL =  EPSILON (XSMALL)
         XHUGE  =  HALF / SQRT (XSMALL)
C********XMAX   =  MIN (XINF, ONE / (XMIN * 1.772454E0))
         XDEN   =  0.8 * (XMIN * XINF) ! 0.8 < .5 sqrt(pi) avoids XMAX overflow
         XMAX   =  MIN (XINF, XDEN / XMIN)

C        Newton iteration for XBIG:  Solve the following for v:
C           f(v) = (exp(-v^2)) (v^2 - 0.5) / sqrt(pi) - XMIN v^3 = 0
C        where multiplying throughout by XINF helps avoid extremes, and we
C        use the existing variable SQRPI = 1/sqrt(pi).

         T1     = XMIN * XINF     ! O(1.0)
         ALPHA  = ZERO            ! For printout purposes only
         FNORM0 = 1.E+30          ! I.e., big to avoid iteration 0 test
         V      = -XNEG           ! Starting guess suggested by the above tables
         TOL    = 1.E-5 * ABS (V) ! No need for high precision, but must scale
         ITMAX  = 10              ! 32-bit arithmetic can hit the limit with
         CONVGD = .FALSE.         !    |f| stuck at ~1.E-3, but result is OK
         FAIL   = .FALSE.         ! For step-halving inner iteration

         DO ITER = 0, ITMAX

C           Evaluate f(v), the magnitude of which should converge to ~0.

            VSQ   = V * V
            T2    = SQRPI * XINF * EXP (-VSQ)
            F     = T2 * (VSQ - HALF) - V * VSQ * T1
            FNORM = ABS (F)

C           Halve the step until |f| is reduced (except first time through).

            DO WHILE (FNORM > FNORM0) 
               IF (ALPHA > TWO * XSMALL) THEN ! Assume XSMALL = machine epsilon
                  ALPHA = HALF * ALPHA
                  V = VLAST - ALPHA * DV
                  VSQ   = V * V
                  T2    = SQRPI * XINF * EXP (-VSQ)
                  F     = T2 * (VSQ - HALF) - V * VSQ * T1
                  FNORM = ABS (F)
               ELSE
                  WRITE (*, '(A)') ' CALERF: Step halving failed.'
                  FAIL = .TRUE. ! To the cease outer iteration
                  EXIT
               END IF
            END DO

C****       WRITE (*, '(A, I3, A, 1P, E9.2, A, E9.2, A, E14.6)')
C****>         ' CALERF:', ITER, '  |f|:', FNORM, '  step:', ALPHA,
C****>         '  v:', V

            IF (FNORM < TOL) CONVGD = .TRUE.

            IF (CONVGD .OR. FAIL) EXIT

C           Calculate step dv = f(v) / f'(v):

            FNORM0 = FNORM
            DV     = F / (T2 * V * (3. - TWO * VSQ) - 3. * VSQ * T1)
            VLAST  = V
            V      = V - DV
            ALPHA  = ONE

         END DO ! Next iteration

         IF (.NOT. CONVGD) THEN
            IF (DV > TOL .OR. FAIL) THEN
               WRITE (*, '(A)') ' CALERF: Iteration failed. XBIG <- 9.'
               V = 9.
            END IF
         END IF

         XBIG = V - V * TWO * XSMALL ! Play safe

C****    WRITE (*, '(A, 1P, E25.15)')
C****>   ' XMIN:  ', XMIN,
C****>   ' XINF:  ', XINF,
C****>   ' XNEG:  ', XNEG,
C****>   ' XSMALL:', XSMALL,
C****>   ' XHUGE: ', XHUGE,
C****>   ' XMAX:  ', XMAX,
C****>   ' XBIG:  ', XBIG

      END IF

      X = ARG
      Y = ABS (X)

      IF (Y <= THRESH) THEN ! Evaluate  erf  for  |X| <= 0.46875:

         YSQ = ZERO
         IF (Y > XSMALL) YSQ = Y * Y
         XNUM = A(5)*YSQ
         XDEN = YSQ
         DO I = 1, 3
            XNUM = (XNUM + A(I)) * YSQ
            XDEN = (XDEN + B(I)) * YSQ
         END DO
         RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
         IF (JINT /= 0) RESULT = ONE - RESULT
         IF (JINT == 2) RESULT = EXP (YSQ) * RESULT
         GO TO 800

      ELSE IF (Y <= FOUR) THEN ! Evaluate  erfc  for 0.46875 <= |X| <= 4.0:

         XNUM = C(9)*Y
         XDEN = Y
         DO I = 1, 7
            XNUM = (XNUM + C(I)) * Y
            XDEN = (XDEN + D(I)) * Y
         END DO
         RESULT = (XNUM + C(8)) / (XDEN + D(8))
         IF (JINT /= 2) THEN
            YSQ = AINT (Y*SIXTEN) / SIXTEN
            DEL = (Y - YSQ)*(Y + YSQ)
            RESULT = EXP (-YSQ*YSQ) * EXP (-DEL) * RESULT
         END IF

      ELSE ! Evaluate  erfc  for |X| > 4.0:

         RESULT = ZERO
         IF (Y >= XBIG) THEN
            IF ((JINT /= 2) .OR. (Y >= XMAX)) GO TO 300
            IF (Y >= XHUGE) THEN
               RESULT = SQRPI / Y
               GO TO 300
            END IF
         END IF
         YSQ = ONE / (Y * Y)
         XNUM = P(6)*YSQ
         XDEN = YSQ
         DO I = 1, 4
            XNUM = (XNUM + P(I)) * YSQ
            XDEN = (XDEN + Q(I)) * YSQ
         END DO
         RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
         RESULT = (SQRPI -  RESULT) / Y
         IF (JINT /= 2) THEN
            YSQ = AINT (Y*SIXTEN) / SIXTEN
            DEL = (Y - YSQ)*(Y + YSQ)
            RESULT = EXP (-YSQ*YSQ) * EXP (-DEL) * RESULT
         END IF

      END IF

C     Fix up for negative argument, erf, etc.:

  300 IF (JINT == 0) THEN
         RESULT = (HALF - RESULT) + HALF
         IF (X < ZERO) RESULT = -RESULT
      ELSE IF (JINT == 1) THEN
         IF (X < ZERO) RESULT = TWO - RESULT
      ELSE
         IF (X < ZERO) THEN
            IF (X < XNEG) THEN
               RESULT = XINF
            ELSE
               YSQ = AINT (X*SIXTEN) / SIXTEN
               DEL = (X - YSQ)*(X + YSQ)
               Y = EXP (YSQ*YSQ) * EXP (DEL)
               RESULT = (Y+Y) - RESULT
            END IF
         END IF
      END IF

  800 RETURN

      END SUBROUTINE CALERF

C-----------------------------------------------------------------------
C
      REAL FUNCTION ERF (X)
C
C     This subprogram computes approximate values for erf(x).
C     (See comments heading CALERF.)
C
C     Author/date: W. J. Cody, January 8, 1985
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

      REAL X
      REAL RESULT

      CALL CALERF (X, 0, RESULT)

      ERF = RESULT

      END FUNCTION ERF

C-----------------------------------------------------------------------
C
      REAL FUNCTION ERFC (X)
C
C     This subprogram computes approximate values for erfc(x).
C     (See comments heading CALERF.)
C
C     Author/date: W. J. Cody, January 8, 1985
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

      REAL X
      REAL RESULT

      CALL CALERF (X, 1, RESULT)

      ERFC = RESULT

      END FUNCTION ERFC

C-----------------------------------------------------------------------
C
      REAL FUNCTION ERFCX (X)
C
C     This subprogram computes approximate values for exp(x*x) * erfc(x).
C     (See comments heading CALERF.)
C
C     Author/date: W. J. Cody, March 30, 1987
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

      REAL X
      REAL RESULT

      CALL CALERF (X, 2, RESULT)

      ERFCX = RESULT

      END FUNCTION ERFCX

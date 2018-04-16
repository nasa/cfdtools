C+------------------------------------------------------------------------------
C
      SUBROUTINE ELLIP3D (IDM, JDM, KDM,  I1, I2, J1, J2, K1, K2,
     >                    X,      ARCS,   PRINT,  ITMAX,
     >                    CONVG,  DMAX,   OMG3D,  BGPHI, BGPSI, BGOMG,
     >                    FGMODE, FGLIMIT,FGRELAX,EXPI1, EXPI2, EXPJ1,
     >                    EXPJ2,  EXPK1,  EXPK2,  NNTIP, ILE,   KTIP)
C
C     DESCRIPTION:
C
C        ELLIP3D smooths the interior X,Y,Z coordinates of a volume (sub)grid.
C     Thomas-Middlecoff-type spacing control is imposed on the grid interior;
C     these terms are referred to as background terms.  Near the boundaries,
C     spacing and orthogonality are controlled by the method of Sorenson
C     (GRAPE and 3DGRAPE); these terms are referred to as the foreground terms.
C     All terms follow Sorenson's formulation of the coupled set of Poisson's
C     elliptic PDEs:
C
C        LHS = A11*Xxx + A22*Xee + A33*Xzz + 2(A12*Xxe + A13*Xxz + A23*Xez)
C        RHS = -J**2 (P*Xx + Q*Xe +R*Xz)
C     where
C        X is the 3D position vector of each point
C        x is a derivative in xsi (one of the computational directions)
C        y is a derivative in eta
C        z is a derivative in zeta
C        J is the determinant of the Jacobian matrix (first row: Xx Xe Xz)
C        Aij is the sum (l=1:3) of GAM(li)*GAM(lj)
C        GAM(li) is the li cofactor of the Jacobian matrix
C        P = PHI which is a blend of the foreground and background terms off
C            xsi faces
C        Q = PSI which is a blend of the foreground and background terms off
C            eta faces
C        R = OMEGA which is a blend of the foreground and background terms off
C            zeta faces
C
C        Solution of the nonlinear set of three PDEs is by SLOR (Successive Line
C     Over-Relaxation) in the xsi (i) direction.  Derivation of the background
C     terms appears in NUMERICAL GRID GENERATION by Thompson, Warsi, Mastin
C     (1985), p. 231.  The code appears slightly different in order to suit
C     Sorenson's formulation.  The foreground terms are discussed in Sorenson
C     (May 1980) NASA Tech. Mem. 81198.
C
C     LIMITATIONS:
C
C     >  The unusually large number of work-space arrays have not been made
C        arguments for practical reasons.  It is expected that Fortran 77
C        applications such as typical CFD flow solvers will have more than
C        enough existing work-space which can be reused here via a suitable
C        common block or include file.  Automatic arrays serve under Fortran 90.
C
C     >  Arguments NNTIP, ILE, KTIP are application-specific (for more than
C        one topology about a wing/body).  Check the code for their usage.
C
C     ARGUMENTS:
C
C        Suggested argument values apply to an Euler-type C grid about a wing,
C     smoothed in two halves.
C
C     PRINT:   .TRUE. shows smoothing iterations on unit 6
C     ITMAX:   maximum number of SLOR smoothing iterations;
C              suggested value ITMAX = 100-200
C     CONVG:   orders of magnitude drop in maximum correction before
C              terminating SLOR iteration;
C              suggested value CONVG = 4.
C     DMAX:    tolerance on maximum correction for terminating SLOR iteration;
C              suggested value DMAX = 0.001 (depends on units)
C     OMG3D:   the over-relaxation parameter for the SLOR iteration;
C              suggested value OMG3D = 1.1
C     BGPHI:   'YYYY' turns on background control of interior spacing from
C              the I faces (Phi effects) at edges J1, J2, K1, K2 resp.
C     BGPSI:   likewise for background control from J faces (Psi effects)
C              at edges I1, I2, K1, K2 resp.
C     BGOMG:   likewise for background control from K faces (Omg effects)
C              at edges I1, I2, J1, J2 resp.
C              suggestions:    BGPHI = 'YYYY'
C                              BGPSI = 'YYYY'
C                              BGOMG = 'YYYY'
C              set them to 'NNNN' or 'YYNN' etc. to turn background effects off
C     FGMODE:  activates foreground control at each face;
C              suggestions:    FGMODE = 'YYYYYY'
C              use N to turn off selected foreground effects
C     FGLIMIT: used to limit growth of P, Q, and R foreground terms;
C              combines with the curvature;
C              suggested value FGLIMIT = 1.0;
C              lower value may help with extreme curvature.
C     FGRELAX: used to under-relax changes in P, Q, and R foreground terms;
C              suggested value FGRELAX = 0.1;
C              lower value may help with extreme curvature.
C     EXPI1:   foreground decay factors at the indicated faces;
C     etc.     the decay is exponential and based on indices, not arc lengths;
C              suggestions:    EXPI1 = 0.5
C                              EXPI2 = 0.5
C                              EXPJ1 = 0.5
C                              EXPJ2 = 0.5
C                              EXPK1 = 0.6
C                              EXPK2 = 0.6
C              lower values hold foreground terms further into the interior;
C              larger values reduce foreground effects.
C
C        One must play with the DECAY knobs because these can have drastic
C     effects on the rate of convergence and the beauty of the grid.  This is
C     especially true if you are trying to force orthogonality in a "corner"
C     that is very skewed.  It is a good idea to decide which edge is most
C     important at the corner.  Tighten it, and loosen the other intersecting
C     edges.
C
C     HISTORY:
C
C     ????????  Lockheed     Original implementation as TTM3D from WBGRID.
C     09/13/94  D.Saunders   Replaced TRIB with TRID3R; allowed for zero
C                            iterations; moved test for itn. 1 residuals.
C     02/07/95    "   "      Removed tests for coincident boundary pts;
C                            replace IFs with MAX (MIN ( to vectorize.
C     04/22/95    "   "      J, K switched to match JL, KL usage.
C     05/01/95    "   "      NNTIP argument allows use on C-H + C-O grid.
C     05/08/95 J.Reuther/DAS Thomas-Middlecoff forcing functions installed
C                            for all six boundaries; interior interpolation
C                            revised along Gridgen3d lines.
C     05/21/95     DAS       3D source term limiters separated from 2D's;
C                            application at boundaries suffices;
C                            force PHI = 0. at I = ILE.
C     06/19/95      "        Arranged for application to subvolumes.
C     7/95-7/96  S.Edwards   ELLIP3D enhancement of TTM3D with Sorenson-type
C                            orthogonality control.
C     8/96         DAS       Some polishing.
C     6/97       SJE/DAS     Jacobian is needed only on the faces; merged
C                            SINTERP function into the GAM* routines; TFI-type
C                            interpolation of edge increments into upper face
C                            interiors had been using lower face coefs.
C     07/27/97      "        (1 - "s") factors added to foreground decay coefs.
C     04/27/98     DAS       Fortran 90 allows automatic arrays for local work-
C                            space; GAM* AND PQR* procedures are now internal.
C     05/18/98      "        I2/J2/K2 controls of background were those of
C                            I1/J1/K1 by mistake.
C ------------------------------------------------------------------------------

C     Arguments:

      INTEGER
     >   IDM, JDM, KDM, I1, I2, J1, J2, K1, K2, ITMAX, NNTIP, ILE, KTIP
      REAL
     >   X(IDM,JDM,KDM,3), ARCS(IDM,JDM,KDM,3),
     >   CONVG, DMAX, OMG3D, FGLIMIT, FGRELAX,
     >   EXPI1, EXPI2, EXPJ1, EXPJ2, EXPK1, EXPK2
      LOGICAL
     >   PRINT
      CHARACTER
     >   BGPHI*4, BGPSI*4, BGOMG*4, FGMODE*6

C-------------------------------------------------------------------------------

C     Local constants:

      CHARACTER, PARAMETER ::
     >   YES*1 = 'Y'
      REAL, PARAMETER ::
     >   EPS = 1.E-16, ! Safeguards 0/0
     >   FOURTH = 0.25, HALF = 0.5, ONE = 1., TWO = 2., ZERO = 0.

C     Automatic arrays:

      REAL, DIMENSION (I1:I2) ::
     >   A, B, C, RHS1, RHS2, RHS3, DECAYI1, DECAYI2

      REAL, DIMENSION (J1:J2) ::
     >   DECAYJ1, DECAYJ2

      REAL, DIMENSION (K1:K2) ::
     >   DECAYK1, DECAYK2

      REAL, DIMENSION (J1:J2,K1:K2,2) ::
     >   CURVI,  F1XSI,  F2XSI,  F3XSI,  FX1, FX2, FX3,
     >   GAMX11, GAMX12, GAMX13, GAMX21, GAMX22, GAMX23,
     >   GAMX31, GAMX32, GAMX33, PSII, OMGI, PX, QX, RX, RJACI

      REAL, DIMENSION (I1:I2,K1:K2,2) ::
     >   CURVJ,  F1ETA,  F2ETA,  F3ETA,  FE1, FE2, FE3,
     >   GAME11, GAME12, GAME13, GAME21, GAME22, GAME23,
     >   GAME31, GAME32, GAME33, OMGJ, PHIJ, PE, QE, RE, RJACJ

      REAL, DIMENSION (I1:I2,J1:J2,2) ::
     >   CURVK,  F1ZETA, F2ZETA, F3ZETA, FZ1, FZ2, FZ3,
     >   GAMZ11, GAMZ12, GAMZ13, GAMZ21, GAMZ22, GAMZ23,
     >   GAMZ31, GAMZ32, GAMZ33, PHIK, PSIK, PZ, QZ, RZ, RJACK

C     Since the grid is ordered X(I,J,K,*), each array above is dimensioned
C     to suit the way Fortran is traditionally more efficient for traversing
C     a corresponding plane of X.

C     Local variables:

      LOGICAL
     >   CO, FOREI, FOREJ, FOREK,
     >   FRCSI1, FRCSI2, FRCPJ1, FRCPJ2, FRCOK1, FRCOK2

C     Execution:

      IF (ITMAX == 0) GO TO 999

      SRCPJ1 = ONE
      IF (BGPHI(1:1) /= YES) SRCPJ1 = ZERO ! P = "PHI"
      SRCPJ2 = ONE
      IF (BGPHI(2:2) /= YES) SRCPJ2 = ZERO
      SRCPK1 = ONE
      IF (BGPHI(3:3) /= YES) SRCPK1 = ZERO
      SRCPK2 = ONE
      IF (BGPHI(4:4) /= YES) SRCPK2 = ZERO

      SRCSI1 = ONE
      IF (BGPSI(1:1) /= YES) SRCSI1 = ZERO ! S = "PSI"
      SRCSI2 = ONE
      IF (BGPSI(2:2) /= YES) SRCSI2 = ZERO
      SRCSK1 = ONE
      IF (BGPSI(3:3) /= YES) SRCSK1 = ZERO
      SRCSK2 = ONE
      IF (BGPSI(4:4) /= YES) SRCSK2 = ZERO

      SRCOI1 = ONE
      IF (BGOMG(1:1) /= YES) SRCOI1 = ZERO ! O = "OMEGA"
      SRCOI2 = ONE
      IF (BGOMG(2:2) /= YES) SRCOI2 = ZERO
      SRCOJ1 = ONE
      IF (BGOMG(3:3) /= YES) SRCOJ1 = ZERO
      SRCOJ2 = ONE
      IF (BGOMG(4:4) /= YES) SRCOJ2 = ZERO

      FRCSI1 = FGMODE(1:1) == YES
      FRCSI2 = FGMODE(2:2) == YES
      FRCPJ1 = FGMODE(3:3) == YES
      FRCPJ2 = FGMODE(4:4) == YES
      FRCOK1 = FGMODE(5:5) == YES
      FRCOK2 = FGMODE(6:6) == YES

C     Normalized arc-lengths for the starting-guess grid:

      CALL PARAMXYZ (1, IDM, 1, JDM, 1, KDM, I1, I2, J1, J2, K1, K2,
     >               X(1,1,1,1), X(1,1,1,2), X(1,1,1,3), ARCS)

      I1P1 = I1 + 1
      J1P1 = J1 + 1
      K1P1 = K1 + 1
      I2M1 = I2 - 1
      J2M1 = J2 - 1
      K2M1 = K2 - 1

C     Set up the background forcing functions at the faces, Jacobians,
C     and all the derivatives used in the foreground terms except
C     second derivatives off the faces.  Also calculate curvatures
C     at each face point to be used as foreground limiters.
C     Use first planes of PX, PE, PZ for SXSI, SETA, SZETA work-space.

      BACKSI = SRCSI1
      BACKOI = SRCOI1
      BACKPJ = SRCPJ1
      BACKOJ = SRCOJ1
      BACKPK = SRCPK1
      BACKSK = SRCSK1

      DO L = 1, 2

C        I face:

         CALL GAMXSI (L,               RJACI(J1,K1,L),
     >                GAMX11(J1,K1,L), GAMX12(J1,K1,L), GAMX13(J1,K1,L),
     >                GAMX21(J1,K1,L), GAMX22(J1,K1,L), GAMX23(J1,K1,L),
     >                GAMX31(J1,K1,L), GAMX32(J1,K1,L), GAMX33(J1,K1,L),
     >                FX1(J1,K1,L),    FX2(J1,K1,L),    FX3(J1,K1,L),
     >                F1XSI(J1,K1,L),  F2XSI(J1,K1,L),  F3XSI(J1,K1,L),
     >                PSII(J1,K1,L),   OMGI(J1,K1,L),   CURVI(J1,K1,L),
     >                PX,              BACKSI,          BACKOI)

C        J face:

         CALL GAMETA (L,               RJACJ(I1,K1,L),
     >                GAME11(I1,K1,L), GAME12(I1,K1,L), GAME13(I1,K1,L),
     >                GAME21(I1,K1,L), GAME22(I1,K1,L), GAME23(I1,K1,L),
     >                GAME31(I1,K1,L), GAME32(I1,K1,L), GAME33(I1,K1,L),
     >                FE1(I1,K1,L),    FE2(I1,K1,L),    FE3(I1,K1,L),
     >                F1ETA(I1,K1,L),  F2ETA(I1,K1,L),  F3ETA(I1,K1,L),
     >                PHIJ(I1,K1,L),   OMGJ(I1,K1,L),   CURVJ(I1,K1,L),
     >                PE,              BACKPJ,          BACKOJ)

C        K face:

         CALL GAMZETA (L,              RJACK(I1,J1,L),
     >                GAMZ11(I1,J1,L), GAMZ12(I1,J1,L), GAMZ13(I1,J1,L),
     >                GAMZ21(I1,J1,L), GAMZ22(I1,J1,L), GAMZ23(I1,J1,L),
     >                GAMZ31(I1,J1,L), GAMZ32(I1,J1,L), GAMZ33(I1,J1,L),
     >                FZ1(I1,J1,L),    FZ2(I1,J1,L),    FZ3(I1,J1,L),
     >                F1ZETA(I1,J1,L), F2ZETA(I1,J1,L), F3ZETA(I1,J1,L),
     >                PHIK(I1,J1,L),   PSIK(I1,J1,L),   CURVK(I1,J1,L),
     >                PZ,              BACKPK,          BACKSK)

         BACKSI = SRCSI2
         BACKOI = SRCOI2
         BACKPJ = SRCPJ2
         BACKOJ = SRCOJ2
         BACKPK = SRCPK2
         BACKSK = SRCSK2

      END DO

      CO = NNTIP /= 0 ! True  -> C-O for K >= KTIP (singular line at I=ILE);
                      ! False -> no special handling (C-H everywhere).

C     Evaluate index-based foreground decay terms off each face:

      I1POWER = NINT (-LOG10 (MIN (ARCS(I1P1,J1,K1,1),
     >                             ARCS(I1P1,J2,K1,1),
     >                             ARCS(I1P1,J1,K2,1),
     >                             ARCS(I1P1,J1,K2,1))))
      I1POWER = MAX ((I1POWER + 1)/ 2, 3)
      I2POWER = NINT (-LOG10 (ONE - MAX (ARCS(I2M1,J1,K1,1),
     >                                   ARCS(I2M1,J2,K1,1),
     >                                   ARCS(I2M1,J1,K2,1),
     >                                   ARCS(I2M1,J2,K2,1))))
      I2POWER = MAX ((I2POWER + 1)/ 2, 3)

      J1POWER = NINT (-LOG10 (MIN (ARCS(I1,J1P1,K1,2),
     >                             ARCS(I2,J1P1,K1,2),
     >                             ARCS(I1,J1P1,K2,2),
     >                             ARCS(I2,J1P1,K2,2))))
      J1POWER = MAX ((J1POWER + 1)/ 2, 3)
      J2POWER = NINT (-LOG10 (ONE - MAX (ARCS(I1,J2M1,K1,2),
     >                                   ARCS(I2,J2M1,K1,2),
     >                                   ARCS(I1,J2M1,K2,2),
     >                                   ARCS(I2,J2M1,K2,2))))
      J2POWER = MAX ((J2POWER + 1)/ 2, 3)

      K1POWER = NINT (-LOG10 (MIN (ARCS(I1,J1,K1P1,3),
     >                             ARCS(I2,J1,K1P1,3),
     >                             ARCS(I1,J2,K1P1,3),
     >                             ARCS(I2,J2,K1P1,3))))
      K1POWER = MAX ((K1POWER + 1)/ 2, 3)
      K2POWER = NINT (-LOG10 (ONE - MAX (ARCS(I1,J1,K2M1,3),
     >                                   ARCS(I2,J1,K2M1,3),
     >                                   ARCS(I1,J2,K2M1,3),
     >                                   ARCS(I2,J2,K2M1,3))))
      K2POWER = MAX ((K2POWER + 1)/ 2, 3)

      RI2MI1 = ONE / REAL (I2 - I1)
      RJ2MJ1 = ONE / REAL (J2 - J1)
      RK2MK1 = ONE / REAL (K2 - K1)

      IF (FRCSI1) THEN
         DO I = I1P1, I2M1
            DECAYI1(I) = (RI2MI1 * REAL (I2 - I)) ** I1POWER
     >             * EXP (-EXPI1 * REAL (I - I1))
         END DO
      ELSE
         DECAYI1 = ZERO
      END IF

      IF (FRCSI2) THEN
         DO I = I1P1, I2M1
            DECAYI2(I) = (RI2MI1 * REAL (I - I1)) ** I2POWER
     >             * EXP (-EXPI2 * REAL (I2 - I))
         END DO
      ELSE
         DECAYI2 = ZERO
      END IF

      IF (FRCPJ1) THEN
         DO J = J1P1, J2M1
            DECAYJ1(J) = (RJ2MJ1 * REAL (J2 - J)) ** J1POWER
     >             * EXP (-EXPJ1 * REAL (J - J1))
         END DO
      ELSE
         DECAYJ1 = ZERO
      END IF

      IF (FRCPJ2) THEN
         DO J = J1P1, J2M1
            DECAYJ2(J) = (RJ2MJ1 * REAL (J - J1)) ** J2POWER
     >             * EXP (-EXPJ2 * REAL (J2 - J))
         END DO
      ELSE
         DECAYJ2 = ZERO
      END IF

      IF (FRCOK1) THEN
         DO K = K1P1, K2M1
            DECAYK1(K) = (RK2MK1 * REAL (K2 - K)) ** K1POWER
     >             * EXP (-EXPK1 * REAL (K - K1))
         END DO
      ELSE
         DECAYK1 = ZERO
      END IF

      IF (FRCOK2) THEN
         DO K = K1P1, K2M1
            DECAYK2(K) = (RK2MK1 * REAL (K - K1)) ** K2POWER
     >             * EXP (-EXPK2 * REAL (K2 - K))
         END DO
      ELSE
         DECAYK2 = ZERO
      END IF

C     Initialize face P, Q, R values (I1:I2,J1:J2,1:2), etc.:

      PZ = ZERO
      QZ = ZERO
      RZ = ZERO

      PX = ZERO
      QX = ZERO
      RX = ZERO

      PE = ZERO
      QE = ZERO
      RE = ZERO

C***********************************************************************
C     Begin SLOR iterations.
C***********************************************************************

      ROMEGA = ONE / OMG3D

      IF (PRINT) THEN
         I = I2M1 - I1P1 + 3
         J = J2M1 - J1P1 + 3
         K = K2M1 - K1P1 + 3
         WRITE (6,1500) I, J, K, ITMAX, CONVG
      END IF

      DO IT = 1, ITMAX ! Or until convergence

C        Set up the P, Q, and R terms for each face. These routines
C        update only the second derivatives off the faces; all other
C        other derivatives are taken from the GAM* subroutines.

         FOREI = FRCSI1 ! As with the GAM* routines, we save some
         FOREJ = FRCPJ1 ! subscripting by treating just 2-D arrays
         FOREK = FRCOK1 ! at the lower level.

         DO L = 1, 2

            IF (FOREI) CALL PQRXSI (L,            RJACI(J1,K1,L),
     >          GAMX11(J1,K1,L), GAMX12(J1,K1,L), GAMX13(J1,K1,L),
     >          GAMX21(J1,K1,L), GAMX22(J1,K1,L), GAMX23(J1,K1,L),
     >          GAMX31(J1,K1,L), GAMX32(J1,K1,L), GAMX33(J1,K1,L),
     >          FX1(J1,K1,L),    FX2(J1,K1,L),    FX3(J1,K1,L),
     >          F1XSI(J1,K1,L),  F2XSI(J1,K1,L),  F3XSI(J1,K1,L),
     >          PX(J1,K1,L),     QX(J1,K1,L),     RX(J1,K1,L),
     >          CURVI(J1,K1,L))

            IF (FOREJ) CALL PQRETA (L,            RJACJ(I1,K1,L),
     >          GAME11(I1,K1,L), GAME12(I1,K1,L), GAME13(I1,K1,L),
     >          GAME21(I1,K1,L), GAME22(I1,K1,L), GAME23(I1,K1,L),
     >          GAME31(I1,K1,L), GAME32(I1,K1,L), GAME33(I1,K1,L),
     >          FE1(I1,K1,L),    FE2(I1,K1,L),    FE3(I1,K1,L),
     >          F1ETA(I1,K1,L),  F2ETA(I1,K1,L),  F3ETA(I1,K1,L),
     >          PE(I1,K1,L),     QE(I1,K1,L),     RE(I1,K1,L),
     >          CURVJ(I1,K1,L))

            IF (FOREK) CALL PQRZETA (L,           RJACK(I1,J1,L),
     >          GAMZ11(I1,J1,L), GAMZ12(I1,J1,L), GAMZ13(I1,J1,L),
     >          GAMZ21(I1,J1,L), GAMZ22(I1,J1,L), GAMZ23(I1,J1,L),
     >          GAMZ31(I1,J1,L), GAMZ32(I1,J1,L), GAMZ33(I1,J1,L),
     >          FZ1(I1,J1,L),    FZ2(I1,J1,L),    FZ3(I1,J1,L),
     >          F1ZETA(I1,J1,L), F2ZETA(I1,J1,L), F3ZETA(I1,J1,L),
     >          PZ(I1,J1,L),     QZ(I1,J1,L),     RZ(I1,J1,L),
     >          CURVK(I1,J1,L))

            FOREI = FRCSI2
            FOREJ = FRCPJ2
            FOREK = FRCOK2

         END DO

         DXMAX = ZERO
         DYMAX = ZERO
         DZMAX = ZERO

         DO 400 K = K1P1, K2M1

            DECAYK = MAX (DECAYK1(K), DECAYK2(K))
            RK = REAL (K-K1) * RK2MK1
            OMRK = ONE - RK

            DO 300 J = J1P1, J2M1

               DECAYJ = MAX (DECAYJ1(J), DECAYJ2(J))
               RJ = REAL (J-J1) * RJ2MJ1
               OMRJ = ONE - RJ

               DO 200 I = I1P1, I2M1

                  DECAYI = MAX (DECAYI1(I), DECAYI2(I))
                  DECAY  = ONE - MAX (DECAYI, DECAYJ, DECAYK)

C                 Set up actual derivatives in current iteration:

                  X1  = (X(I+1,J,K,1) - X(I-1,J,K,1)) * HALF
                  X2  = (X(I,J+1,K,1) - X(I,J-1,K,1)) * HALF
                  X3  = (X(I,J,K+1,1) - X(I,J,K-1,1)) * HALF
                  X11 = X(I+1,J,K,1) - TWO * X(I,J,K,1) + X(I-1,J,K,1)
                  X22 = X(I,J+1,K,1) - TWO * X(I,J,K,1) + X(I,J-1,K,1)
                  X33 = X(I,J,K+1,1) - TWO * X(I,J,K,1) + X(I,J,K-1,1)
                  X12 = (X(I+1,J+1,K,1) - X(I+1,J-1,K,1) -
     >                   X(I-1,J+1,K,1) + X(I-1,J-1,K,1)) * HALF ! 0.5, not 0.25
                  X23 = (X(I,J+1,K+1,1) - X(I,J+1,K-1,1) -       ! avoids factor
     >                   X(I,J-1,K+1,1) + X(I,J-1,K-1,1)) * HALF ! of 2. below
                  X13 = (X(I+1,J,K+1,1) - X(I+1,J,K-1,1) -
     >                   X(I-1,J,K+1,1) + X(I-1,J,K-1,1)) * HALF
                  Y1  = (X(I+1,J,K,2) - X(I-1,J,K,2)) * HALF
                  Y2  = (X(I,J+1,K,2) - X(I,J-1,K,2)) * HALF
                  Y3  = (X(I,J,K+1,2) - X(I,J,K-1,2)) * HALF
                  Y11 = X(I+1,J,K,2) - TWO * X(I,J,K,2) + X(I-1,J,K,2)
                  Y22 = X(I,J+1,K,2) - TWO * X(I,J,K,2) + X(I,J-1,K,2)
                  Y33 = X(I,J,K+1,2) - TWO * X(I,J,K,2) + X(I,J,K-1,2)
                  Y12 = (X(I+1,J+1,K,2) - X(I+1,J-1,K,2) -
     >                   X(I-1,J+1,K,2) + X(I-1,J-1,K,2)) * HALF
                  Y23 = (X(I,J+1,K+1,2) - X(I,J+1,K-1,2) -
     >                   X(I,J-1,K+1,2) + X(I,J-1,K-1,2)) * HALF
                  Y13 = (X(I+1,J,K+1,2) - X(I+1,J,K-1,2) -
     >                   X(I-1,J,K+1,2) + X(I-1,J,K-1,2)) * HALF
                  Z1  = (X(I+1,J,K,3) - X(I-1,J,K,3)) * HALF
                  Z2  = (X(I,J+1,K,3) - X(I,J-1,K,3)) * HALF
                  Z3  = (X(I,J,K+1,3) - X(I,J,K-1,3)) * HALF
                  Z11 = X(I+1,J,K,3) - TWO * X(I,J,K,3) + X(I-1,J,K,3)
                  Z22 = X(I,J+1,K,3) - TWO * X(I,J,K,3) + X(I,J-1,K,3)
                  Z33 = X(I,J,K+1,3) - TWO * X(I,J,K,3) + X(I,J,K-1,3)
                  Z12 = (X(I+1,J+1,K,3) - X(I+1,J-1,K,3) -
     >                   X(I-1,J+1,K,3) + X(I-1,J-1,K,3)) * HALF
                  Z23 = (X(I,J+1,K+1,3) - X(I,J+1,K-1,3) -
     >                   X(I,J-1,K+1,3) + X(I,J-1,K-1,3)) * HALF
                  Z13 = (X(I+1,J,K+1,3) - X(I+1,J,K-1,3) -
     >                   X(I-1,J,K+1,3) + X(I-1,J,K-1,3)) * HALF

                  GAM11 = Y2 * Z3 - Y3 * Z2
                  GAM21 = X3 * Z2 - X2 * Z3
                  GAM31 = X2 * Y3 - X3 * Y2
                  GAM12 = Y3 * Z1 - Z3 * Y1
                  GAM22 = X1 * Z3 - Z1 * X3
                  GAM32 = X3 * Y1 - Y3 * X1
                  GAM13 = Y1 * Z2 - Z1 * Y2
                  GAM23 = X2 * Z1 - Z2 * X1
                  GAM33 = X1 * Y2 - Y1 * X2
                  ALP11 = GAM11 * GAM11 + GAM21 * GAM21 + GAM31 * GAM31
                  ALP22 = GAM12 * GAM12 + GAM22 * GAM22 + GAM32 * GAM32
                  ALP33 = GAM13 * GAM13 + GAM23 * GAM23 + GAM33 * GAM33
                  ALP12 = GAM11 * GAM12 + GAM21 * GAM22 + GAM31 * GAM32
                  ALP13 = GAM11 * GAM13 + GAM21 * GAM23 + GAM31 * GAM33
                  ALP23 = GAM12 * GAM13 + GAM22 * GAM23 + GAM32 * GAM33
                  AJ2   = (X1 * GAM11 + X2 * GAM12 + X3 * GAM13) ** 2

                  RI   = REAL (I-I1) * RI2MI1
                  OMRI = ONE - RI
                  PHI1 = (OMRJ*PHIJ(I,K,1) + RJ*PHIJ(I,K,2) +
     >                    OMRK*PHIK(I,J,1) + RK*PHIK(I,J,2) -
     >                    OMRK*(OMRJ*PHIJ(I,K1,1) + RJ*PHIJ(I,K1,2)) -
     >                      RK*(OMRJ*PHIJ(I,K2,1) + RJ*PHIJ(I,K2,2))) *
     >                   ALP11

                  DENOM  = ONE / (DECAYI + DECAYJ + DECAYK + EPS)

                  PHIFOR = ((PX(J,K,1)*DECAYI1(I) +
     >                       PX(J,K,2)*DECAYI2(I))*DECAYI +
     >                      (PE(I,K,1)*DECAYJ1(J) +
     >                       PE(I,K,2)*DECAYJ2(J))*DECAYJ +
     >                      (PZ(I,J,1)*DECAYK1(K) +
     >                       PZ(I,J,2)*DECAYK2(K))*DECAYK) * DENOM

                  PHI    = DECAY*PHI1 + PHIFOR*AJ2 ! Blend fore- and background
                  SGNPHI = SIGN (ONE, PHI)

                  PSI1 = (OMRI*PSII(J,K,1) + RI*PSII(J,K,2) +
     >                    OMRK*PSIK(I,J,1) + RK*PSIK(I,J,2) -
     >                    OMRK*(OMRI*PSII(J,K1,1) + RI*PSII(J,K1,2)) -
     >                      RK*(OMRI*PSII(J,K2,1) + RI*PSII(J,K2,2))) *
     >                   ALP22
                  PSIFOR = ((QX(J,K,1)*DECAYI1(I) +
     >                       QX(J,K,2)*DECAYI2(I))*DECAYI +
     >                      (QE(I,K,1)*DECAYJ1(J) +
     >                       QE(I,K,2)*DECAYJ2(J))*DECAYJ +
     >                      (QZ(I,J,1)*DECAYK1(K) +
     >                       QZ(I,J,2)*DECAYK2(K))*DECAYK) * DENOM

                  PSI    = DECAY*PSI1 + PSIFOR*AJ2 ! Blend ...
                  SGNPSI = SIGN (ONE, PSI)

                  OMG1 = (OMRI*OMGI(J,K,1) + RI*OMGI(J,K,2) +
     >                    OMRJ*OMGJ(I,K,1) + RJ*OMGJ(I,K,2) -
     >                    OMRJ*(OMRI*OMGI(J1,K,1) + RI*OMGI(J1,K,2)) -
     >                      RJ*(OMRI*OMGI(J2,K,1) + RI*OMGI(J2,K,2))) *
     >                   ALP33
                  OMGFOR = ((RX(J,K,1)*DECAYI1(I) +
     >                       RX(J,K,2)*DECAYI2(I))*DECAYI +
     >                      (RE(I,K,1)*DECAYJ1(J) +
     >                       RE(I,K,2)*DECAYJ2(J))*DECAYJ +
     >                      (RZ(I,J,1)*DECAYK1(K) +
     >                       RZ(I,J,2)*DECAYK2(K))*DECAYK) * DENOM
                  OMG    = DECAY*OMG1 + OMGFOR*AJ2
                  SGNOMG = SIGN (ONE, OMG)

                  A(I) = (ALP11 + HALF*PHI*(SGNPHI-ONE))*ROMEGA
                  C(I) = (ALP11 + HALF*PHI*(SGNPHI+ONE))*ROMEGA
                  B(I) = (-TWO * (ALP11 + ALP22 + ALP33) -
     >                    SGNPHI*PHI - SGNPSI*PSI - SGNOMG*OMG)*ROMEGA

                  RHS1(I) =-(ALP11*X11 + ALP22*X22 + ALP33*X33) -
     >                      (ALP12*X12 + ALP13*X13 + ALP23*X23) - HALF *
     >                      (PHI*((SGNPHI+ONE)*X(I+1,J,K,1) -TWO*SGNPHI*
     >                       X(I,J,K,1) + (SGNPHI-ONE)*X(I-1,J,K,1)) +
     >                       PSI*((SGNPSI+ONE)*X(I,J+1,K,1) -TWO*SGNPSI*
     >                       X(I,J,K,1) + (SGNPSI-ONE)*X(I,J-1,K,1)) +
     >                       OMG*((SGNOMG+ONE)*X(I,J,K+1,1) -TWO*SGNOMG*
     >                       X(I,J,K,1) + (SGNOMG-ONE)*X(I,J,K-1,1)))
                  RHS2(I) =-(ALP11*Y11 + ALP22*Y22 + ALP33*Y33) -
     >                      (ALP12*Y12 + ALP13*Y13 + ALP23*Y23) - HALF *
     >                      (PHI*((SGNPHI+ONE)*X(I+1,J,K,2) -TWO*SGNPHI*
     >                       X(I,J,K,2) + (SGNPHI-ONE)*X(I-1,J,K,2)) +
     >                       PSI*((SGNPSI+ONE)*X(I,J+1,K,2) -TWO*SGNPSI*
     >                       X(I,J,K,2) + (SGNPSI-ONE)*X(I,J-1,K,2)) +
     >                       OMG*((SGNOMG+ONE)*X(I,J,K+1,2) -TWO*SGNOMG*
     >                       X(I,J,K,2) + (SGNOMG-ONE)*X(I,J,K-1,2)))
                  RHS3(I) =-(ALP11*Z11 + ALP22*Z22 + ALP33*Z33) -
     >                      (ALP12*Z12 + ALP13*Z13 + ALP23*Z23) - HALF *
     >                      (PHI*((SGNPHI+ONE)*X(I+1,J,K,3) -TWO*SGNPHI*
     >                       X(I,J,K,3) + (SGNPHI-ONE)*X(I-1,J,K,3)) +
     >                       PSI*((SGNPSI+ONE)*X(I,J+1,K,3) -TWO*SGNPSI*
     >                       X(I,J,K,3) + (SGNPSI-ONE)*X(I,J-1,K,3)) +
     >                       OMG*((SGNOMG+ONE)*X(I,J,K+1,3) -TWO*SGNOMG*
     >                       X(I,J,K,3) + (SGNOMG-ONE)*X(I,J,K-1,3)))
  200          CONTINUE

               IF (CO) THEN
                  IF (K >= KTIP) THEN ! Singular line stays fixed
                     A(ILE) = ZERO
                     C(ILE) = ZERO
                     RHS1(ILE) = ZERO
                     RHS2(ILE) = ZERO
                     RHS3(ILE) = ZERO
                  END IF
               END IF


C              Solve three tridiagonal systems with the same LHS.
C              Note that the tridiagonal solver is in line for efficiency.

               W = ONE / B(I1P1)
               RHS1(I1P1) = W * RHS1(I1P1)
               RHS2(I1P1) = W * RHS2(I1P1)
               RHS3(I1P1) = W * RHS3(I1P1)

               DO I = I1P1+1, I2M1
                  B(I-1) = C(I-1) * W
                  W =  ONE / (B(I) - A(I) * B(I-1))
                  RHS1(I) = W * (RHS1(I) - A(I) * RHS1(I-1))
                  RHS2(I) = W * (RHS2(I) - A(I) * RHS2(I-1))
                  RHS3(I) = W * (RHS3(I) - A(I) * RHS3(I-1))
               END DO

               DO I = I2M1-1, I1P1, -1
                  RHS1(I) = RHS1(I) - B(I) * RHS1(I+1)
                  RHS2(I) = RHS2(I) - B(I) * RHS2(I+1)
                  RHS3(I) = RHS3(I) - B(I) * RHS3(I+1)
               END DO

               DO I = I1P1, I2M1 ! Doesn't vectorize

                  IF (DXMAX < ABS (RHS1(I))) THEN
                      DXMAX = ABS (RHS1(I))
                      IX = I
                      JX = J
                      KX = K
                  END IF

                  IF (DYMAX < ABS (RHS2(I))) THEN
                      DYMAX = ABS (RHS2(I))
                      IY = I
                      JY = J
                      KY = K
                  END IF

                  IF (DZMAX < ABS (RHS3(I))) THEN
                      DZMAX = ABS (RHS3(I))
                      IZ = I
                      JZ = J
                      KZ = K
                  END IF

               END DO

               DO I = I1P1, I2M1
                  X(I,J,K,1) = X(I,J,K,1) + RHS1(I)
                  X(I,J,K,2) = X(I,J,K,2) + RHS2(I)
                  X(I,J,K,3) = X(I,J,K,3) + RHS3(I)
               END DO

  300       CONTINUE ! Next J

  400    CONTINUE ! Next K

         IF (PRINT) WRITE (6,1510)
     >      IT, DXMAX, IX, JX, KX, DYMAX, IY, JY, KY, DZMAX, IZ, JZ, KZ

         IF (IT <= 1) THEN
            DXM1 = DXMAX
            DYM1 = DYMAX
            DZM1 = DZMAX
         END IF

         IF (LOG10 (MAX (DXM1/DXMAX, DYM1/DYMAX, DZM1/DZMAX)) > CONVG
     >      .OR.    MAX (DXMAX, DYMAX, DZMAX) < DMAX) EXIT

      END DO ! Next SLOR iteration

  999 RETURN

C     Formats:

 1500 FORMAT (/, ' Mesh size:', I5, '  x', I4, '  x', I4,
     >   '   Smoothing iterations:', I4,
     >   '   Orders of magnitude:', F5.2, //, ' ITER',
     >   15X, 'DXMAX    I    J    K',
     >   15X, 'DYMAX    I    J    K',
     >   15X, 'DZMAX    I    J    K', /)
 1510 FORMAT (I4, 1X, 3(1P, E20.7, 3I5))
 1600 FORMAT (2I4, 1P, E14.6)

C     ELLIP3D internal procedures:

      CONTAINS

!***********************************************************************
!
      SUBROUTINE GAMXSI (IFACE, RJACI, GAMX11, GAMX12, GAMX13, GAMX21,
     >                   GAMX22, GAMX23, GAMX31, GAMX32, GAMX33,
     >                   FX1, FX2, FX3, F1XSI, F2XSI, F3XSI, PSII, OMGI,
     >                   CURVI, SXSI, PSION, OMGON)

!     GAMXSI calculates the gamma values, etc., for an I face.
!     It also returns, at each interior face point, the larger of the
!     two curvatures in the eta and zeta directions.
!
!***********************************************************************

!     Arguments:

      INTEGER
     >   IFACE
      REAL, DIMENSION (J1:J2,K1:K2) ::
     >   RJACI,  GAMX11, GAMX12, GAMX13, GAMX21, GAMX22, GAMX23,
     >   GAMX31, GAMX32, GAMX33, FX1,    FX2,    FX3,
     >   F1XSI,  F2XSI,  F3XSI,  PSII,   OMGI,   CURVI,  SXSI
      REAL
     >   PSION, OMGON

!     Execution:

      IF (IFACE == 1) THEN
         I  = I1
         IP = I1P1
         IM = I
      ELSE ! IFACE = 2
         I  = I2
         IP = I
         IM = I2M1
      END IF

!     Off-edge increments:

      DO K = K1, K2, K2 - K1
         DO J = J1, J2
            SXSI(J,K) = SQRT ((X(IP,J,K,1) - X(IM,J,K,1))**2 +
     >                        (X(IP,J,K,2) - X(IM,J,K,2))**2 +
     >                        (X(IP,J,K,3) - X(IM,J,K,3))**2)
         END DO
      END DO

      DO J = J1, J2, J2 - J1
         DO K = K1, K2
            SXSI(J,K) = SQRT ((X(IP,J,K,1) - X(IM,J,K,1))**2 +
     >                        (X(IP,J,K,2) - X(IM,J,K,2))**2 +
     >                        (X(IP,J,K,3) - X(IM,J,K,3))**2)
         END DO
      END DO

!     Set up the first derivative off the I plane on the J1 and J2 edges.
!     Also set up the background Omega:

      DO K = K1P1, K2M1
         XZ = (X(I,J1,K+1,1) - X(I,J1,K-1,1)) * HALF
         YZ = (X(I,J1,K+1,2) - X(I,J1,K-1,2)) * HALF
         ZZ = (X(I,J1,K+1,3) - X(I,J1,K-1,3)) * HALF

         XZZ = X(I,J1,K+1,1) - TWO * X(I,J1,K,1) + X(I,J1,K-1,1)
         YZZ = X(I,J1,K+1,2) - TWO * X(I,J1,K,2) + X(I,J1,K-1,2)
         ZZZ = X(I,J1,K+1,3) - TWO * X(I,J1,K,3) + X(I,J1,K-1,3)

!        Since face values do not get smoothed inside ELLIP3D,
!        at edges a forward or backward difference may be used to
!        calculate the first derivative in the xsi direction.
!        Note these are first-order derivatives.

         F1XSI(J1,K) = X(IP,J1,K,1) - X(IM,J1,K,1)
         F2XSI(J1,K) = X(IP,J1,K,2) - X(IM,J1,K,2)
         F3XSI(J1,K) = X(IP,J1,K,3) - X(IM,J1,K,3)

         OMG = -(XZ*XZZ + YZ*YZZ + ZZ*ZZZ) /
     >          (XZ*XZ  + YZ*YZ  + ZZ*ZZ + EPS)
         OMGI(J1,K) = MIN (MAX (OMG, -OMGON), OMGON)

         XZ = (X(I,J2,K+1,1) - X(I,J2,K-1,1)) * HALF
         YZ = (X(I,J2,K+1,2) - X(I,J2,K-1,2)) * HALF
         ZZ = (X(I,J2,K+1,3) - X(I,J2,K-1,3)) * HALF

         XZZ = X(I,J2,K+1,1) - TWO * X(I,J2,K,1) + X(I,J2,K-1,1)
         YZZ = X(I,J2,K+1,2) - TWO * X(I,J2,K,2) + X(I,J2,K-1,2)
         ZZZ = X(I,J2,K+1,3) - TWO * X(I,J2,K,3) + X(I,J2,K-1,3)

         F1XSI(J2,K) = X(IP,J2,K,1) - X(IM,J2,K,1)
         F2XSI(J2,K) = X(IP,J2,K,2) - X(IM,J2,K,2)
         F3XSI(J2,K) = X(IP,J2,K,3) - X(IM,J2,K,3)

         OMG = -(XZ*XZZ + YZ*YZZ + ZZ*ZZZ) /
     >          (XZ*XZ  + YZ*YZ  + ZZ*ZZ + EPS)
         OMGI(J2,K) = MIN (MAX (OMG, -OMGON), OMGON)
      END DO

!     Set up derivatives on the K1 and K2 edges of the face.
!     Also set up the background PSI:

      DO J = J1P1, J2M1
         XE = (X(I,J+1,K1,1) - X(I,J-1,K1,1)) * HALF
         YE = (X(I,J+1,K1,2) - X(I,J-1,K1,2)) * HALF
         ZE = (X(I,J+1,K1,3) - X(I,J-1,K1,3)) * HALF

         XEE = X(I,J+1,K1,1) - TWO * X(I,J,K1,1) + X(I,J-1,K1,1)
         YEE = X(I,J+1,K1,2) - TWO * X(I,J,K1,2) + X(I,J-1,K1,2)
         ZEE = X(I,J+1,K1,3) - TWO * X(I,J,K1,3) + X(I,J-1,K1,3)

         F1XSI(J,K1) = X(IP,J,K1,1) - X(IM,J,K1,1)
         F2XSI(J,K1) = X(IP,J,K1,2) - X(IM,J,K1,2)
         F3XSI(J,K1) = X(IP,J,K1,3) - X(IM,J,K1,3)

         PSI = -(XE*XEE + YE*YEE + ZE*ZEE) /
     >          (XE*XE  + YE*YE  + ZE*ZE + EPS)
         PSII(J,K1) = MIN (MAX (PSI, -PSION), PSION)

         XE = (X(I,J+1,K2,1) - X(I,J-1,K2,1)) * HALF
         YE = (X(I,J+1,K2,2) - X(I,J-1,K2,2)) * HALF
         ZE = (X(I,J+1,K2,3) - X(I,J-1,K2,3)) * HALF

         XEE = X(I,J+1,K2,1) - TWO * X(I,J,K2,1) + X(I,J-1,K2,1)
         YEE = X(I,J+1,K2,2) - TWO * X(I,J,K2,2) + X(I,J-1,K2,2)
         ZEE = X(I,J+1,K2,3) - TWO * X(I,J,K2,3) + X(I,J-1,K2,3)

         F1XSI(J,K2) = X(IP,J,K2,1) - X(IM,J,K2,1)
         F2XSI(J,K2) = X(IP,J,K2,2) - X(IM,J,K2,2)
         F3XSI(J,K2) = X(IP,J,K2,3) - X(IM,J,K2,3)

         PSI = -(XE*XEE + YE*YEE + ZE*ZEE) /
     >          (XE*XE  + YE*YE  + ZE*ZE + EPS)
         PSII(J,K2) = MIN (MAX (PSI, -PSION), PSION)
      END DO

!     Set up gammas and first derivatives at all interior face points:

      DO K = K1P1, K2M1
         DO J = J1P1, J2M1
            XE = (X(I,J+1,K,1) - X(I,J-1,K,1)) * HALF
            YE = (X(I,J+1,K,2) - X(I,J-1,K,2)) * HALF
            ZE = (X(I,J+1,K,3) - X(I,J-1,K,3)) * HALF

            XEE = X(I,J+1,K,1) - TWO * X(I,J,K,1) + X(I,J-1,K,1)
            YEE = X(I,J+1,K,2) - TWO * X(I,J,K,2) + X(I,J-1,K,2)
            ZEE = X(I,J+1,K,3) - TWO * X(I,J,K,3) + X(I,J-1,K,3)

            XZ = (X(I,J,K+1,1) - X(I,J,K-1,1)) * HALF
            YZ = (X(I,J,K+1,2) - X(I,J,K-1,2)) * HALF
            ZZ = (X(I,J,K+1,3) - X(I,J,K-1,3)) * HALF

            XZZ = X(I,J,K+1,1) - TWO * X(I,J,K,1) + X(I,J,K-1,1)
            YZZ = X(I,J,K+1,2) - TWO * X(I,J,K,2) + X(I,J,K-1,2)
            ZZZ = X(I,J,K+1,3) - TWO * X(I,J,K,3) + X(I,J,K-1,3)

            GAMX11(J,K) = YE*ZZ - YZ*ZE
            GAMX21(J,K) = XZ*ZE - XE*ZZ
            GAMX31(J,K) = XE*YZ - XZ*YE

!           Transfinite interpolation of edge increments into the face
!           interior, using interior arc lengths, not just edge arc lengths:

            P   = ARCS(I,J,K,2)
            PM1 = ONE - ARCS(I,J,K,2)
            Q   = ARCS(I,J,K,3)
            QM1 = ONE - ARCS(I,J,K,3)

            SXSIJK = SXSI(J2,K) * P + SXSI(J1,K) * PM1     +
     >               SXSI(J,K2) * Q + SXSI(J,K1) * QM1     -
     >               P  *(SXSI(J2,K2)*Q + SXSI(J2,K1)*QM1) -
     >               PM1*(SXSI(J1,K2)*Q + SXSI(J1,K1)*QM1)

            SXSIJK = SXSIJK / (SQRT (GAMX11(J,K)*GAMX11(J,K) +
     >                               GAMX21(J,K)*GAMX21(J,K) +
     >                               GAMX31(J,K)*GAMX31(J,K) + EPS))
            XX = SXSIJK*GAMX11(J,K)
            YX = SXSIJK*GAMX21(J,K)
            ZX = SXSIJK*GAMX31(J,K)

            F1XSI(J,K) = XX
            F2XSI(J,K) = YX
            F3XSI(J,K) = ZX

            GAMX12(J,K) = YZ*ZX - YX*ZZ
            GAMX13(J,K) = YX*ZE - YE*ZX
            GAMX22(J,K) = XX*ZZ - XZ*ZX
            GAMX23(J,K) = XE*ZX - XX*ZE
            GAMX32(J,K) = XZ*YX - XX*YZ
            GAMX33(J,K) = XX*YE - XE*YX

            RJACI(J,K) = ONE /
     >         (XX*GAMX11(J,K) +XE*GAMX12(J,K) +XZ*GAMX13(J,K) + EPS)
         END DO
      END DO

!     Set up cross derivatives and save foreground terms with all but
!     the second derivatives off the face:

      DO K = K1P1, K2M1
         DO J = J1P1, J2M1
            XE = (X(I,J+1,K,1) - X(I,J-1,K,1)) * HALF
            YE = (X(I,J+1,K,2) - X(I,J-1,K,2)) * HALF
            ZE = (X(I,J+1,K,3) - X(I,J-1,K,3)) * HALF

            XEE = X(I,J+1,K,1) - TWO * X(I,J,K,1) + X(I,J-1,K,1)
            YEE = X(I,J+1,K,2) - TWO * X(I,J,K,2) + X(I,J-1,K,2)
            ZEE = X(I,J+1,K,3) - TWO * X(I,J,K,3) + X(I,J-1,K,3)

            XZ = (X(I,J,K+1,1) - X(I,J,K-1,1)) * HALF
            YZ = (X(I,J,K+1,2) - X(I,J,K-1,2)) * HALF
            ZZ = (X(I,J,K+1,3) - X(I,J,K-1,3)) * HALF

            XZZ = X(I,J,K+1,1) - TWO * X(I,J,K,1) + X(I,J,K-1,1)
            YZZ = X(I,J,K+1,2) - TWO * X(I,J,K,2) + X(I,J,K-1,2)
            ZZZ = X(I,J,K+1,3) - TWO * X(I,J,K,3) + X(I,J,K-1,3)

            XEZ = (X(I,J+1,K+1,1) - X(I,J-1,K+1,1) -
     >             X(I,J+1,K-1,1) + X(I,J-1,K-1,1)) * FOURTH
            YEZ = (X(I,J+1,K+1,2) - X(I,J-1,K+1,2) -
     >             X(I,J+1,K-1,2) + X(I,J-1,K-1,2)) * FOURTH
            ZEZ = (X(I,J+1,K+1,3) - X(I,J-1,K+1,3) -
     >             X(I,J+1,K-1,3) + X(I,J-1,K-1,3)) * FOURTH

            XXZ = (F1XSI(J,K+1) - F1XSI(J,K-1)) * HALF
            YXZ = (F2XSI(J,K+1) - F2XSI(J,K-1)) * HALF
            ZXZ = (F3XSI(J,K+1) - F3XSI(J,K-1)) * HALF

            XXE = (F1XSI(J+1,K) - F1XSI(J-1,K)) * HALF
            YXE = (F2XSI(J+1,K) - F2XSI(J-1,K)) * HALF
            ZXE = (F3XSI(J+1,K) - F3XSI(J-1,K)) * HALF

            ALP22 = GAMX12(J,K)*GAMX12(J,K) + GAMX22(J,K)*GAMX22(J,K) +
     >              GAMX32(J,K)*GAMX32(J,K)
            ALP33 = GAMX13(J,K)*GAMX13(J,K) + GAMX23(J,K)*GAMX23(J,K) +
     >              GAMX33(J,K)*GAMX33(J,K)

            ALP12 = GAMX11(J,K)*GAMX12(J,K) + GAMX21(J,K)*GAMX22(J,K) +
     >              GAMX31(J,K)*GAMX32(J,K)
            ALP13 = GAMX11(J,K)*GAMX13(J,K) + GAMX21(J,K)*GAMX23(J,K) +
     >              GAMX31(J,K)*GAMX33(J,K)
            ALP23 = GAMX12(J,K)*GAMX13(J,K) + GAMX22(J,K)*GAMX23(J,K) +
     >              GAMX32(J,K)*GAMX33(J,K)

            RAJ2  = -RJACI(J,K)**2

            FX1(J,K) = (ALP22*XEE + ALP33*XZZ + TWO * (ALP12*XXE +
     >                  ALP13*XXZ + ALP23*XEZ)) * RAJ2
            FX2(J,K) = (ALP22*YEE + ALP33*YZZ + TWO * (ALP12*YXE +
     >                  ALP13*YXZ + ALP23*YEZ)) * RAJ2
            FX3(J,K) = (ALP22*ZEE + ALP33*ZZZ + TWO * (ALP12*ZXE +
     >                  ALP13*ZXZ + ALP23*ZEZ)) * RAJ2

            DENOM = ONE / (XE**2 + YE**2 + ZE**2 + EPS)
            PSI = -(XE*XEE + YE*YEE + ZE*ZEE) * DENOM
            PSII(J,K) = MIN (MAX (PSI, -PSION), PSION)
            CURVE1 = ((XE*YEE - XEE*YE)**2 +
     >                (YE*ZEE - YEE*ZE)**2 +
     >                (ZE*XEE - ZEE*XE)**2) * DENOM

            DENOM = ONE / (XZ**2 + YZ**2 + ZZ**2 + EPS)
            OMG = -(XZ*XZZ + YZ*YZZ + ZZ*ZZZ) * DENOM
            OMGI(J,K) = MIN (MAX (OMG, -OMGON), OMGON)
            CURVE2 = ((XZ*YZZ - XZZ*YZ)**2 +
     >                (YZ*ZZZ - YZZ*ZZ)**2 +
     >                (ZZ*XZZ - ZZZ*XZ)**2) * DENOM

            CURVI(J,K) = SQRT (MAX (CURVE1, CURVE2))
         END DO
      END DO

      END SUBROUTINE GAMXSI

!***********************************************************************
!
      SUBROUTINE GAMETA (JFACE, RJACJ, GAME11, GAME12, GAME13, GAME21,
     >                   GAME22, GAME23, GAME31, GAME32, GAME33,
     >                   FE1, FE2, FE3, F1ETA, F2ETA, F3ETA,
     >                   PHIJ, OMGJ, CURVJ, SETA, PHION, OMGON)
!
!     GAMETA calculates the gamma values, etc., for a J face.
!     It also returns, at each interior face point, the larger of the
!     two curvatures in the xsi and zeta directions.
!
!***********************************************************************

!     Arguments:

      INTEGER
     >   JFACE
      REAL, DIMENSION (I1:I2,K1:K2) ::
     >   RJACJ,  GAME11, GAME12, GAME13, GAME21, GAME22, GAME23,
     >   GAME31, GAME32, GAME33, FE1,    FE2,    FE3,
     >   F1ETA,  F2ETA,  F3ETA,  PHIJ,   OMGJ,   CURVJ,  SETA
      REAL
     >   PHION, OMGON

!     Execution:

      IF (JFACE == 1) THEN
         J  = J1
         JP = J1P1
         JM = J
      ELSE ! JFACE = 2
         J  = J2
         JP = J
         JM = J2M1
      END IF

!     Edge increments:

      DO K = K1, K2, K2 - K1
         DO I = I1, I2
            SETA(I,K) = SQRT ((X(I,JP,K,1) - X(I,JM,K,1))**2 +
     >                        (X(I,JP,K,2) - X(I,JM,K,2))**2 +
     >                        (X(I,JP,K,3) - X(I,JM,K,3))**2)
         END DO
      END DO

      DO I = I1, I2, I2 - I1
         DO K = K1, K2
            SETA(I,K) = SQRT ((X(I,JP,K,1) - X(I,JM,K,1))**2 +
     >                        (X(I,JP,K,2) - X(I,JM,K,2))**2 +
     >                        (X(I,JP,K,3) - X(I,JM,K,3))**2)
         END DO
      END DO

!     Set up the first derivative off the J plane on the I1 and I2 edges.
!     Also set up the background Omega:

      DO K = K1P1, K2M1
         XZ = (X(I1,J,K+1,1) - X(I1,J,K-1,1)) * HALF
         YZ = (X(I1,J,K+1,2) - X(I1,J,K-1,2)) * HALF
         ZZ = (X(I1,J,K+1,3) - X(I1,J,K-1,3)) * HALF

         XZZ = X(I1,J,K+1,1) - TWO * X(I1,J,K,1) + X(I1,J,K-1,1)
         YZZ = X(I1,J,K+1,2) - TWO * X(I1,J,K,2) + X(I1,J,K-1,2)
         ZZZ = X(I1,J,K+1,3) - TWO * X(I1,J,K,3) + X(I1,J,K-1,3)

!        Since face values do not get smoothed inside ELLIP3D,
!        at edges a forward or backward difference may be used to
!        calculate the first derivative in the xsi direction.
!        Note these are first-order derivatives.

         F1ETA(I1,K) = X(I1,JP,K,1) - X(I1,JM,K,1)
         F2ETA(I1,K) = X(I1,JP,K,2) - X(I1,JM,K,2)
         F3ETA(I1,K) = X(I1,JP,K,3) - X(I1,JM,K,3)

         OMG = -(XZ*XZZ + YZ*YZZ + ZZ*ZZZ) /
     >          (XZ**2  + YZ**2  + ZZ**2 + EPS)
         OMGJ(I1,K) = MIN (MAX (OMG, -OMGON), OMGON)

         XZ = (X(I2,J,K+1,1) - X(I2,J,K-1,1)) * HALF
         YZ = (X(I2,J,K+1,2) - X(I2,J,K-1,2)) * HALF
         ZZ = (X(I2,J,K+1,3) - X(I2,J,K-1,3)) * HALF

         XZZ = X(I2,J,K+1,1) - TWO * X(I2,J,K,1) + X(I2,J,K-1,1)
         YZZ = X(I2,J,K+1,2) - TWO * X(I2,J,K,2) + X(I2,J,K-1,2)
         ZZZ = X(I2,J,K+1,3) - TWO * X(I2,J,K,3) + X(I2,J,K-1,3)

         F1ETA(I2,K) = X(I2,JP,K,1) - X(I2,JM,K,1)
         F2ETA(I2,K) = X(I2,JP,K,2) - X(I2,JM,K,2)
         F3ETA(I2,K) = X(I2,JP,K,3) - X(I2,JM,K,3)

         OMG = -(XZ*XZZ + YZ*YZZ + ZZ*ZZZ) /
     >          (XZ**2  + YZ**2  + ZZ**2 + EPS)
         OMGJ(I2,K) = MIN (MAX (OMG, -OMGON), OMGON)
      END DO

!     Set up the first derivative off the J plane on the K1 and K2 edges.
!     Also set up the background Phi:

      DO I = I1P1, I2M1
         XX = (X(I+1,J,K1,1) - X(I-1,J,K1,1)) * HALF
         YX = (X(I+1,J,K1,2) - X(I-1,J,K1,2)) * HALF
         ZX = (X(I+1,J,K1,3) - X(I-1,J,K1,3)) * HALF

         XXX = X(I+1,J,K1,1) - TWO * X(I,J,K1,1) + X(I-1,J,K1,1)
         YXX = X(I+1,J,K1,2) - TWO * X(I,J,K1,2) + X(I-1,J,K1,2)
         ZXX = X(I+1,J,K1,3) - TWO * X(I,J,K1,3) + X(I-1,J,K1,3)

         F1ETA(I,K1) = X(I,JP,K1,1) - X(I,JM,K1,1)
         F2ETA(I,K1) = X(I,JP,K1,2) - X(I,JM,K1,2)
         F3ETA(I,K1) = X(I,JP,K1,3) - X(I,JM,K1,3)

         PHI = -(XX*XXX + YX*YXX + ZX*ZXX) /
     >          (XX**2  + YX**2  + ZX**2 + EPS)
         PHIJ(I,K1) = MIN (MAX (PHI, -PHION), PHION)

         XX = (X(I+1,J,K2,1) - X(I-1,J,K2,1)) * HALF
         YX = (X(I+1,J,K2,2) - X(I-1,J,K2,2)) * HALF
         ZX = (X(I+1,J,K2,3) - X(I-1,J,K2,3)) * HALF

         XXX = X(I+1,J,K2,1) - TWO * X(I,J,K2,1) + X(I-1,J,K2,1)
         YXX = X(I+1,J,K2,2) - TWO * X(I,J,K2,2) + X(I-1,J,K2,2)
         ZXX = X(I+1,J,K2,3) - TWO * X(I,J,K2,3) + X(I-1,J,K2,3)

         F1ETA(I,K2) = X(I,JP,K2,1) - X(I,JM,K2,1)
         F2ETA(I,K2) = X(I,JP,K2,2) - X(I,JM,K2,2)
         F3ETA(I,K2) = X(I,JP,K2,3) - X(I,JM,K2,3)

         PHI = -(XX*XXX + YX*YXX + ZX*ZXX) /
     >          (XX**2  + YX**2  + ZX**2 + EPS)
         PHIJ(I,K2) = MIN (MAX (PHI, -PHION), PHION)
      END DO

!     Set up gammas and first derivatives at all interior face points:

      DO K = K1P1, K2M1
         DO I = I1P1, I2M1
            XX = (X(I+1,J,K,1) - X(I-1,J,K,1)) * HALF
            YX = (X(I+1,J,K,2) - X(I-1,J,K,2)) * HALF
            ZX = (X(I+1,J,K,3) - X(I-1,J,K,3)) * HALF

            XZ = (X(I,J,K+1,1) - X(I,J,K-1,1)) * HALF
            YZ = (X(I,J,K+1,2) - X(I,J,K-1,2)) * HALF
            ZZ = (X(I,J,K+1,3) - X(I,J,K-1,3)) * HALF

            GAME12(I,K) = YZ*ZX - YX*ZZ
            GAME22(I,K) = XX*ZZ - XZ*ZX
            GAME32(I,K) = XZ*YX - XX*YZ

!           Transfinite interpolation of edge increments into the face interiors
!           using normalized arc-lengths for the whole face (not just edges):

            P   = ARCS(I,J,K,1)
            PM1 = ONE - ARCS(I,J,K,1)
            Q   = ARCS(I,J,K,3)
            QM1 = ONE - ARCS(I,J,K,3)

            SETAIK = SETA(I2,K) * P + SETA(I1,K) * PM1     +
     >               SETA(I,K2) * Q + SETA(I,K1) * QM1     -
     >               P  *(SETA(I2,K2)*Q + SETA(I2,K1)*QM1) -
     >               PM1*(SETA(I1,K2)*Q + SETA(I1,K1)*QM1)

            SETAIK = SETAIK / SQRT (GAME12(I,K)**2 +
     >                              GAME22(I,K)**2 +
     >                              GAME32(I,K)**2 + EPS)
            XE = SETAIK*GAME12(I,K)
            YE = SETAIK*GAME22(I,K)
            ZE = SETAIK*GAME32(I,K)

            F1ETA(I,K) = XE
            F2ETA(I,K) = YE
            F3ETA(I,K) = ZE

            GAME11(I,K) = YE*ZZ - YZ*ZE
            GAME13(I,K) = YX*ZE - YE*ZX
            GAME21(I,K) = XZ*ZE - XE*ZZ
            GAME23(I,K) = XE*ZX - XX*ZE
            GAME31(I,K) = XE*YZ - XZ*YE
            GAME33(I,K) = XX*YE - XE*YX

            RJACJ(I,K)  = ONE /
     >         (XX*GAME11(I,K) +XE*GAME12(I,K) +XZ*GAME13(I,K) + EPS)
         END DO
      END DO

!     Cross derivatives & foreground terms except for 2nd derivs. off the face:

      DO K = K1P1, K2M1
         DO I = I1P1, I2M1
            XX = (X(I+1,J,K,1) - X(I-1,J,K,1)) * HALF
            YX = (X(I+1,J,K,2) - X(I-1,J,K,2)) * HALF
            ZX = (X(I+1,J,K,3) - X(I-1,J,K,3)) * HALF

            XXX = X(I+1,J,K,1) - TWO * X(I,J,K,1) + X(I-1,J,K,1)
            YXX = X(I+1,J,K,2) - TWO * X(I,J,K,2) + X(I-1,J,K,2)
            ZXX = X(I+1,J,K,3) - TWO * X(I,J,K,3) + X(I-1,J,K,3)

            XZ = (X(I,J,K+1,1) - X(I,J,K-1,1)) * HALF
            YZ = (X(I,J,K+1,2) - X(I,J,K-1,2)) * HALF
            ZZ = (X(I,J,K+1,3) - X(I,J,K-1,3)) * HALF

            XZZ = X(I,J,K+1,1) - TWO * X(I,J,K,1) + X(I,J,K-1,1)
            YZZ = X(I,J,K+1,2) - TWO * X(I,J,K,2) + X(I,J,K-1,2)
            ZZZ = X(I,J,K+1,3) - TWO * X(I,J,K,3) + X(I,J,K-1,3)

            XXZ = (X(I+1,J,K+1,1) - X(I+1,J,K-1,1) -
     >             X(I-1,J,K+1,1) + X(I-1,J,K-1,1)) * FOURTH
            YXZ = (X(I+1,J,K+1,2) - X(I+1,J,K-1,2) -
     >             X(I-1,J,K+1,2) + X(I-1,J,K-1,2)) * FOURTH
            ZXZ = (X(I+1,J,K+1,3) - X(I+1,J,K-1,3) -
     >             X(I-1,J,K+1,3) + X(I-1,J,K-1,3)) * FOURTH

            XEZ = (F1ETA(I,K+1) - F1ETA(I,K-1)) * HALF
            YEZ = (F2ETA(I,K+1) - F2ETA(I,K-1)) * HALF
            ZEZ = (F3ETA(I,K+1) - F3ETA(I,K-1)) * HALF

            XXE = (F1ETA(I+1,K) - F1ETA(I-1,K)) * HALF
            YXE = (F2ETA(I+1,K) - F2ETA(I-1,K)) * HALF
            ZXE = (F3ETA(I+1,K) - F3ETA(I-1,K)) * HALF

            ALP11 = GAME11(I,K)**2 + GAME21(I,K)**2 + GAME31(I,K)**2
            ALP33 = GAME13(I,K)**2 + GAME23(I,K)**2 + GAME33(I,K)**2

            ALP12 = GAME11(I,K)*GAME12(I,K) + GAME21(I,K)*GAME22(I,K) +
     >              GAME31(I,K)*GAME32(I,K)
            ALP13 = GAME11(I,K)*GAME13(I,K) + GAME21(I,K)*GAME23(I,K) +
     >              GAME31(I,K)*GAME33(I,K)
            ALP23 = GAME12(I,K)*GAME13(I,K) + GAME22(I,K)*GAME23(I,K) +
     >              GAME32(I,K)*GAME33(I,K)

            RAJ2  = -RJACJ(I,K)**2

            FE1(I,K) = (ALP11*XXX + ALP33*XZZ + TWO * (ALP12*XXE +
     >                  ALP13*XXZ + ALP23*XEZ)) * RAJ2
            FE2(I,K) = (ALP11*YXX + ALP33*YZZ + TWO * (ALP12*YXE +
     >                  ALP13*YXZ + ALP23*YEZ)) * RAJ2
            FE3(I,K) = (ALP11*ZXX + ALP33*ZZZ + TWO * (ALP12*ZXE +
     >                  ALP13*ZXZ + ALP23*ZEZ)) * RAJ2

            DENOM = ONE / (XX**2 + YX**2 + ZX**2 + EPS)
            PHI = -(XX*XXX + YX*YXX + ZX*ZXX) * DENOM
            PHIJ(I,K) = MIN (MAX (PHI, -PHION), PHION)
            CURVE1 = ((XX*YXX - XXX*YX)**2 +
     >                (YX*ZXX - YXX*ZX)**2 +
     >                (ZX*XXX - ZXX*XX)**2) * DENOM

            DENOM = ONE / (XZ**2 + YZ**2 + ZZ**2 + EPS)
            OMG = -(XZ*XZZ + YZ*YZZ + ZZ*ZZZ) * DENOM
            OMGJ(I,K) = MIN (MAX (OMG, -OMGON), OMGON)
            CURVE2 = ((XZ*YZZ - XZZ*YZ)**2 +
     >                (YZ*ZZZ - YZZ*ZZ)**2 +
     >                (ZZ*XZZ - ZZZ*XZ)**2) * DENOM

            CURVJ(I,K) = SQRT (MAX (CURVE1, CURVE2))
         END DO
      END DO

      END SUBROUTINE GAMETA

!***********************************************************************
!
      SUBROUTINE GAMZETA (KFACE, RJACK, GAMZ11, GAMZ12, GAMZ13, GAMZ21,
     >                    GAMZ22, GAMZ23, GAMZ31, GAMZ32, GAMZ33,
     >                    FZ1, FZ2, FZ3, F1ZETA, F2ZETA, F3ZETA,
     >                    PHIK, PSIK, CURVK, SZETA, PHION, PSION)
!
!     GAMZETA calculates the gamma values, etc., for an I face.
!     It also returns, at each interior face point, the larger of the
!     two curvatures in the xsi and eta directions.
!
!***********************************************************************

!     Arguments:

      INTEGER
     >   KFACE
      REAL, DIMENSION (I1:I2,J1:J2) ::
     >   RJACK,  GAMZ11, GAMZ12, GAMZ13, GAMZ21, GAMZ22, GAMZ23,
     >   GAMZ31, GAMZ32, GAMZ33, FZ1,    FZ2,    FZ3,
     >   F1ZETA, F2ZETA, F3ZETA, PHIK,   PSIK,   CURVK,  SZETA
      REAL
     >   PHION, PSION

!     Execution:

      IF (KFACE == 1) THEN
         K  = K1
         KP = K1P1
         KM = K
      ELSE ! KFACE = 2
         K  = K2
         KP = K
         KM = K2M1
      END IF

!     Off-edge increments:

      DO J = J1, J2, J2 - J1
         DO I = I1, I2
            SZETA(I,J) = SQRT ((X(I,J,KP,1) - X(I,J,KM,1))**2 +
     >                         (X(I,J,KP,2) - X(I,J,KM,2))**2 +
     >                         (X(I,J,KP,3) - X(I,J,KM,3))**2)
         END DO
      END DO

      DO I = I1, I2, I2 - I1
          DO J = J1, J2
            SZETA(I,J) = SQRT ((X(I,J,KP,1) - X(I,J,KM,1))**2 +
     >                         (X(I,J,KP,2) - X(I,J,KM,2))**2 +
     >                         (X(I,J,KP,3) - X(I,J,KM,3))**2)
         END DO
      END DO

!     Set up the first derivative off the K plane on the I1 and I2 edges.
!     Also set up the background Psi:

      DO J = J1P1, J2M1
         XE = (X(I1,J+1,K,1) - X(I1,J-1,K,1)) * HALF
         YE = (X(I1,J+1,K,2) - X(I1,J-1,K,2)) * HALF
         ZE = (X(I1,J+1,K,3) - X(I1,J-1,K,3)) * HALF

         XEE = X(I1,J+1,K,1) - TWO * X(I1,J,K,1) + X(I1,J-1,K,1)
         YEE = X(I1,J+1,K,2) - TWO * X(I1,J,K,2) + X(I1,J-1,K,2)
         ZEE = X(I1,J+1,K,3) - TWO * X(I1,J,K,3) + X(I1,J-1,K,3)

!        Since face value do not get smoothed inside ELLIP3D,
!        at edges a forward or backward difference may be used to
!        calculate the first derivative in the zeta direction.
!        Note these are first-order derivatives.

         F1ZETA(I1,J) = X(I1,J,KP,1) - X(I1,J,KM,1)
         F2ZETA(I1,J) = X(I1,J,KP,2) - X(I1,J,KM,2)
         F3ZETA(I1,J) = X(I1,J,KP,3) - X(I1,J,KM,3)

         PSI = -(XE*XEE + YE*YEE + ZE*ZEE) /
     >          (XE*XE  + YE*YE  + ZE*ZE + EPS)
         PSIK(I1,J) = MIN (MAX (PSI, -PSION), PSION)

         XE = (X(I2,J+1,K,1) - X(I2,J-1,K,1)) * HALF
         YE = (X(I2,J+1,K,2) - X(I2,J-1,K,2)) * HALF
         ZE = (X(I2,J+1,K,3) - X(I2,J-1,K,3)) * HALF

         XEE = X(I2,J+1,K,1) - TWO * X(I2,J,K,1) + X(I2,J-1,K,1)
         YEE = X(I2,J+1,K,2) - TWO * X(I2,J,K,2) + X(I2,J-1,K,2)
         ZEE = X(I2,J+1,K,3) - TWO * X(I2,J,K,3) + X(I2,J-1,K,3)

         F1ZETA(I2,J) = X(I2,J,KP,1) - X(I2,J,KM,1)
         F2ZETA(I2,J) = X(I2,J,KP,2) - X(I2,J,KM,2)
         F3ZETA(I2,J) = X(I2,J,KP,3) - X(I2,J,KM,3)

         PSI = -(XE*XEE + YE*YEE + ZE*ZEE) /
     >          (XE*XE  + YE*YE  + ZE*ZE + EPS)
         PSIK(I2,J) = MIN (MAX (PSI, -PSION), PSION)
      END DO

!     Set up first derivatives on the J1 and J2 edges of the face.
!     Also set up background Phi terms:

      DO I = I1P1, I2M1
         XX = (X(I+1,J1,K,1) - X(I-1,J1,K,1)) * HALF
         YX = (X(I+1,J1,K,2) - X(I-1,J1,K,2)) * HALF
         ZX = (X(I+1,J1,K,3) - X(I-1,J1,K,3)) * HALF

         XXX = X(I+1,J1,K,1) - TWO * X(I,J1,K,1) + X(I-1,J1,K,1)
         YXX = X(I+1,J1,K,2) - TWO * X(I,J1,K,2) + X(I-1,J1,K,2)
         ZXX = X(I+1,J1,K,3) - TWO * X(I,J1,K,3) + X(I-1,J1,K,3)

         F1ZETA(I,J1) = X(I,J1,KP,1) - X(I,J1,KM,1)
         F2ZETA(I,J1) = X(I,J1,KP,2) - X(I,J1,KM,2)
         F3ZETA(I,J1) = X(I,J1,KP,3) - X(I,J1,KM,3)

         PHI = -(XX*XXX + YX*YXX + ZX*ZXX) /
     >          (XX*XX  + YX*YX  + ZX*ZX + EPS)
         PHIK(I,J1) = MIN (MAX (PHI, -PHION), PHION)

         XX = (X(I+1,J2,K,1) - X(I-1,J2,K,1)) * HALF
         YX = (X(I+1,J2,K,2) - X(I-1,J2,K,2)) * HALF
         ZX = (X(I+1,J2,K,3) - X(I-1,J2,K,3)) * HALF

         XXX = X(I+1,J2,K,1) - TWO * X(I,J2,K,1) + X(I-1,J2,K,1)
         YXX = X(I+1,J2,K,2) - TWO * X(I,J2,K,2) + X(I-1,J2,K,2)
         ZXX = X(I+1,J2,K,3) - TWO * X(I,J2,K,3) + X(I-1,J2,K,3)

         F1ZETA(I,J2) = X(I,J2,KP,1) - X(I,J2,KM,1)
         F2ZETA(I,J2) = X(I,J2,KP,2) - X(I,J2,KM,2)
         F3ZETA(I,J2) = X(I,J2,KP,3) - X(I,J2,KM,3)

         PHI = -(XX*XXX + YX*YXX + ZX*ZXX) /
     >          (XX*XX  + YX*YX  + ZX*ZX + EPS)
         PHIK(I,J2) = MIN (MAX (PHI, -PHION), PHION)
      END DO

!     Set up gammas and first derivatives at all interior face points:

      DO J = J1P1, J2M1
         DO I = I1P1, I2M1
            XX = (X(I+1,J,K,1) - X(I-1,J,K,1)) * HALF
            YX = (X(I+1,J,K,2) - X(I-1,J,K,2)) * HALF
            ZX = (X(I+1,J,K,3) - X(I-1,J,K,3)) * HALF

            XE = (X(I,J+1,K,1) - X(I,J-1,K,1)) * HALF
            YE = (X(I,J+1,K,2) - X(I,J-1,K,2)) * HALF
            ZE = (X(I,J+1,K,3) - X(I,J-1,K,3)) * HALF

            GAMZ13(I,J) = YX*ZE - YE*ZX
            GAMZ23(I,J) = XE*ZX - XX*ZE
            GAMZ33(I,J) = XX*YE - XE*YX

!           Transfinite interpolation of edge increments into the face interiors
!           using normalized arc-lengths for the whole face (not just edges):

            P   = ARCS(I,J,K,1)
            PM1 = ONE - ARCS(I,J,K,1)
            Q   = ARCS(I,J,K,2)
            QM1 = ONE - ARCS(I,J,K,2)

            SZETAIJ = SZETA(I2,J) * P + SZETA(I1,J) * PM1      +
     >                SZETA(I,J2) * Q + SZETA(I,J1) * QM1      -
     >                P   *(SZETA(I2,J2)*Q + SZETA(I2,J1)*QM1) -
     >                PM1 *(SZETA(I1,J2)*Q + SZETA(I1,J1)*QM1)

            SZETAIJ = SZETAIJ / SQRT (GAMZ13(I,J)**2 +
     >                                GAMZ23(I,J)**2 +
     >                                GAMZ33(I,J)**2 + EPS)
            XZ = SZETAIJ*GAMZ13(I,J)
            YZ = SZETAIJ*GAMZ23(I,J)
            ZZ = SZETAIJ*GAMZ33(I,J)

            F1ZETA(I,J) = XZ
            F2ZETA(I,J) = YZ
            F3ZETA(I,J) = ZZ

            GAMZ11(I,J) = YE*ZZ - YZ*ZE
            GAMZ21(I,J) = XZ*ZE - XE*ZZ
            GAMZ31(I,J) = XE*YZ - XZ*YE

            GAMZ12(I,J) = YZ*ZX - YX*ZZ
            GAMZ22(I,J) = XX*ZZ - XZ*ZX
            GAMZ32(I,J) = XZ*YX - XX*YZ

            RJACK(I,J) = ONE /
     >         (XX*GAMZ11(I,J) +XE*GAMZ12(I,J) +XZ*GAMZ13(I,J) + EPS)
         END DO
      END DO

!     Set up cross derivatives and save foreground terms with all but
!     the second derivatives off the face:

      DO J = J1P1, J2M1
         DO I = I1P1, I2M1
            XX = (X(I+1,J,K,1) - X(I-1,J,K,1)) * HALF
            YX = (X(I+1,J,K,2) - X(I-1,J,K,2)) * HALF
            ZX = (X(I+1,J,K,3) - X(I-1,J,K,3)) * HALF

            XXX = X(I+1,J,K,1) - TWO * X(I,J,K,1) + X(I-1,J,K,1)
            YXX = X(I+1,J,K,2) - TWO * X(I,J,K,2) + X(I-1,J,K,2)
            ZXX = X(I+1,J,K,3) - TWO * X(I,J,K,3) + X(I-1,J,K,3)

            XE = (X(I,J+1,K,1) - X(I,J-1,K,1)) * HALF
            YE = (X(I,J+1,K,2) - X(I,J-1,K,2)) * HALF
            ZE = (X(I,J+1,K,3) - X(I,J-1,K,3)) * HALF

            XEE = X(I,J+1,K,1) - TWO * X(I,J,K,1) + X(I,J-1,K,1)
            YEE = X(I,J+1,K,2) - TWO * X(I,J,K,2) + X(I,J-1,K,2)
            ZEE = X(I,J+1,K,3) - TWO * X(I,J,K,3) + X(I,J-1,K,3)

            XXE = (X(I+1,J+1,K,1) - X(I+1,J-1,K,1) -
     >             X(I-1,J+1,K,1) + X(I-1,J-1,K,1)) * FOURTH
            YXE = (X(I+1,J+1,K,2) - X(I+1,J-1,K,2) -
     >             X(I-1,J+1,K,2) + X(I-1,J-1,K,2)) * FOURTH
            ZXE = (X(I+1,J+1,K,3) - X(I+1,J-1,K,3) -
     >             X(I-1,J+1,K,3) + X(I-1,J-1,K,3)) * FOURTH

            XEZ = (F1ZETA(I,J+1) - F1ZETA(I,J-1)) * HALF
            YEZ = (F2ZETA(I,J+1) - F2ZETA(I,J-1)) * HALF
            ZEZ = (F3ZETA(I,J+1) - F3ZETA(I,J-1)) * HALF

            XXZ = (F1ZETA(I+1,J) - F1ZETA(I-1,J)) * HALF
            YXZ = (F2ZETA(I+1,J) - F2ZETA(I-1,J)) * HALF
            ZXZ = (F3ZETA(I+1,J) - F3ZETA(I-1,J)) * HALF

            ALP11 = GAMZ11(I,J)**2 + GAMZ21(I,J)**2 + GAMZ31(I,J)**2
            ALP22 = GAMZ12(I,J)**2 + GAMZ22(I,J)**2 + GAMZ32(I,J)**2

            ALP12 = GAMZ11(I,J)*GAMZ12(I,J) + GAMZ21(I,J)*GAMZ22(I,J) +
     >              GAMZ31(I,J)*GAMZ32(I,J)
            ALP13 = GAMZ11(I,J)*GAMZ13(I,J) + GAMZ21(I,J)*GAMZ23(I,J) +
     >              GAMZ31(I,J)*GAMZ33(I,J)
            ALP23 = GAMZ12(I,J)*GAMZ13(I,J) + GAMZ22(I,J)*GAMZ23(I,J) +
     >              GAMZ32(I,J)*GAMZ33(I,J)

            RAJ2  = -RJACK(I,J)**2

            FZ1(I,J) = (ALP11*XXX + ALP22*XEE + TWO * (ALP12*XXE +
     >                  ALP13*XXZ + ALP23*XEZ)) * RAJ2
            FZ2(I,J) = (ALP11*YXX + ALP22*YEE + TWO * (ALP12*YXE +
     >                  ALP13*YXZ + ALP23*YEZ)) * RAJ2
            FZ3(I,J) = (ALP11*ZXX + ALP22*ZEE + TWO * (ALP12*ZXE +
     >                  ALP13*ZXZ + ALP23*ZEZ)) * RAJ2

            DENOM = ONE / (XX**2 + YX**2 + ZX**2 + EPS)
            PHI = -(XX*XXX + YX*YXX + ZX*ZXX) * DENOM
            PHIK(I,J) = MIN (MAX (PHI, -PHION), PHION)
            CURVE1 = ((XX*YXX - XXX*YX)**2 +
     >                (YX*ZXX - YXX*ZX)**2 +
     >                (ZX*XXX - ZXX*XX)**2) * DENOM

            DENOM = ONE / (XE**2 + YE**2 + ZE**2 + EPS)
            PSI = -(XE*XEE + YE*YEE + ZE*ZEE) * DENOM
            PSIK(I,J) = MIN (MAX (PSI, -PSION), PSION)
            CURVE2 = ((XE*YEE - XEE*YE)**2 +
     >                (YE*ZEE - YEE*ZE)**2 +
     >                (ZE*XEE - ZEE*XE)**2) * DENOM

            CURVK(I,J) = MAX (CURVE1, CURVE2)
         END DO
      END DO

      END SUBROUTINE GAMZETA

!***********************************************************************
!
      SUBROUTINE PQRXSI (IFACE, RJACI, GAMX11, GAMX12, GAMX13, GAMX21,
     >                   GAMX22, GAMX23, GAMX31, GAMX32, GAMX33,
     >                   FX1, FX2, FX3, F1XSI, F2XSI, F3XSI,
     >                   PX, QX, RX, CURVI)
!
!     PQRXSI calculates the second derivatives normal to a constant xsi
!     surface and returns Sorenson's P, Q and R terms.
!
!***********************************************************************

!     Arguments:

      INTEGER
     >   IFACE
      REAL, DIMENSION (J1:J2,K1:K2) ::
     >   RJACI,  GAMX11, GAMX12, GAMX13, GAMX21, GAMX22, GAMX23,
     >   GAMX31, GAMX32, GAMX33, FX1,    FX2,    FX3,
     >   F1XSI,  F2XSI,  F3XSI,  PX,     QX,     RX,     CURVI

!     Execution:

      IF (IFACE == 1) THEN
         I  = I1
         IP = I1P1
         IQ = IP + 1
         SIGN1 = 3.
      ELSE ! IFACE = 2
         I  = I2
         IP = I2M1
         IQ = IP - 1
         SIGN1 = -3.
      END IF

!     Set up the second derivatives and the P, Q, R terms:

      DO K = K1 + 1, K2 - 1
         DO J = J1 + 1, J2 - 1
            XXX = -3.5*X(I,J,K,1) + 4.*X(IP,J,K,1) - HALF*X(IQ,J,K,1) -
     >             SIGN1*F1XSI(J,K)
            YXX = -3.5*X(I,J,K,2) + 4.*X(IP,J,K,2) - HALF*X(IQ,J,K,2) -
     >             SIGN1*F2XSI(J,K)
            ZXX = -3.5*X(I,J,K,3) + 4.*X(IP,J,K,3) - HALF*X(IQ,J,K,3) -
     >             SIGN1*F3XSI(J,K)

            RAJ  = RJACI(J,K)
            RAJ2 = RAJ**2

            ALPHA = (GAMX11(J,K)**2 + GAMX21(J,K)**2 + GAMX31(J,K)**2) *
     >               RAJ2

            H1 = FX1(J,K) - ALPHA*XXX
            H2 = FX2(J,K) - ALPHA*YXX
            H3 = FX3(J,K) - ALPHA*ZXX

            TEMPLIM = MIN (ONE, FGLIMIT / (ONE + CURVI(J,K)))

            DIF = (H1*GAMX11(J,K) + H2*GAMX21(J,K) + H3*GAMX31(J,K)) *
     >             RAJ - PX(J,K)
            TEMP1 = TEMPLIM * MAX (ABS (PX(J,K)), ONE) ! Limiting term
            TEMP2 = MIN (FGRELAX * ABS (DIF), TEMP1)   ! Relaxation term
            PX(J,K) = PX(J,K) + SIGN (TEMP2, DIF)

            DIF = (H1*GAMX12(J,K) + H2*GAMX22(J,K) + H3*GAMX32(J,K)) *
     >             RAJ - QX(J,K)
            TEMP1 = TEMPLIM * MAX (ABS (QX(J,K)), ONE)
            TEMP2 = MIN (FGRELAX * ABS (DIF), TEMP1)
            QX(J,K) = QX(J,K) + SIGN (TEMP2, DIF)

            DIF = (H1*GAMX13(J,K) + H2*GAMX23(J,K) + H3*GAMX33(J,K)) *
     >             RAJ - RX(J,K)
            TEMP1 = TEMPLIM * MAX (ABS (RX(J,K)), ONE)
            TEMP2 = MIN (FGRELAX * ABS (DIF), TEMP1)
            RX(J,K) = RX(J,K) + SIGN (TEMP2, DIF)
         END DO
      END DO

      END SUBROUTINE PQRXSI

!***********************************************************************
!
      SUBROUTINE PQRETA (JFACE, RJACJ, GAME11, GAME12, GAME13, GAME21,
     >                   GAME22, GAME23, GAME31, GAME32, GAME33,
     >                   FE1, FE2, FE3, F1ETA, F2ETA, F3ETA,
     >                   PE, QE, RE, CURVJ)
!
!     PQRETA calculates the second derivatives normal to a constant eta
!     surface and returns Sorenson's P, Q and R terms.
!
!***********************************************************************

!     Arguments:

      INTEGER
     >   JFACE
      REAL, DIMENSION (I1:I2,K1:K2) ::
     >   RJACJ,  GAME11, GAME12, GAME13, GAME21, GAME22, GAME23,
     >   GAME31, GAME32, GAME33, FE1,    FE2,    FE3,
     >   F1ETA,  F2ETA,  F3ETA,  PE,     QE,     RE,     CURVJ

!     Execution:

      IF (JFACE == 1) THEN
         J  = J1
         JP = J1P1
         JQ = JP + 1
         SIGN1 = 3.
      ELSE ! JFACE = 2
         J  = J2
         JP = J2M1
         JQ = JP - 1
         SIGN1 = -3.
      END IF

!     Set up the second derivatives and the P, Q, R terms:

      DO K = K1 + 1, K2 - 1
         DO I = I1 + 1, I2 - 1
            XEE = -3.5*X(I,J,K,1) + 4.*X(I,JP,K,1) - HALF*X(I,JQ,K,1) -
     >             SIGN1*F1ETA(I,K)
            YEE = -3.5*X(I,J,K,2) + 4.*X(I,JP,K,2) - HALF*X(I,JQ,K,2) -
     >             SIGN1*F2ETA(I,K)
            ZEE = -3.5*X(I,J,K,3) + 4.*X(I,JP,K,3) - HALF*X(I,JQ,K,3) -
     >             SIGN1*F3ETA(I,K)

            RAJ  = RJACJ(I,K)
            RAJ2 = RAJ**2

            ALPHA = (GAME12(I,K)**2 + GAME22(I,K)**2 + GAME32(I,K)**2) *
     >               RAJ2

            H1 = FE1(I,K) - ALPHA*XEE
            H2 = FE2(I,K) - ALPHA*YEE
            H3 = FE3(I,K) - ALPHA*ZEE

            TEMPLIM = MIN (ONE, FGLIMIT / (ONE + CURVJ(I,K)))

            DIF = (H1*GAME11(I,K) + H2*GAME21(I,K) + H3*GAME31(I,K)) *
     >             RAJ - PE(I,K)
            TEMP1 = TEMPLIM * MAX (ABS (PE(I,K)), ONE)
            TEMP2 = MIN (FGRELAX * ABS (DIF), TEMP1)
            PE(I,K) = PE(I,K) + SIGN (TEMP2, DIF)

            DIF = (H1*GAME12(I,K) + H2*GAME22(I,K) + H3*GAME32(I,K)) *
     >             RAJ - QE(I,K)
            TEMP1 = TEMPLIM * MAX (ABS (QE(I,K)), ONE)
            TEMP2 = MIN (FGRELAX * ABS (DIF), TEMP1)
            QE(I,K) = QE(I,K) + SIGN (TEMP2, DIF)

            DIF = (H1*GAME13(I,K) + H2*GAME23(I,K) + H3*GAME33(I,K)) *
     >             RAJ - RE(I,K)
            TEMP1 = TEMPLIM * MAX (ABS (RE(I,K)), ONE)
            TEMP2 = MIN (FGRELAX * ABS (DIF), TEMP1)
            RE(I,K) = RE(I,K) + SIGN (TEMP2, DIF)
         END DO
      END DO

      END SUBROUTINE PQRETA

!***********************************************************************
!
      SUBROUTINE PQRZETA (KFACE, RJACK, GAMZ11, GAMZ12, GAMZ13, GAMZ21,
     >                    GAMZ22, GAMZ23, GAMZ31, GAMZ32, GAMZ33,
     >                    FZ1, FZ2, FZ3, F1ZETA, F2ZETA, F3ZETA,
     >                    PZ, QZ, RZ, CURVK)
!
!     PQRZETA calculates the second derivatives normal to a constant zeta
!     surface and returns Sorenson's P, Q and R terms.
!
!***********************************************************************

!     Arguments:

      INTEGER
     >   KFACE
      REAL, DIMENSION (I1:I2,J1:J2) ::
     >   RJACK,  GAMZ11, GAMZ12, GAMZ13, GAMZ21, GAMZ22, GAMZ23,
     >   GAMZ31, GAMZ32, GAMZ33, FZ1,    FZ2,    FZ3,
     >   F1ZETA, F2ZETA, F3ZETA, PZ,     QZ,     RZ,     CURVK

!     Execution:

      IF (KFACE == 1) THEN
         K  = K1
         KP = K1P1
         KQ = KP + 1
         SIGN1 = 3.
      ELSE ! KFACE = 2
         K  = K2
         KP = K2M1
         KQ = KP - 1
         SIGN1 = -3.
      END IF

!     Set up the second derivatives and the P, Q, R terms:

      DO J = J1 + 1, J2 - 1
         DO I = I1 + 1, I2 - 1
            XZZ = -3.5*X(I,J,K,1) + 4.*X(I,J,KP,1) - HALF*X(I,J,KQ,1) -
     >             SIGN1*F1ZETA(I,J)
            YZZ = -3.5*X(I,J,K,2) + 4.*X(I,J,KP,2) - HALF*X(I,J,KQ,2) -
     >             SIGN1*F2ZETA(I,J)
            ZZZ = -3.5*X(I,J,K,3) + 4.*X(I,J,KP,3) - HALF*X(I,J,KQ,3) -
     >             SIGN1*F3ZETA(I,J)

            RAJ  = RJACK(I,J)
            RAJ2 = RAJ**2

            ALPHA = (GAMZ13(I,J)**2 + GAMZ23(I,J)**2 + GAMZ33(I,J)**2) *
     >               RAJ2

            H1 = FZ1(I,J) - ALPHA*XZZ
            H2 = FZ2(I,J) - ALPHA*YZZ
            H3 = FZ3(I,J) - ALPHA*ZZZ

            TEMPLIM = MIN (ONE, FGLIMIT / (ONE + CURVK(I,J)))

            DIF = (H1*GAMZ11(I,J) + H2*GAMZ21(I,J) + H3*GAMZ31(I,J)) *
     >             RAJ - PZ(I,J)

            TEMP1 = TEMPLIM * MAX (ABS (PZ(I,J)), ONE)
            TEMP2 = MIN (FGRELAX * ABS (DIF), TEMP1)
            PZ(I,J) = PZ(I,J) + SIGN (TEMP2, DIF)

            DIF = (H1*GAMZ12(I,J) + H2*GAMZ22(I,J) + H3*GAMZ32(I,J)) *
     >             RAJ - QZ(I,J)

            TEMP1 = TEMPLIM * MAX (ABS (QZ(I,J)), ONE)
            TEMP2 = MIN (FGRELAX * ABS (DIF), TEMP1)
            QZ(I,J) = QZ(I,J) + SIGN (TEMP2, DIF)

            DIF = (H1*GAMZ13(I,J) + H2*GAMZ23(I,J) + H3*GAMZ33(I,J)) *
     >             RAJ - RZ(I,J)

            TEMP1 = TEMPLIM * MAX (ABS (RZ(I,J)), ONE)
            TEMP2 = MIN (FGRELAX * ABS (DIF), TEMP1)
            RZ(I,J) = RZ(I,J) + SIGN (TEMP2, DIF)
         END DO
      END DO

      END SUBROUTINE PQRZETA

      END SUBROUTINE ELLIP3D

*     Modules for program CONVERTQ:

      MODULE AERO ! Aerodynamic data

      INTEGER           :: NMA, NAL, NPR 
      REAL, ALLOCATABLE :: MACH(:), ALPHA(:), PRES(:), RE(:,:,:),
     >                     CL(:,:,:), CD(:,:,:), CM(:,:,:), ALT(:,:,:)

      END MODULE AERO
      
      MODULE ATMOS ! Atmospheric data

      INTEGER           :: MAXH
      REAL, ALLOCATABLE :: H(:), PH(:), TEMPH(:)

      END MODULE ATMOS

********************************************************************************
*
      PROGRAM CONVERTQ
*
*     Description:
*
*        CONVERTQ is a variant of CALC_RE, q.v.  Instead of reading HAVOC-
*        type aerodynamic data, it reads CL, CD, CM data that have first
*        been combined from CBAERO single-function tables into the format
*        produced by COMBINE_TABLES, with dynamic pressure in atmospheres,
*        not psf, and reference area and length inserted as shown below.
*
*        Like CALC_RE, it reads CL, CD, CM as functions of Mach, Alpha,
*        and Q and converts Q into log (Reynolds number) plus altitude,
*        saving results in the "Aerodyn Ma_Alpha_logRe" custom aerodynamic
*        model format handled by Traj/Traj_opt.
*
*     Input format:
*
*        HL-20 aerodynamic database
*         16  ! # Machs
*         16  ! # Alphas
*          8  ! # Qs (bars)
*        34.00251 ! Sref, m^2
*        8.607552 ! Lref, m
*          1  ! # multi-function tables
*        Multi-function table  1
*        Mach  Alpha  Q    CL       CD       CM
*         0.3  0.0 0.005  -0.0244   0.0654   0.00642
*         0.6  0.0 0.005  -0.0306   0.0822   0.00715
*         0.9  0.0 0.005  -0.0205   0.1112   0.00089
*         1.2  0.0 0.005   0.0054   0.1473  -0.00045
*         2.0  0.0 0.005   0.0010   0.1158  -0.00007
*         4.0  0.0 0.005  -0.0216   0.0996   0.00322
*         6.0  0.0 0.005  -0.0305   0.0936   0.00451
*         9.0  0.0 0.005  -0.0381   0.0944   0.00528
*          :    :   :       :        :        :
*        23.0  0.0 0.005  -0.0380   0.1122   0.00551
*        26.0  0.0 0.005  -0.0383   0.1164   0.00574
*          :    :   :       :        :        :
*          :    :   :       :        :        :
*         6.0 50.0 1.000   0.7307   0.9973  -0.06070
*         9.0 50.0 1.000   0.7367   1.0051  -0.06023
*        11.0 50.0 1.000   0.7411   1.0116  -0.06040
*        13.0 50.0 1.000   0.7417   1.0119  -0.06027
*        15.0 50.0 1.000   0.7427   1.0140  -0.06034
*        17.0 50.0 1.000   0.7438   1.0154  -0.06036
*        19.0 50.0 1.000   0.7455   1.0175  -0.06046
*        21.0 50.0 1.000   0.7470   1.0207  -0.06065
*        23.0 50.0 1.000   0.7480   1.0228  -0.06078
*        26.0 50.0 1.000   0.7489   1.0238  -0.06081
*
*     History:
*
*        Early 2000?  James Reuther   Original implementation (of CALC_RE).
*        07/07/00     David Saunders  Streamlined & modified to extrapolate
*                                     the coefficients to give rectangularity
*                                     in terms of Reynolds number (lost
*                                     during the conversion from Q).
*        07/09/00       "      "      Save results in additional PLOT3D form.
*        08/08/02       "      "      CONVERTQ adapted from CALC_RE.
*        11/05/02       "      "      Save original data in PLOT3D form too.
*        03/11/05       "      "      Leave the "Multi-function table  1" line
*                                     from COMBINE_TABLES in there - it is not
*                                     worth deleting it when these utilities
*                                     are expected to run in the ModelCenter
*                                     environment.
*        05/09/05       "      "      CBAERO switched from atmospheres to bars,
*                                     so do likewise here.
*        10/25/10       "      "      The atmosphere table being read has
*                                     pressures in millibars, so it should be
*                                     making use of the psf/bar conversion.
*                                     Also, it appears that Reynolds numbers
*                                     have been too small all along: GET_RE
*                                     should be given feet, not meters, because
*                                     it works with English units.
*
*     Origination:    NASA Ames Research Center, Moffett Field, CA
*
********************************************************************************

*     Global variables:

      USE AERO
      USE ATMOS

      IMPLICIT NONE

*     Local constants:

      INTEGER, PARAMETER ::
     >   LATMS = 1, LREAD = 2, LUNKBD = 5,
     >   LOUT1 = 7, LOUT2 = 8, LUNCRT = 6

      REAL, PARAMETER ::
     >   GAMMA = 1.4,
     >   HALF  = 0.5,
     >   PSF   = 2088.54342, ! psf per bar  (2116.4,  ! psf per atmosphere)
     >   R     = 1716.5      ! Gas constant for Sutherland's viscosity calc.

*     Local variables:

      INTEGER   :: I, J, K, N, NUNI, IOUT1, IOUT2, LEFT
      REAL      :: LENGTH, REFA, REFL, REMAX, REMIN, REYN
      CHARACTER :: FILENAME * 64

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   ORIGINALRE, ORIGINALCF, UNIFORMRE, UNIFORMCF

      REAl, ALLOCATABLE, DIMENSION (:,:,:) ::
     >   CLU, CDU, CMU, ALTU, CLOVERCD

*     Execution:

      WRITE (LUNCRT, '(/, A)', ADVANCE='NO')
     >   ' Input aero. database file: '
      READ (LUNKBD, '(A)') FILENAME

      OPEN (UNIT = LREAD, FILE = FILENAME,    STATUS='OLD')
      OPEN (UNIT = LATMS, FILE = 'atmos.inp', STATUS='OLD')

*     Read the atmospheric table.
*     Note that its pressures are labeled psf, but are actually millibars.

      READ (LATMS, *)
      READ (LATMS, *)
      READ (LATMS, *) MAXH

      ALLOCATE (H(MAXH), PH(MAXH), TEMPH(MAXH))

      READ (LATMS, *)

      DO I = 1, MAXH
         READ (LATMS,*) H(I), PH(I), TEMPH(I)
***      PH(I) = PH(I) * 2.0885468           ! Explain
         PH(I) = PH(I) * (PSF * 0.001)       ! Millibars -> psf
      END DO

      CLOSE (LATMS)

C     Read the aerodynamic database:

      READ (LREAD, *) ! Title
      READ (LREAD, *) NMA
      READ (LREAD, *) NAL
      READ (LREAD, *) NPR
      READ (LREAD, *) REFA  ! Square meters
      READ (LREAD, *) REFL  ! Meters

      ALLOCATE (MACH(NMA), ALPHA(NAL), PRES(NPR),
     >          CL(NMA,NAL,NPR), CD(NMA,NAL,NPR), CM(NMA,NAL,NPR),
     >          RE(NMA,NAL,NPR), ALT(NMA,NAL,NPR))

      READ (LREAD, *) ! "Multi-function table  1"
      READ (LREAD, *) ! "Mach  Alpha  Q  CL  CD  CM"
      READ (LREAD, *)
     >   (((MACH(I), ALPHA(J), PRES(K), CL(I,J,K), CD(I,J,K), CM(I,J,K),
     >      I = 1, NMA), J = 1, NAL), K = 1, NPR)

*     Convert Q to Reynolds # then to LOG (RE) and determine the range.
*     GET_RE works in English units, so pass the length in feet.

      LENGTH = REFL / 0.3048
      REMAX  = -30.
      REMIN  =  30.

      DO K = 1, NPR
         DO J = 1, NAL
            DO I = 1, NMA

               CALL GET_RE (PRES(K) * PSF, MACH(I), LENGTH,
     >                      REYN, GAMMA, R, ALT(I,J,K))

               RE(I,J,K) = LOG10 (REYN)
               REMAX = MAX (REMAX, RE(I,J,K))
               REMIN = MIN (REMIN, RE(I,J,K))

               write (20, '(3f10.6, 1p, 2e16.7)')
     >            mach(i), alpha(j), pres(k), reyn, re(i,j,k)
            END DO 
         END DO
      END DO

*     Round off the extremes to the nearest 0.5:

      REMIN = HALF * REAL (INT (REMIN / HALF))
      REMAX = HALF * REAL (INT (REMAX / HALF) + 1)

*     Set up log (Re) at uniform intervals of 0.5:

      NUNI = NINT ((REMAX - REMIN) / HALF) + 1

      ALLOCATE (UNIFORMRE(NUNI), UNIFORMCF(NUNI))

      DO I = 1, NUNI
         UNIFORMRE(I) = REMIN + REAL (I - 1) * HALF
      END DO

      ALLOCATE (CLU(NMA,NAL,NUNI),  CDU(NMA,NAL,NUNI),
     >          CMU(NMA,NAL,NUNI), ALTU(NMA,NAL,NUNI),
     >          CLOVERCD(NMA,NAL,NUNI),
     >          ORIGINALRE(NPR),   ORIGINALCF(NPR))

*     Interpolate CL, CD, CM, and Altitude at uniform log (Re) points:

      DO J = 1, NAL

         DO I = 1, NMA

            ORIGINALRE = RE(I,J,1:NPR)
            ORIGINALCF = CL(I,J,1:NPR)
            LEFT = NPR ! Since Re is decreasing

            CALL QINTERP (NPR, ORIGINALRE, ORIGINALCF,
     >                    NUNI, UNIFORMRE,  UNIFORMCF, LEFT)

            CLU(I,J,1:NUNI) = UNIFORMCF

            ORIGINALCF = CD(I,J,1:NPR)
            LEFT = NPR

            CALL QINTERP (NPR, ORIGINALRE, ORIGINALCF,
     >                    NUNI, UNIFORMRE,  UNIFORMCF, LEFT)

            CDU(I,J,1:NUNI) = UNIFORMCF

            DO K = 1, NUNI ! For PLOT3D purposes
               CLOVERCD(I,J,K) = CLU(I,J,K) / CDU(I,J,K)
            END DO

            ORIGINALCF = CM(I,J,1:NPR)
            LEFT = NPR

            CALL QINTERP (NPR, ORIGINALRE, ORIGINALCF,
     >                    NUNI, UNIFORMRE,  UNIFORMCF, LEFT)

            CMU(I,J,1:NUNI) = UNIFORMCF

            ORIGINALCF = ALT(I,J,1:NPR)
            LEFT = NPR

            CALL QINTERP (NPR, ORIGINALRE, ORIGINALCF,
     >                    NUNI, UNIFORMRE,  UNIFORMCF, LEFT)

            ALTU(I,J,1:NUNI) = MAX (UNIFORMCF, -150000.)

         END DO

      END DO

*     Save results in Traj format:

      OPEN (UNIT = LOUT1, FILE = 'out1.dat',  STATUS='UNKNOWN')

      WRITE (LOUT1, '(2F8.3, F6.2, 2F10.6, F11.6, F12.2)')
     >   (((MACH(I), ALPHA(J), UNIFORMRE(K),
     >      CLU(I,J,K), CDU(I,J,K), CMU(I,J,K), ALTU(I,J,K),
     >      I = 1, NMA), J = 1, NAL), K = 1, NUNI)

      CLOSE (LOUT1)

*     Save results in more readable form:

      OPEN (UNIT = LOUT2, FILE = 'out2.dat',  STATUS='UNKNOWN')

      WRITE (LOUT2, '(3I5, 2F8.3, F6.2, 2F10.6, F11.6, F12.2)')
     >   (((I, J, K,
     >      MACH(I), ALPHA(J), UNIFORMRE(K),
     >      CLU(I,J,K), CDU(I,J,K), CMU(I,J,K), ALTU(I,J,K),
     >      I = 1, NMA), J = 1, NAL), K = 1, NUNI)

      CLOSE (LOUT2)

*     Save results in PLOT3D form:

      OPEN (UNIT = LOUT1, FILE = 'aero-coefs.grid',    STATUS='UNKNOWN')
      OPEN (UNIT = LOUT2, FILE = 'aero-coefs.plot3d',  STATUS='UNKNOWN')

      WRITE (LOUT1, '(I1, /, 3I5)') 1, NMA, NAL, NUNI
      WRITE (LOUT1, '(1P, 5E15.7)')
     >   (( MACH,                      J = 1, NAL), K = 1, NUNI),
     >   (((ALPHA(J),     I = 1, NMA), J = 1, NAL), K = 1, NUNI),
     >   (((UNIFORMRE(K), I = 1, NMA), J = 1, NAL), K = 1, NUNI)

      CLOSE (LOUT1)

      WRITE (LOUT2, '(I1, /, 3I5)') 1, NMA, NAL, NUNI
      WRITE (LOUT2, '(4F3.0)')      1., 1., 1., 1.
      WRITE (LOUT2, '(1P, 5E15.7)')
     >   (((CLU(I,J,K),       I = 1, NMA), J = 1, NAL), K = 1, NUNI),
     >   (((CDU(I,J,K),       I = 1, NMA), J = 1, NAL), K = 1, NUNI),
     >   (((CMU(I,J,K),       I = 1, NMA), J = 1, NAL), K = 1, NUNI),
     >   (((CLOVERCD(I,J,K),  I = 1, NMA), J = 1, NAL), K = 1, NUNI),
     >   (((ALTU(I,J,K),      I = 1, NMA), J = 1, NAL), K = 1, NUNI)

      CLOSE (LOUT2)

*     Save original data in PLOT3D form also:

      OPEN (LOUT1, FILE = 'MAlphaQ-coefs.grid',   STATUS='UNKNOWN')
      OPEN (LOUT2, FILE = 'MAlphaQ-coefs.plot3d', STATUS='UNKNOWN')

      WRITE (LOUT1, '(I1, /, 3I5)') 1, NMA, NAL, NPR
      WRITE (LOUT1, '(1P, 5E15.7)')
     >   (( MACH,                  J = 1, NAL), K = 1, NPR),
     >   (((ALPHA(J), I = 1, NMA), J = 1, NAL), K = 1, NPR),
     >   (((PRES(K),  I = 1, NMA), J = 1, NAL), K = 1, NPR)

      WRITE (LOUT2, '(I1, /, 3I5)') 1, NMA, NAL, NPR
      WRITE (LOUT2, '(4F3.0)')      1., 1., 1., 1.
      WRITE (LOUT2, '(1P, 5E15.7)')
     >   (((CL(I,J,K),           I = 1, NMA), J = 1, NAL), K = 1, NPR),
     >   (((CD(I,J,K),           I = 1, NMA), J = 1, NAL), K = 1, NPR),
     >   (((CM(I,J,K),           I = 1, NMA), J = 1, NAL), K = 1, NPR),
     >   (((CL(I,J,K)/CD(I,J,K), I = 1, NMA), J = 1, NAL), K = 1, NPR),
     >   (((RE(I,J,K),           I = 1, NMA), J = 1, NAL), K = 1, NPR)

      END PROGRAM CONVERTQ

************************************************************************
*
      SUBROUTINE GET_RE (Q, MACH, LENGTH, RE, GAMMA, R, ALT)
*
*     Calculate Reynolds number from dynamic pressure, etc.
*     All quantities are in English units, except that LENGTH has been
*     passed as meters in Traj_opt usage from c. 2000 to October 2010.
*     LENGTH is now passed in feet to give a consistent nondimensional
*     Reynolds number, RE.
*     Q is in psf.
*     R = 1716.5 is the gas constant for air (perfect gas).
*     Temperatures are in degrees Rankine.
*     ALT is output as an altitude in feet.
*
************************************************************************

*     Global variables:

      USE ATMOS

      IMPLICIT NONE

*     Arguments:

      REAL    :: Q, MACH, LENGTH, RE, GAMMA, R, ALT

*     Local variables:

      INTEGER :: IT, MODE
      REAL    :: PRES, DHDP, DTDH, TEMP, V, MU, RHO

*     Execution:

      PRES = (Q + Q) / (GAMMA * MACH * MACH)

      MODE = 0

      CALL INTRP1 (MODE, MAXH, IT, PH, PRES, H, ALT, DHDP)

      CALL INTRP1 (MODE, MAXH, IT, H, ALT, TEMPH, TEMP, DTDH)

      RHO = PRES / (R * TEMP)

*     Use Sutherland's formula to compute viscosity.

      MU  = 2.270E-8 * ((TEMP ** 1.5) / (TEMP + 198.6))
      V   = MACH * SQRT (R * GAMMA * TEMP)
      RE  = RHO * V * LENGTH / MU

      END SUBROUTINE GET_RE

C***********************************************************************

      SUBROUTINE INTRP1 (MODE, N, IT, X, XT, F, FT, FACT)

C     INTRP1 performs one dimensional linear interpolation of a given
C     discrete function F(X) at a given target location XT.
C
C     ARG     DIM    TYPE  I/O/S     DESCRIPTION
C     ----    ---    ----  -----     -----------
C     MODE     -      I    I         MODE = 0:  interval not known;
C                                    MODE = 1:  interval known.
C     N        -      I    I         Dimension of X and F arrays.
C     IT       -      I    I/O       The upper end of the interval used
C                                    to interpolate the function.
C     X       (N)     R    I         Independent variable.
C     F       (N)     R    I         Dependent     "
C     XT       -      R    I         Target X.
C     FT       -      R      O       Linearly interpolated F.
C
C     Author:  James ???
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER  MODE, N, IT
      REAL     X(N), XT, F(N), FT, FACT

C     Local variables:

      INTEGER  I, IL, IU
      REAL     XA, XB, XC

C     Execution:

      IF (MODE == 0) THEN ! Find the relevant interval:

         XA = XT   - X(1)
         XB = XT   - X(N)
         XC = X(N) - X(1)

         IF (XA * XB < 0.) THEN

            DO I = 2, N

               XA = XT - X(I-1)
               XB = XT - X(I)

               IF (XA * XB < 0.) THEN
                  IT = I
                  EXIT
               END IF

            END DO

         ELSE

            IF (XC > 0.) THEN

               IF (XB > 0.) THEN
                  IT = N
               ELSE
                  IT = 2
               END IF

            ELSE

               IF (XA > 0.) THEN
                  IT = 2
               ELSE
                  IT = N
               END IF

            END IF

         END IF

      END IF

      IL = IT - 1
      IU = IT

      FACT = (F(IU) - F(IL)) / (X(IU) - X(IL))
      FT   = F(IL) + (XT - X(IL)) * FACT

      END SUBROUTINE INTRP1

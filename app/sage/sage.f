C*******************************************************************************
C
C     This is version 4 of SAGE. It started as sagev3 on April 2000. 
C     New features:
C       + Cosmetics by David Saunders in preparation for possible substantive
C         changes later.
C       + No extra I/O occurs when grid size changes.
C       + A new routine called DXSP computes a variable spacing for wall and 
C         outer bound expansions.  It can be a function of x, not s, and takes
C         two new inputs (xmax and xmin).
C       + Reynolds number can be used for wall spacing.
C       + SHOCKLE, BLCLUST and PUSHIT all work on a portion of the zone;
C         RECLUST is now 0 or 1 and ST and END parameters control the range.
C       + Formatted write statements are no longer free form.
C         They are appropriate for input to GASP and GRIDPRO, etc.
C       + WRITMULT no longer exists: MULTIO performs all the tasks.
C       + REWIND, using FORMOO, to correct changes between formatted and
C         unformatted files in multiple runs.
C         Note: don't mix formats in export. 5/31/01
C       + If Q file is a function file, input FUN = .TRUE.;
C         FUN = T need only be used on the first adaption.
C       + For non-standard Q file (#5), input NQDIM = number of variables.
C       + TWOD Q data are no longer shifted, so INDQ and IQ are consistent with
C         location.
C       + ALL COMMON BLOCKS HAVE BEEN REMOVED: ALLOCATABLE VARIABLES ARE USED.
C         THE MINIMUM SIZE OF EACH VARIABLE IS USED TO ALLOW LARGER GRIDS.  USE 
C         'unlimit' if there is still insufficient memory.
C       + MVC1, MVC2 added to allow 'constant' outbound add.
C
C*******************************************************************************
C
C     AUTHOR:  Carol Davies, ELORET/ASA Branch, NASA Ames Research Center
C              cdavies@mail.arc.nasa.gov
C
C     And, now for some legalese:
C
C*********************** SOFTWARE DISCLAIMER ***********************************
C
C     This software is provided "as is" without warranty of any kind, either
C     express, implied, or statutory, including, but not limited to, any
C     guarantee that the software will conform to specifications, any implied
C     warranties of merchantability, fitness for a particular purpose, and
C     freedom from infringement, and any warranty that the documentation will
C     conform to the program, or any warranty that the software will be error
C     free.  In no event shall NASA be liable for any damages, including, but
C     not limited to direct, indirect, special or consequential damages, arising
C     out of, resulting from, or in any way connected with this software,
C     whether or not based upon warranty, contract, tort or otherwise, and
C     whether or not loss was sustained from, or arose out of the results of,
C     or the use of, the software or services provided hereunder.
C
C     Copyright 1999, United States Government as represented by the 
C     Administrator of the National Aeronautics and Space Administration.
C
C     The end. Phew!     
C
C*******************************************************************************

      MODULE SAGEMOD

      INTEGER, PARAMETER ::
     >   NDIM = 5,        ! Explain
     >   NA   = 200       !    "

      REAL, PARAMETER ::
     >   ONE = 1.0, ZERO = 0.0

      LOGICAL 
     >   ISTEP, JSTEP, KSTEP, IINVERSE, JINVERSE, KINVERSE, TWOD,
     >   IJPLANE, IKPLANE, JKPLANE, MARCH, MARCHPL, SAVE, OK,
     >   FORMI, FORMO, FORMOO, FUN, GEOM, QFUN, NOQ, NOUP, FV, SINGRD,
     >   EXPORT, IMPORT, FORMCK, MVIO

      LOGICAL, DIMENSION(2) :: 
     >   ORTHS, ORTHE

      INTEGER 
     >   IMAX, JMAX, KMAX, IST, IEND, JST, JEND, KST, KEND, NIPTS, NADS,
     >   NITGI, NITQI, NITGO, NITQO, MGRID, NFILT, INTER, MGSTEPS,
     >   MGPLS, NQDIM, ADD, LSTADD, LENDADD, SUB, LSTSUB, LENDSUB,
     >   MAXITS, NFLAG, LNSING, PLSING, NOSESM, REMOVE, RECLUST, IMS,
     >   JME, NGRID, HIOFLAG, INDQ, NEDGE, MG1, MG2, MXFG, MPLANE,
     >   IS, IE, JS, JE, KS, KE, ID, JD, KD, IMX

      INTEGER ::
     >   IERA(9), IQ(50)

      INTEGER, ALLOCATABLE, DIMENSION (:) ::
     >   IMQ,  JMQ,  KMQ,  MGCT, IM, JM, KM, IMO, JMO, KMO,
     >   IMQO, JMQO, KMQO, NQFUN

      REAL ::
     >   FSMACH, ALPHA, RE, TIME, CONV, BK1, MVBOUND, MVP1, MVP2, DSW,
     >   DSE, DSW1, DSW2, A, B, WDS, WDE, DSMAX, DSMIN, DLENGS, DLENGE,
     >   WDSS, WDES, RDSMIN, RDSMAX, FGW, FQW, XMIN, XMAX, XMINP, XMAXP,
     >   MVC1, MVC2, MVC

      REAL, DIMENSION(2) ::
     >   CLAM, CT, TN, CNM, CLAMW, CTM, TNM

      REAL, ALLOCATABLE, DIMENSION (:) :: 
     >   F, FB, WEIGHT, SS, SMS, DS, SN, SNM, SNMK, WT, SP, SPPL, DAP,
     >   DAPPL, FQ, FG, FGS, SMSS

      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   COSE, COSU, COSB, COSEK, COSUK, COSBK, XP, YP, ZP

      REAL, ALLOCATABLE, DIMENSION (:,:,:) ::
     >   X, Y, Z, XJ, YJ, ZJ, QP

      REAL, ALLOCATABLE, DIMENSION (:,:,:,:) ::
     >   Q

C     Local work-space:

      REAL, ALLOCATABLE, DIMENSION (:,:,:) ::
     >   T1

      END MODULE SAGEMOD

C*******************************************************************************
C
      PROGRAM SAGE

C     SAGE (Self-Adaptive Grid codE) adjusts a 2-D or 3-D multiblock grid in a
C     choice of ways, most notably adapting it to a preliminary flow solution.
C
C     SAGE Version 2 is documented in NASA/TM-1995-110350.
C     SAGE Version 3 is documented in NASA/TM-1999-208792.
C
C     This driving program loops over j lines within k planes.
C     Subroutines perform all the action.
C
C*******************************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Local variables:

      INTEGER
     >   I, J, K

      LOGICAL
     >   NOMORE, MTCH

C     Procedures:

      EXTERNAL
     >   BLOCK, DEALLO, FBAR, INITIAL, LINE1, MARCHJ, NOADAPT, OUTPUT,
     >   SETUPJ, SETUPK, SINGPLN, SOLUT, TORCOF, TORSION, UPDATE, WTEDGE

C     Execution:

C     Start adaptions:  NADS keeps the count of NAMEL input passes.

      DO 200 NADS = 1, NA

C        Read input controls and initialize.
C        Note that all reclustering is handled by INITIAL.

         CALL INITIAL (NOMORE, MTCH)

         IF (NOMORE) GO TO 900 ! Last NAMEL has been read.
         IF (MTCH)   GO TO 200
         IF (NOUP)   GO TO 150

C        Start main adaption loop.
C        Adapt each line j in plane k, then step to the next j plane.

         DO 125 K = KST, KEND

            IF (K /= KST) CALL TORCOF (2, K, KST, KEND, MGPLS, MARCHPL)

            DO 100 J = JST, JEND

               IF ((J == JST .AND. MGSTEPS /= 0) .OR.
     >             (K /= KST .AND. LNSING  == J) .OR.
     >             (K == KST .AND. MGPLS   /= 0)) THEN

                  CALL NOADAPT (J, K)
                  CYCLE

               END IF

               CALL BLOCK (J, K)
               CALL FBAR  (J, K)

               IF (.NOT. OK) GO TO 900

               IF (J == JST .AND. K == KST) THEN
                  CALL LINE1
                  IF (PLSING == KST) THEN
                     CALL SINGPLN
                     GO TO 125
                  END IF
               ELSE
                  IF (J /= JST) CALL SETUPJ (J, K)
                  IF (K /= KST) CALL SETUPK (J, K)
                  CALL TORCOF  (1, J, JST, JEND, MGSTEPS, MARCH)
                  CALL TORSION (J, K)
                  CALL SOLUT   (J, K)
               END IF

               CALL UPDATE (J, K)

               IF (.NOT. OK) THEN
                  CALL OUTPUT
                  GO TO 900
               END IF

               IF (NEDGE /= 0) CALL WTEDGE (J, K)

  100       CONTINUE ! Next line

            IF (MARCH .AND. J /= JMAX) CALL MARCHJ (K)

  125    CONTINUE ! Next plane

  150    CONTINUE

         CALL OUTPUT
         CALL DEALLO

  200 CONTINUE ! Next adaption

  900 CONTINUE ! Done

C     Inform user of unit numbers for final output files:

      IF (NITGO /= 0) THEN
         IF (NOQ) THEN 
            WRITE (6, '(/, A, I3)') ' Output file on unit', NITGO
         ELSE
            WRITE (6, '(/, A, I3, A, I3)')
     >         ' Output files on units', NITGO, ' and', NITQO
         END IF
      ELSE
         WRITE (6, '(/, A)') ' No output files.'
      END IF

      END PROGRAM SAGE

C***********************************************************************
C
      SUBROUTINE ADDPTS
C
C     ADDPTS increases the number of points along the adaption line.
C
C     CALLED BY: INITIAL
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Local variables:

      INTEGER
     >   I, IMAXB, INTO, J, K, L, L1, M, MEND, N

      REAL
     >   DENOM, XS, YS, ZS,
     >   XT(IMX), YT(IMX), ZT(IMX), QT(IMX,NDIM)

C     Execution:

C     Initialize start & end limits to be consistent with adaption domain:

      IF (LSTADD  == 0) LSTADD  = IST
      IF (LENDADD == 0) LENDADD = IEND
      IF (LSTADD > LENDADD) THEN
         L1 = LENDADD
         LENDADD = LSTADD
         LSTADD = L1
      END IF
      IF (LSTADD == LENDADD) THEN
         WRITE (6, '(/, A)') ' No points added.'
         GO TO 999
      END IF
      IMAXB = IMAX
      IMAX  = IMAX + (LENDADD - LSTADD)*ADD

C     Increase points in requested region:

      MEND = (LENDADD - LSTADD - 1)*(ADD + 1) + LSTADD
      DENOM = REAL (ADD + 1)

      DO 600 K = 1, KMAX

         DO 500 J = 1, JMAX

            M = LSTADD - (ADD + 1)

            DO I = LSTADD, LENDADD - 1

               M = M + ADD + 1
               XT(M) = X(I,J,K)
               YT(M) = Y(I,J,K)
               ZT(M) = Z(I,J,K)

               IF (.NOT. NOQ) QT(M,1:NQDIM) = Q(I,J,K,1:NQDIM)

               XS = (X(I+1,J,K) - X(I,J,K)) / DENOM
               YS = (Y(I+1,J,K) - Y(I,J,K)) / DENOM
               ZS = (Z(I+1,J,K) - Z(I,J,K)) / DENOM

               DO L = 1, ADD
                  XT(M+L) = X(I,J,K) + XS*L
                  YT(M+L) = Y(I,J,K) + YS*L
                  ZT(M+L) = Z(I,J,K) + ZS*L
                  IF (.NOT. NOQ) THEN
                     DO N = 1, NQDIM
                        QT(M+L,N) = Q(I,J,K,N) + (REAL (L) / DENOM) *
     >                             (Q(I+1,J,K,N) - Q(I,J,K,N))
                     END DO
                  END IF
               END DO

            END DO

C           Put remaining data into temporary arrays:

            DO I = LENDADD, IMAXB
               INTO = M + ADD + 1 + (I - LENDADD)
               XT(INTO) = X(I,J,K)
               YT(INTO) = Y(I,J,K)
               ZT(INTO) = Z(I,J,K)
               IF (.NOT. NOQ) QT(INTO,1:NQDIM) = Q(I,J,K,1:NQDIM)
            END DO

C           Return data in original arrays:

            DO I = LSTADD, IMAX
               X(I,J,K) = XT(I)
               Y(I,J,K) = YT(I)
               Z(I,J,K) = ZT(I)
               IF (.NOT. NOQ) Q(I,J,K,1:NQDIM) = QT(I,1:NQDIM)
            END DO

  500    CONTINUE ! Next J

  600 CONTINUE ! Next K

      IF (IEND > IST) IEND = IEND + (LENDADD - LSTADD)*ADD
      IF (IST  > IEND) IST = IST  + (LENDADD - LSTADD)*ADD

      WRITE (6, '(/, A, I4, A, I4)')
     >   ' Number of points increased from', IMAXB, ' to', IMAX

  999 CONTINUE

      END SUBROUTINE ADDPTS

C***********************************************************************
C
      SUBROUTINE ADDV (COSX1, COSY1, COSZ1, A1, COSX2, COSY2, COSZ2, A2,
     >                 COSX,  COSY,  COSZ)
C
C     ADDV calculates the direction cosines of a linear combination of
C     two 3-space vectors.
C
C     CALLED BY: NORMPT, SETUPJ, SETUPK, TORSION, VMERGE
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      REAL, INTENT (IN) ::
     >   COSX1, COSY1, COSZ1,      ! Elements of input vectors
     >   COSX2, COSY2, COSZ2,
     >   A1, A2                    ! Multipliers

      REAL, INTENT (OUT) ::
     >   COSX, COSY, COSZ          ! Desired direction cosines

C     Local variables:

      REAL
     >   VMOD

      COSX = A1*COSX1 + A2*COSX2
      COSY = A1*COSY1 + A2*COSY2
      COSZ = A1*COSZ1 + A2*COSZ2
      VMOD = SQRT (COSX**2 + COSY**2 + COSZ**2)

      IF (VMOD /= 0.0) THEN
         COSX = COSX/VMOD
         COSY = COSY/VMOD
         COSZ = COSZ/VMOD
      END IF

      END SUBROUTINE ADDV

C***********************************************************************
C
      SUBROUTINE ALOALO
C
C     Allocate (almost) all arrays that are dependent on the grid size.
C
C     CALLED BY: HEADIO, which first allocates the dimension arrays.
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Local variables:

      INTEGER
     >   IDIOM, JDIOM, KDIOM, IDZM,  JDZM,  KDZM,
     >   IMAXT, JMAXT, KMAXT, MSI2,  MSJ2,  MSK2,
     >   IMX2, LADD, LEND, N, NEWPTS

C     Execution:

C     (1) Determine the size to enable all zones to be read.
C         This is needed even if it is not a requested zone. 

      MSI2 = IM(1)
      MSJ2 = JM(1)
      MSK2 = KM(1)

      DO N = 2, NGRID
         MSI2 = MAX (MSI2, IM(N))
         MSJ2 = MAX (MSJ2, JM(N))
         MSK2 = MAX (MSK2, KM(N))
      END DO

      IMX2  = MAX (MSI2, MSJ2, MSK2)
      IDIOM = MSI2
      JDIOM = MSJ2
      KDIOM = MSK2

C     (2) For this adaption pass, does the zone increase in size?
C         Note that we don't care if it decreases.

      IMAXT = IM(MGRID)
      JMAXT = JM(MGRID)
      KMAXT = KM(MGRID)

C     Carol:  What is the precise definition of MGRID?

      IF (ADD /= 0) THEN
         IF (IEND == 0) IEND = IM(MGRID)
         IF (JEND == 0) JEND = JM(MGRID)
         IF (KEND == 0) KEND = KM(MGRID)
         LADD = LSTADD
         IF (LSTADD == 0) THEN
            IF (IJPLANE .AND. JSTEP .OR. IKPLANE .AND. KSTEP) LADD = IST
            IF (IJPLANE .AND. ISTEP .OR. JKPLANE .AND. KSTEP) LADD = JST
            IF (IKPLANE .AND. ISTEP .OR. JKPLANE .AND. JSTEP) LADD = KST
         END IF
         LEND = LENDADD
         IF (LENDADD == 0) THEN
           IF (IJPLANE .AND. JSTEP .OR. IKPLANE .AND. KSTEP) LEND = IEND
           IF (IJPLANE .AND. ISTEP .OR. JKPLANE .AND. KSTEP) LEND = JEND
           IF (IKPLANE .AND. ISTEP .OR. JKPLANE .AND. JSTEP) LEND = KEND
         END IF
         NEWPTS = ABS ((LEND - LADD)*ADD)
         IF (IJPLANE .AND. JSTEP .OR. IKPLANE .AND. KSTEP) THEN
            IMAXT = NEWPTS + IMAXT
         END IF
         IF (IJPLANE .AND. ISTEP .OR. JKPLANE .AND. KSTEP) THEN
            JMAXT = NEWPTS + JMAXT
         END IF
         IF (IKPLANE .AND. ISTEP .OR. JKPLANE .AND. JSTEP) THEN
            KMAXT = NEWPTS + KMAXT
         END IF
      END IF

C     (3) For this pass, is there any swapping? 

      IF (TWOD) KMAXT = 1
      IDZM = IMAXT
      JDZM = JMAXT
      KDZM = KMAXT

      IF (IJPLANE .AND. ISTEP) THEN
         IDZM = MAX (IMAXT, JMAXT)
         JDZM = MAX (IMAXT, JMAXT)
      END IF
      IF (JKPLANE .AND. KSTEP) THEN
         IDZM = MAX (IMAXT, JMAXT, KMAXT)
         JDZM = MAX (IMAXT, JMAXT, KMAXT)
         KDZM = MAX (IMAXT, JMAXT, KMAXT)
      END IF  
      IF (IKPLANE .AND. ISTEP) THEN
         IDZM = MAX (IMAXT, JMAXT, KMAXT)
         JDZM = MAX (IMAXT, JMAXT, KMAXT)
         KDZM = MAX (IMAXT, JMAXT, KMAXT)
      END IF
      IF (JKPLANE .AND. JSTEP) THEN
         IDZM = MAX (IMAXT, KMAXT)
         KDZM = MAX (IMAXT, KMAXT)
      END IF
      IF (IKPLANE .AND. KSTEP) THEN
         JDZM = MAX (JMAXT, KMAXT)
         KDZM = MAX (JMAXT, KMAXT)
      END IF

C     Now find the maximum of the three dimensions:

      ID  = MAX (IDIOM, IDZM)
      JD  = MAX (JDIOM, JDZM)
      KD  = MAX (KDIOM, KDZM)
      IMX = MAX (IDZM, JDZM, KDZM)

C     Allocate space:

      ALLOCATE (F(IMX), FB(IMX), WEIGHT(IMX), STAT=IERA(1))
      ALLOCATE (SS(IMX), SMS(IMX), DS(IMX), SN(IMX), 
     >          SNM(IMX), SNMK(IMX), WT(IMX), STAT=IERA(2))
      ALLOCATE (X(ID,JD,KD), Y(ID,JD,KD), Z(ID,JD,KD), STAT=IERA(3))
      IF (.NOT. NOQ) ALLOCATE (Q(ID,JD,KD,NQDIM), STAT=IERA(4))
      ALLOCATE (COSE(IMX, 3), COSU(IMX, 3), COSB(IMX, 3), 
     >          COSEK(IMX,3), COSUK(IMX,3), COSBK(IMX,3), STAT=IERA(5))
      ALLOCATE (XJ(IMX,3,3), YJ(IMX,3,3), ZJ(IMX,3,3), STAT=IERA(6))
      ALLOCATE (SP(IMX), SPPL(IMX), DAP(IMX), DAPPL(IMX), STAT=IERA(7))
      ALLOCATE (FQ(IMX), FG(IMX), FGS(IMX), SMSS(IMX), STAT=IERA(8))
      ALLOCATE (T1(ID,JD,KD), STAT=IERA(9))

      DO N = 1, 9
         IF (IERA(N) /= 0) THEN
            WRITE (6, '(/, A, I2, A)')
     >         ' ***** Trouble allocating space. Case:', N, ' *****'
            STOP
         END IF
      END DO

      END SUBROUTINE ALOALO

C***********************************************************************
C
      SUBROUTINE BLCLUST
C
C     BLCLUST reclusters at the surface when there is no outer boundary
C     movement.
C     Carol:  Is this true?  Is reclustering too difficult if the outer
C             boundary does move?
C
C     CALLED BY: INITIAL
C
C     CALLS: GETDS0, LAGCOF, VINOKUR
C
C***********************************************************************
C
      USE SAGEMOD

      IMPLICIT NONE

C     Local variables:

      INTEGER
     >   I, IER, J, K, L, M, N, NPTS

      REAL
     >   DS0, DS1, P1, P2, P3,
     >   SNW(IMX), DSC(IMX), SC(IMX), SW(IMX),
     >   XM(IMX),  YM(IMX),  ZM(IMX), QM(IMX,NDIM)

C     Execution:

      NPTS = IEND - IST + 1

      DO 900 K = KST, KEND

         DO 800 J = JST, JEND

C           If this is a common line, copy over from the first plane:

            IF (K /= KST .AND. J == LNSING) THEN
               DO I = 1, IMAX
                  X(I,J,K) = X(I,J,1)
                  Y(I,J,K) = Y(I,J,1)
                  Z(I,J,K) = Z(I,J,1)
               END DO
               IF (.NOT. NOQ) THEN
                  DO I = 1, IMAX
                     Q(I,J,K,1:NQDIM) = Q(I,J,1,1:NQDIM)
                  END DO
               END IF
               CYCLE ! Next J
            END IF

C           Cumulative arc lengths for this line:

            SC(1) = ZERO
            L = IST
            DO I = 1, NPTS - 1
               DSC(I) = SQRT ((X(L+1,J,K) - X(L,J,K))**2 +
     >                        (Y(L+1,J,K) - Y(L,J,K))**2 +
     >                        (Z(L+1,J,K) - Z(L,J,K))**2)
               IF (DSC(I) == ZERO) THEN
                  WRITE (6, '(A, 3(I3, 1X), A, I4, A)')
     >               ' WARNING: Zero edge at ', I, J, K,
     >               ', MGRID:', MGRID, '  SAGE will try to correct.'
                  IF (I < NPTS - 1) THEN
                     DSC(I) = (SQRT ((X(L+2,J,K) - X(L+1,J,K))**2 +
     >                               (Y(L+2,J,K) - Y(L+1,J,K))**2 +
     >                               (Z(L+2,J,K) - Z(L+1,J,K))**2)) *0.5
                  END IF
               END IF
               SC(I+1) = SC(I) + DSC(I)
               L = L + 1
            END DO

C           Find delta s at b.l. edge (DS1 defined by RECLUST) & wall (DS0):

C           Carol: Are DS0 and DS1 described correctly?
C           Carol: Is I = 1 always at a wall?
C           Carol: Why is DSE = 5.0 (default) treated differently for DS1?

            DS1 = DSC(NPTS-1)
            IF (DSE /= 5.0) DS1 = DSE*SC(NPTS) / REAL (NPTS-1)
            IF (DSW == ZERO) THEN
               DS0 = DSC(1)
            ELSE
               CALL GETDS0 (DS0, J, K)
            END IF

C           Vinokur distribution along the arc:

            SNW(1) = ZERO
            SNW(NPTS) = SC(NPTS)

CCCCCC      CALL CLUST2 (SNW, DS0, DS1, NPTS, IMX)
            CALL VINOKUR (1, NPTS, DS0, DS1, SNW, 6, IER)

            IF (IER /= 0) THEN ! ! Highly unlikely
               WRITE (6, '(/, A, /, 4I5, /, A, 1P, 2E16.7)')
     >         ' Vinokur distribution did not converge in BLCLUST.',
     >         ' IER, J, K, NPTS:', IER, J, K, NPTS,
     >         ' DS0, DS1: ', DS0, DS1
               STOP
            END IF

C           Interpolate at new arc length locations:

            DO I = 1, NPTS - 1

               CALL LAGCOF (SNW(I), SC, NPTS, M, P1, P2, P3)

               M = IST - 1 + M
               XM(I) = P1*X(M,J,K) + P2*X(M+1,J,K) + P3*X(M+2,J,K)
               YM(I) = P1*Y(M,J,K) + P2*Y(M+1,J,K) + P3*Y(M+2,J,K)
               ZM(I) = P1*Z(M,J,K) + P2*Z(M+1,J,K) + P3*Z(M+2,J,K)

               IF (.NOT. NOQ) THEN
                  DO N = 1, NQDIM
                     QM(I,N) = P1*Q(M,J,K,N) + P2*Q(M+1,J,K,N) +
     >                         P3*Q(M+2,J,K,N)
                  END DO
               END IF

            END DO

C           Transfer results back to X, Y, Z, Q:

            L = IST
            DO I = 1, NPTS - 1
               X(L,J,K) = XM(I)
               Y(L,J,K) = YM(I)
               Z(L,J,K) = ZM(I)
               IF (.NOT. NOQ) Q(L,J,K,1:NQDIM) = QM(I,1:NQDIM)
               L = L + 1
            END DO

  800    CONTINUE ! Next J

  900 CONTINUE ! Next K

      END SUBROUTINE BLCLUST

C***********************************************************************
C
      SUBROUTINE BLOCK (J, K)
C
C     This routine extracts the "block" of the grid around the current
C     J line, i.e., planes K-1 : K+1 and lines J-1 : J+1 for all I.
C     This block is used for all the major calculations, particularly
C     for the normals.  It prevents the original grid from being
C     contaminated when projections are required.
C
C     CALLED BY: Main progam
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: J, K

C     Local variables:

      INTEGER  I, J1, K1, J1ST, J1END, K1ST, K1END, LI, LJ, LK

C     Execution:

      K1ST  = 1
      K1END = 3
      J1ST  = 1
      J1END = 3

      IF (J == 1) THEN
         J1ST = 2
         XJ(1:NIPTS,1,1:3) = 999.0
      END IF

      IF (J == JMAX) THEN
         J1END = 2
         XJ(1:NIPTS,3,1:3) = 999.0
      END IF

      IF (K == 1) THEN
         K1ST = 2
         XJ(1:NIPTS,1:3,1) = 999.0
      END IF

      IF (K == KMAX) THEN
         K1END = 2
         XJ(1:NIPTS,1:3,3) = 999.0
      END IF

      DO K1 = K1ST, K1END
         LK = K - 2 + K1
         DO J1 = J1ST, J1END
            LJ = J - 2 + J1
            LI = IST
            DO I = 1, NIPTS
               XJ(I,J1,K1) = X(LI,LJ,LK)
               YJ(I,J1,K1) = Y(LI,LJ,LK)
               ZJ(I,J1,K1) = Z(LI,LJ,LK)
               LI = LI + 1
            END DO
         END DO
      END DO

      END SUBROUTINE BLOCK

C***********************************************************************
C
      SUBROUTINE CHSMALL (MG)
C
C     KLUDGE ALERT:  With 64-bit arithmetic (-r8 or equivalent switch,
C     always recommended for CFD), very small numbers with 3-digit
C     exponents have been observed.  This conflicts with the 1P, E19.11
C     choice of format in MULTIO.  (List-directed output evidently has
C     its own drawbacks.)  Therefore, fudge the tiny numbers.  This
C     version preserves true zeroes and the signs of the numbers.
C
C     CALLED BY: MULTIO
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: MG ! Grid block number

      REAL, PARAMETER :: EPS = 1.0E-17

      INTEGER I, J, K, N

      DO K = 1, KM(MG)
         DO J = 1, JM(MG)
            DO I = 1, IM(MG)
               IF (ABS (X(I,J,K)) < EPS) THEN
                  IF (X(I,J,K) /= ZERO) X(I,J,K) = SIGN (EPS, X(I,J,K))
               END IF
               IF (ABS (Y(I,J,K)) < EPS) THEN
                  IF (Y(I,J,K) /= ZERO) Y(I,J,K) = SIGN (EPS, Y(I,J,K))
               END IF
               IF (ABS (Z(I,J,K)) < EPS) THEN
                  IF (Z(I,J,K) /= ZERO) Z(I,J,K) = SIGN (EPS, Z(I,J,K))
               END IF

               IF (.NOT. NOQ) THEN
                  DO N = 1, NQDIM
                     IF (ABS (Q(I,J,K,N)) < EPS) THEN
                          IF (Q(I,J,K,N) /= ZERO)
     >                        Q(I,J,K,N)  = SIGN (EPS, Q(I,J,K,N))
                     END IF
                  END DO 
               END IF

            END DO 
         END DO 
      END DO 

      END SUBROUTINE CHSMALL

C***********************************************************************

      SUBROUTINE CROSSV (XT, YT, ZT, XT1, YT1, ZT1, DST, COSV, AAP,
     >                   DAPT, ICROSS, J)
C
C     CROSSV finds the intersection of vector T (from line J-1) and
C     J line segment.  COSV is the vector cosine of T.  <Direction cosines?>
C     DAP, AAP and ICROSS are computed.  <Needs clearer description.>
C
C     CALLED BY: TORSION
C
C     CALLS: DETERM, PURPLE, UNITV
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

      REAL, INTENT (IN), DIMENSION (IMX) ::
     >   XT, YT, ZT, XT1, YT1, ZT1, DST

      REAL, INTENT (IN) ::
     >   COSV(IMX,3)

      REAL, INTENT (OUT) ::
     >   AAP(IMX), DAPT(IMX)

      INTEGER, INTENT (OUT) ::
     >   ICROSS(IMX)

      INTEGER, INTENT (IN) ::
     >   J

C     Local variables:

      INTEGER
     >   I, L

      REAL
     >   DAPX, DAPY, DAPZ,
     >   ADX,  ADY,  ADZ,
     >   DAX,  DAY,  DAZ,                ! Direction cosines
     >   DNX,  DNY,  DNZ,                ! Unit normal
     >   V1,   V2,   V3,
     >   AAD,  BIGD,
     >   SSX(IMX), SSY(IMX), SSZ(IMX)

C     Execution:

C     Calculate s segments (unit vectors at each I) along the indicated J line:
C     Carol: J doesn't get used here.  Can you clarify what's going on?

      DO I = 1, NIPTS - 1

         CALL UNITV (XT(I), YT(I), ZT(I), XT(I+1), YT(I+1), ZT(I+1),
     >               SSX(I), SSY(I), SSZ(I))
      END DO

C     Find DAPT, AAP lengths that define s':

      ICROSS(1) = 1

      DO 200 I = 2, NIPTS - 1

         DO 120 L = ICROSS(I-1), NIPTS - 1

C           Find ad vector and n, normal to plane:

            CALL UNITV (XT1(I), YT1(I), ZT1(I), XT(L), YT(L), ZT(L),
     >                  DAX, DAY, DAZ)
            CALL PURPLE (SSX(L), SSY(L), SSZ(L), DAX, DAY, DAZ,
     >                   DNX, DNY, DNZ, NFLAG)

            V1 = -COSV(I,1)
            V2 = -COSV(I,2)
            V3 = -COSV(I,3)
            ADX = XT1(I) - XT(L)
            ADY = YT1(I) - YT(L)
            ADZ = ZT1(I) - ZT(L)

            CALL DETERM (SSX(L), V1, DNX, SSY(L), V2, DNY, SSZ(L), V3,
     >                   DNZ, BIGD)
            IF (BIGD == ZERO) CYCLE

            CALL DETERM (ADX, V1, DNX, ADY, V2, DNY, ADZ, V3, DNZ, AAD)

            AAP(I) = AAD/BIGD
            IF (AAP(I)  <  ZERO) AAP(I) = ZERO

            IF (AAP(I)  <  DST(L)) THEN ! AAP found: now find DAPT
               DAPX = AAP(I)*SSX(L) - ADX
               DAPY = AAP(I)*SSY(L) - ADY
               DAPZ = AAP(I)*SSZ(L) - ADZ
               DAPT(I) = SQRT (DAPX**2 + DAPY**2 + DAPZ**2)
               ICROSS(I) = L
               GO TO 199
            END IF   

  120    CONTINUE ! Next L

         ICROSS(I) = ICROSS(I-1)
         DAPT(I) = DAPT(I-1)

  199    CONTINUE
         IF (DAPT(I) <= ZERO) DAPT(I) = DAPT(I-1)

  200 CONTINUE ! Next I

      END SUBROUTINE CROSSV

C***********************************************************************
C
      SUBROUTINE DEALLO
C
C     Deallocate zonal work-space so it can be allocated for each pass.
C
C     CALLED BY: Main program
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

      INTEGER N

      DEALLOCATE (F, FB, WEIGHT, STAT=IERA(1))
      DEALLOCATE (SS, SMS, DS, SN, SNM, SNMK, WT, STAT=IERA(2))
      DEALLOCATE (X, Y, Z, STAT=IERA(3))
      IF (.NOT. NOQ) DEALLOCATE (Q, STAT=IERA(4))
      DEALLOCATE (COSE, COSU, COSB, COSEK, COSUK, COSBK, STAT=IERA(5))
      DEALLOCATE (XJ, YJ, ZJ, STAT=IERA(6))
      DEALLOCATE (SP, SPPL, DAP, DAPPL, STAT=IERA(7))
      DEALLOCATE (FQ, FG, FGS, SMSS, STAT=IERA(8))
      DEALLOCATE (T1, STAT=IERA(9))

      DO N = 1, 9
         IF (IERA(N) /= 0) THEN
            WRITE (6, '(/, A, I2, A)')
     >         ' ***** Trouble deallocating space. Case:', N, ' *****'
            STOP
         END IF
      END DO

      END SUBROUTINE DEALLO

C***********************************************************************
C
      SUBROUTINE DETERM (A1, B1, C1, A2, B2, C2, A3, B3, C3, DET)
C
C     Compute the indicated 3rd-order determinant.
C
C     CALLED BY: CROSSV
C
C***********************************************************************

      REAL  A1, B1, C1, A2, B2, C2, A3, B3, C3, DET

      DET = A1*(B2*C3 - C2*B3) - B1*(A2*C3 - C2*A3) + C1*(A2*B3 - B2*A3)

      END SUBROUTINE DETERM

C***********************************************************************
C
      SUBROUTINE DLENG (JL, K)
C
C     Finds the multipliers for edge treatment when NEDGE /= 0.
C
C     CALLED BY: WTEDGE
C
C***********************************************************************
C
      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN ):: JL, K

C     Local variables:

      INTEGER  I, JL1, JL2

C     Execution:

C     Calculations depend on whether data exists external to adaption domain.

C     At first edge point:

      JL1 = JL-1
      IF (JL1 == 0) JL1 = 2
      JL2 = JL+1
      IF (JL2 > JMAX) JL2 = JMAX

C     Carol:   Does the following assume preset to zero for undefined data?
C     Carol:   Why is X being tested for zero, and not Y or Z?
C     Carol:   Also, can we be sure it won't access X(0,*,*) if IST = 1?

      IF (IST == 1 .OR. (X(IST-1,JL,K) == ZERO  .AND. 
     >    (X(IST-1,JL1,K) == ZERO .OR. X(IST-1,JL2,K) == ZERO))) THEN
         I = IST
      ELSE
         I = IST - 1
      END IF

      DLENGS = SQRT ((X(I+1,JL,K) - X(I,JL,K))**2 +
     >               (Y(I+1,JL,K) - Y(I,JL,K))**2 +
     >               (Z(I+1,JL,K) - Z(I,JL,K))**2)

C     At last edge point:

      IF (IEND == IMAX .OR. (X(IEND+1,JL,K) == ZERO .AND. 
     >    (X(IEND+1,JL1,K) == ZERO .OR. X(IEND+1,JL2,K) == ZERO))) THEN
         I = IEND
      ELSE
         I = IEND + 1
      END IF

      DLENGE = SQRT ((X(I,JL,K) - X(I-1,JL,K))**2 +
     >               (Y(I,JL,K) - Y(I-1,JL,K))**2 +
     >               (Z(I,JL,K) - Z(I-1,JL,K))**2)

      END SUBROUTINE DLENG

C***********************************************************************
C
      SUBROUTINE DXSP (DS0, II, PST, PEND, XAT1, XAT2, J, K)
C
C     This routine computes a value based on the current X or S using
C     the input first and last values (in PST and PEND).
C     If DSW = 999., PST and PEND contain first and last wall spacing
C     either at input XAT1 and XAT2, or at first and last S for this zone.
C     If MVBOUND = 900., then PST and PEND contain first and last % for 
C     a boundary move, again using XAT1 and XAT2 if provided, otherwise
C     using the first and last S of the current zone.
C
C     CALLED BY: GETDS0 if DSW = 999. (varying wall spacing for Vinokur
C                                      algorithm);
C                SHOCKLE if MVBOUND = 900. (varying % boundary move)
C                        or if a varying constant move.
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)  :: II, J, K
      REAL,    INTENT (IN)  :: PST, PEND, XAT1, XAT2
      REAL,    INTENT (OUT) :: DS0

C     Local variables:

      INTEGER  L
      REAL     DSCW, DXMAX, XT(IMX)

C     Execution:

C     For a boundary move, look at distance along outer edge.
C     For a body reclustering, use the surface point.

      IF (XAT1 == ZERO .AND. XAT2 == ZERO) THEN

C        If no XMIN and XMAX input, compute s for this zone as X:

         XT(1) = ZERO     ! XT(II) looked wrong here
         DO L = 1, JMAX - 1
            DSCW = SQRT ((X(II,L+1,K) - X(II,L,K))**2 +
     >                   (Y(II,L+1,K) - Y(II,L,K))**2 +
     >                   (Z(II,L+1,K) - Z(II,L,K))**2)
            XT(L+1) = XT(L) + DSCW
         END DO
         DXMAX = XT(JMAX)

      ELSE

C        Need X at all J, but it could be stored in Y or Z:
C        Carol: Can the second two of the following three IFs be ELSE IFs?
C        Carol: Actually, it seems only XT(J) is needed here, not XT(1:KMX).

         DO L = 1, JMAX
            IF (IJPLANE .AND. JSTEP .OR.
     >          IKPLANE .AND. KSTEP) XT(L) = X(II,L,K)
            IF (IJPLANE .AND. ISTEP .OR.
     >          IKPLANE .AND. ISTEP) XT(L) = Y(II,L,K)
            IF (JKPLANE) XT(L) = Z(II,L,K)
         END DO
         DXMAX = XAT2 - XAT1
      END IF

C     Compute local DS or % based on current J location:

      DS0 = PST + (XT(J) - XAT1) * ((PEND - PST) / DXMAX)

      IF (XAT1 /= ZERO .OR. XAT2 /= ZERO) THEN
         IF (XT(J) < XAT1) DS0 = PST
         IF (XT(J) > XAT2) DS0 = PEND
      END IF

      IF (DS0 < MIN (PST, PEND)) DS0 = PST
      IF (DS0 > MAX (PST, PEND)) DS0 = PEND

      END SUBROUTINE DXSP

C***********************************************************************
C
      SUBROUTINE EDGEMG (VAR)
C
C     EDGEMG merges the end values of VAR into the next 'MG' mesh pts.
C     L1 is the start edge (often 1); L2 is the end edge (often NIPTS).
C     Carol:  L1 and L2 aren't used here.  What was meant?
C     If either is zero, then no merging is done at that end.
C
C     CALLED BY: LINE1, SOLUT
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

      REAL, INTENT (INOUT) :: VAR(IMX)

      INTEGER  I, L

C     Execution:

      IF (MG1 /= 0) THEN
         DO I = 2, MG1
            VAR(I) = (VAR(I)*REAL (I-1) + VAR(1)*REAL (MG1+1-I)) /
     >               REAL (MG1)
         END DO
      END IF

      IF (MG2 /= 0) THEN
         DO I = 2, MG2
            L = NIPTS - I
            VAR(L) = (VAR(L)*REAL (I-1) + VAR(NIPTS-1)*REAL (MG2+1-I)) /
     >               REAL (MG2)
         END DO
      END IF

      END SUBROUTINE EDGEMG

C***********************************************************************
C
      SUBROUTINE FBAR (J, K)
C
C     FBAR computes the flow field gradients (FQ), and the geometry
C     gradients (FG).  It normalizes (to FBAR), finds B by calling GETB,
C     computes S from X and Y, finds the mid-points of S, and
C     reproportions S for the initial guess within the data block.
C     Finally, it interpolates for X, Y, FB, etc., at the new S.
C
C     CALLED BY: Main program
C
C     CALLS: FILTER, GETB, INTF, MGWALLS, NORM, PROPS, WALLS
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: J, K

      REAL, PARAMETER :: GAM = 1.4

C     Local variables:

      INTEGER
     >   I, L, N

      REAL
     >   AVEDS, DSNMK, FDQAM, FDQP, FDQT, P1,
     >   AMACH(IMX), PRES(IMX), TRAT(IMX), FDQ(IMX)

C     Execution:

C     Calculate S array.  SS and SMS are evaluated at input points
C     and stay fixed for each line - used for interpolation (?).
C     SN is the updated S array; SNMK is the S array at K-1.

      SS(1)   = ZERO
      SNMK(1) = ZERO

      L = IST
      DO I = 1, NIPTS - 1

         DS(I) = SQRT ((X(L+1,J,K) - X(L,J,K))**2 +
     >                 (Y(L+1,J,K) - Y(L,J,K))**2 +
     >                 (Z(L+1,J,K) - Z(L,J,K))**2)

         IF (DS(I) == ZERO) THEN
            OK = .FALSE.
            WRITE (6, '(/, 4(A, I4))')
     >         ' *** Grid has identical points at', L, ' and', L+1,
     >         ' on line', J, ', plane', K,
     >         ' Try RECLUST option to remedy this.'
            GO TO 900
         END IF

         SS(I+1) = SS(I) + DS(I)

         IF (K > KST) THEN
            DSNMK = SQRT ((X(L+1,J,K-1) - X(L,J,K-1))**2 +
     >                    (Y(L+1,J,K-1) - Y(L,J,K-1))**2 +
     >                    (Z(L+1,J,K-1) - Z(L,J,K-1))**2)
            SNMK(I+1) = SNMK(I) + DSNMK
         ELSE
            SNMK(I+1) = ZERO
         END IF

         L = L + 1
      END DO

      DO I = 1, NIPTS - 1
         SMS(I) = (SS(I+1) + SS(I))*0.5
      END DO

      AVEDS = SS(NIPTS) / REAL (NIPTS-1)
      DSMAX = RDSMAX*AVEDS
      DSMIN = RDSMIN*AVEDS

C     Find DQ/DS:

      IF (.NOT. NOQ) THEN

         DO N = 1, NQDIM
            L = IST
            DO I = 1, NIPTS - 1
               FDQ(I) = (ABS(Q(L+1,J,K,N) - Q(L,J,K,N))) / DS(I)
               L = L + 1
            END DO

            CALL NORM (NIPTS-1, FDQ)

            DO I = 1, NIPTS - 1
               IF (N == 1) FQ(I) = ZERO
               FQ(I) = FQ(I) + IQ(N)*FDQ(I)
            END DO
         END DO

C        Evaluate pressure, Mach number or temp ratio if used:

         IF (IQ(NQDIM+1) /= 0 .OR. IQ(NQDIM+2) /= 0 .OR.
     >       IQ(NQDIM+3) /= 0) THEN
            L = IST
            DO I = 1, NIPTS 
               IF (TWOD) THEN
                  P1 = (Q(L,J,K,2)**2 + Q(L,J,K,3)**2) / Q(L,J,K,1)
                  PRES(L) = (GAM-1.0)*(Q(L,J,K,4) - 0.5*P1)
               ELSE 
                  P1 = (Q(L,J,K,2)**2 + Q(L,J,K,3)**2 + Q(L,J,K,4)**2) /
     >                  Q(L,J,K,1)
                  PRES(L) = (GAM-1.0)*(Q(L,J,K,5) - 0.5*P1)
               END IF
               AMACH(L) = SQRT (P1 / (GAM*PRES(L)))
               TRAT(L)  = GAM*PRES(L) / Q(L,J,K,1)
               L = L + 1
            END DO 

            L = IST
            DO I = 1, NIPTS - 1
               FDQP  = ABS (PRES(L+1)  - PRES(L))
               FDQAM = ABS (AMACH(L+1) - AMACH(L))
               FDQT  = ABS (TRAT(L+1)  - TRAT(L))
               FQ(I) = FQ(I) + (IQ(NQDIM+1)*FDQP + IQ(NQDIM+2)*FDQAM
     >                       +  IQ(NQDIM+3)*FDQT) / DS(I)
               L = L + 1
            END DO
         END IF

      END IF

C     Calculate the geometry gradients, DG/DS:

      IF (GEOM) THEN

         IF (J == JST .OR. .NOT.QFUN) CALL WALLS (J, K)

         IF (QFUN .AND. J == ((JEND-JST)/2+1)) CALL WALLS (JEND, K)

C        Find FGW, coefficient of geometry function:

         IF (.NOT. QFUN) THEN
            FQW = ZERO
            FGW = 1.0
         ELSE
            CALL MGWALLS (J)
            IF (FGW == ZERO) GO TO 460
         END IF

C        Interpolate FG (using FGS computed in MGWALLS)  
C        and normalize FG and FQ:

         CALL INTF (INTER, NIPTS-1, SMSS, SS, FGS, FG)
         CALL NORM (NIPTS-1, FG)
         CALL NORM (NIPTS-1, FQ)

      END IF

  460 CONTINUE

C     Compute F:

      DO I = 1, NIPTS-1
         F(I) = FQW*FQ(I) + FGW*FG(I)
      END DO

C     If required, filter F:

      IF (NFILT > 0) CALL FILTER (NIPTS-1, NFILT, F)

C     F is now finalized and evaluated at SS.
C     Store F in FB before normalizing.

      CALL INTF (INTER, NIPTS-1, SMS, SS, F, FB)
      CALL NORM (NIPTS-1, FB)

C     Calculate B using data evaluated at input nodes:

      CALL GETB (J, K)

C     Initialise SN array - if 1st time, = SS, else proportion from
C     converged solution at J-1 line and re-evaluate F, DS:

      IF (J == JST .AND. K == KST) THEN
         DO I = 1, NIPTS
            SN(I) = SS(I)   
         END DO
      ELSE
         CALL PROPS(J,K)
      END IF

      DO I = 1, NIPTS - 1
         DS(I) = SN(I+1) - SN(I)
      END DO

C     Interpolate for F at new S, and compute FB:

      CALL INTF (INTER, NIPTS-1, SMS, SN, F, FB)
      CALL NORM (NIPTS-1, FB)

  900 CONTINUE

      END SUBROUTINE FBAR

C***********************************************************************

      SUBROUTINE FILTER (NPOINTS, NFILT, VAR)
C
C     This subroutine filters (smooths) the given array.
C
C     CALLED BY: FBAR, LINE1, SOLUT
C
C***********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT (IN)    :: NPOINTS, NFILT
      REAL,    INTENT (INOUT) :: VAR(NPOINTS)

      INTEGER  I, L
      REAL     VT(NPOINTS)

      DO L = 1, NFILT
         DO I = 2, NPOINTS - 1
            VT(I) = 0.75*VAR(I) + 0.125*(VAR(I+1) + VAR(I-1))
         END DO
         DO I = 2, NPOINTS - 1
            VAR(I) = VT(I)
         END DO
      END DO

      END SUBROUTINE FILTER

C***********************************************************************
C
      SUBROUTINE FVORG (IND)
C
C     This routine rearranges the finite volume data. 
C     IND = 1 means transfer cell-centered Q to grid points;
C     IND = 2 means output cell-centered values.
C
C     CALLED BY: INITIAL, OUTPUT 
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: IND

C     Local variables:

      INTEGER
     >   I, J, K, KPL, N

      REAL
     >   Q1, Q2, QS(IMX)

C     Execution:

      IF (IND == 1) THEN

C        Transform Q (given at cell ceneters) to grid point locations;
C        at borders, use same points as one point in.

         KPL = KMAX - 2
         IF (KMAX == 1) KPL = 1

         DO 500 N = 1, NQDIM

            DO 100 K = 1, KPL

               DO J = 2, JMAX - 1

                  DO I = 2, IMAX - 1
                     Q1 = Q(I-1,J,K,N) + Q(I,J-1,K,N) +
     >                    Q(I,J,K,N) + Q(I-1,J-1,K,N)
                     IF (KPL == 1) THEN
                       QS(I) = 0.25*Q1
                     ELSE
                       Q2 = Q(I-1,J,K+1,N) + Q(I,J-1,K+1,N) +
     >                      Q(I,J,K+1,N) + Q(I-1,J-1,K+1,N)
                       QS(I) = 0.125*(Q1 + Q2)
                     END IF
                  END DO
                  DO I = 2, IMAX - 1
                     Q(I,J-1,K,N) = QS(I)
                  END DO

               END DO

               DO J = JMAX - 1, 2, -1
                  DO I = 2, IMAX - 1
                     Q(I,J,K,N) = Q(I,J-1,K,N)
                  END DO
               END DO
               DO I = 2, IMAX - 1
                  Q(I,JMAX,K,N) = Q(I,JMAX-1,K,N)
               END DO
               DO J = 1, JMAX
                  Q(1,J,K,N) = Q(2,J,K,N)
                  Q(IMAX,J,K,N) = Q(IMAX-1,J,K,N)
               END DO

  100       CONTINUE ! Next K

C           Shift first plane into 2nd, etc.:

            IF (KMAX > 1) THEN
               DO K = KMAX - 1, 2, -1
                  DO J = 1, JMAX
                     Q(1:IMAX,J,K,N) = Q(1:IMAX,J,K-1,N)
                  END DO
               END DO

C              Set KMAX-1 & KMAX data the same:

               DO J = 1, JMAX
                  Q(1:IMAX,J,KMAX,N) = Q(1:IMAX,J,KMAX-1,N)
               END DO
            END IF

  500    CONTINUE ! Next N

      ELSE ! IND == 2; transform Q back to center

         KPL = KMAX - 1
         IF (KMAX == 1) KPL = 1

         DO 700 N = 1, NQDIM

            DO K = 1, KPL

               DO I = 1, IMAX - 1

                  DO J = 1, JMAX - 1
                     Q1 = Q(I,J,K,N) + Q(I,J+1,K,N) + Q(I+1,J,K,N) +
     >                    Q(I+1,J+1,K,N)
                     IF (KMAX == 1 .OR. K == KEND) THEN
                        Q(I,J,1,N) = .25*Q1
                     ELSE
                        Q2 = Q(I,J,K+1,N)   + Q(I,J+1,K+1,N) +
     >                       Q(I+1,J,K+1,N) + Q(I+1,J+1,K+1,N)
                        Q(I,J,K,N) = 0.125*(Q1 + Q2)
                     END IF
                  END DO

               END DO ! Next I

            END DO ! Next K

C           Ensure that all 1,IMAX-1,JMAX-1 points have been mirrored:

            IF (KMAX > 1) THEN
               DO J = 2, JMAX - 2
                  DO I = 2, IMAX - 2
                     Q(I,J,1,N) = Q(I,J,2,N)
                     Q(I,J,KMAX-1,N) = Q(I,J,KMAX-2,N)
                  END DO
               END DO
               DO J = 2, JMAX - 2
                  DO K = 2, KMAX - 2
                     Q(1,J,K,N) = Q(2,J,K,N)
                     Q(IMAX-1,J,K,N) = Q(IMAX-2,J,K,N)
                  END DO
               END DO
               DO I = 2, IMAX - 2
                  DO K = 2, KMAX - 2
                     Q(I,1,K,N) = Q(I,2,K,N)
                     Q(I,JMAX-1,K,N) = Q(I,JMAX-2,K,N)
                  END DO
               END DO
            END IF

C           Now reflect outer edges:

            DO I = 2, IMAX - 2
               Q(I,1,1,N) = Q(I,2,1,N)
               Q(I,JMAX-2,1,N) = Q(I,JMAX-1,1,N)
               IF (KMAX > 1) THEN
                  Q(I,1,KMAX-1,N) = Q(I,2,KMAX-1,N)
                  Q(I,JMAX-2,KMAX-1,N) = Q(I,JMAX-1,KMAX-1,N)
               END IF
            END DO
            DO J = 1, JMAX - 1
               Q(1,J,1,N) = Q(2,J,1,N)
               Q(IMAX-1,J,1,N) = Q(IMAX-2,J,1,N)
               IF (KMAX > 1) THEN
                  Q(1,J,KMAX-1,N) = Q(2,J,KMAX-1,N)
                  Q(IMAX-1,J,KMAX-1,N) = Q(IMAX-2,J,KMAX-1,N)
               END IF
            END DO
            IF (KMAX > 1) THEN
               DO K = 2, KMAX - 2
                  Q(1,1,K,N) = Q(2,1,K,N)
                  Q(1,JMAX-1,K,N) = Q(2,JMAX-1,K,N)
                  Q(IMAX-1,1,K,N) = Q(IMAX-2,1,K,N)
                  Q(IMAX-1,JMAX-1,K,N) = Q(IMAX-2,JMAX-1,K,N)
               END DO
            END IF

  700    CONTINUE ! Next N

         WRITE (6, '(/, 1X, A)')
     >      'Solution file is finite volume; check your boundary cells!'

      END IF

      END SUBROUTINE FVORG

C***********************************************************************
C
      SUBROUTINE GETB (J, K)
C
C     GETB calculates B iteratively, using the initial grid spacing.
C     B converges when input DSmin = computed DSmin.
C     B controls the minimum allowed S mesh spacing.
C
C     CALLED BY: FBAR
C
C     CALLS: INTF
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: J, K

C     Local variables:

      INTEGER
     >   I, L, ITER, MAXBITS

      REAL
     >   BCONV, BJ1, DB, DMINSDB, DSIMIN, SUM, TSTCONV, WTSUM,
     >   ALNF(IMX), FINT(IMX), WIT(IMX), SNB(IMX), DSB(IMX)

C     Execution:

C     If first line, set initial guess of B = 1.0 and permit more
C     iterations for convergence, otherwise first guess on B is the
C     current value from the J-1 line (or if JST, then K-1 line).

      MAXBITS = MAXITS
      IF (J == JST .OR. (J == JST+1 .AND. MGSTEPS /= 0)) THEN
         IF (K == KST .OR. (K == KST+1 .AND. MGPLS /= 0)) THEN
            B = 1.0
            MAXBITS = MAXITS + 20
         ELSE
            B = BK1
         END IF
      END IF
      IF (J == JST+1) BK1 = B
      BJ1 = B
      IF (A == ZERO) GO TO 999

      DB = ZERO
      SNB(1) = ZERO

C     Convergence of DSMIN must be better than requested DSMIN,
C     hence BCONV is a function of DSMIN.

      BCONV = DSMIN*0.02

C     Iterate to find B:

      DO 500 ITER = 1, MAXBITS

C        Compute weight based on current B using initial grid spacing:

         DO L = 1, NIPTS - 1
            WEIGHT(L) = 1.0 + A*FB(L)**B
         END DO

C        Find new DS, using this B:

         WTSUM = ZERO
         DO I = 1, NIPTS - 1
            WTSUM = WTSUM + 1.0/WEIGHT(I)
         END DO
         DO I = 1, NIPTS - 1
            DSB(I) = SS(NIPTS)/(WEIGHT(I)*WTSUM)
            SNB(I+1) = SNB(I) + DSB(I)
         END DO

C        Find min. value of DS:

         DSIMIN = DSB(1)
         DO L = 2, NIPTS -1
            IF (DSB(L) < DSIMIN) DSIMIN = DSB(L)
         END DO

C        Test for convergence:

         TSTCONV = ABS (DSMIN - DSIMIN)
         IF (TSTCONV <= BCONV) GO TO 999       

C        No convergence; compute DB for next iteration.
C        (1) Find new FB and W at this new spacing:

         CALL INTF (INTER, NIPTS-1, SMS, SNB, WEIGHT, WIT)
         CALL INTF (INTER, NIPTS-1, SMS, SNB, FB, FINT)

C        (2) Find derivative of DXIMIN w.r.t. B:

         SUM = ZERO
         DO L = 1, NIPTS - 1
            ALNF(L) = LOG (FINT(L))
            SUM = SUM + FINT(L)**B * ALNF(L) / WIT(L)**2
         END DO
         DMINSDB = A*(1.0 + A)*DSIMIN**2 * SUM / SS(NIPTS)
         IF (DMINSDB == ZERO) GO TO 998

         DB = (DSMIN - DSIMIN) / DMINSDB
         B = B + DB

C        Maintain reasonable limit for DB.  First iteration may be extreme.

         IF (B > 5.0)  B = 3.0
         IF (B < ZERO) B = 0.1

  500 CONTINUE ! Next iteration

C     No convergence on B: use old value.

  998 CONTINUE

      IF (J == JST) THEN
         B = 1.0
      ELSE
         B = BJ1
      END IF

  999 CONTINUE

      END SUBROUTINE GETB

C***********************************************************************
C
      SUBROUTINE GETDS0 (DS0, J, K)
C
C     GETDS0 computes the wall spacing for the Vinokur algorithm when
C     DSW > 0 (i.e., the original first spacing is not used).
C
C     CALLED BY: BLCLUST, PUSHIT
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

      INTEGER, INTENT (IN)  :: J, K
      REAL,    INTENT (OUT) :: DS0

C     Local variables:

      REAL, DIMENSION (IMX) ::
     >   SNW, DSC, SC, SW, XM, YM, ZM

      REAL
     >   REYC, QM(IMX,NDIM)

C     Execution:

      IF (DSW  > ZERO .AND. DSW  <  997.) DS0 = DSW

      IF (DSW == 999.) CALL DXSP (DS0, 1, DSW1, DSW2, XMIN, XMAX, J, K)

      IF (DSW == 997. .AND. .NOT. NOQ) THEN
         REYC = 2.0
         IF (DSW1 >= 1.0 .AND. DSW1 < 20.) REYC = DSW1
         DS0 = REYC*Q(IST,J,K,3) / (Q(IST,J,K,1)*Q(IST,J,K,2))
         IF (DSW2 /= ZERO .AND. DS0 > DSW2) DS0 = DSW2
      END IF

      IF (DSW == 997. .AND. NOQ) THEN
         WRITE (6, '(/, A)')
     >      ' *** Need Q file for Rey wall spacing *** '
         STOP
      END IF

      END SUBROUTINE GETDS0

C***********************************************************************
C
      SUBROUTINE GETWT
C
C     GETWT determines the multiplier of the weight constants, WT.
C     These are an additional control on the max. and min. allowable DSs.
C     AM is a variable that prevents a "flip-flop" condition.
C
C     CALLED BY: LINE1, SOLUT
C
C***********************************************************************

      USE SAGEMOD   

      IMPLICIT NONE

      REAL, PARAMETER :: AM = 0.5
      INTEGER I, I1, I2
      REAL    DT

      WT(1:NIPTS-1) = 1.0

C     If NEDGE is set, permit DS to be too small.

      I1 = 1
      IF (MG1 /= 0) I1 = MG1
      I2 = NIPTS - 1
      IF (MG2 /= 0) I2 = NIPTS - MG2

      DO I = 1, NIPTS - 1
         IF (DS(I) > DSMAX) THEN
            DT = DS(I) / DSMAX
            WT(I) = (DT - 1.0)*AM + 1.0
         END IF
      END DO
      DO I = I1, I2
         IF (DS(I) < DSMIN) THEN
            DT = DS(I) / DSMIN
            WT(I) = (DT - 1.0)*AM + 1.0
         END IF
      END DO

      END SUBROUTINE GETWT

C***********************************************************************
C
      SUBROUTINE HEADIO (IC, IER)
C
C     HEADIO reads and writes headers for grids. 
C     IC = 1 means READ ONLY; IC = 2 means WRITE ONLY.
C
C     CALLED BY: MATCH, READMULT
C
C     CALLS: ALOALO
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

      INTEGER, INTENT (IN)  :: IC
      INTEGER, INTENT (OUT) :: IER

      INTEGER  IERB, N, NGRIDQ

      IER = 0

      IF (IC == 1) THEN ! Read the number of zones (aka MGRID):

         IF (SINGRD) THEN 
            NGRID = 1
         ELSE
            IF (FORMI) THEN
               READ (NITGI, *) NGRID
               IF (.NOT. NOQ) READ (NITQI, *) NGRIDQ
            ELSE
               READ(NITGI) NGRID
               IF (.NOT. NOQ) READ (NITQI) NGRIDQ
            END IF
            IF (.NOT. NOQ) THEN
               IF (NGRID /= NGRIDQ) THEN
                  WRITE (6, '(/, A, 2I6)')
     >               ' Number of zones in grid & Q file do no match: ',
     >               NGRID, NGRIDQ
                  IER = 1
               END IF
            END IF
         END IF

C        If first time through, allocate zone dimensions:

         IF (HIOFLAG == 1) THEN
            ALLOCATE (IM(NGRID),  JM(NGRID),  KM(NGRID),
     >                IMO(NGRID), JMO(NGRID), KMO(NGRID),
     >                IMQO(NGRID),JMQO(NGRID),KMQO(NGRID),
     >                IMQ(NGRID), JMQ(NGRID), KMQ(NGRID),
     >                NQFUN(NGRID), MGCT(NGRID), STAT=IERB)

            MGCT = 0
         END IF

C        Read zone dimensions and, if Function file, # of functions.
C        Read grid file and Q/Fun file if it exists.

         IF (TWOD) THEN

            IF (FORMI) THEN
               READ (NITGI, *) (IM(N), JM(N), N = 1, NGRID)
               IF (.NOT. NOQ) THEN
                  IF (FUN) THEN
                     READ (NITQI, *)
     >                  (IMQ(N), JMQ(N), NQFUN(N), N = 1, NGRID)
                  ELSE
                     READ (NITQI, *) (IMQ(N), JMQ(N), N = 1, NGRID)
                  END IF
               END IF
            ELSE
               READ (NITGI) (IM(N), JM(N), N = 1, NGRID)
               IF (.NOT. NOQ) THEN
                  IF (FUN) THEN
                     READ (NITQI)
     >                  (IMQ(N), JMQ(N), NQFUN(N), N = 1, NGRID)
                  ELSE
                     READ (NITQI) (IMQ(N), JMQ(N), N = 1, NGRID)
                  END IF
               END IF
            END IF
            
            KM = 1 ! 1:NGRID

         ELSE ! 3-D

            IF (FORMI) THEN
               READ (NITGI, *) (IM(N), JM(N), KM(N), N = 1, NGRID)
               IF (.NOT. NOQ) THEN
                  IF (FUN) THEN
                     READ (NITQI, *)
     >                  (IMQ(N), JMQ(N), KMQ(N), NQFUN(N), N = 1, NGRID)
                  ELSE
                     READ (NITQI, *)
     >                  (IMQ(N), JMQ(N), KMQ(N), N = 1, NGRID)
                  END IF
               END IF
            ELSE
               READ (NITGI) (IM(N), JM(N), KM(N), N = 1, NGRID)
               IF (.NOT. NOQ) THEN
                  IF (FUN) THEN
                     READ (NITQI)
     >                  (IMQ(N), JMQ(N), KMQ(N), NQFUN(N), N = 1, NGRID)
                  ELSE
                     READ (NITQI) (IMQ(N), JMQ(N), KMQ(N), N = 1, NGRID)
                  END IF
               END IF
            END IF

         END IF

C        Zone dimensions read: now allocate zone arrays.

         CALL ALOALO

         IF (HIOFLAG == 1) THEN
            IF (FUN) NQDIM = NQFUN(1)
            HIOFLAG = 0
         END IF

      ELSE ! IC == 2: write header records

         IF (.NOT. SINGRD) THEN
            IF (FORMO) THEN
               WRITE (NITGO, '(I4)') NGRID
            ELSE
               WRITE (NITGO) NGRID
            END IF
         END IF

         IF (TWOD) THEN
            IF (FORMO) THEN
               WRITE (NITGO, '(I5, I6)') (IMO(N), JMO(N), N = 1, NGRID)
            ELSE
               WRITE (NITGO) (IMO(N), JMO(N), N = 1, NGRID)
            END IF
         ELSE
            IF (FORMO) THEN
               WRITE (NITGO, '(I5, 2I6)')
     >            (IMO(N), JMO(N), KMO(N), N = 1, NGRID)
            ELSE
               WRITE (NITGO) (IMO(N), JMO(N), KMO(N), N = 1, NGRID)
            END IF
         END IF

         IF (.NOT. NOQ) THEN
            IF (.NOT. SINGRD) THEN
               IF (FORMO) THEN
                  WRITE (NITQO, '(I4)') NGRID
               ELSE
                  WRITE (NITQO) NGRID
               END IF
            END IF
            IF (TWOD) THEN
               IF (FORMO) THEN
                  IF (FUN) THEN
                     WRITE (NITQO, '(I5, 2I6)')
     >                  (IMQO(N), JMQO(N), NQFUN(N), N = 1, NGRID)
                  ELSE
                     WRITE (NITQO, '(I5, I6)')
     >                  (IMQO(N), JMQO(N), N = N, NGRID)
                  END IF
               ELSE
                  IF (FUN) THEN
                     WRITE (NITQO)
     >                  (IMQO(N), JMQO(N), NQFUN(N), N = 1, NGRID)
                  ELSE
                     WRITE (NITQO) (IMQO(N), JMQO(N), N = 1, NGRID)
                  END IF
               END IF
            ELSE
               IF (FORMO) THEN
                  IF (FUN) THEN
                     WRITE (NITQO, '(I5, 3I6)')
     >               (IMQO(N), JMQO(N), KMQO(N), NQFUN(N), N = 1, NGRID)
                  ELSE
                     WRITE (NITQO, '(I5, 2I6)')
     >                  (IMQO(N), JMQO(N), KMQO(N), N = 1, NGRID)
                  END IF
               ELSE
                  IF (FUN) THEN
                     WRITE (NITQO)
     >               (IMQO(N), JMQO(N), KMQO(N), NQFUN(N), N = 1, NGRID)
                  ELSE
                     WRITE (NITQO)
     >                  (IMQO(N), JMQO(N), KMQO(N), N = 1, NGRID)
                  END IF
               END IF
            END IF
         END IF

      END IF

      END SUBROUTINE HEADIO

C***********************************************************************
C
      SUBROUTINE INITIAL (NOMORE, MTCH)
C
C     Read the control variables, initialize defaults, and swap data in
C     X, Y arrays if necessary.
C
C     CALLED BY: Main program
C
C     CALLS: ADDPTS, BLCLUST, FVORG, MATCH, PUSHIT,
C            READMULT, SWAP2D, SWAPINV, SWAPXYZ, SUBPTS
C
C***********************************************************************
C                                                                      I
C IST,IEND = FIRST AND LAST ADAPTION POINTS IN I DIRECTION             I
C JST,JEND =  """""""""""""""""""""""""""""""" J """""""""             I
C KST,KEND = """"""""""""""""""""""""""""""""" K """""""""             I
C ISTEP,JSTEP,KSTEP = TRUE FOR  STEPPING DIRECTION                     I
C IJPLANE,IKPLANE,JKPLANE = TRUE FOR MARCHING PLANE                    I
C RDSMIN,RDSMAX = MIN AND MAX ALLOWABLE GRID SPACINGS                  I
C CLAM(2) = MAGNITUDE OF LINE TORSION AND PLANE TORSION                I
C NEDGE = EDGE SPACING CONTROL,1=BOTH,2=ST,3=END                       I
C CT(2) = PROP. OF STRAIGHTNESS TO ORTHOGONALITY,1=STRAIGHT,0=NORMAL   I
C MGSTEPS = NO OF STEPS BEFORE FULL ADAPTION REQUIRED                  I
C MGPLS = NO OF PLANES BEFORE FULL ADAPTION OCCURS                     I
C GEOM = TRUE, INCLUDE GEOMETRY GRADIENTS CLOSE TO WALL                I
C QFUN = FALSE-USED ONLY WHEN GEOM=TRUE AND ONLY ADAPT TO GEOMETRY     I
C INDQ = INDEX ON FLOWFIELD VARIABLE                                   I
C IQ(NDIM+3) = IF INDQ=0, PROPORTION OF FLOWFIELD VARIABLES            I
C NOUP = TRUE IF NO ADAPTION REQUIRED (EG IF ADDING PTS ONLY)          I
C NOQ = FALSE: SET TO TRUE IF Q FILE DOES NOT EXIST                    I
C MARCH = TRUE FOR EXTRAPOLATING LAST ADAPTED LINE (IF JST LT JMAX)    I
C MARCHPL = TRUE FOR " LAST ADAPTED PLANE TO REMAINING PLANES          I
C MGRID = INDEX OF MULTIPLE GRID TO ADAPT                              I
C IMPORT/EXPORT = true TO PERFORM MULTIGRID BOUNDARY MATCHING          I
C MPLANE = PLANE NO FOR EXPORT/IMPORT MATCHING                         I
C IS,IE,JS,JE,KS,KE = DOMAINS FOR EXPORT/IMPORT PLANES                 I
C ADD = INTEGER FOR ADDING POINTS IN ADAPTION DIRECTION                I
C LSTADD,LENDADD = RANGE WITHIN TO ADD PTS                             I
C SUB = INTEGER FOR DELETING POINTS                                    I
C LSTSUB,LENDSUB = RANGE FOR DELETING POINTS			       I
C REMOVE = NO OF PTS TO REMOVE FROM OUTER BOUNDARY OF ADAPTION LINE    I
C*** FOR INNER BOUNDARY REMOVAL, EXCHANGE ST AND END PARAMETERS        I
C SAVE = FALSE TO SUPPRESS O/P OF ADAPTED FILES for single grids only  I
C ORTHS,ORTHE = FALSE TO SUPPRESS ORTHOGONALITY AT MARCHING BOUNDARIES I
C NFILT = INDEX ON AMOUNT OF FILTERING OF RAW DATA                     I
C INTER = ORDER OF INTERPOLATION                                       I
C LNSING = N, SAME LINE IN EACH PLANE: ADAPT ONCE ONLY                 I
C PLSING = N, PLANE COLLAPSED TO A LINE                                I
C MG1,MG2 = CHANGES NO. OF POINTS IN NEDGE CONTROL                     I
C TWOD = .true. FOR 2-D DATASETS				       I
C FV = .true. FOR FINITE VOLUME Q file                                 I
C MVBOUND = MOVE OUTER BOUNDARY:% or if > 998,fits to data             I
C MVC1,MVC2 = constant added to MVBOUND result                         I
C XMINP,XMAXP,MVP1,MVP2 relates to % outer boundary move               I
C DSN: USED WITH MVBOUND. RELATES TO # OF CELLS OUTSIDE OF SHOCK       I
C NSM,NSM1 = NUMBER OF PASSES TO SMOOTH 'MOVED' OUTER BOUNDARY         I
C NOSESM to average kink at nose with boundary move                    I
C RECLUST  =  1 TO RECLUSTER POINTS WRT MARCEL V. ALGORITHM            I
C DSW = SPACING AT WALL FOR RECLUSTERING ALGORITHM                     I
C      = 999 or 997, enter DSW1 and DSW2                               I
C     and maybe XMIN AND XMAX                                          I
C DSE = MULTIPLE OF DS TO COMPUTE OUTER SPACING FOR M-V ALGORITHM      I
C FORMI = TRUE IF INPUT FILES ARE FORMATTED. UNFORMATTED IS ASSUMED    I
C FORMO = TRUE IF OUTPUT FILES ARE TO BE FORMATTED                     I
C FUN = TRUE IF Q FILE IS A FUNCTION FILE                              I
C                                                                      I
CIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

      USE SAGEMOD

      IMPLICIT NONE

      INTEGER  IER, IMAXQ, L, NSM, NSM1
      REAL     DSN
      LOGICAL  MTCH, NOMORE, RSWAP

      NAMELIST /NAMEL/
     >   IST, IEND, JST, JEND, KST, KEND, RDSMAX, RDSMIN, ISTEP, JSTEP,
     >   KSTEP, IJPLANE, IKPLANE, JKPLANE, MGRID, XMIN, XMAX, CLAM,
     >   NEDGE, CT, MGSTEPS, MGPLS, INDQ, IQ, NOUP, TWOD, REMOVE, NQDIM,
     >   MARCH, MARCHPL, ADD, LSTADD, LENDADD, SUB, LSTSUB, LENDSUB,
     >   SAVE, ORTHS, ORTHE, NFILT, INTER, GEOM, QFUN, NOQ, LNSING,
     >   PLSING, MG1, MG2, MPLANE, EXPORT, IMPORT, IS, IE, JS, JE, KS,
     >   KE, FV, FORMI, FORMO, FUN, MVBOUND, XMINP, XMAXP, MVP1, MVP2,
     >   DSN, NSM, NSM1, IMS, JME, NOSESM, RECLUST, DSW, DSW1, DSW2,
     >   DSE, MVC1, MVC2

C     Set input default values:

      IST      = 1
      IEND     = 0
      JST      = 1
      JEND     = 0
      KST      = 1
      KEND     = 0
      ISTEP    = .FALSE.
      JSTEP    = .TRUE.
      KSTEP    = .FALSE.
      IJPLANE  = .TRUE.
      JKPLANE  = .FALSE.
      IKPLANE  = .FALSE.
      RDSMAX   = 2.0
      RDSMIN   = .5
      CLAM(1)  = .01
      CLAM(2)  = .0001
      NEDGE    = 0
      CT(1)    = .5
      CT(2)    = .5
      MGSTEPS  = 0
      MGRID    = 0
      MGPLS    = 0
      LNSING   = 0
      PLSING   = 0
      INDQ     = 1
      NOUP     = .FALSE.
      MTCH     = .FALSE.
      MARCH    = .FALSE.
      MARCHPL  = .FALSE.
      ADD      = 0
      LSTADD   = 0
      LENDADD  = 0
      SUB      = 0
      LSTSUB   = 0
      LENDSUB  = 0
      SAVE     = .TRUE.
      ORTHS(1) = .TRUE.
      ORTHS(2) = .TRUE.
      ORTHE(1) = .TRUE.
      ORTHE(2) = .TRUE.
      NFILT    = 2
      GEOM     = .FALSE.
      QFUN     = .TRUE.
      INTER    = 2
      MPLANE   = 1
      EXPORT   = .FALSE.
      IMPORT   = .FALSE.
      IS       = 0
      IE       = 0
      JS       = 0
      JE       = 0
      KS       = 0
      KE       = 0
      MVBOUND  = ZERO
      RECLUST  = 0
      IMS      = 0
      JME      = 0
      DSW      = ZERO
      DSW1     = ZERO
      DSW2     = ZERO
      MVP1     = ZERO
      MVP2     = ZERO
      MVC1     = ZERO
      MVC2     = ZERO
      MVIO     = .FALSE.
      DSE      = 5.0
      DSN      = 10.0
      NSM      = 10
      NSM1     = 10
      NOSESM   = 0
      FORMI    = .FALSE.
      FORMO    = .FALSE.
      XMIN     = ZERO
      XMAX     = ZERO
      XMINP    = ZERO
      XMAXP    = ZERO

C     Set other start values:

      OK       = .TRUE.
      FQW      = 1.0
      FGW      = ZERO
      NOMORE   = .FALSE.
      IINVERSE = .FALSE.
      JINVERSE = .FALSE.
      KINVERSE = .FALSE.
      TN(1)    = .5
      TN(2)    = .5
      CONV     = .001
      MG1      = 0
      MG2      = 0
      MAXITS   = 20
      RSWAP    = .FALSE.
      NFLAG    = 1
      WDS      = ZERO
      WDE      = ZERO
      REMOVE   = 0

      IF (NADS == 1) THEN
         TWOD    = .FALSE.
         FV      = .FALSE.
         NOQ     = .FALSE.
         NITGI   = 7
         NITQI   = 8
         NITGO   = 0
         NITQO   = 0
         FUN     = .FALSE.
         FORMOO  = .FALSE.
         NQDIM   = NDIM
         HIOFLAG = 1
      END IF

      IQ(1 : NDIM + 3) = 0
      IERA(1 : 9) = 0

C     Read user control file:

      READ (5, NAMEL, END=230)

      IF (NOQ) NQDIM = 1
      IF (IKPLANE .OR. JKPLANE) IJPLANE = .FALSE.
      IF (MGRID == 0) THEN
         SINGRD = .TRUE.
         MGRID = 1
      ELSE
         SINGRD = .FALSE.
      END IF

      IF (EXPORT) THEN
         CALL MATCH (NOMORE, IER)
         IF (IER == 1) GO TO 230
         MTCH = .TRUE.
         GO TO 250
      END IF

      IF (TWOD .AND. .NOT. JSTEP) ISTEP = .TRUE.
      IF (ISTEP .OR. KSTEP) JSTEP = .FALSE.
      IF (IJPLANE .AND. KSTEP .OR. IKPLANE .AND. JSTEP .OR.
     >    JKPLANE .AND. ISTEP) THEN
         WRITE (6, '(/, A)') ' *** INCONSISTENT PLANE and STEP ***'
         GO TO 230
      END IF
      IF (NQDIM == 5 .AND. TWOD .AND. .NOT. FUN) NQDIM = 4
      IF (DSW == 997 .AND. NOQ) THEN
         WRITE (6, '(/, A)') ' *** Need Q file for Rey wall spacing ***'
         GO TO 230
      END IF

C     Read grid and function files:

      CALL READMULT (IER)
      IF (IER == 1) GO TO 230

      MGCT(MGRID) = MGCT(MGRID) + 1
      SAVE = .TRUE.

      IF (.NOT. NOQ) THEN
         IMAXQ = IM(MGRID)
         IF (IMAX /= IMAXQ .AND. (IMAX /= IMAXQ+1 .OR. .NOT. FV)) THEN
            WRITE (6, '(/, A)') ' ** INPUT FILE SIZES DO NOT MATCH ** '
            IF (IMAX == IMAXQ + 1) THEN
               WRITE (6, '(/, A)')
     >           ' ** Is this finite volume? If so, set FV = .TRUE. **'
            END IF
            GO TO 230
         END IF
      END IF

      IF (IER  == 1) GO TO 230
      IF (IEND == 0) IEND = IMAX
      IF (JEND == 0) JEND = JMAX
      IF (KEND == 0) KEND = KMAX
      A = RDSMAX/RDSMIN - 1.0

C     Initialize to adaption variable requested:

      IF (INDQ /= 0) THEN
         DO L = 1, NDIM + 3
            IF (L == INDQ) THEN
               IQ(L) = 1
            ELSE
               IQ(L) = 0
            END IF
         END DO
      END IF

C     Set up indices for edge treatment, based on NEDGE:

      IF ((NEDGE == 1 .OR. NEDGE == 2) .AND. MG1 == 0) MG1 = 4
      IF ((NEDGE == 1 .OR. NEDGE == 3) .AND. MG2 == 0) MG2 = 4

C     Swap coordinates if necessary.
C     (The code always adapts 'I' stepping in 'J' lines and 'K' planes.)

      IF (TWOD) CALL SWAP2D

      IF (FV) THEN
         IF (KMAX == 2 .OR. KMAX == 3) THEN
            WRITE (6, '(/, A)')
     >         ' ** FINITE VOLUME method needs 1 or 4+ planes. **'
            GO TO 230
         ELSE
            IF (.NOT. NOQ) CALL FVORG (1)
         END IF
      END IF
      IF (.NOT. (IJPLANE .AND. JSTEP)) CALL SWAPXYZ (RSWAP)

C     IF ADD or SUB, change number of points in the adapted direction:

      IF (ADD /= 0) CALL ADDPTS
      IF (SUB /= 0) CALL SUBPTS

C     If IST > IEND or JST > JEND, swap data:
C     (The analysis assumes i and j indices always increase.)

      IF (IST  > IEND) IINVERSE = .TRUE. 
      IF (JST  > JEND) JINVERSE = .TRUE.
      IF (KST  > KEND) KINVERSE = .TRUE.
      IF (JST == JEND .AND. JST == JMAX) JINVERSE = .TRUE.

      IF (IINVERSE .OR. JINVERSE .OR. KINVERSE) CALL SWAPINV

C     Remove any points from outer boundary:

      IF (REMOVE /= 0) THEN
         IMAX = IMAX - REMOVE
         WRITE (6, '(/, I4, A, I4)')
     >      REMOVE, ' point(s) removed; new number of points:', IMAX
      END IF
      IF (IEND > IMAX) IEND = IMAX

C     Move outer boundary if requested:

      IF (MVBOUND /= ZERO .OR. MVC1 /= ZERO .OR. MVC2 /= ZERO)
     >    MVIO = .TRUE.
      IF (IMS == 0) IMS = JMAX
      IF (JME == 0) JME = KMAX
      IF (RECLUST == 1 .OR. MVIO) NOUP = .TRUE.

      IF (MVIO) CALL PUSHIT(DSN,NSM,NSM1)

      IF (.NOT. MVIO .AND. RECLUST == 1) CALL BLCLUST

C     Compute the number of points in the adaptive section:

      NIPTS = IEND - IST + 1
      IF (NIPTS <= 9 .AND. NEDGE == 1) THEN
         WRITE (6, '(/, A)')
     >      ' ** TOO FEW POINTS for adaption with NEDGE = 1 **'
         GO TO 230
      END IF

C     Correct the merge numbers to be the actual values:

      IF (MGSTEPS /= 0) MGSTEPS = MGSTEPS + JST
      IF (MGPLS   /= 0) MGPLS   = MGPLS   + KST
      GO TO 250

  230 NOMORE = .TRUE.
  250 CONTINUE

      END SUBROUTINE INITIAL

C***********************************************************************
C
      SUBROUTINE INTF (INTER, NPTS, SMID, S1, F1, F2)
C
C     F1 is a variable associated with the mid-point (SMID) of the  grid
C     element (e.g., derivatives).  INTF computes the mid-points of the
C     new S array (S1) and interpolates for F1, returning results in F2.
C     F2 is forced to be at least EPS > 0.
C
C     CALLED BY: FBAR, GETB, LINE1, SOLUT
C
C     CALLS: CSFIT, CSEVAL, LAGCOF
C
C***********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT (IN)  ::
     >   INTER,              ! INTER = 4 means spline interpolation
     >   NPTS                !           else parabolic
      REAL,    INTENT (IN)  ::
     >   SMID(NPTS),         ! Data abscissas
     >   S1(NPTS+1),         ! Target abscissas
     >   F1(NPTS)            ! Data ordinates
      REAL,    INTENT (OUT) ::
     >   F2(NPTS)            ! Interpolated values of F1 at S1
                             ! (confusing nomenclature!)
C     Local constants:

      REAL, PARAMETER ::
     >    EPS = 1.E-05,   ! Explain why F2 must be > 0
     >    HALF = 0.5,
     >    ZERO = 0.0
 
C     Local variables:

      INTEGER
     >    I, IER, M
      REAL
     >    P1, P2, P3
      REAL, DIMENSION (NPTS) ::
     >    BCOF, CCOF, DCOF, SM

C     Execution:

      DO I = 1, NPTS
         SM(I) = (S1(I+1) + S1(I)) * HALF
      END DO

      IF (INTER == 4) THEN

         CALL CSFIT (NPTS, SMID, F1, 0, ZERO, 0, ZERO,
     >               BCOF, CCOF, DCOF, IER)

         IF (IER /= 0) THEN ! Highly unlikely
            WRITE (6, '(/, A, I2)') ' INTF: IER =', IER
            STOP
         END IF

         CALL CSEVAL (NPTS, SMID, F1, NPTS, SM, BCOF, CCOF, DCOF, F2)

      ELSE

         DO I = 1, NPTS

            CALL LAGCOF (SM(I), SMID, NPTS, M, P1, P2, P3)

            F2(I) = P1*F1(M) + P2*F1(M+1) + P3*F1(M+2)

         END DO

      END IF

      DO I = 1, NPTS ! Previously, F2 could be positive but less than EPS.
         IF (F2(I) < EPS) F2(I) = EPS
      END DO

      END SUBROUTINE INTF

C***********************************************************************
C
      SUBROUTINE INTXYZQ (J, K, J1, K1, SS9, SN9, QJ)
C
C     Given X, Y, Z, Q at SS, find the same at SN.
C
C     Carol:  Why is QJ an argument but not XJ, YJ, ZJ?
C
C     CALLED BY: MARCHJ, MARCHK, PROPS, UPDATE
C
C     CALLS: CSFIT, CSEVAL, LAGCOF
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   J,  K,
     >   J1, K1             ! Explain

      REAL, INTENT (IN) ::
     >   SS9(IMX),          ! Explain
     >   SN9(IMX)           !  "

      REAL, INTENT (OUT) ::
     >   QJ(IMX,NDIM)       !  "

C     Local variables:

      INTEGER
     >   I, IER, M, MS, N, NEVAL

      REAL
     >   P1, P2, P3, XJ1, XJN, YJ1, YJN, ZJ1, ZJN

      REAL, DIMENSION (IMX) ::
     >   BCOF, CCOF, DCOF

C     Execution:

      IF (INTER == 4) THEN ! Spline interpolation

         NEVAL = NIPTS - 2 ! Explain why boundary (x,y,z)s aren't needed

         CALL CSFIT  (NIPTS, SS9, X(IST,J,K), 0, ZERO, 0, ZERO,
     >                BCOF, CCOF, DCOF, IER)
         CALL CSEVAL (NIPTS, SS9, X(IST,J,K), NEVAL, SN9(2),
     >                BCOF, CCOF, DCOF, XJ(2,J1,K1))

         CALL CSFIT  (NIPTS, SS9, Y(IST,J,K), 0, ZERO, 0, ZERO,
     >                BCOF, CCOF, DCOF, IER)
         CALL CSEVAL (NIPTS, SS9, Y(IST,J,K), NEVAL, SN9(2),
     >                BCOF, CCOF, DCOF, YJ(2,J1,K1))

         CALL CSFIT  (NIPTS, SS9, Z(IST,J,K), 0, ZERO, 0, ZERO,
     >                BCOF, CCOF, DCOF, IER)
         CALL CSEVAL (NIPTS, SS9, Z(IST,J,K), NEVAL, SN9(2),
     >                BCOF, CCOF, DCOF, ZJ(2,J1,K1))

         IF (.NOT. NOQ) THEN ! Carol: Is IST /= 1 handled at the higher level?

            DO N = 1, NQDIM
               CALL CSFIT  (NIPTS, SS9, Q(IST,J,K,N), 0, ZERO, 0, ZERO,
     >                      BCOF, CCOF, DCOF, IER)
               CALL CSEVAL (NIPTS, SS9, Q(IST,J,K,N), NIPTS, SN9,
     >                      BCOF, CCOF, DCOF, QJ(1,N))
            END DO

         END IF

      ELSE ! Lagrange polynomial interpolation

C        Carol: Why are the boundary (x,y,z)s preserved but not the Qs?

         XJ1 = XJ(1,J1,K1);  XJN = XJ(NIPTS,J1,K1)
         YJ1 = YJ(1,J1,K1);  YJN = YJ(NIPTS,J1,K1)
         ZJ1 = ZJ(1,J1,K1);  ZJN = ZJ(NIPTS,J1,K1)

         DO I = 1, NIPTS

            CALL LAGCOF (SN9(I), SS9, NIPTS, MS, P1, P2, P3)
 
            M = MS + IST - 1
            XJ(I,J1,K1) = P1*X(M,J,K) + P2*X(M+1,J,K) + P3*X(M+2,J,K)
            YJ(I,J1,K1) = P1*Y(M,J,K) + P2*Y(M+1,J,K) + P3*Y(M+2,J,K)
            ZJ(I,J1,K1) = P1*Z(M,J,K) + P2*Z(M+1,J,K) + P3*Z(M+2,J,K)

            IF (.NOT. NOQ) THEN
               DO N = 1, NQDIM
                  QJ(I,N) =
     >               P1*Q(M,J,K,N) + P2*Q(M+1,J,K,N) + P3*Q(M+2,J,K,N)
               END DO
            END IF

         END DO

C        Carol:  Can  we dispense with this?

         XJ(1,J1,K1) = XJ1;  XJ(NIPTS,J1,K1) = XJN
         YJ(1,J1,K1) = YJ1;  YJ(NIPTS,J1,K1) = YJN
         ZJ(1,J1,K1) = ZJ1;  ZJ(NIPTS,J1,K1) = ZJN

      END IF

      END SUBROUTINE INTXYZQ

C***********************************************************************
C
      SUBROUTINE LAGCOF (SNEW, SARR, NPTS, M, P1, P2, P3)
C
C     Calculate 2-pt. or 3-pt. Lagrange polynomial coefficients
C     P1, P2, P3 needed for interpolating at SNEW in the SARR array.
C     M is the first index to use in the interpolated function array. 
C
C     CALLED BY: INTF, INTXYZQ if INTER = 2 or 3
C
C***********************************************************************

      USE SAGEMOD  ! Needed only for INTER, which should be an argument

      IMPLICIT NONE

      INTEGER, INTENT (IN)  ::  NPTS
      INTEGER, INTENT (OUT) ::  M
      REAL,    INTENT (IN)  ::  SNEW, SARR(NPTS)
      REAL,    INTENT (OUT) ::  P1, P2, P3

      INTEGER  I
      REAL     S1, S2, S3

C     Execution:

C     If data outside range, extrapolate linearly:

      IF (SNEW <= SARR(1)) THEN
         M = 1
         IF (SNEW == SARR(1)) THEN
            P1 = ONE
            P2 = ZERO
            P3 = ZERO
         ELSE
            P1 = (SARR(2) - SNEW) / (SARR(2) - SARR(1))
            P2 = ONE - P1
            P3 = ZERO
         END IF
         GO TO 999  ! Carol:   This GO TO had been above the END IF
      END IF

      IF (SNEW >= SARR(NPTS)) THEN
         M  = NPTS - 1
         P2 = (SNEW - SARR(NPTS-1))/(SARR(NPTS)-SARR(NPTS-1))       
         P1 = ONE - P2
         P3 = ZERO
         GO TO 999
      END IF

      IF (INTER == 2) THEN ! Linear interpolation

         DO I = 1, NPTS
            IF (SNEW > SARR(I) .AND. SNEW <= SARR(I+1)) THEN
               M  = I
               P1 = (SARR(I+1) - SNEW) / (SARR(I+1) - SARR(I))
               P2 = ONE - P1
               P3 = ZERO
               GO TO 999
            END IF
         END DO

      END IF

C     3-pt. interpolation:

C     Carol:  This is linear in the first or last interval, but need not be.

      IF (SNEW <= SARR(2)) THEN
         M  = 1
         P1 = (SARR(2) - SNEW) / (SARR(2) - SARR(1))
         P2 = ONE - P1
         P3 = ZERO
         GO TO 999
      END IF        

      IF (SNEW >= SARR(NPTS-1)) THEN
         M = NPTS - 1
         P2 = (SNEW - SARR(NPTS-1)) / (SARR(NPTS) - SARR(NPTS-1))
         P1 = ONE - P2
         P3 = ZERO
         GO TO 999
      END IF

      DO I = 2, NPTS - 2
         IF (SNEW >= SARR(I) .AND. SNEW < SARR(I+1)) THEN
            M = I            ! Is backward or forward differencing is best?
            IF (SNEW <= 0.5*(SARR(I) + SARR(I+1))) M = I - 1

            S1 = SARR(M)
            S2 = SARR(M+1)
            S3 = SARR(M+2)
            P1 = (SNEW - S2)*(SNEW - S3)/((S1 - S2)*(S1 - S3))
            P2 = (SNEW - S1)*(SNEW - S3)/((S2 - S1)*(S2 - S3))
            P3 = (SNEW - S2)*(SNEW - S1)/((S3 - S1)*(S3 - S2))
            GO TO 999
         END IF
      END DO

  999 CONTINUE

      END SUBROUTINE LAGCOF

C*********************************************************************** 
C
      SUBROUTINE LINE1
C
C     LINE1 adapts the first line of the first plane using a 1-D technique.
C
C     CALLED BY: Main progam
C
C     CALLS: EDGEMG, FILTER, GETWT, INTF, NORM, WTEDGE
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Local variables:

      INTEGER
     >   I, ITER, L, MAXLITS

      REAL
     >   ERR, ERRMIN, WTSUM
      REAL, DIMENSION (IMX) ::
     >   SNT, SNBEST, WTBEST, DSBEST

C     Execution:

      ERRMIN = ZERO
      MAXLITS = MAXITS + 30

      DO 500 ITER = 1, MAXLITS

         DO I = 1, NIPTS - 1             ! Explain
            WEIGHT(I) = ONE + A*FB(I)**B
         END DO

C        If edge problems (?), override the influence of FB: 

         IF (NEDGE /= 0) THEN
            IF (ITER == 1) CALL WTEDGE (JST-1, KST)
            IF (MG1 /= 0)  WEIGHT(1) = WDS
            IF (MG2 /= 0)  WEIGHT(NIPTS-1) = WDE
            CALL EDGEMG (WEIGHT)
         END IF

         IF (ITER /= 1) THEN
            CALL GETWT
            DO I = 1, NIPTS - 1
               WEIGHT(I) = WEIGHT(I)*WT(I)
            END DO
         END IF

         IF (NFILT /= 0) CALL FILTER (NIPTS-1, NFILT, WEIGHT)

         WTSUM = ZERO
         DO L = 1, NIPTS - 1
            WTSUM = WTSUM + ONE / WEIGHT(L)
         END DO

C        New values of S, DS:

         SNT(1) = ZERO
         DO I = 1, NIPTS - 1
            DS(I) = SS(NIPTS) / (WEIGHT(I)*WTSUM)
            SNT(I+1) = SNT(I) + DS(I)
         END DO

C        Convergence test:

         ERR = ZERO
         DO L = 1, NIPTS - 1
            ERR = ERR + ABS ((SNT(L) - SN(L)) / SS(NIPTS))
            SN(L) = SNT(L)
         END DO
         IF (ERR < CONV) GO TO 999

         IF (ITER == 1 .OR. (ITER > 1 .AND. ERR < ERRMIN)) THEN
            ERRMIN = ERR
            DO L = 1, NIPTS
               SNBEST(L) = SNT(L)
               WTBEST(L) = WEIGHT(L)
               DSBEST(L) = DS(L)
            END DO
         END IF

C        If no convergence, interpolate for F at new values of S,
C        and re-evaluate FB:

         CALL INTF (INTER, NIPTS-1, SMS, SN, F, FB)
         CALL NORM (NIPTS-1, FB)

  500 CONTINUE ! Next iteration

      DO I = 1, NIPTS
         SN(I)     = SNBEST(I)
         DS(I)     = DSBEST(I)
         WEIGHT(I) = WTBEST(I)
      END DO

      WRITE (6, '(/, A, 1P, E10.3)')
     >   ' No convergence along initial line. ERRMIN:', ERRMIN

999   CONTINUE

      END SUBROUTINE LINE1

C***********************************************************************
C
      SUBROUTINE MARCHJ (K)
C
C     Extend the proportions of the last adapted J line throughout the
C     rest of the current plane up to JMAX.
C
C     CALLED BY: Main program
C
C     CALLS: INTXYZQ
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER K

C     Local variables:

      INTEGER I, J, L

      REAL    DSL, QJ(IMX,NDIM)

C     Execution:

C     Proportion remaining lines (explain):

      SS(1) = ZERO

      DO 200 J = JEND + 1, JMAX

         L = IST
         DO I = 1, NIPTS - 1
            DSL = SQRT ((X(L+1,J,K) - X(L,J,K))**2 +
     >                  (Y(L+1,J,K) - Y(L,J,K))**2 +
     >                  (Z(L+1,J,K) - Z(L,J,K))**2)
            SS(I+1) = SS(I) + DSL
            L = L + 1
         END DO

         DO I = 1, NIPTS
            SN(I) = SNM(I) * SS(NIPTS) / SNM(NIPTS)
         END DO

C        Interpolate X, Y Z and Q at new S:

         CALL INTXYZQ (J, K, 2, 2, SS, SN, QJ)

         L = IST
         DO I = 2, NIPTS - 1
            L = L + 1
            X(L,J,K) = XJ(I,2,2)
            Y(L,J,K) = YJ(I,2,2)
            Z(L,J,K) = ZJ(I,2,2)
            IF (.NOT. NOQ) Q(L,J,K,1:NQDIM) = QJ(I,1:NQDIM)
         END DO

  200 CONTINUE ! Next J

      END SUBROUTINE MARCHJ

C***********************************************************************
C
      SUBROUTINE MARCHK
C
C     Extend the proportions of the last adapted K plane throughout the
C     remaining grid, up to KMAX.
C
C     CALLED BY: OUTPUT
C
C     CALLS: INTXYZQ
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Local variables:

      INTEGER I, J, K, L, JSTOP

      REAL    DSL, SNK(IMX), QJ(IMX,NDIM)

C     Execution:

      IF (MARCH) THEN 
         JSTOP = JMAX
      ELSE
         JSTOP = JEND
      END IF

C     Calculate S on last adapted plane, for each line:

      SNK(1) = ZERO

      DO 200 J = 1, JSTOP

         L = IST
         DO I = 1, NIPTS - 1
            DSL = SQRT ((X(L+1,J,KEND) - X(L,J,KEND))**2 +
     >                  (Y(L+1,J,KEND) - Y(L,J,KEND))**2 +
     >                  (Z(L+1,J,KEND) - Z(L,J,KEND))**2)
            SNK(I+1) = SNK(I) + DSL
            L = L + 1
         END DO 

C        For each K plane with this line, proportion new S:

         DO K = KEND + 1, KMAX

            SS(1) = ZERO
            L = IST
            DO I = 1, NIPTS - 1
               DSL = SQRT ((X(L+1,J,K) - X(L,J,K))**2 +
     >                     (Y(L+1,J,K) - Y(L,J,K))**2 +
     >                     (Z(L+1,J,K) - Z(L,J,K))**2)
               SS(I+1) = SS(I) + DSL    
               L = L + 1
            END DO

            DO I = 1, NIPTS
               SN(I) = SNK(I) * SS(NIPTS) / SNK(NIPTS)
            END DO

C           Calculate X, Y, Z and Q at new S:

            CALL INTXYZQ (J, K, 2, 2, SS, SN, QJ)

            L = IST
            DO I = 2, NIPTS - 1
               L = L + 1
               X(L,J,K) = XJ(I,2,2)
               Y(L,J,K) = YJ(I,2,2)
               Z(L,J,K) = ZJ(I,2,2)
               IF (.NOT. NOQ) Q(L,J,K,1:NQDIM) = QJ(I,1:NQDIM)
            END DO
       
         END DO ! Next K

  200 CONTINUE ! Next J

      END SUBROUTINE MARCHK

C***********************************************************************
C
      SUBROUTINE MATCH (NOMORE, IER)
C
C     This routine controls the transfer of data for matching domains.
C     No adaption is performed.  Export control data have already been
C     read.
C
C     CALLED BY: INITIAL
C
C     CALLS: HEADIO, MULTIO, READMULT, REWND, STOREX, STORIM
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      LOGICAL, INTENT (OUT) :: NOMORE
      INTEGER, INTENT (OUT) :: IER

C     Local variables:

      INTEGER  MGE, MPLE, NEX1, NEX2, NSM, NSM1
      REAL     DSN

      NAMELIST /NAMEL/
     >   IST, IEND, JST, JEND, KST, KEND, RDSMAX, RDSMIN, ISTEP, JSTEP,
     >   KSTEP, IJPLANE, IKPLANE, JKPLANE, MGRID, XMIN, XMAX, CLAM,
     >   NEDGE, CT, MGSTEPS, MGPLS, INDQ, IQ, NOUP, TWOD, REMOVE, NQDIM,
     >   MARCH, MARCHPL, ADD, LSTADD, LENDADD, SUB, LSTSUB, LENDSUB,
     >   SAVE, ORTHS, ORTHE, NFILT, INTER, GEOM, QFUN, NOQ, LNSING,
     >   PLSING, MG1, MG2, MPLANE, EXPORT, IMPORT, IS, IE, JS, JE, KS,
     >   KE, FV, FORMI, FORMO, FUN, MVBOUND, XMINP, XMAXP, MVP1, MVP2,
     >   DSN, NSM, NSM1, IMS, JME, NOSESM, RECLUST, DSW, DSW1, DSW2,
     >   DSE, MVC1, MVC2

C     Execution:

C     Test that export card has consistent formats:

      FORMO = FORMI

C     Read file header records:

      IF (NADS > 1) CALL REWND

      CALL HEADIO (1, IER)
      IF (IER == 1) GO TO 999

C     Allocate storage for a plane:

      ALLOCATE (XP(IMX,IMX), YP(IMX,IMX), ZP(IMX,IMX), STAT=IERA(1))

      IF (NOQ) THEN
         IERA(2) = 0
      ELSE
         ALLOCATE (QP(IMX,IMX,NQDIM), STAT=IERA(2))
      END IF

      IF (IERA(1) /= 0 .AND. IERA(2) /= 0) THEN
         WRITE (6, '(/, A)') ' ** Export/import allocation problem **'
         STOP
      END IF

C     Read grid and Q files up to mgrid (the export grid):

      CALL MULTIO (1, MGRID, 0, 0)

      IMAX = IM(MGRID)
      JMAX = JM(MGRID)
      KMAX = KM(MGRID)
      MGE  = MGRID
      MPLE = MPLANE

C     Grid containing export grid is now in memory; save plane:

      CALL STOREX (NEX1, NEX2)

C     Restore plane defaults:

      IJPLANE = .TRUE.
      JKPLANE = .FALSE.
      IKPLANE = .FALSE.
      FORMCK  = FORMI

C     Read import control parameter file:

      READ (5, NAMEL)

      IF (IKPLANE .OR. JKPLANE) IJPLANE = .FALSE.
      FORMI = FORMCK
      FORMO = FORMI
      IF (.NOT. IMPORT) THEN
         WRITE  (6, '(/, A)') ' ** Import card expected; not found **'
         NOMORE = .TRUE.
         GO TO 999
      END IF

C     Reread input file up to import grid:

      CALL DEALLO
      CLOSE(NITGI)
      IF (FORMI) OPEN (NITGI, FORM='FORMATTED')
      IF (.NOT. FORMI) OPEN (NITGI, FORM='UNFORMATTED')
      IF (.NOT. NOQ) THEN
         CLOSE (NITQI)
         IF (.NOT. FORMI) OPEN (NITQI, FORM='UNFORMATTED')
         IF (FORMI) OPEN (NITQI, FORM='FORMATTED')
      END IF
      CALL READMULT (IER)
      IF (IER == 1) GO TO 999

C     The import grid is now in memory. 
C     Place the stored plane into the correct import plane:

      CALL STORIM (NEX1, NEX2)
      CALL MULTIO (MGRID+1, NGRID, MGRID, NGRID)

      WRITE (6, '(/, 4(A, I4))')
     >   ' Plane', MPLE, ', grid', MGE,
     >   ' copied to plane', MPLANE, ', grid', MGRID

  999 CONTINUE

      FORMOO = FORMI
      DEALLOCATE (XP, YP, ZP, STAT=IERA(1))
      IF (.NOT. NOQ) DEALLOCATE (QP, STAT=IERA(2))
      IF (IERA(1) /= 0 .AND. IERA(2) /= 0) THEN
         WRITE (6, '(/, A)') ' ** Export/import deallocation problem **'
         STOP
      END IF
      CALL DEALLO

      END SUBROUTINE MATCH

C***********************************************************************
C
      SUBROUTINE MGWALLS (J)
C
C     When wall geometry has large gradients, an input parameter (GEOM)
C     requests the inclusion of geometry gradients (F(G)) into the
C     adaption variable.  This function should be gradually decreased
C     away from both wall boundaries by computing the coefficient FGW,
C     which is a function of aspect ratio. 
C
C     CALLED BY: FBAR
C
C***********************************************************************
C
      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER J

C     Local variables:

      REAL    DSA1, DSA2, FGASP

C     Execution:

C     Compute aspect ratio using average of new and old delta S.
C     No geometry effect is included if aspect ratio > 1/4.

      IF (J <= JST + 1 .OR. J >= JEND - 1) THEN
         FGW = ONE
      ELSE
         DSA1  = SS(MXFG+1) - SS(MXFG)
         DSA2  = SN(MXFG+1) - SN(MXFG)
         FGASP = (8.0 * ABS (DAP(MXFG))) / (DSA1 + DSA2)

C        If aspect ratio indicates no geom effect, ensure smooth turn off:

         IF (FGASP > ONE) THEN
            IF (J > JST + 4 .AND. J < JEND - 4) THEN
               FGW = ZERO
            ELSE
               IF (J <= JST  + 4) FGW = FGW * (JST - REAL (J - 2) *0.25)
               IF (J >= JEND - 4) FGW = REAL ((J - JEND + 5)) * 0.25
            END IF
         ELSE
            FGW = ONE - FGASP
         END IF
      END IF

      END SUBROUTINE MGWALLS

C***********************************************************************
C
      SUBROUTINE MULTIO (IN1, IN2, IOUT1, IOUT2)
C
C     MULTIO reads and writes the indicated grid blocks.
C     Header records are treated elsewhere.
C
C     CALLED BY: MATCH, READMULT, OUTPUT
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: IN1, IN2, IOUT1, IOUT2 ! Explain

C     Local variables:

      INTEGER I, J, K, N, I1, I2, MG

C     Execution:

      IF (IN1   == 0) I1 = IOUT1
      IF (IOUT1 == 0) I1 = IN1
      IF (IN1   /= 0 .AND. IOUT1 /= 0) I1 = MIN (IN1, IOUT1)
      I2 = MAX (IN2, IOUT2)

      DO 100 MG = I1, I2

       IF (MG >= IN1 .AND. MG <= IN2) THEN

        IF (.NOT. TWOD) THEN

         IF (.NOT. FORMI) THEN
          READ (NITGI)  (((X(I,J,K),I=1,IM(MG)),J=1,JM(MG)),K=1,KM(MG)),
     >                  (((Y(I,J,K),I=1,IM(MG)),J=1,JM(MG)),K=1,KM(MG)),
     >                  (((Z(I,J,K),I=1,IM(MG)),J=1,JM(MG)),K=1,KM(MG))
         ELSE
          READ(NITGI,*) (((X(I,J,K),I=1,IM(MG)),J=1,JM(MG)),K=1,KM(MG)),
     >                  (((Y(I,J,K),I=1,IM(MG)),J=1,JM(MG)),K=1,KM(MG)),
     >                  (((Z(I,J,K),I=1,IM(MG)),J=1,JM(MG)),K=1,KM(MG))
         END IF

         IF (.NOT. NOQ) THEN
          IF (.NOT. FORMI) THEN
           IF (.NOT. FUN) READ (NITQI) FSMACH,ALPHA,RE,TIME
           READ  (NITQI)  ((((Q(I,J,K,N), I=1,IMQ(MG)),J=1,JMQ(MG)),
     >                                    K=1,KMQ(MG)),N=1,NQDIM)
          ELSE
           IF (.NOT. FUN) READ (NITQI,*) FSMACH,ALPHA,RE,TIME
           READ (NITQI,*) ((((Q(I,J,K,N), I=1,IMQ(MG)),J=1,JMQ(MG)),
     >                                    K=1,KMQ(MG)),N=1,NQDIM)
          END IF
         END IF

        ELSE ! 2-D

         IF (.NOT. FORMI) THEN
          READ  (NITGI)  ((X(I,J,1),I=1,IM(MG)),J=1,JM(MG)),
     >                   ((Y(I,J,1),I=1,IM(MG)),J=1,JM(MG))
         ELSE
          READ (NITGI,*) ((X(I,J,1),I=1,IM(MG)),J=1,JM(MG)),
     >                   ((Y(I,J,1),I=1,IM(MG)),J=1,JM(MG))
         END IF

         IF (.NOT. NOQ) THEN
          IF (.NOT. FORMI) THEN
           IF (.NOT. FUN) READ (NITQI) FSMACH,ALPHA,RE,TIME
           READ  (NITQI)  (((Q(I,J,1,N), I=1,IMQ(MG)),J=1,JMQ(MG)),
     >                                   N=1,NQDIM)
          ELSE
           IF (.NOT. FUN) READ (NITQI,*) FSMACH,ALPHA,RE,TIME
           READ (NITQI,*) (((Q(I,J,1,N), I=1,IMQ(MG)),J=1,JMQ(MG)),
     >                                   N=1,NQDIM)
          END IF
         END IF

        END IF ! 2-D/3-D block

       END IF ! Read block

C      Write block:

       IF (MG >= IOUT1 .AND. MG <= IOUT2) THEN

        IF (FORMO) CALL CHSMALL (MG) ! Guard against 3-digit exponents

        IF (.NOT. TWOD) THEN

         IF (.NOT. FORMO) THEN
          WRITE (NITGO) (((X(I,J,K),I=1,IM(MG)),J=1,JM(MG)),K=1,KM(MG)),
     >                  (((Y(I,J,K),I=1,IM(MG)),J=1,JM(MG)),K=1,KM(MG)),
     >                  (((Z(I,J,K),I=1,IM(MG)),J=1,JM(MG)),K=1,KM(MG))
         ELSE
          WRITE (NITGO, 2000)
     >                (((X(I,J,K),I = 1,IM(MG)),J=1,JM(MG)),K=1,KM(MG)),
     >                (((Y(I,J,K),I = 1,IM(MG)),J=1,JM(MG)),K=1,KM(MG)),
     >                (((Z(I,J,K),I = 1,IM(MG)),J=1,JM(MG)),K=1,KM(MG))
         END IF

         IF (.NOT. NOQ) THEN
          IF (.NOT. FORMO) THEN
           IF (.NOT. FUN) WRITE (NITQO) FSMACH,ALPHA,RE,TIME
           WRITE (NITQO)      ((((Q(I,J,K,N),I=1,IMQ(MG)),J=1,JMQ(MG)),
     >                                       K=1,KMQ(MG)),N=1,NQDIM)
          ELSE
           IF (.NOT. FUN) WRITE (NITQO, 2000) FSMACH,ALPHA,RE,TIME
           WRITE (NITQO, 2000) ((((Q(I,J,K,N),I=1,IMQ(MG)),J=1,JMQ(MG)),
     >                                        K=1,KMQ(MG)),N=1,NQDIM)
          END IF
         END IF

        ELSE ! 2-D

         IF (.NOT. FORMO) THEN
          WRITE (NITGO)      ((X(I,J,1),I=1,IM(MG)),J=1,JM(MG)),
     >                       ((Y(I,J,1),I=1,IM(MG)),J=1,JM(MG))
         ELSE
          WRITE (NITGO,2000) ((X(I,J,1),I=1,IM(MG)),J=1,JM(MG)),
     >                       ((Y(I,J,1),I=1,IM(MG)),J=1,JM(MG))
         END IF

         IF (.NOT. NOQ) THEN
          IF (.NOT. FORMO) THEN
           IF (.NOT. FUN) WRITE (NITQO) FSMACH,ALPHA,RE,TIME
           WRITE (NITQO)      (((Q(I,J,1,N),I=1,IMQ(MG)),J=1,JMQ(MG)),
     >                                      N=1,NQDIM)
          ELSE
           IF (.NOT. FUN) WRITE (NITQO, 2000) FSMACH,ALPHA,RE,TIME
           WRITE (NITQO,2000) (((Q(I,J,1,N),I=1,IMQ(MG)),J=1,JMQ(MG)),
     >                                      N=1,NQDIM)
          END IF
         END IF

        END IF ! 2-D/3-D block

       END IF ! Write block

  100 CONTINUE ! Next grid block

 2000 FORMAT (1P, 4E19.11)

      END SUBROUTINE MULTIO

C***********************************************************************
C
      SUBROUTINE NOADAPT (J, K)
C
C     This routine updates variables needed for stepping to the next
C     line when the current line is not adapted.
C
C     CALLED BY: Main program
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER J, K

C     Local variables:

      INTEGER I, L
      REAL    DSS

C     Execution:

C     If singular line, store data from line 1, plane 1:

      IF (J == LNSING .AND. K /= KST) THEN
         L = IST
         DO I = 1, NIPTS
            X(L,J,K) = X(L,J,KST)
            Y(L,J,K) = Y(L,J,KST)
            Z(L,J,K) = Z(L,J,KST)
            IF (.NOT. NOQ) Q(L,J,K,1:NQDIM) = Q(L,J,KST,1:NQDIM)
            L = L + 1
         END DO
      END IF

      SNM(1) = ZERO
      L = IST
      DO I = 1, NIPTS - 1
         DSS = SQRT ((X(L+1,J,K) - X(L,J,K))**2 +
     >               (Y(L+1,J,K) - Y(L,J,K))**2 +
     >               (Z(L+1,J,K) - Z(L,J,K))**2)
         SNM(I+1) = SNM(I) + DSS
         L = L + 1
      END DO

      END SUBROUTINE NOADAPT

C***********************************************************************
C
      SUBROUTINE NORM (NPTS, F1)
C
C     Normalize input flow-field vector F1 unless it is essentially
C     constant, in which case, reset all elements to a small number.
C     F1 is updated in-place.
C
C     CALLED BY: FBAR, LINE1, SOLUT
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: NPTS
      REAL, INTENT (INOUT) :: F1(NPTS)

C     Local constants:

      REAL, PARAMETER :: EPS = 1.E-4, SMALL = 1.E-5

C     Local variables:

      INTEGER I
      REAL    DENOM, FMAX, FMIN

C     Execution:

      FMIN = F1(1)
      FMAX = F1(1)
      DO I = 2, NPTS
         IF (F1(I) > FMAX) THEN
            FMAX = F1(I)
         ELSE IF (F1(I) < FMIN) THEN
            FMIN = F1(I)
         END IF
      END DO

C     Check for ~constant flow-field:

      IF (ABS (FMAX - FMIN) < EPS) THEN
         F1 = SMALL ! 1:NPTS
      ELSE
         DENOM = 1.0 / (FMAX - FMIN)
         DO I = 1, NPTS
            F1(I) = MAX (SMALL, (F1(I) - FMIN) * DENOM)
         END DO
      END IF

      END SUBROUTINE NORM

C***********************************************************************
C
      SUBROUTINE NORMPT (IP, JP, KP, INDPL, PA, PB, PC, SING)
C
C     Calculate the components PA, PB, PC of the unit vector normal to
C     the plane defined by INDPL = 1 (IKPLANE) or INDPL = 2 (IJPLANE)
C     at grid point A (IP, JP, KP).  This vector is the average of the
C     four normals defined by the four pairs of cell edges at point A
C     in the indicated index plane.  The points surrounding A in this
C     plane are B, C, D, E.  SING = 1. on return for a degenerate case.
C
C     CALLED BY: SETUPJ, SETUPK
C
C     CALLS: ADDV, PURPLE, UNITV
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: IP, JP, KP, INDPL
      REAL,   INTENT (OUT) :: PA, PB, PC, SING

C     Local variables:

      REAL XA, YA, ZA, XB, YB, ZB, XC, YC, ZC, XD, YD, ZD, XE, YE, ZE,
     >     XAB, YAB, ZAB, XAC, YAC, ZAC, XDA, YDA, ZDA, XEA, YEA, ZEA,
     >     S1A, S1B, S1C, S2A, S2B, S2C, S3A, S3B, S3C, S4A, S4B, S4C,
     >     S12A, S12B, S12C, S34A, S34B, S34C

C     Execution:

      XA = XJ(IP,JP,KP)
      YA = YJ(IP,JP,KP)
      ZA = ZJ(IP,JP,KP)

      IF (INDPL == 1) THEN ! (IKPLANE)
         XB = XJ(IP+1,JP,KP)
         YB = YJ(IP+1,JP,KP)
         ZB = ZJ(IP+1,JP,KP)
         XC = XJ(IP,JP,KP+1)
         YC = YJ(IP,JP,KP+1)
         ZC = ZJ(IP,JP,KP+1)
         XD = XJ(IP-1,JP,KP)
         YD = YJ(IP-1,JP,KP)
         ZD = ZJ(IP-1,JP,KP)
         XE = XJ(IP,JP,KP-1)
         YE = YJ(IP,JP,KP-1)
         ZE = ZJ(IP,JP,KP-1)
      ELSE      ! INDPL == 2 (IJPLANE)
         XB = XJ(IP,JP+1,KP)
         YB = YJ(IP,JP+1,KP)
         ZB = ZJ(IP,JP+1,KP)
         XC = XJ(IP+1,JP,KP)
         YC = YJ(IP+1,JP,KP)
         ZC = ZJ(IP+1,JP,KP)
         XD = XJ(IP,JP-1,KP)
         YD = YJ(IP,JP-1,KP)
         ZD = ZJ(IP,JP-1,KP)
         XE = XJ(IP-1,JP,KP)
         YE = YJ(IP-1,JP,KP)
         ZE = ZJ(IP-1,JP,KP)
      END IF

C     Calculate the unit vectors between A and its four neighbors in the plane:

      IF (XB == 999.0) THEN ! Explain the 999.0
         XAB = ZERO
         YAB = ZERO
         ZAB = ZERO
      ELSE
         CALL UNITV (XA, YA, ZA, XB, YB, ZB, XAB, YAB, ZAB)
      END IF

      IF (XC == 999.0) THEN
         XAC = ZERO
         YAC = ZERO
         ZAC = ZERO
      ELSE
         CALL UNITV (XA, YA, ZA, XC, YC, ZC, XAC, YAC, ZAC)
      END IF

      IF (XD == 999.0) THEN
         XDA = ZERO
         YDA = ZERO
         ZDA = ZERO
      ELSE
         CALL UNITV (XD, YD, ZD, XA, YA, ZA, XDA, YDA, ZDA)
      END IF

      IF (XE == 999.0) THEN
         XEA = ZERO
         YEA = ZERO
         ZEA = ZERO
      ELSE
         CALL UNITV (XE, YE, ZE, XA, YA, ZA, XEA, YEA, ZEA)
      END IF

C     Find the normal to each of the planes defined by four pairs of vectors:

      CALL PURPLE (XAC, YAC, ZAC, XAB, YAB, ZAB, S1A, S1B, S1C, NFLAG)
      CALL PURPLE (XAC, YAC, ZAC, XDA, YDA, ZDA, S2A, S2B, S2C, NFLAG)
      CALL PURPLE (XEA, YEA, ZEA, XDA, YDA, ZDA, S3A, S3B, S3C, NFLAG)
      CALL PURPLE (XEA, YEA, ZEA, XAB, YAB, ZAB, S4A, S4B, S4C, NFLAG)

C     Average the normals:

      CALL ADDV (S1A, S1B, S1C, ONE, S2A, S2B, S2C, ONE, S12A,S12B,S12C)
      CALL ADDV (S3A, S3B, S3C, ONE, S4A, S4B, S4C, ONE, S34A,S34B,S34C)
      CALL ADDV (S12A, S12B, S12C, ONE, S34A, S34B, S34C, ONE, PA,PB,PC)

      IF (PA**2 + PB**2 + PC**2 < 0.9) THEN
         SING = ONE  ! Normally this is 1.0 by construction in ADDV
      ELSE
         SING = ZERO ! Not singular (the normal case)
      END IF

      END SUBROUTINE NORMPT

C***********************************************************************
C
      SUBROUTINE OUTPUT
C
C     Write the adapted grid and flow field as PLOT3D files.
C     First, update remaining K planes if continuity is requested at
C     non-adapted points (MARCH = T).
C     If necessary, reinvert arrays and/or swap X, Y, Z, Q back to the
C     original I, J, K order.
C     Unit numbers for the output files depend on the number of passes.
C
C     CALLED BY: Main program
C
C     CALLS: FVORG, MARCHK, SWAPINV, SWAPXYZ, MULTIO
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Local variables:

      LOGICAL RSWAP

C     Execution:

C     If continuity is required at non-adapted points (MARCH = T),
C     interpolate for new X, Y, Z, Q at remaining K planes:

      IF (MARCHPL .AND. KEND < KMAX) CALL MARCHK

C     If necessary, re-invert arrays and swap data back to original order:

      IF (IINVERSE .OR. JINVERSE .OR. KINVERSE) CALL SWAPINV       

      RSWAP = .TRUE.
      IF (.NOT. (IJPLANE .AND. JSTEP)) CALL SWAPXYZ (RSWAP)

      IF (FV .AND. .NOT. NOQ) CALL FVORG (2)

C     Output adapted grid and flow field:

      IF (OK) WRITE (6, '(/, A, I3, A, I3)')
     >      ' Adaption', MGCT(MGRID), ' complete for grid', MGRID

      IF (.NOT. OK .OR. SAVE) THEN
         IF (ADD /= 0 .OR. SUB /= 0 .OR. REMOVE /= 0) THEN
            IM(MGRID)  = IMAX
            JM(MGRID)  = JMAX
            KM(MGRID)  = KMAX
            IMQ(MGRID) = IMAX
            JMQ(MGRID) = JMAX
            KMQ(MGRID) = KMAX
            IF (FV) THEN
               IMQ(MGRID) = IMAX - 1
               JMQ(MGRID) = JMAX - 1
               IF (KMAX /= 1) KMQ(MGRID) = KMAX - 1
            END IF
         END IF
         CALL MULTIO (MGRID+1, NGRID, MGRID, NGRID)
      END IF

      FORMOO = FORMO

      END SUBROUTINE OUTPUT

C***********************************************************************
C
      SUBROUTINE PROPS (J, K)
C
C     This routine takes the converged value of S on the J-1, K line
C     and proportions it to the remaining block of XJ, YJ, ZJ.
C
C     CALLED BY: FBAR
C
C     CALLS: INTXYZQ
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER J, K

C     Local variables:

      INTEGER I, L
      REAL    DSSK, RATIO,
     >        SSNK(IMX), SSK(IMX), SNK(IMX), SSJK(IMX), SNJK(IMX),
     >        DUM(IMX,NDIM)

C     Execution:

C     Proportion J line into SN:

      SN(1) = ZERO

      IF (J == JST) THEN
         RATIO = SS(NIPTS) / SNMK(NIPTS)
         DO I = 2, NIPTS
            SN(I) = SNMK(I) * RATIO
         END DO
      ELSE
         RATIO = SS(NIPTS) / SNM(NIPTS)
         DO I = 2, NIPTS
            SN(I) = SNM(I) * RATIO
         END DO
      END IF

      CALL INTXYZQ (J, K, 2, 2, SS, SN, DUM)

C     At J+1, K (for COSBK):  ! Explain

      IF (J /= JEND) THEN

         SSK(1) = ZERO
         L = IST
         DO I = 1, NIPTS - 1
            DSSK = SQRT ((X(L+1,J+1,K) - X(L,J+1,K))**2 +
     >                   (Y(L+1,J+1,K) - Y(L,J+1,K))**2 +
     >                   (Z(L+1,J+1,K) - Z(L,J+1,K))**2)
            SSK(I+1) = SSK(I) + DSSK
            L = L + 1
         END DO

         IF (J == JST) THEN
            RATIO = SSK(NIPTS) / SNMK(NIPTS)
            DO I = 1, NIPTS
               SSNK(I) = SNMK(I) * RATIO
            END DO
         ELSE
            RATIO = SSK(NIPTS) / SNM(NIPTS)
            DO I = 1, NIPTS
               SSNK(I) = SNM(I) * RATIO
            END DO
         END IF

         CALL INTXYZQ (J+1, K, 3, 2, SSK, SSNK, DUM)

      END IF

C     K+1 plane, both J-1 (for normal U) and J (for normal B):

      IF (K /= KMAX .AND. J /= JST) THEN

C        At J-1, K+1:

         SSK(1) = ZERO
         L = IST
         DO I = 1, NIPTS - 1
            DSSK = SQRT ((X(L+1,J-1,K+1) - X(L,J-1,K+1))**2 +
     >                   (Y(L+1,J-1,K+1) - Y(L,J-1,K+1))**2 +
     >                   (Z(L+1,J-1,K+1) - Z(L,J-1,K+1))**2)
            SSK(I+1) = SSK(I) + DSSK
            L = L + 1
         END DO

         RATIO = SSK(NIPTS) / SNM(NIPTS)
         DO I = 1, NIPTS
            SNK(I) = SNM(I) * RATIO
         END DO

         CALL INTXYZQ (J-1, K+1, 1, 3, SSK, SNK, DUM)

C        At J, K+1:

         SSJK(1) = ZERO
         L = IST
         DO I = 1, NIPTS - 1
            DSSK = SQRT ((X(L+1,J,K+1) - X(L,J,K+1))**2 +
     >                   (Y(L+1,J,K+1) - Y(L,J,K+1))**2 + 
     >                   (Z(L+1,J,K+1) - Z(L,J,K+1))**2)
            SSJK(I+1) = SSJK(I) + DSSK
            L = L + 1
         END DO

         RATIO = SSJK(NIPTS) / SNM(NIPTS)
         DO I = 1, NIPTS
            SNJK(I) = SNM(I) * RATIO
         END DO

         CALL INTXYZQ (J, K+1, 2, 3, SSJK, SNJK, DUM)

      END IF

      END SUBROUTINE PROPS

C***********************************************************************
C
      SUBROUTINE PURPLE (A1, A2, A3, B1, B2, B3, V1, V2, V3, NFL)
C
C     (PURPLE == "perpendicular".)
C     Calculate the cross product of unit vectors A and B and return it
C     as a unit vector V perpendicular to the plane defined by A and B.
C     Flag NFL = 1 or -1 permits reversing the direction.
C
C     CALLED BY: CROSSV, NORMPT
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      REAL,    INTENT (IN)  :: A1, A2, A3, B1, B2, B3
      INTEGER, INTENT (IN)  :: NFL
      REAL,    INTENT (OUT) :: V1, V2, V3

C     Local variables:

      REAL     VMOD

C     Execution:

      V1 = A2*B3 - A3*B2
      V2 = A3*B1 - A1*B3
      V3 = A1*B2 - A2*B1

      VMOD = SQRT (V1**2 + V2**2 + V3**2)

      IF (VMOD /= 0.0) THEN
         VMOD = REAL (NFL) / VMOD
         V1   = V1 * VMOD
         V2   = V2 * VMOD
         V3   = V3 * VMOD
      END IF

      END SUBROUTINE PURPLE

C***********************************************************************
C
      SUBROUTINE PUSHIT (DSN, NSM, NSM1)
C
C     This routine controls the outer boundary movement.
C
C     CALLED BY: INITIAL
C
C     CALLS: SHOCKLE, LAGCOF, GETDS0, VINOKUR
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      REAL,    INTENT (IN) :: DSN
      INTEGER, INTENT (IN) :: NSM, NSM1

C     Local variables:

      INTEGER  I, IER, J, K, L, M, N, NPTS
      REAL     DS0, DS1, P1, P2, P3, RATIO, SMAX,
     >         DSC(IMX), SC(IMX), SNW(IMX), SW(IMX),
     >         XM(IMX), YM(IMX), ZM(IMX),
     >         QM(IMX,NDIM), SNMAX(IMX,IMX)

C     Execution:

C     Calculate new SMAX (based on shock or %) for each J:

      CALL SHOCKLE (SNMAX, DSN, NSM, NSM1)

      NPTS = IEND - IST + 1

C     With a boundary move, the defined domain must be the total domain:

      IF (IMAX /= NPTS .AND. RECLUST == 1) WRITE (6, '(/, A)')
     >    ' Moving outer boundary?  Entire domain must be reclustered.'

      DO 900 K = KST, KEND

         DO 800 J = JST, JEND

            IF (J == LNSING .AND. K /= 1) THEN

               DO I = 1, IMAX
                  X(I,J,K) = X(I,J,1)
                  Y(I,J,K) = Y(I,J,1)
                  Z(I,J,K) = Z(I,J,1)
                  IF (.NOT. NOQ) Q(I,J,K,1:NQDIM) = Q(I,J,1,1:NQDIM)
               END DO

               CYCLE ! Next J

            END IF

C           Arc lengths for this J line:

            SC(1) = ZERO
            L = IST
            DO I = 1, NPTS - 1
               DSC(I) = SQRT ((X(L+1,J,K) - X(L,J,K))**2 +
     >                        (Y(L+1,J,K) - Y(L,J,K))**2 +
     >                        (Z(L+1,J,K) - Z(L,J,K))**2)
               SC(I+1) = SC(I) + DSC(I)
               L = L + 1
            END DO
            SMAX = SC(NPTS)

C           Evaluate X,Y,Z at new SMAX:

            IF (SNMAX(J,K) >= SMAX) THEN
               RATIO = (SNMAX(J,K) - SMAX) / (SMAX - SC(NPTS-1))
               XM(NPTS) = RATIO*(X(IMAX,J,K)-X(IMAX-1,J,K)) +X(IMAX,J,K)
               YM(NPTS) = RATIO*(Y(IMAX,J,K)-Y(IMAX-1,J,K)) +Y(IMAX,J,K)
               ZM(NPTS) = RATIO*(Z(IMAX,J,K)-Z(IMAX-1,J,K)) +Z(IMAX,J,K)
            ELSE
               CALL LAGCOF (SNMAX(J,K), SC, NPTS, M, P1, P2, P3)
               M = IST + M - 1
               XM(NPTS) = P1*X(M,J,K) + P2*X(M+1,J,K) + P3*X(M+2,J,K)
               YM(NPTS) = P1*Y(M,J,K) + P2*Y(M+1,J,K) + P3*Y(M+2,J,K)
               ZM(NPTS) = P1*Z(M,J,K) + P2*Z(M+1,J,K) + P3*Z(M+2,J,K)
            END IF

            X(IMAX,J,K) = XM(NPTS)
            Y(IMAX,J,K) = YM(NPTS)
            Z(IMAX,J,K) = ZM(NPTS)

C           Compute new S:

            IF (RECLUST == 0) THEN
               RATIO = SNMAX(J,K) / SMAX
               DO I = 1, NPTS
                  SNW(I) = SC(I) * RATIO
               END DO
            ELSE
               DS1 = DSE * SNMAX(J,K) / REAL (NPTS-1)
               IF (DSW == ZERO) THEN
                  DS0 = DSC(1)
               ELSE
                  CALL GETDS0 (DS0, J, K)
               END IF
               SNW(1) = ZERO
               SNW(NPTS) = SNMAX(J,K)

CCCCCC         CALL CLUST2 (SNW, DS0, DS1, NPTS, IMX)
               CALL VINOKUR (1, NPTS, DS0, DS1, SNW, 6, IER)

               IF (IER /= 0) THEN ! ! Highly unlikely
                  WRITE (6, '(/, A, /, 4I5, /, A, 1P, 2E16.7)')
     >            ' Vinokur distribution did not converge in PUSHIT.',
     >            ' IER, J, K, NPTS:', IER, J, K, NPTS,
     >            ' DS0, DS1: ', DS0, DS1
                  STOP
               END IF

            END IF

C           Put new SMAX into old S array MAX location.
C           Now the new Xs and SC(IMAX) are consistent for LAGCOF.

            SC(NPTS) = SNMAX(J,K)

C           Interpolate at new S:

            DO I = 1, NPTS - 1
               CALL LAGCOF (SNW(I), SC, NPTS, M, P1, P2, P3)
               M = IST + M - 1
               XM(I) = P1*X(M,J,K) + P2*X(M+1,J,K) + P3*X(M+2,J,K)
               YM(I) = P1*Y(M,J,K) + P2*Y(M+1,J,K) + P3*Y(M+2,J,K)
               ZM(I) = P1*Z(M,J,K) + P2*Z(M+1,J,K) + P3*Z(M+2,J,K)
               IF (.NOT. NOQ) THEN
                  DO N = 1, NQDIM
                     QM(I,N) = P1*Q(M,J,K,N) + P2*Q(M+1,J,K,N) +
     >                         P3*Q(M+2,J,K,N)
                  END DO
               END IF
            END DO

C           Store new values back in X, Y, Z, Q:

            L = IST
            DO I = 1, NPTS - 1
               X(L,J,K) = XM(I)
               Y(L,J,K) = YM(I)
               Z(L,J,K) = ZM(I)
               IF (.NOT. NOQ) Q(L,J,K,1:NQDIM) = QM(I,1:NQDIM)
               L = L + 1
            END DO

  800    CONTINUE ! Next J line

  900 CONTINUE ! Next K plane

      END SUBROUTINE PUSHIT

C***********************************************************************
C
      SUBROUTINE READMULT (IER)
C
C     Read grid and functions files, all PLOT3D formats.
C
C     CALLED BY: INITIAL, MATCH
C
C     CALLS: HEADIO, MULTIO, REWND
C 
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (OUT) :: IER

C     Local variables:

      INTEGER  IG, LADD, LEND, LSUB, NEWPTS

C     Execution:

      IF (NITGO == 0) THEN
         NITGO = 10
         NITQO = 11
      ELSE
         CALL REWND
         NITGO = NITGO + 2
         NITQO = NITQO + 2
      END IF

      CALL HEADIO (1, IER)
      IF (IER == 1) GO TO 999

      IMAX = IM(MGRID)
      JMAX = JM(MGRID)
      KMAX = KM(MGRID)

      DO IG = 1, NGRID
         IMO(IG)  = IM(IG)
         IMQO(IG) = IMQ(IG)
         JMO(IG)  = JM(IG)
         JMQO(IG) = JMQ(IG)
         KMO(IG)  = KM(IG)
         KMQO(IG) = KMQ(IG)
      END DO

C     Will this zone change size?  If so, amend header record accordingly:

      IF (ADD /= 0 .OR. SUB /= 0 .OR. REMOVE /= 0) THEN
         NEWPTS = 0
         IF (IEND == 0) IEND = IM(MGRID)
         IF (JEND == 0) JEND = JM(MGRID)
         IF (KEND == 0) KEND = KM(MGRID)
         IF (ADD  /= 0) THEN
            LADD = LSTADD
            IF (LADD == 0) THEN
               IF (IJPLANE .AND. JSTEP .OR.
     >             IKPLANE .AND. KSTEP) LADD = IST
               IF (IJPLANE .AND. ISTEP .OR.
     >             JKPLANE .AND. KSTEP) LADD = JST
               IF (IKPLANE .AND. ISTEP .OR.
     >             JKPLANE .AND. JSTEP) LADD = KST
            END IF
            LEND = LENDADD
            IF (LEND == 0) THEN
               IF (IJPLANE .AND. JSTEP .OR.
     >             IKPLANE .AND. KSTEP) LEND = IEND
               IF (IJPLANE .AND. ISTEP .OR.
     >             JKPLANE .AND. KSTEP) LEND = JEND
               IF (IKPLANE .AND. ISTEP .OR.
     >             JKPLANE .AND. JSTEP) LEND = KEND
            END IF
            NEWPTS = NEWPTS + ABS ((LEND - LADD)*ADD)
         END IF
         IF (SUB /= 0) THEN
            LSUB = LSTSUB
            IF (LSUB == 0) THEN
               IF (IJPLANE .AND. JSTEP .OR.
     >             IKPLANE .AND. KSTEP) LADD = IST
               IF (IJPLANE .AND. ISTEP .OR.
     >             JKPLANE .AND. KSTEP) LADD = JST
               IF (IKPLANE .AND. ISTEP .OR.
     >             JKPLANE .AND. JSTEP) LADD = KST
            END IF
            LEND = LENDSUB
            IF (LENDSUB == 0) THEN
               IF (IJPLANE .AND. JSTEP .OR.
     >             IKPLANE .AND. KSTEP) LEND = IEND
               IF (IJPLANE .AND. ISTEP .OR.
     >             JKPLANE .AND. KSTEP) LEND = JEND
               IF (IKPLANE .AND. ISTEP .OR.
     >             JKPLANE .AND. JSTEP) LEND = KEND
            END IF
            NEWPTS = NEWPTS - ABS ((LEND - LADD)*SUB/(SUB+1))
         END IF
         IF (REMOVE /= 0) THEN
            NEWPTS = NEWPTS - REMOVE
         END IF
         IF (IJPLANE .AND. JSTEP .OR. IKPLANE .AND. KSTEP) THEN
            IMO(MGRID)  = NEWPTS + IM(MGRID)
            IMQO(MGRID) = NEWPTS + IMQ(MGRID)
         END IF
         IF (IJPLANE .AND. ISTEP .OR. JKPLANE .AND. KSTEP) THEN
            JMO(MGRID)  = NEWPTS + JM(MGRID)
            JMQO(MGRID) = NEWPTS + JMQ(MGRID)
         END IF
         IF (IKPLANE .AND. ISTEP .OR. JKPLANE .AND. JSTEP) THEN
            KMO(MGRID)  = NEWPTS + KM(MGRID)
            KMQO(MGRID) = NEWPTS + KMQ(MGRID)
         END IF
      END IF

C     Write header records and transcribe X,Y,Z and Q grids to
C     output files for any grids preceding the specified grid:

      CALL HEADIO (2, IER)
      CALL MULTIO (1, MGRID, 1, MGRID-1)

C     Requested grid MGRID is now in memory.

  999 CONTINUE

      END SUBROUTINE READMULT

C***********************************************************************
C
      SUBROUTINE REWND
C
C     Rewind the output file to become the input file and update unit #s.
C
C     CALLED BY: MATCH, READMULT
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

      NITGI = NITGO
      NITQI = NITQO
      CLOSE (NITGI)

      IF (FORMOO) THEN
         OPEN (NITGI, FORM='FORMATTED')
      ELSE
         OPEN (NITGI, FORM='UNFORMATTED')
      END IF

      IF (.NOT. NOQ) THEN
         CLOSE (NITQI)
         IF (FORMOO) THEN
            OPEN (NITQI, FORM='FORMATTED')
         ELSE
            OPEN (NITQI, FORM='UNFORMATTED')
         END IF
      END IF
      FORMI = FORMOO

      END SUBROUTINE REWND

C***********************************************************************
C
      SUBROUTINE SETUPJ (J, K)
C
C     Calculate the direction cosines of the vectors used in SOLUT that
C     are not functions of the iterations:  normal vector U from J - 1
C     line, normal vector B from J line, and straightness vector E.
C
C     CALLED BY: Main program
C
C     CALLS: ADDV, NORMPT, UNITV, VMERGE
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: J, K

C     Local variables:

      INTEGER  I, L, M
      REAL     CFX, CFY, CFZ, SING, COSD(IMX,3)

C     Execution:

C     Find U vector as average of normals to 4 planes at (I,J-1,K):

      DO I = 2, NIPTS - 1
         CALL NORMPT (I, 1, 2, 1, COSU(I,1), COSU(I,2), COSU(I,3), SING)
         IF (SING == ONE)
     >      CALL UNITV (XJ(I,1,2), YJ(I,1,2), ZJ(I,1,2),
     >                  XJ(I,2,2), YJ(I,2,2), ZJ(I,2,2),
     >                  COSU(I,1), COSU(I,2), COSU(I,3))
      END DO

C     Find B vector normal to J line similarly:

      DO I = 2, NIPTS - 1
         CALL NORMPT (I, 2, 2, 1, COSB(I,1), COSB(I,2), COSB(I,3), SING)
         IF (SING == ONE)
     >      CALL UNITV (XJ(I,1,2), YJ(I,1,2), ZJ(I,1,2),
     >                  XJ(I,2,2), YJ(I,2,2), ZJ(I,2,2),
     >                  COSB(I,1), COSB(I,2), COSB(I,3))
      END DO

C     Find D vector (used for straightness vector E). 
C     If 2nd J line and JST = 1, then no J-2 line exists.

      IF (J == 2) THEN
         DO I = 2, NIPTS - 1
            CALL UNITV (XJ(I,1,2), YJ(I,1,2), ZJ(I,1,2),
     >                  XJ(I,2,2), YJ(I,2,2), ZJ(I,2,2),
     >                  COSD(I,1), COSD(I,2), COSD(I,3))
         END DO
      ELSE
         L = IST
         DO I = 2, NIPTS - 1
            L = L + 1
            CALL UNITV (X(L,J-2,K), Y(L,J-2,K), Z(L,J-2,K),
     >                  X(L,J-1,K), Y(L,J-1,K), Z(L,J-1,K),
     >                  COSD(I,1),  COSD(I,2),  COSD(I,3))
            IF (COSD(I,1) == ZERO .AND. COSD(I,2) == ZERO .AND.
     >          COSD(I,3) == ZERO)
     >         CALL UNITV (XJ(I,1,2), YJ(I,1,2), ZJ(I,1,2),
     >                     XJ(I,2,2), YJ(1,2,2), ZJ(I,2,2),
     >                     COSD(I,1), COSD(I,2), COSD(I,3))
         END DO
      END IF

C     Find actual boundary line for edge conditions:

      CALL UNITV (X(IST,J-1,K),  Y(IST,J-1,K),  Z(IST,J-1,K),
     >            X(IST,J,K),    Y(IST,J,K),    Z(IST,J,K),
     >            COSD(1,1),     COSD(1,2),     COSD(1,3))
      CALL UNITV (X(IEND,J-1,K), Y(IEND,J-1,K), Z(IEND,J-1,K),
     >            X(IEND,J,K),   Y(IEND,J,K),   Z(IEND,J,K),
     >            COSD(NIPTS,1), COSD(NIPTS,2), COSD(NIPTS,3))
      CALL VMERGE (COSD, 1, NIPTS, IMX)

C     Before continuing, check that signs for normal agrees with D vector
C     (a possible check on right-handedness):

      IF (J == JST+1 .AND. K == KST) THEN
         DO M = 1, 3
            IF (ABS (COSU(2,M)) > .01 .AND. ABS (COSD(2,M)) > .01) THEN
               IF (COSU(2,M)*COSD(2,M) < ZERO) THEN
                  WRITE (6, '(/, 2A)') ' WARNING: Direction cosines ',
     >             'indicate coordinate system may not be right-handed.'
                  EXIT
               END IF
            END IF
         END DO
      END IF

C     Merge actual ends into normal vectors too:

      DO M = 1, 3
         COSU(1,M) = COSD(1,M)
         COSB(1,M) = COSD(1,M)
         COSE(1,M) = COSD(1,M)
         COSU(NIPTS,M) = COSD(NIPTS,M)
         COSB(NIPTS,M) = COSD(NIPTS,M)
         COSE(NIPTS,M) = COSD(NIPTS,M)
      END DO
      CALL VMERGE (COSU, 1, NIPTS, IMX)
      CALL VMERGE (COSB, 1, NIPTS, IMX)

C     Find E as average D (I+1,I,I-1) (at ends, E = D):

      DO I = 2, NIPTS - 1
         CALL ADDV (COSD(I-1,1), COSD(I-1,2), COSD(I-1,3), ONE,
     >              COSD(I,1),   COSD(I,2),   COSD(I,3),   ONE,
     >              CFX,         CFY,         CFZ)
         CALL ADDV (CFX,         CFY,         CFZ,         ONE,
     >              COSD(I+1,1), COSD(I+1,2), COSD(I+1,3), ONE,
     >              COSE(I,1),   COSE(I,2),   COSE(I,3))
      END DO

      END SUBROUTINE SETUPJ

C***********************************************************************
C
      SUBROUTINE SETUPK (J, K)
C
C     Calculate the directions cosines of the vectors used to compute
C     the torsion vector from the K-1 plane, normal vector U at (J,K-1),
C     normal vector B from (J,K), and straightness vector E.
C
C     CALLED BY: Main program
C
C     CALLS: ADDV, NORMPT, UNITV, VMERGE
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: J, K

C     Local variables:

      INTEGER  I, L, M
      REAL     CFX, CFY, CFZ, SING, COSD(IMX,3)

C     Execution:

C     Find U as average of normals to 4 planes at (I,J,K-1):

      DO I = 2, NIPTS - 1
         CALL NORMPT (I, 2, 1, 2, COSUK(I,1), COSUK(I,2), COSUK(I,3),
     >                SING)
      END DO

C     Find B vector normal to J line:

      DO I = 2, NIPTS - 1
         CALL NORMPT (I, 2, 2, 2, COSBK(I,1), COSBK(I,2), COSBK(I,3),
     >                SING)
      END DO

C     Carol:  I notice SING is not checked here.

C     Find D vector (used for straightness  vector E).
C     If 2nd K plane and KST = 1, then no K-2 line exists.

      IF (K == 2) THEN
         DO I = 2, NIPTS - 1
            CALL UNITV (XJ(I,2,1), YJ(I,2,1), ZJ(I,2,1),
     >                  XJ(I,2,2), YJ(I,2,2), ZJ(I,2,2),
     >                  COSD(I,1), COSD(I,2), COSD(I,3))
         END DO
      ELSE
         L = IST
         DO I = 2, NIPTS - 1
            L = L + 1
            CALL UNITV (X(L,J,K-2), Y(L,J,K-2), Z(L,J,K-2),
     >                  X(L,J,K-1), Y(L,J,K-1), Z(L,J,K-1),
     >                  COSD(I,1),  COSD(I,2),  COSD(I,3))
            IF (COSD(I,1) == ZERO .AND. COSD(I,2) == ZERO .AND.
     >          COSD(I,3) == ZERO) THEN
               CALL UNITV (XJ(I,2,1), YJ(I,2,1), ZJ(I,2,1),
     >                     XJ(I,2,2), YJ(I,2,2), ZJ(I,2,2),
     >                     COSD(I,1), COSD(I,2), COSD(I,3))
            END IF
         END DO
      END IF

C     At end, D is actual boundary line:

      CALL UNITV (X(IST,J,K-1),  Y(IST,J,K-1),  Z(IST,J,K-1),
     >            X(IST,J,K),    Y(IST,J,K),    Z(IST,J,K),
     >            COSD(1,1),     COSD(1,2),     COSD(1,3))
      CALL UNITV (X(IEND,J,K-1), Y(IEND,J,K-1), Z(IEND,J,K-1),
     >            X(IEND,J,K),   Y(IEND,J,K),   Z(IEND,J,K),
     >            COSD(NIPTS,1), COSD(NIPTS,2), COSD(NIPTS,3))
      CALL VMERGE (COSD, 1, NIPTS, IMX)       

C     Merge into normal vectors too:

      DO M = 1, 3
         COSUK(1,M) = COSD(1,M)
         COSBK(1,M) = COSD(1,M)
         COSEK(1,M) = COSD(1,M)
         COSUK(NIPTS,M) = COSD(NIPTS,M)
         COSBK(NIPTS,M) = COSD(NIPTS,M)
         COSEK(NIPTS,M) = COSD(NIPTS,M)
      END DO
      CALL VMERGE (COSUK, 1, NIPTS, IMX)
      CALL VMERGE (COSBK, 1, NIPTS, IMX)

C     Find E as average D (I+1,I,I-1) (at ends, E = D):

      DO I = 2, NIPTS - 1
         CALL ADDV (COSD(I-1,1), COSD(I-1,2), COSD(I-1,3), ONE,
     >              COSD(I,1),   COSD(I,2),   COSD(I-1,3), ONE,
     >              CFX,         CFY,         CFZ)
         CALL ADDV (CFX,         CFY,         CFZ,         ONE,
     >              COSD(I+1,1), COSD(I+1,2), COSD(I+1,3), ONE,
     >              COSEK(I,1),  COSEK(I,2),  COSEK(I,3))
      END DO

      END SUBROUTINE SETUPK

C***********************************************************************
C
      SUBROUTINE SHOCKLE (SNMAX, DSN, NSM, NSM1)
C
C     Find the new Smax for every plane based on shock location or user
C     request.  Smooth in both directions.
C
C     CALLED BY: PUSHIT
C
C     CALLS: FILTER, DXSP
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      REAL,    INTENT (INOUT) :: SNMAX(IMX,IMX), DSN
      INTEGER, INTENT (IN)    :: NSM, NSM1

C     Local variables:

      INTEGER  I, J, K, L, L0, L1, L2, NPTS, NWL
      REAL     AVESN1, DEL, DQMAX, DSMX, PERCENT, QVAL, RATIO, SMAX

      INTEGER, DIMENSION (IMX)     :: ISHST
      REAL,    DIMENSION (IMX)     :: DSC, DQDS, SC, SC1
      REAL,    DIMENSION (IMX+2)   :: STMP
      REAL,    DIMENSION (IMX,IMX) :: SOMAX

C     Execution:

      NPTS = IEND - IST + 1

      DO 600 K = KST, KEND

         DO 500 J = JST, JEND

C           Compute Si for this J, store in SC:

            SC(1) = ZERO
            L = IST
            DO I = 1, NPTS - 1
               DSC(I) = SQRT ((X(L+1,J,K) - X(L,J,K))**2 +
     >                        (Y(L+1,J,K) - Y(L,J,K))**2 +
     >                        (Z(L+1,J,K) - Z(L,J,K))**2)
               IF (DSC(I) == ZERO) THEN
                  WRITE (6, '(/, A, 3I4, A)')
     >               ' WARNING: Zero cell at', I, J, K,
     >               '.  SAGE will try to correct it.'
                  IF (I < NPTS - 1) THEN
                     DSC(I) = (SQRT ((X(L+2,J,K) - X(L+1,J,K))**2 +
     >                               (Y(L+2,J,K) - Y(L+1,J,K))**2 +
     >                               (Z(L+2,J,K) - Z(L+1,J,K))**2))*0.5
                  END IF
               END IF
               SC(I+1) = SC(I) + DSC(I)
               L = L + 1
            END DO
            SMAX = SC(NPTS)
            SOMAX(J,K) = SMAX

            IF (MVC1 /= MVC2) THEN
               CALL DXSP (MVC, IMAX, MVC1, MVC2, XMINP, XMAXP, J, K)
            ELSE
               MVC = MVC1
            END IF

            IF (MVBOUND == ZERO .AND. MVC /= ZERO) THEN
               SNMAX(J,K) = SOMAX(J,K)
               GO TO 450
            END IF

            IF (MVBOUND <= 998.) THEN

C              MVBOUND is either +/- p, where p is a percentage of line length,
C              or MVBOUND = 900., allowing linear variation of percentage:

               PERCENT = MVBOUND * 0.01

               IF (MVBOUND == 900.) THEN

                  IF (MVP1 == ZERO .AND. MVP2 == ZERO) THEN
                     WRITE (6, '(/, 2A)') ' MVBOUND = 900. but ',
     >                  'end percentages not input in MVP1 and MVP2.'
                  END IF

                  CALL DXSP (PERCENT, IMAX, MVP1, MVP2, XMINP, XMAXP,
     >                       J, K)
                  PERCENT = PERCENT * 0.01

               END IF

               DSMX = SMAX * PERCENT
               SNMAX(J,K) = SMAX + DSMX

            ELSE IF (MVBOUND >= 999.) THEN

C              Move is based on shock location.

               IF (NOQ) WRITE (6, '(/, A)')
     >            ' Cannot use MVBOUND >= 999. without a Q file.'
               IF (INDQ > NQDIM) WRITE (6, '(/, A)')
     >            ' Cannot use INDQ > number of Q in solution file.'

               IF (MVBOUND == 999.) THEN ! Need shock location of INDQ

C                 Find gradient of Qindq, DQDS:

                  L = IST
                  DO I = 1, NPTS - 1
                     DQDS(I) = (ABS (Q(L+1,J,K,INDQ) - Q(L,J,K,INDQ))) /
     >                         DSC(I)
                     L = K + 1
                  END DO
                  DQMAX = DQDS(1)
                  DO I = 2, NPTS - 1
                     IF (DQDS(I) > DQMAX) DQMAX = DQDS(I)
                  END DO
                  DO I = 1, NPTS - 1
                     DQDS(I) = DQDS(I) / DQMAX
                  END DO

C                 Where is the shock?

                  DO I = NPTS - 1, 2, -1
                     IF (ABS (DQDS(I) - DQDS(I-1)) > 0.0001) THEN
                        ISHST(J) = I
                        GO TO 405
                     END IF
                  END DO
                  ISHST(J) = IMAX

  405             CONTINUE

C*****            SNMAX(J,K) = SC(ISHST(J)) + SOMAX(J,K)*DSN/100.
C*****            SNMAX(J,K) = (1.0 + DSN/(IMAX-1))*SC(ISHST(J))
                  SNMAX(J,K) = (1.0 + DSN * 0.01)  *SC(ISHST(J))

C                 ISHST smooth?  Look at previous two max. locations. 
C                 If not, improve the new Smax for J-1.

                  IF (J > 2) THEN
                     L0 = ISHST(J)
                     L1 = ISHST(J-1)
                     L2 = ISHST(J-2)
                     IF ((L1 < L2 .AND. L1 < L0) .OR.
     >                   (L1 > L2 .AND. L1 > L0)) THEN
                        NWL = (L0 + L2) / 2
C****                   SNMAX(J-1,K) = SC1(NWL) + SOMAX(J-1,K)*DSN/100.
C****                   SNMAX(J-1,K) = (1.0 + DSN/(IMAX-1))*SC1(NWL)
                        SNMAX(J-1,K) = (1.0 + DSN * 0.01)  *SC1(NWL)
                     END IF
                  END IF
                  DO I = 1, NPTS
                     SC1(I) = SC(I)
                  END DO

               ELSE ! MVBOUND > 999.

C                 Boundary parallels a Q file contour defined by INDQ.

                  QVAL = MVBOUND - 1000. ! Contour value

                  L = IST + NPTS - 1
                  DO I = NPTS - 1, 1, -1
                     L = L - 1
                     IF ((QVAL <= Q(L+1,J,K,INDQ) .AND.
     >                    QVAL >  Q(L,J,K,INDQ))  .OR.
     >                   (QVAL >= Q(L+1,J,K,INDQ) .AND.
     >                    QVAL <  Q(L,J,K,INDQ))) THEN
                        RATIO = (QVAL - Q(L,J,K,INDQ)) /
     >                          (Q(L+1,J,K,INDQ) - Q(L,J,K,INDQ))
                        SNMAX(J,K) = SC(I) + (SC(I+1) - SC(I))*RATIO
C****                   SNMAX(J,K) = (ONE + DSN/(IMAX-1))*SNMAX(J,K)
                        SNMAX(J,K) = (ONE + DSN * 0.01)  *SNMAX(J,K)
C****                   SNMAX(J,K) = SNMAX(J,K) + SOMAX(J,K)*DSN*0.01
                        GO TO 430
                     END IF
                  END DO
                  IF (J == 1 .AND. K == 1) WRITE (6, '(/, A)')
     >               ' Input contour value not found on initial line.'
                  IF (J /= 1) SNMAX(J,K) = SNMAX(J-1,K)
                  IF (J == 1 .AND. K /= 1) SNMAX(J,K) = SNMAX(J,K-1)

  430             CONTINUE

               END IF ! End of MVBOUND == 999. or > 1000. branch

            END IF ! End of MVBOUND options

  450       CONTINUE
            IF (MVC /= ZERO) SNMAX(J,K) = SNMAX(J,K) + MVC

  500    CONTINUE ! Next J line

  600 CONTINUE ! Next K plane

C     Smooth SNMAX in both directions:

      IF ((.NOT. SINGRD .AND. NGRID > 1 .AND. NSM  > 0) .OR. 
     >    (.NOT. SINGRD .AND. NGRID > 1 .AND. NSM1 > 0)) 
     >   WRITE (6, '(/, A)')
     >      ' Check zonal boundaries: outer smoothing may mismatch.'

      IF (NSM /= 0) THEN
         DO K = KST, KEND
            DO J = JST, JEND
               STMP(J) = SNMAX(J,K)
            END DO

            CALL FILTER (JMAX, NSM, STMP)

            DO J = JST, JEND
               SNMAX(J,K) = STMP(J)
            END DO
         END DO
      END IF

C     If NOSESM input, remove the discontinuity caused by smoothing:

      IF (NOSESM /= 0) THEN
         AVESN1 = 0.5*(SNMAX(NOSESM+1,1) + SNMAX(NOSESM+1,KMAX))
         DO K = KST, KEND
            SNMAX(1,K) = AVESN1
            IF (NOSESM > 1) THEN
               DEL = (SNMAX(NOSESM+1,K) - SNMAX(1,K)) / REAL (NOSESM)
               DO L = 2, NOSESM
                  SNMAX(L,K) = SNMAX(L-1,K) + DEL
               END DO
            END IF
         END DO
      END IF

      IF (NSM1 /= 0) THEN
         DO J = JST, JEND
            DO K = KST, KEND
               STMP(K) = SNMAX(J,K)
            END DO

            CALL FILTER (KMAX, NSM1, STMP)

            DO K = KST, KEND
               SNMAX(J,K) = STMP(K)
            END DO
         END DO
      END IF

      END SUBROUTINE SHOCKLE

C***********************************************************************
C
      SUBROUTINE SINGPLN
C
C     When a plane is only a single line (but stored as a plane), each
C     line must have the same results.
C
C     CALLED BY: Main program
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

      INTEGER I, J

      DO J = JST + 1, JEND
         DO I = IST, IEND
C           Carol:   M = IST+I-1 looked wrong if (say) IST = 5.  Use I, not M.
            X(I,J,1) = X(I,JST,1)
            Y(I,J,1) = Y(I,JST,1)
            Z(I,J,1) = Z(I,JST,1)
            IF (.NOT. NOQ) Q(I,J,1,1:NQDIM) = Q(I,JST,1,1:NQDIM)
         END DO
      END DO

      END SUBROUTINE SINGPLN

C***********************************************************************
C
      SUBROUTINE SOLUT (J, K)
C
C     Solve the tridiagonal system of nonlinear equations that underly
C     the adaption technique, producing new S(i) for the current J line.
C
C     CALLED BY: Main program
C
C     CALLS: EDGEMG, FILTER, GETWT, INTF, NORM, WTEDGE
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: J, K

C     Local variables:

      INTEGER
     >   I, JERR, KERR, L, LOOP

      REAL
     >   ASPECT, DSL, ERR, SPC, SPCPL, TEMP

      REAL, DIMENSION (IMX) ::
     >   AA, BB, CC, FF, SOLD, TAU, TAUPL

C     Execution:

      DO I = 1, NIPTS          ! Initialize previous iterate
         SOLD(I) = SN(I)
      END DO

C     Calculate torsion spring constants, TAU, including aspect ratio of DAP:

      SPC   = (ONE + A)*CLAMW(1)
      SPCPL = (ONE + A)*CLAMW(2)

      DO I = 2, NIPTS - 1
         IF (J > JST) THEN
            ASPECT = 0.5*(SNM(I+1) - SNM(I-1)) / ABS (DAP(I))
            TAU(I) = SPC*ASPECT / ABS (DAP(I))
         ELSE
            TAU(I) = ZERO
         END IF
         IF (K > KST) THEN
            ASPECT = 0.5*(SNMK(I+1) - SNMK(I-1)) / ABS (DAPPL(I))
            TAUPL(I) = SPCPL*ASPECT / ABS (DAPPL(I))
         ELSE
            TAUPL(I) = ZERO
         END IF
      END DO

C     Newton-type iteration for new S along line J:

      DO 600 LOOP = 1, MAXITS

C        Calculate tension spring constant W:

         DO I = 1, NIPTS - 1
            WEIGHT(I) = ONE + A*FB(I)**B
         END DO

C        Override edge weights, if necessary:

         IF (NEDGE /= 0) THEN
            IF (WDS == ZERO .AND. WDE == ZERO) CALL WTEDGE (J-1, K)
            IF (MG1 /= 0) WEIGHT(1) = WDS
            IF (MG2 /= 0) WEIGHT(NIPTS-1) = WDE

            CALL EDGEMG (WEIGHT)
         END IF

C        Find mutiplier of WEIGHT for spacing control:

         IF (LOOP > 1) THEN
            CALL GETWT
            DO I = 1, NIPTS - 1
               WEIGHT(I) = WEIGHT(I)*WT(I)
            END DO
         END IF

C        Filter if specified:

         IF (NFILT > 0) CALL FILTER (NIPTS-1, NFILT, WEIGHT)

C        Set up coeffcients of S-1, S, S+1 and RHS:

         DO I = 2, NIPTS - 1
            AA(I) =   WEIGHT(I-1)
            BB(I) = -(WEIGHT(I) + WEIGHT(I-1) + TAU(I) + TAUPL(I))
            CC(I) =   WEIGHT(I)
            FF(I) = -TAU(I)*SP(I) - TAUPL(I)*SPPL(I)
         END DO

C        Correct first and last equations to reflect known end S'S:

         FF(2) = FF(2) - AA(2)*SN(1)
         AA(2) = ZERO
         FF(NIPTS-1) = FF(NIPTS-1) - CC(NIPTS-1)*SN(NIPTS)
         CC(NIPTS-1) = ZERO

C        Tridiagonal solution by LU factorization:

         DO I = 3, NIPTS - 1
            TEMP  = AA(I) / BB(I-1)
            BB(I) = BB(I) - CC(I-1)*TEMP
            FF(I) = FF(I) - FF(I-1)*TEMP
         END DO

C        Back substitute:

         SN(NIPTS-1) = FF(NIPTS-1) / BB(NIPTS-1)
         DO I = NIPTS - 2, 2, -1
            SN(I) = (FF(I) - CC(I)*SN(I+1)) / BB(I)
         END DO

C        Check SN for possible numeric error:

         DO I = NIPTS - 1, 2, -1
            IF (SN(I) < 1.0E-6) THEN
               DSL = SN(I+1) / REAL (I)
               DO L = 2, I
                  SN(L) = SN(L-1) + DSL
               END DO
               EXIT
            END IF
         END DO

         DO I = 1, NIPTS - 1
            DS(I) = SN(I+1) - SN(I)
         END DO

C        Re-evaluate FB at these SN:

         CALL INTF (INTER, NIPTS-1, SMS, SN, F, FB)
         CALL NORM (NIPTS-1, FB)

C        Convergence test:

         ERR = ZERO
         DO I = 2, NIPTS - 1
            ERR = ERR + ABS (SN(I) - SOLD(I))
         END DO
         ERR = ERR / SS(NIPTS)
         IF (ERR < CONV) GO TO 700 ! Converged

         DO I = 2, NIPTS
           SOLD(I) = SN(I)
         END DO

  600 CONTINUE ! Next iteration

      KERR = K
      IF (KINVERSE) KERR = KMAX - K + 1
      JERR = J
      IF (JINVERSE) JERR = JMAX - J + 1
      IF (ERR > CONV*100.0) WRITE (6, '(/, A, I4, A, I4, A, E11.3)')
     >   ' No convergence on line J =', JERR, ', plane K =', KERR,
     >   '.  ERR:', ERR

  700 CONTINUE

C     Proceed, converged or not.
C     Test that S is monotonically increasing.  If not, signal save & quit:

      DO I = 2, NIPTS
         IF (SN(I) <= SN(I-1)) THEN
            JERR = J
            KERR = K
            IF (JINVERSE) JERR = JMAX - J + 1
            IF (KINVERSE) KERR = KMAX - K + 1
            WRITE (6, '(/, A, I4, A, I4)')
     >         ' ** S is not monotonic on line', JERR, ', plane', KERR
            OK = .FALSE.
            EXIT
         END IF
      END DO

      END SUBROUTINE SOLUT

C***********************************************************************
C
      SUBROUTINE STOREX (NEX1, NEX2)
C
C     Transfer an "export" plane from the current grid to XP, YP, ZP, QP.
C
C     CALLED BY: MATCH
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (OUT) :: NEX1, NEX2 ! Explain

C     Local variables:

      INTEGER  I, J, K, IB, JB, KB, II, JJ, KK, MPLANEQ

C     Execution:

      IF (TWOD) MPLANE = 1
      MPLANEQ = MPLANE

      IF (IS == 0 .AND. TWOD) IS = IMAX
      IF (IS == 0 .AND. .NOT. TWOD) IS = 1
      IF (IE == 0) IE = IMAX
      IF (JS == 0) JS = 1
      IF (JE == 0) JE = JMAX
      IF (KS == 0) KS = 1
      IF (KE == 0) KE = KMAX

C     Allow for inverse transfers:

      IB = 1
      IF (IS > IE) IB = -1
      JB = 1
      IF (JS > JE) JB = -1
      KB = 1
      IF (KS > KE) KB = -1

      IF (JKPLANE) THEN
         IF (FV .AND. MPLANE == IMAX .AND. IMAX /= 1) MPLANEQ = MPLANE-1
         NEX1 = ABS (JE - JS) + 1
         NEX2 = ABS (KE - KS) + 1
         DO K = KS, KE, KB
            KK = (K - KS)*KB + 1
            DO J = JS, JE, JB
               JJ = (J - JS)*JB + 1
               XP(JJ,KK) = X(MPLANE,J,K)
               YP(JJ,KK) = Y(MPLANE,J,K)
               ZP(JJ,KK) = Z(MPLANE,J,K)
               IF (.NOT. NOQ) QP(JJ,KK,1:NQDIM) = Q(MPLANEQ,J,K,1:NQDIM)
            END DO
         END DO
      END IF

      IF (IJPLANE) THEN
         IF (FV .AND. MPLANE == KMAX .AND. KMAX /= 1) MPLANEQ = MPLANE-1
         NEX1 = ABS (IE - IS) + 1
         NEX2 = ABS (JE - JS) + 1
         DO J = JS, JE, JB
            JJ = (J - JS)*JB + 1
            DO I = IS, IE, IB
               II = (I - IS)*IB + 1
               XP(II,JJ) = X(I,J,MPLANE)
               YP(II,JJ) = Y(I,J,MPLANE)
               ZP(II,JJ) = Z(I,J,MPLANE)
               IF (.NOT. NOQ) QP(II,JJ,1:NQDIM) = Q(I,J,MPLANEQ,1:NQDIM)
            END DO
         END DO
      END IF

      IF (IKPLANE) THEN
         IF (FV .AND. MPLANE == JMAX .AND. JMAX /= 1) MPLANEQ = MPLANE-1
         NEX1 = ABS (IE - IS) + 1
         NEX2 = ABS (KE - KS) + 1
         DO K = KS, KE, KB
            KK = (K - KS)*KB + 1
            DO I = IS, IE, IB
               II = (I - IS)*IB + 1
               XP(II,KK) = X(I,MPLANE,K)
               YP(II,KK) = Y(I,MPLANE,K)
               ZP(II,KK) = Z(I,MPLANE,K)
               IF (.NOT. NOQ) QP(II,KK,1:NQDIM) = Q(I,MPLANEQ,K,1:NQDIM)
            END DO
         END DO
      END IF

      IS = 0
      IE = 0
      JS = 0
      JE = 0
      KS = 0
      KE = 0

      END SUBROUTINE STOREX

C***********************************************************************
C
      SUBROUTINE STORIM (NEX1, NEX2)
C
C     Transfer "export" data to the "import" plane of the current grid.
C
C     CALLED BY: MATCH
C
C***********************************************************************
C
      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: NEX1, NEX2 ! Explain

C     Local variables:

      INTEGER  I, J, K, IB, JB, KB, IDUM, IIMP, JIMP, KIMP, MPLANEQ,
     >         N1IM, N2IM, N3IM

C     Execution:

      IF (IS == 0) IS = 1
      IF (IE == 0 .AND. .NOT. TWOD) IE = IMAX
      IF (IE == 0 .AND. TWOD) IE = 1
      IF (JS == 0) JS = 1
      IF (JE == 0) JE = JMAX
      IF (KS == 0) KS = 1
      IF (KE == 0) KE = KMAX

      N1IM = ABS (IE - IS) + 1
      N2IM = ABS (JE - JS) + 1
      N3IM = ABS (KE - KS) + 1

C     Allow for inverse transfers:

      IB = 1
      IF (IS > IE) IB = -1
      JB = 1
      IF (JS > JE) JB = -1
      KB = 1
      IF (KS > KE) KB = -1

      MPLANEQ = MPLANE

      IF (JKPLANE) THEN
         IF (FV .AND. MPLANE == IMAX .AND. IMAX /= 1) MPLANEQ = MPLANE-1
         DO K = KS, KE, KB
            KIMP = (K - KS)*KB + 1
            DO J = JS, JE, JB
               JIMP = (J - JS)*JB + 1
               IF (N2IM /= NEX1) THEN
                  IF (N2IM == NEX2 .AND. N3IM == NEX1) THEN
                     IDUM = JIMP
                     JIMP = KIMP
                     KIMP = IDUM
                  ELSE
                     WRITE (6, 1000)
                     GO TO 999
                  END IF
               END IF
               X(MPLANE,J,K) = XP(JIMP,KIMP)
               Y(MPLANE,J,K) = YP(JIMP,KIMP)
               Z(MPLANE,J,K) = ZP(JIMP,KIMP)
               IF (.NOT. NOQ)
     >            Q(MPLANEQ,J,K,1:NQDIM) = QP(JIMP,KIMP,1:NQDIM)
            END DO
         END DO
      END IF

      IF (IJPLANE) THEN
         IF (FV .AND. MPLANE == KMAX .AND. KMAX /= 1) MPLANEQ = MPLANE-1
         DO J = JS, JE, JB
            JIMP = (J - JS)*JB + 1
            DO I = IS, IE, IB
               IIMP = (I - IS)*IB + 1
               IF (N1IM /= NEX1) THEN
                  IF (N1IM == NEX2 .AND. N2IM == NEX1) THEN
                     IDUM = IIMP
                     IIMP = JIMP
                     JIMP = IDUM
                  ELSE
                     WRITE (6, 1000)
                     GO TO 999
                  END IF
               END IF
               X(I,J,MPLANE) = XP(IIMP,JIMP)
               Y(I,J,MPLANE) = YP(IIMP,JIMP)
               Z(I,J,MPLANE) = ZP(IIMP,JIMP)
               IF (.NOT. NOQ)
     >            Q(I,J,MPLANEQ,1:NQDIM) = QP(IIMP,JIMP,1:NQDIM)
            END DO
         END DO
      END IF

      IF (IKPLANE) THEN
         IF (FV .AND. MPLANE == JMAX .AND. JMAX /= 1) MPLANEQ = MPLANE-1
         DO K = KS, KE, KB
            KIMP = (K - KS)*KB + 1
            DO I = IS, IE, IB
               IIMP = (I - IS)*IB + 1
               IF (N1IM /= NEX1) THEN
                  IF (N1IM == NEX2 .AND. N3IM == NEX1) THEN
                     IDUM = IIMP
                     IIMP = KIMP
                     KIMP = IDUM
                  ELSE
                     WRITE (6, 1000)
                     GO TO 999
                  END IF
               END IF
               X(I,MPLANE,K) = XP(IIMP,KIMP)
               Y(I,MPLANE,K) = YP(IIMP,KIMP)
               Z(I,MPLANE,K) = ZP(IIMP,KIMP)
               IF (.NOT. NOQ)
     >            Q(I,MPLANEQ,K,1:NQDIM) = QP(IIMP,KIMP,1:NQDIM)
            END DO
         END DO
      END IF

  999 CONTINUE

 1000 FORMAT (/, ' ** STORIM:  Error in matching directions. **')

      END SUBROUTINE STORIM

C***********************************************************************
C
      SUBROUTINE SUBPTS
C
C     SUBPTS decreases the number of points along the adaption line.
C
C     CALLED BY: INITIAL
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Local variables:

      INTEGER  I, J, K, L, M, IMAXB, L1, NUMSUB
      REAL     XT(IMX), YT(IMX), ZT(IMX), QT(IMX,NDIM)

C     Execution:

C     Initialize start and end limits to be consistent with adaption domain:

      IF (LSTSUB  == 0) LSTSUB = IST
      IF (LENDSUB == 0) LENDSUB = IEND
      IF (LSTSUB > LENDSUB) THEN
         L1 = LENDSUB
         LENDSUB = LSTSUB
         LSTSUB = L1
      END IF
      IF (LSTSUB > LENDSUB) THEN
         WRITE (6, '(/, A)') ' No points subtracted.'
         GO TO 999
      END IF
      NUMSUB = (LENDSUB - LSTSUB) * SUB / (SUB + 1)
      IMAXB  = IMAX
      IMAX   = IMAX - NUMSUB
      IF (IMAX < 10) WRITE (6, '(/, A, I2)')
     >   ' SUB option produces less than 10 points:', IMAX

      DO 600 K = 1, KMAX

         DO 500 J = 1, JMAX

C           Move over first region:

            DO I = 1, LSTSUB
               XT(I) = X(I,J,K)
               YT(I) = Y(I,J,K)
               ZT(I) = Z(I,J,K)
               IF (.NOT. NOQ) QT(I,1:NQDIM) = Q(I,J,K,1:NQDIM)
            END DO

C           Subtract points in requested region:

            DO I = 1, NUMSUB
               M = LSTSUB + I*(SUB + 1)
               L = LSTSUB + I
               IF (M >= LENDSUB) CYCLE
               XT(L) = X(M,J,K)
               YT(L) = Y(M,J,K)
               ZT(L) = Z(M,J,K)
               IF (.NOT. NOQ) QT(L,1:NQDIM) = Q(M,J,K,1:NQDIM)
            END DO

C           Move over last region:

            DO I = LENDSUB, IMAXB
               L = I - NUMSUB
               XT(L) = X(I,J,K)
               YT(L) = Y(I,J,K)
               ZT(L) = Z(I,J,K)
               IF (.NOT. NOQ) QT(L,1:NQDIM) = Q(I,J,K,1:NQDIM)
            END DO

C           Return data in original arrays:

            DO I = 1, IMAX
               X(I,J,K) = XT(I)
               Y(I,J,K) = YT(I)
               Z(I,J,K) = ZT(I)
               IF (.NOT. NOQ) Q(I,J,K,1:NQDIM) = QT(I,1:NQDIM)
            END DO

  500    CONTINUE ! Next J line

  600 CONTINUE ! Next K plane

      IEND = IEND - NUMSUB
      WRITE (6, '(/, A, I4, A, I4)')
     >   ' Number of points decreased from', IMAXB, ' to', IMAX

999   CONTINUE

      END SUBROUTINE SUBPTS

C***********************************************************************
C
      SUBROUTINE SWAPINV
C
C     SWAPINV swaps data around to provide consistent ordering within
C     the code, regardless of input demands.  The code requires that
C     IST < IEND, JST < JEND, and KST < KEND.  Hence, data need to be
C     swapped at the start and reswapped at the end of the run.
C
C     Carol:  Can you compare and contrast with SWAPXYZ?
C
C     CALLED BY: INITIAL, OUTPUT
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Local variables:

      INTEGER I, J, K, II, JJ, KK
      REAL    XT(IMX), YT(IMX), ZT(IMX), QT(IMX,NDIM)

C     Execution:

      IF (IINVERSE) THEN ! I swap
         NFLAG = -NFLAG
         DO K = 1, KMAX
            DO J = 1, JMAX
               II = IMAX
               DO I = 1, IMAX
                  XT(I) = X(II,J,K)
                  YT(I) = Y(II,J,K)
                  ZT(I) = Z(II,J,K)
                  IF (.NOT. NOQ) QT(I,1:NQDIM) = Q(II,J,K,1:NQDIM)
                  II = II - 1
               END DO
               DO I = 1, IMAX
                  X(I,J,K) = XT(I)
                  Y(I,J,K) = YT(I)
                  Z(I,J,K) = ZT(I)
                  IF (.NOT. NOQ) Q(I,J,K,1:NQDIM) = QT(I,1:NQDIM)
               END DO
            END DO
         END DO
         IST  = IMAX - IST  + 1
         IEND = IMAX - IEND + 1
      END IF

      IF (JINVERSE) THEN ! J swap
         NFLAG = -NFLAG
         DO K = 1, KMAX
            DO I = 1, IMAX
               JJ = JMAX
               DO J = 1, JMAX
                  XT(J) = X(I,JJ,K)
                  YT(J) = Y(I,JJ,K)
                  ZT(J) = Z(I,JJ,K)
                  IF (.NOT. NOQ) QT(J,1:NQDIM) = Q(I,JJ,K,1:NQDIM)
                  JJ = J - 1
               END DO
               DO J = 1, JMAX
                  X(I,J,K) = XT(J)
                  Y(I,J,K) = YT(J)
                  Z(I,J,K) = ZT(J)
                  IF (.NOT. NOQ) Q(I,J,K,1:NQDIM) = QT(J,1:NQDIM)
               END DO
            END DO
         END DO
         JST  = JMAX - JST  + 1
         JEND = JMAX - JEND + 1
         IF (LNSING /= 0) LNSING = JMAX - LNSING + 1
      END IF

      IF (KINVERSE) THEN ! K swap
         NFLAG = -NFLAG
         DO J = 1, JMAX
            DO I = 1, IMAX
               KK = KMAX
               DO K = 1, KMAX
                  XT(K) = X(I,J,KK)
                  YT(K) = Y(I,J,KK)
                  ZT(K) = Z(I,J,KK)
                  IF (.NOT. NOQ) QT(K,1:NQDIM) = Q(I,J,KK,1:NQDIM)
               END DO
               DO K = 1, KMAX
                  X(I,J,K) = XT(K)
                  Y(I,J,K) = YT(K)
                  Z(I,J,K) = ZT(K)
                  IF (.NOT. NOQ) Q(I,J,K,1:NQDIM) = QT(K,1:NQDIM)
               END DO
            END DO
         END DO
         KST  = KMAX - KST  + 1
         KEND = KMAX - KEND + 1
      END IF                     

      END SUBROUTINE SWAPINV

C***********************************************************************
C
      SUBROUTINE SWAPXYZ (RSWAP)
C
C     SAGE requires adapting in the I direction, stepping along J lines
C     and marching up K planes.  Hence, SWAPXYZ permutes indices for
C     X, Y, Z, Q (inefficiently, to save temporary storage).
C
C     Carol:  How come it also swaps X,Y,Z coordinates, not just indices???
C
C     CALLED BY: INITIAL, OUTPUT
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      LOGICAL, INTENT (IN) :: RSWAP ! Explain

C     Local variables:

      INTEGER  I, J, K, IMAXT, JMAXT, KMAXT, ISTT, JSTT, KSTT,
     >         IENDT, JENDT, KENDT, MDIR, N

C     Execution:

      IF (IJPLANE .AND. ISTEP) MDIR = 1
      IF (IKPLANE .AND. KSTEP) MDIR = 2
      IF (JKPLANE .AND. JSTEP) MDIR = 3
      IF (JKPLANE .AND. KSTEP .AND. RSWAP .OR.
     >    IKPLANE .AND. ISTEP .AND. .NOT. RSWAP) MDIR = 4
      IF (IKPLANE .AND. ISTEP .AND. RSWAP .OR.
     >    JKPLANE .AND. KSTEP .AND. .NOT. RSWAP) MDIR = 5

C     Store grid size into temporaries:

      IMAXT = IMAX
      JMAXT = JMAX
      KMAXT = KMAX
      ISTT  = IST
      IENDT = IEND
      JSTT  = JST
      JENDT = JEND
      KSTT  = KST
      KENDT = KEND

C     Swap indices (and coordinates), starting with Q:

      IF (.NOT. NOQ) THEN

         DO N = 1, NQDIM

            DO K = 1, KMAX
               DO J = 1, JMAX
                  DO I = 1, IMAX
                     T1(I,J,K) = Q(I,J,K,N)
                  END DO
               END DO
            END DO

            DO K = 1, KMAX
               DO J = 1, JMAX
                  DO I = 1, IMAX
                     SELECT CASE (MDIR)
                     CASE (1)
                        Q(J,I,K,N) = T1(I,J,K)
                     CASE (2)
                        Q(I,K,J,N) = T1(I,J,K)
                     CASE (3)
                        Q(K,J,I,N) = T1(I,J,K)
                     CASE (4)
                        Q(K,I,J,N) = T1(I,J,K)
                     CASE (5)
                        Q(J,K,I,N) = T1(I,J,K)
                     END SELECT
                  END DO
               END DO
            END DO

         END DO ! Next N

      END IF

      SELECT CASE (MDIR)

      CASE (1)

         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  T1(I,J,K) = Z(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  Z(J,I,K) = T1(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  T1(I,J,K) = X(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  X(J,I,K) = Y(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  Y(J,I,K) = T1(I,J,K)
               END DO
            END DO
         END DO

      CASE (2)

         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  T1(I,J,K) = X(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  X(I,K,J) = T1(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  T1(I,J,K) = Y(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  Y(I,K,J) = Z(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  Z(I,K,J) = T1(I,J,K)
               END DO
            END DO
         END DO

      CASE (3)

         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  T1(I,J,K) = Y(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  Y(K,J,I) = T1(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  T1(I,J,K) = X(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  X(K,J,I) = Z(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  Z(K,J,I) = T1(I,J,K)
               END DO
            END DO
         END DO

      CASE (4)

         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  T1(I,J,K) = X(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  X(K,I,J) = Z(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  Z(K,I,J) = Y(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  Y(K,I,J) = T1(I,J,K)
               END DO
            END DO
         END DO

      CASE (5)

         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  T1(I,J,K) = X(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  X(J,K,I) = Y(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  Y(J,K,I) = Z(I,J,K)
               END DO
            END DO
         END DO
         DO K = 1, KMAX
            DO J = 1, JMAX
               DO I = 1, IMAX
                  Z(J,K,I) = T1(I,J,K)
               END DO
            END DO
         END DO

      END SELECT

      IF (MDIR == 1 .OR. MDIR == 5) THEN
         IMAX = JMAXT
         IST  = JSTT
         IEND = JENDT
      END IF
      IF (MDIR == 3 .OR. MDIR == 4) THEN
         IMAX = KMAXT
         IST  = KSTT
         IEND = KENDT
      END IF
      IF (MDIR == 1 .OR. MDIR == 4) THEN
         JMAX = IMAXT
         JST  = ISTT
         JEND = IENDT
      END IF
      IF (MDIR == 2 .OR. MDIR == 5) THEN
         JMAX = KMAXT
         JST  = KSTT
         JEND = KENDT
      END IF
      IF (MDIR == 2 .OR. MDIR == 4) THEN
         KMAX = JMAXT
         KST  = JSTT
         KEND = JENDT
      END IF
      IF (MDIR == 3 .OR. MDIR == 5) THEN
         KMAX = IMAXT
         KST  = ISTT
         KEND = IENDT
      END IF

      END SUBROUTINE SWAPXYZ

C***********************************************************************
C
      SUBROUTINE SWAP2D
C
C     Turn a 2-D case into a single-plane 3-D case at Z = 0.
C
C     CALLED BY: INITIAL, OUTPUT
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

      IKPLANE = .FALSE.
      JKPLANE = .FALSE.
      IJPLANE = .TRUE.
      KMAX = 1
      KEND = 1
      IF (.NOT. JSTEP) ISTEP = .TRUE.
      Z(1:IMAX, 1:JMAX, 1) = ZERO

      END SUBROUTINE SWAP2D

C***********************************************************************
C
      SUBROUTINE TORCOF (L, JK, JKST, JKEND, MGNOS, MARCHJK)
C
C     Compute the coefficients that are used to calculate the two
C     torsion vectors.  L = 1 implies from J-1; L = 2 implies from K-1.
C
C     CALLED BY: Main program
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: L, JK, JKST, JKEND, MGNOS
      LOGICAL, INTENT (IN) :: MARCHJK

C     Local variables:

C     Execution:

C     Set up constants for the indicated torsion coeffcient.
C     If merging region, compute proportional coeff (CNM).
C     Enforce more orthogonality at walls (can be overridden by input):
C        if at either wall boundary, increase CLAM;
C        if at a start wall, increase U;
C        if at end boundary, increase B.

      CLAMW(L) = CLAM(L)
      CTM(L)   = CT(L)
      TNM(L)   = TN(L)
      CNM(L)   = ZERO

      IF (JK < MGNOS) THEN
         CNM(L)   = REAL (MGNOS - JK) / REAL (MGNOS - JKST)
         CLAMW(L) = CLAM(L)*(ONE + 5.0*CNM(L))
      ELSE
         IF (ORTHS(L)) THEN
            IF (JK == 2) THEN
               TNM(L)   = ZERO
               CTM(L)   = 0.25*CT(L)
               CLAMW(L) = 10.0*CLAM(L)
            ELSE IF (JK == 3) THEN
               TNM(L)   = 0.5*TN(L)
               CTM(L)   = 0.5*CT(L)
               CLAMW(L) = 3.0*CLAM(L)
            END IF 
         END IF
      END IF

C     Increase CLAM at outer wall:

      IF (.NOT. MARCHJK) THEN
         IF (ORTHE(L)) THEN
            IF (JK == JKEND-2) THEN
               CLAMW(L) = 2.0*CLAM(L)
            ELSE IF (JK == JKEND-1) THEN
               CLAMW(L) = 5.0*CLAM(L)
               TNM(L)   = TN(L) + (ONE - TN(L))/3.0
               CTM(L)   = 0.5*CT(L)
            ELSE IF (JK == JKEND) THEN
               CLAMW(L) = 10.0*CLAM(L)
               TNM(L)   = TN(L) + (ONE - TN(L))*(2.0/3.0)
               CTM(L)   = 0.25*CT(L)
            END IF
         END IF
      END IF

      END SUBROUTINE TORCOF

C***********************************************************************
C
      SUBROUTINE TORSION (J, K)
C
C     Calculate bothe torsion vectors T as functions of vectors U, B, E,
C     and hence S' for (J-1,K) and S* for (J,K-1) required in the main
C     adaption equation.
C
C     CALLED BY: Main program
C
C     CALLS: ADDV, CROSSV
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: J, K

C     Local variables:

      INTEGER  I, L, M
      REAL     ACT, ACT1, ACT2, ACT21, SPT, SWFL, TNM1, TNM2

      INTEGER, DIMENSION (IMX)   :: ICROSS
      REAL,    DIMENSION (IMX)   :: AAP, DUM2
      REAL,    DIMENSION (IMX,3) :: COSB1K, COST, COSB1, COSN

C     Execution:

      DO I = 1, NIPTS - 1
         DO M = 1, 3
            COSB1(I,M)  = COSB(I,M)
            COSB1K(I,M) = COSBK(I,M)
         END DO
      END DO

      IF (J == JST) THEN

         SP(1:NIPTS) = ZERO

      ELSE

C        Find the index ICROSS(*) of the DS interval intersecting each COSB:

         CALL CROSSV (XJ(1,2,2), YJ(1,2,2), ZJ(1,2,2),
     >                XJ(1,1,2), YJ(1,1,2), ZJ(1,1,2),
     >                DS, COSB1, AAP, DUM2, ICROSS, J)

C        Choose correct unit normal vector to use for each I point at J-1:

         DO I = 2, NIPTS - 1
            L = ICROSS(I)
            IF ((AAP(I) < DS(L)*0.5 .AND. L /= 1) .OR. L == NIPTS-1)
     >         L = L-1
            CALL ADDV (COSB1(L,1),   COSB1(L,2),   COSB1(L,3),   0.5,
     >                 COSB1(L+1,1), COSB1(L+1,2), COSB1(L+1,3), 0.5,
     >                 COSB(I,1),    COSB(I,2),    COSB(I,3))
         END DO

C        Add U and B to give normal vector N:

         TNM1 = ONE - TNM(1)
         DO I = 2, NIPTS - 1
            CALL ADDV (COSU(I,1), COSU(I,2), COSU(I,3), TNM1,
     >                 COSB(I,1), COSB(I,2), COSB(I,3), TNM(1),
     >                 COSN(I,1), COSN(I,2), COSN(I,3))
         END DO

C        Find T by adding N and E; proportion by CT:

         ACT  = CTM(1)*(ONE - CNM(1)) + CNM(1)
         ACT1 = ONE - ACT
         DO I = 2, NIPTS - 1
            CALL ADDV (COSE(I,1), COSE(I,2), COSE(I,3), ACT,
     >                 COSN(I,1), COSN(I,2), COSN(I,3), ACT1,
     >                 COST(I,1), COST(I,2), COST(I,3))
         END DO

C        Compute S':

         CALL CROSSV (XJ(1,2,2), YJ(1,2,2), ZJ(1,2,2),
     >                XJ(1,1,2), YJ(1,1,2), ZJ(1,1,2),
     >                DS, COST, AAP, DAP, ICROSS, J) 
         SP(1) = ZERO
         SP(NIPTS) = SN(NIPTS)
         DO I = 2, NIPTS - 1
            SP(I) = SN(ICROSS(I)) + AAP(I)
         END DO

C        Check for any cross-over of S'.
C        Try to correct for minor errors, else hope for the best!

         DO L = 1, 10
            SWFL = ZERO
            DO I = 1, NIPTS - 1
               IF (SP(I+1) < SP(I)) THEN
                  SPT = SP(I+1)
                  SP(I+1) = SP(I)
                  SP(I) = SPT
                  SWFL = ONE
               END IF
            END DO
            IF (SWFL == ZERO) EXIT
         END DO

      END IF ! End dealing with the J direction

C     Likewise for SPPL in the K direction:

      IF (K == KST) THEN

         SPPL(1:NIPTS) = ZERO

      ELSE     

         CALL CROSSV (XJ(1,2,2), YJ(1,2,2), ZJ(1,2,2),
     >                XJ(1,2,1), YJ(1,2,1), ZJ(1,2,1),
     >                DS, COSB1K, AAP, DUM2, ICROSS, J)

         DO I = 2, NIPTS - 1
            L = ICROSS(I)
            IF ((AAP(I) < DS(L)*0.5 .AND. L /= 1) .OR. L == NIPTS - 1)
     >          L = L - 1
            CALL ADDV (COSB1K(L,1),   COSB1K(L,2),   COSB1K(L,3),   0.5,
     >                 COSB1K(L+1,1), COSB1K(L+1,2), COSB1K(L+1,3), 0.5,
     >                 COSBK(I,1),    COSBK(I,2),    COSBK(I,3))
         END DO

         TNM2 = ONE - TNM(2)
         DO I = 2, NIPTS - 1
            CALL ADDV (COSUK(I,1), COSUK(I,2), COSUK(I,3), TNM2,
     >                 COSBK(I,1), COSBK(I,2), COSBK(I,3), TNM(2),
     >                 COSN(I,1),  COSN(I,2),  COSN(I,3))   
         END DO

         ACT2  = CTM(2)*(ONE - CNM(2)) + CNM(2)
         ACT21 = ONE - ACT2
         DO I = 2, NIPTS - 1
            CALL ADDV (COSEK(I,1), COSEK(I,2), COSEK(I,3), ACT2,
     >                 COSN(I,1),  COSN(I,2),  COSN(I,3),  ACT21,
     >                 COST(I,1),  COST(I,2),  COST(I,3))
         END DO

C        Compute S* for the plane:

         CALL CROSSV (XJ(1,2,2), YJ(1,2,2), ZJ(1,2,2),
     >                XJ(1,2,1), YJ(1,2,1), ZJ(1,2,1),
     >                DS, COST, AAP, DAPPL, ICROSS, J) 
         SPPL(1) = ZERO
         SPPL(NIPTS) = SN(NIPTS)
         DO I = 2, NIPTS - 1
            SPPL(I) = SN(ICROSS(I)) + AAP(I)
         END DO

C        Similarly check for any cross-over of S* plane:

         DO L = 1, 10
            SWFL = ZERO
            DO I = 1, NIPTS - 1
               IF (SPPL(I+1) < SPPL(I)) THEN
                  SPT = SPPL(I+1)
                  SPPL(I+1) = SPPL(I)
                  SPPL(I) = SPT
                  SWFL = ONE
               END IF
            END DO
            IF (SWFL == ZERO) EXIT
         END DO

      END IF

      END SUBROUTINE TORSION

C***********************************************************************
C
      SUBROUTINE UNITV (X1, Y1, Z1, X2, Y2, Z2, DIRCX, DIRCY, DIRCZ)
C
C     Calculate the direction cosines of the vector joinint two points.
C
C     CALLED BY: CROSSV, NORMPT, SETUPJ, SETUPK
C
C***********************************************************************

      IMPLICIT NONE

      REAL, INTENT (IN)  :: X1, Y1, Z1, X2, Y2, Z2
      REAL, INTENT (OUT) :: DIRCX, DIRCY, DIRCZ

      REAL  VMOD

C     Execution:

      DIRCX = X2 - X1
      DIRCY = Y2 - Y1
      DIRCZ = Z2 - Z1
      VMOD  = SQRT (DIRCX**2 + DIRCY**2 + DIRCZ**2)
      IF (VMOD /= 0.0) THEN
         DIRCX = DIRCX / VMOD
         DIRCY = DIRCY / VMOD
         DIRCZ = DIRCZ / VMOD
      END IF

      END SUBROUTINE UNITV

C***********************************************************************
C
      SUBROUTINE UPDATE (J, K)
C
C     This routine is called after a J line has converged.
C     It re-evaluates X, Y, Z, Q at the converged arc lengths.
C
C     CALLED BY: Main program
C
C     CALLS: INTXYZQ
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: J, K

C     Local variables:

      INTEGER  I, L, M
      REAL     QJ(IMX,NDIM)

C     Execution:

C     Do the interpolation:

      CALL INTXYZQ (J, K, 2, 2, SS, SN, QJ)

C     Update the original grid:

      L = IST
      DO I = 2, NIPTS - 1
         L = L + 1
         X(M,J,K) = XJ(I,2,2)
         Y(M,J,K) = YJ(I,2,2)
         Z(M,J,K) = ZJ(I,2,2)
         IF (.NOT. NOQ) Q(M,J,K,1:NQDIM) = QJ(I,1:NQDIM)
      END DO

      SNM(1:NIPTS) = SN(1:NIPTS)

      END SUBROUTINE UPDATE

C***********************************************************************
C
      SUBROUTINE VMERGE (DIRV, LST, LEND, IMX)
C
C     At each end of a J line, merge towards the end vector (direction
C     cosines) by blending it into the three nearest vectors.
C
C     CALLED BY: SETUPJ, SETUPK
C
C     CALLS: ADDV
C
C***********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: LST, LEND, IMX
      REAL, INTENT (INOUT) :: DIRV(IMX,3)

      REAL, PARAMETER :: PT25 = 0.25, PT50 = 0.50, PT75 = 0.75

      INTEGER L

      IF (LST /= 0) THEN
         L = LST
         CALL ADDV (DIRV(L,1),   DIRV(L,2),   DIRV(L,3),   PT25,
     >              DIRV(L+3,1), DIRV(L+3,2), DIRV(L+3,3), PT75,
     >              DIRV(L+3,1), DIRV(L+3,2), DIRV(L+3,3))
         CALL ADDV (DIRV(L,1),   DIRV(L,2),   DIRV(L,3),   PT50,
     >              DIRV(L+2,1), DIRV(L+2,2), DIRV(L+2,3), PT50,
     >              DIRV(L+2,1), DIRV(L+2,2), DIRV(L+2,3))
         CALL ADDV (DIRV(L,1),   DIRV(L,2),   DIRV(L,3),   PT75,
     >              DIRV(L+1,1), DIRV(L+1,2), DIRV(L+1,3), PT25,
     >              DIRV(L+1,1), DIRV(L+1,2), DIRV(L+1,3))
      END IF

      IF (LEND /= 0) THEN
         L = LEND
         CALL ADDV (DIRV(L,1),   DIRV(L,2),   DIRV(L,3),   PT75,
     >              DIRV(L-1,1), DIRV(L-1,2), DIRV(L-1,3), PT25,
     >              DIRV(L-1,1), DIRV(L-1,2), DIRV(L-1,3))
         CALL ADDV (DIRV(L,1),   DIRV(L,2),   DIRV(L,3),   PT50,
     >              DIRV(L-2,1), DIRV(L-2,2), DIRV(L-2,3), PT50,
     >              DIRV(L-2,1), DIRV(L-2,2), DIRV(L-2,3))
         CALL ADDV (DIRV(L,1),   DIRV(L,2),   DIRV(L,3),   PT25,
     >              DIRV(L-3,1), DIRV(L-3,2), DIRV(L-3,3), PT75,
     >              DIRV(L-3,1), DIRV(L-3,2), DIRV(L-3,3))
      END IF

      END SUBROUTINE VMERGE

C***********************************************************************
C
      SUBROUTINE WALLS (JW, K)
C
C     Calculate geometry functions (radii of curvature) at the boundary
C     line indicated by JW and K.
C
C     CALLED BY: FBAR
C
C     CALLS: CSFIT, CSDVAL
C
C***********************************************************************

      USE SAGEMOD

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: JW, K

C     Local variables:

      INTEGER  I, IER, L, NEVAL
      REAL     DSS
      REAL, DIMENSION (IMX) ::
     >   BCOF, CCOF, DCOF, SSS,
     >   XI, XPI, XPPI, YI, YPI, YPPI, ZI, ZPI, ZPPI

C     Execution:

C     Arc-lengths along the line:

      L      = IST
      NEVAL  = NIPTS - 1
      SSS(1) = ZERO

      DO I = 1, NEVAL
         DSS = SQRT ((X(L+1,JW,K) - X(L,JW,K))**2 +
     >               (Y(L+1,JW,K) - Y(L,JW,K))**2 +
     >               (Z(L+1,JW,K) - Z(L,JW,K))**2)
         SSS(I+1) = SSS(I) + DSS
         SMSS(I) = (SSS(I+1) + SSS(I))*0.5       
         L = L + 1
      END DO

C     Interpolate X, Y, Z arc-length derivatives at mid-points:

      CALL CSFIT  (NIPTS, SSS, X(IST,JW,K), 0, ZERO, 0, ZERO,
     >             BCOF, CCOF, DCOF, IER)
      CALL CSDVAL (NIPTS, SSS, X(IST,JW,K), NEVAL, SMSS,
     >             BCOF, CCOF, DCOF, XI, XPI, XPPI)

      CALL CSFIT  (NIPTS, SSS, Y(IST,JW,K), 0, ZERO, 0, ZERO,
     >             BCOF, CCOF, DCOF, IER)
      CALL CSDVAL (NIPTS, SSS, Y(IST,JW,K), NEVAL, SMSS,
     >             BCOF, CCOF, DCOF, YI, YPI, YPPI)

      CALL CSFIT  (NIPTS, SSS, Z(IST,JW,K), 0, ZERO, 0, ZERO,
     >             BCOF, CCOF, DCOF, IER)
      CALL CSDVAL (NIPTS, SSS, Z(IST,JW,K), NEVAL, SMSS,
     >             BCOF, CCOF, DCOF, ZI, ZPI, ZPPI)

C     Carol:  SMSS was not being used originally.
C             We now include Z derivatives in the curvature formula.

C     Radii of curvature, safeguarded against zero curvature:

      DO I = 1, NEVAL
         FGS(I) = ONE / MAX (1.E-10,
     >                  SQRT (XPPI(I)**2 + YPPI(I)**2 + ZPPI(I)**2))
      END DO

C     Locate the maximum radius of curvature for the aspect ratio calculation:

      MXFG = 2      ! Carol:  Why is I = 1 avoided here?
      DO I = 3, NEVAL
         IF (FGS(I) > FGS(MXFG)) MXFG = I
      END DO

      END SUBROUTINE WALLS

C***********************************************************************
C
      SUBROUTINE WTEDGE (J, K)
C
C     Force the edge weighting function to be a function of DLENGS and
C     DLENGE to provide better continuity by overriding FB.
C     Watch for a singular line - WDE and WDS need to be the same on the
C     next line.
C
C     CALLED BY: MAIN, LINE1, SOLUT
C
C     CALLS: DLENG
C
C***********************************************************************

      USE SAGEMOD    

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: J, K

C     Local variables:

      INTEGER I
      REAL    WDMAX, WTSUM

C     Execution:

      IF (J == JEND) THEN
         WDS = ZERO
         WDE = ZERO
         GO TO 900
      END IF

      IF (J /= 0 .AND. LNSING == J .AND. K /= KST) THEN
         WDE = WDES
         WDS = WDSS
         GO TO 900
      END IF

C     Determine function of DLENG for next J line:

      CALL DLENG (J+1, K)

C     Find average DS*W for interior points:

      WTSUM = ZERO
      DO I = 5, NIPTS - 5
         WTSUM = WTSUM + DS(I)*WEIGHT(I)
      END DO
      WTSUM = WTSUM / REAL (NIPTS - 9)

C     At IST, prevent WDS from being too large when DLENGS is very small:

      WDS = WTSUM / DLENGS       
      WDMAX = (ONE + A)*10.0
      IF (WDS > WDMAX) WDS = WDMAX

C     At IEND:

      WDE = WTSUM / DLENGE
      IF (WDE > WDMAX) WDE = WDMAX

      IF (LNSING == J .AND. K == KST) THEN
         WDES = WDE
         WDSS = WDS
      END IF

  900 CONTINUE

      END SUBROUTINE WTEDGE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               SUBROUTINE CURV1 (N,X,Y,SLP1,SLPN,YP,TEMP,SIGMA)
C
C
C       DIMENSION OF            X(N),Y(N),YF(N),TEMP(N)
C       ARGUMENTS
C
C       LATEST REVISION         FEBRUARY 5, 1974
C
C       PURPOSE                 CURV1 DETERMINES THE PARAMETERS NECESSARY TO
C                               COMPUTE AN INTERPOLATORY SPLINE UNDER TENSION
C                               THROUGH A SEQUENCE OF FUNCTIONAL VALUES.  THE
C                               SLOPES AT THE TWO ENDS OF THE CURVE MAY BE
C                               SPECIFIED OR OMITTED.  FOR ACTUAL COMPUTATION
C                               OF POINTS ON THE CURVE IT IS NECESSARY TO CALL
C                               THE FUNCTION CURV2.
C
C       USAGE                   CALL CURV1 (N,X,Y,SLP1,SLPN,YP,TEMP,SIGMA)
C
C       ARGUMENTS
C
C       ON INPUT                N
C                                 IS THE NUMBER OF VALUES TO BE INTERPOLATED
C                                 (N .GE. 2).
C
C                               X
C                                 IS AN ARRAY OF THE N INCREASING ABSCISSAE OF
C                                 THE FUNCTIONAL VALUES.
C
C                               Y
C                                 IS AN ARRAY OF THE N ORDINATES OF THE VALUES,
C                                 (I.E., Y(K) IS THE FUNCTIONAL VALUE
C                                 CORRESPONDING TO X(K)).
C
C                               SLP1 AND SLPN
C                                 CONTAIN THE DESIRED VALUES FOR THE FIRST
C                                 DERIVATIVE OF THE CURVE AT X(1) AND X(N),
C                                 RESPECTIVELY.  IF THE QUANTITY SIGMA IS
C                                 NEGATIVE THESE VALUES WILL BE DETERMINED
C                                 INTERNALLY AND THE USER NEED ONLY FURNISH
C                                 PLACE-HOLDING PARAMETERS FOR SLP1 AND SLPN.
C                                 SUCH PLACE-HOLDING (DUMMY) PARAMETERS WILL BE
C                                 IGNORED BUT NOT DESTROYED.
C
C                               YP
C                                 IS AN ARRAY OF LENGTH AT LEAST N.
C
C                               TEMP
C                                 IS AN ARRAY OF LENGTH AT LEAST N WHICH IS
C                                 USED FOR SCRATCH STORAGE.
C
C                               SIGMA
C                                 CONTAINS THE TENSION FACTOR.  THIS IS
C                                 NON-ZERO AND INDICATES THE CURVINESS DESIRED.
C                                 IF ABS(SIGMA) IS NEARLY ZERO (E.G., .001) THE
C                                 RESULTING CURVE IS APPROXIMATELY A CUBIC
C                                 SPLINE.  IF ABS(SIGMA) IS LARGE (E.G., 50.)
C                                 THE RESULTING CURVE IS NEARLY A POLYGONAL
C                                 LINE.  THE SIGN OF SIGMA INDICATES WHETHER
C                                 THE DERIVATIVE INFORMATION HAS BEEN INPUT OR
C                                 NOT.  IF SIGMA IS NEGATIVE THE ENDPOINT
C                                 DERIVATIVES WILL BE DETERMINED INTERNALLY.
C                                 A STANDARD VALUE FOR SIGMA IS APPROXIMATELY
C                                 1. IN ABSOLUTE VALUE.
C
C       ON OUTPUT               YP
C                                 CONTAINS VALUES PROPORTIONAL TO THE SECOND
C                                 DERIVATIVE OF THE CURVE AT THE GIVEN NODES.
C
C                               N, X, Y, SLP1, SLPN AND SIGMA ARE UNALTERED.
C
C       ENTRY POINTS            CURV1
C
C       SPECIAL CONDITIONS      NONE
C
C       COMMON BLOCKS           NONE
C
C       I/O                     NONE
C
C       PRECISION               SINGLE
C
C       REQUIRED ULIB           NONE
C       ROUTINES
C
C       SPECIALIST              RUSSELL K. REW, NCAR, BOULDER, COLORADO  80302
C
C       LANGUAGE                FORTRAN
C
C       HISTORY                 WRITTEN BY A. K. CLINE, MARCH 1972.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CURV1(N,X,Y,SLP1,SLPN,YP,TEMP,SIGMA)
C     INTEGER*4 N
C     REAL*4 X(N),Y(N),SLP1,SLPN,YP(N),TEMP(N),SIGMA
      DIMENSION X(N),Y(N),YP(N),TEMP(N)
      NN = N
      NM1 = NN-1
      NP1 = NN+1
      DELX1 = X(2) - X(1)
      DX1 = (Y(2)-Y(1)) / DELX1
C
C     DETERMINE THE SLOPES, IF NECESSARY.
C
      IF (SIGMA .LT. 0.0) GO TO 100
      SLPP1 = SLP1
      SLPPN = SLPN
C
C     DENORMALISE THE TENSION FACTOR.
C
   20 SIGMAP = (ABS(SIGMA)*FLOAT(NN-1)) / (X(NN)-X(1))
C
C     SET UP THE RHS AND THE TRIDIAGONAL SYSTEM FOR YP AND
C     PERFORM FORWARD ELIMINATION.
C
      DELS = SIGMAP*DELX1
      EXPS = EXP(DELS)
      SINHS = 0.5 * (EXPS - 1.0/EXPS)
      SINHIN = 1.0  / (DELX1*SINHS)
      DIAG1 = SINHIN * (DELS * .5 * (EXPS + 1.0/EXPS) - SINHS)
      DIAGIN = 1.0/DIAG1
      YP(1) = DIAGIN * (DX1-SLPP1)
      SPDIAG = SINHIN * (SINHS-DELS)
      TEMP(1) = DIAGIN*SPDIAG
      IF (NN .EQ. 2) GO TO 60
      DO 40 I=2,NM1
        DELX2 = X(I+1)-X(I)
        DX2 = (Y(I+1)-Y(I)) / DELX2
        DELS = SIGMAP*DELX2
        EXPS = EXP(DELS)
        SINHS = 0.5 * (EXPS - 1.0/EXPS)
        SINHIN = 1.0 / (DELX2*SINHS)
        DIAG2 = SINHIN * (DELS * .5 * (EXPS + 1.0/EXPS) - SINHS)
        DIAGIN = 1.0 / (DIAG1 + DIAG2 - SPDIAG*TEMP(I-1))
        YP(I) = DIAGIN * (DX2 - DX1 - SPDIAG*YP(I-1))
        SPDIAG = SINHIN * (SINHS-DELS)
        TEMP(I) = DIAGIN*SPDIAG
        DX1 = DX2
        DIAG1 = DIAG2
   40 CONTINUE
   60 DIAGIN = 1.0 / (DIAG1 - SPDIAG*TEMP(NM1))
      YP(NN) = DIAGIN * (SLPPN - DX2 - SPDIAG*YP(NM1))
C
C     PERFORM BACK SUBSTITUTION.
C
      DO 80 I=2,NN
        IBAK = NP1 - I
        YP(IBAK) = YP(IBAK) - TEMP(IBAK)*YP(IBAK+1)
   80 CONTINUE
      RETURN
  100 IF (NN .EQ. 2) GO TO 120
C
C     IF NO DERIVATIVES ARE GIVEN, USE SECOND ORDER POLYNOMIALS
C     INTERPOLATION ON INPUT DATA FOR VALUES AT END-POINTS.
C
      DELX2 = X(3) - X(2)
      DELX12 = X(3)-X(1)
      C1 = -(DELX12 + DELX1) / (DELX12*DELX1)
      C2 = DELX12 / (DELX1*DELX2)
      C3 = -DELX1 / (DELX12*DELX2)
      SLPP1 = C1*Y(1) + C2*Y(2) + C3*Y(3)
      DELN = X(NN) - X(NM1)
      DELNM1 = X(NM1) - X(NN-2)
      DELNN = X(NN) - X(NN-2)
      C1 = (DELNN + DELN) / (DELNN*DELN)
      C2 = -DELNN / (DELN*DELNM1)
      C3 = DELN / (DELNN*DELNM1)
      SLPPN = C3*Y(NN-2) + C2*Y(NM1) + C1*Y(NN)
      GO TO 20
C
C     IF ONLY TWO POINTS AND NO DERIVATIVES ARE GIVEN -
C     USE A STRAIGHT LINE.
C
  120 YP(1) = 0.0
      YP(2) = 0.0
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FUNCTION CURV2 (T,N,X,Y,YP,SIGMA,IT)
C
C
C       DIMENSION OF            X(N),Y(N),YP(N)
C       ARGUMENTS
C
C       LATEST REVISION         OCTOBER 22, 1973
C
C       PURPOSE                 CURV2 INTERPOLATES A CURVE AT A GIVEN POINT
C                               USING A SPLINE UNDER TENSION.  THE SUBROUTINE
C                               CURV1 SHOULD BE CALLED EARLIER TO DETERMINE
C                               CERTAIN NECESSARY PARAMETERS.
C
C       USAGE                   YT = CURV2 (T,N,X,Y,YP,SIGMA,IT)
C
C       ARGUMENTS
C
C       ON INPUT                T
C                                 CONTAINS A REAL VALUE TO BE MAPPED ONTO THE
C                                 INTERPOLATING CURVE.
C
C                               N
C                                 CONTAINS THE NUMBER OF POINTS WHICH WERE
C                                 INTERPOLATED TO DETERMINE THE CURVE.
C
C                               X AND Y
C                                 ARE ARRAYS CONTAINING THE ABSCISSAE AND
C                                 COORDINATES OF THE INTERPOLATED POINTS.
C
C                               YP
C                                 IS AN ARRAY WITH VALUES PROPORTIONAL TO THE
C                                 SECOND DERIVATIVE OF THE CURVE AT THE NODES.
C
C                               SIGMA
C                                 CONTAINS THE TENSION FACTOR (ITS SIGN IS
C                                 IGNORED).
C
C                               IT
C                                 IS AN INTEGER SWITCH.  IF IT IS NOT 1, THIS
C                                 INDICATES THAT THE FUNCTION HAS BEEN CALLED
C                                 PREVIOUSLY (WITH N, X, Y, YP AND SIGMA
C                                 UNALTERED) AND THAT THIS VALUE OF T EXCEEDS
C                                 THE PREVIOUS VALUE.  WITH SUCH INFORMATION
C                                 THE FUNCTION IS ABLE TO PERFORM THE
C                                 INTERPOLATION MUCH MORE RAPIDLY.  IF A USER
C                                 SEEKS TO INTERPOLATE AT A SEQUENCE OF POINTS
C                                 EFFICIENCY IS GAINED BY ORDERING THE VALUES
C                                 INCREASING AND SETTING IT TO THE INDEX OF THE
C                                 CALL.  IF IT IS 1 THE SEARCH FOR THE INTERVAL
C                                 (X(K),X(K+1)) CONTAINING T STARTS WITH K = 1
C
C                               THE PARAMETERS N, X, Y, YP AND SIGMA SHOULD BE
C                               INPUT UNALTERED FROM THE OUTPUT OF CURV1.
C
C       ON OUTPUT               CURV2
C                                 CONTAINS THE INTERPOLATED VALUE.  FOR T LESS
C                                 THAN X(1), CURV2 = Y(1).  FOR T GREATER THAN
C                                 X(N), CURV2 = Y(N).
C
C                               NONE OF THE INPUT PARAMETERS ARE ALTERED.
C
C       ENTRY POINTS            CURV2
C
C       SPECIAL CONDITIONS      NONE
C
C       COMMON BLOCKS           NONE
C
C       I/O                     NONE
C
C       PRECISION               SINGLE
C
C       REQUIRED ULIB           NONE
C       ROUTINES
C
C       SPECIALIST              RUSSELL K. REW, NCAR, BOULDER, COLORADO  80302
C
C       LANGUAGE                FORTRAN
C
C       HISTORY                 ORIGINALLY WRITTEN BY A. K. CLINE, MARCH, 1972
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION CURV2(T,N,X,Y,YP,SIGMA,IT)
C     INTEGER*4 N,IT
C     REAL*4 T,X(N),Y(N),YP(N),SIGMA
      DIMENSION X(N),Y(N),YP(N)
      S = X(N) - X(1)
C
C     DENORMALISE SIGMA
C
      SIGMAP = (ABS(SIGMA)*FLOAT(N-1)) / S
C
C     IF 'IT' IS NOT EQUAL TO 1, START THE SEARCH WHERE
C      PREVIOUSLY TERMINATED.
C     OTHERWISE START FROM THE BEGINNING.
C
      IF (IT .EQ. 1) I1=2
C
C     SEARCH FOR THE INTERVAL.
C
   20 DO 40 I=I1,N
        IF (X(I)-T) 40,40,60
   40 CONTINUE
      I=N
C
C     CHECK TO ENSURE THAT THE INTERVAL IS CORRECT.
C
   60 IF (X(I-1).LE.T .OR. T.LE.X(1)) GO TO 80
C
C     RESTART SEARCH AND RESET I1,
C      (INPUT 'IT' WAS INCORRECT).
C
      I1=2
      GO TO 20
C
C     SET-UP AND PERFORM INTERPOLATION.
C
   80 DEL1 = T - X(I-1)
      DEL2 = X(I) - T
      DELS = X(I) - X(I-1)
      EXPS1 = EXP(SIGMAP*DEL1)
      SINHD1 = 0.5 * (EXPS1 - 1.0/EXPS1)
      EXPS = EXP(SIGMAP*DEL2)
      SINHD2 = 0.5 * (EXPS - 1.0/EXPS)
      EXPS = EXPS1*EXPS
      SINHS = 0.5 * (EXPS - 1.0/EXPS)
      CURV2 = (YP(I)*SINHD1 + YP(I-1)*SINHD2) / SINHS
     *     + ( (Y(I)-YP(I))*DEL1 + (Y(I-1)-YP(I-1))*DEL2 ) / DELS
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FUNCTION CURVD (T,N,X,Y,YP,SIGMA,IT)
C
C
C       DIMENSION OF            X(N),Y(N),YP(N)
C       ARGUMENTS
C
C       LATEST REVISION         OCTOBER 22, 1973
C
C       PURPOSE                 CURVD DIFFERENTIATES A CURVE AT A GIVEN POINT
C                               USING AN INTERPOLATORY SPLINE UNDER TENSION.
C                               THE SUBROUTINE CURV1 SHOULD BE CALLED EARLIER
C                               TO DETERMINE CERTAIN NECESSARY PARAMETERS.
C
C       USAGE                   YTD = CURVD (T,N,X,Y,YP,SIGMA,IT)
C
C       ARGUMENTS
C
C       ON INPUT                T
C                                 CONTAINS A REAL VALUE AT WHICH THE DERIVATIVE
C                                 IS TO BE DETERMINED.
C
C                               N
C                                 CONTAINS THE NUMBER OF POINTS WHICH WERE
C                                 INTERPOLATED TO DETERMINE THE CURVE.
C
C                               X AND Y
C                                 ARRAYS CONTAINING THE ABSCISSAE AND ORDINATES
C                                 OF THE INTERPOLATED POINTS.
C
C                               YP
C                                 AN ARRAY WITH VALUES PROPORTIONAL TO THE
C                                 SECOND DERIVATIVE OF THE CURVE AT THE NODES.
C
C                               SIGMA
C                                 CONTAINS THE TENSION FACTOR (ITS SIGN IS
C                                 IGNORED).
C
C                               IT
C                                 AN INTEGER SWITCH.  IF IT IS NOT 1 THIS
C                                 INDICATES THAT THE FUNCTION HAS BEEN CALLED
C                                 PREVIOUSLY (WITH N, X, Y, YP AND SIGMA
C                                 UNALTERED) AND THAT THIS VALUE OF T EXCEEDS
C                                 THE PREVIOUS VALUE.  WITH SUCH INFORMATION
C                                 THE FUNCTION IS ABLE TO PERFORM THE
C                                 DIFFERENTIATION MUCH MORE RAPIDLY.  IF A USER
C                                 SEEKS TO DIFFERENTIATE AT A SEQUENCE OF
C                                 POINTS, EFFICIENCY IS GAINED BY ORDERING THE
C                                 VALUES INCREASING AND SETTING IT TO THE INDEX
C                                 OF THE CALL.  IF IT IS 1 THE SEARCH FOR THE
C                                 INTERVAL (X(K),X(K+1)) CONTAINING T STARTS
C                                 WITH K = 1.
C
C                               THE PARAMETERS N, X, Y, YP AND SIGMA SHOULD BE
C                               INPUT UNALTERED FROM THE OUTPUT OF CURV1.
C
C       ON OUTPUT               CURVD
C                                 CONTAINS THE DERIVATIVE VALUE.
C
C                               NONE OF THE INPUT PARAMETERS ARE ALTERED.
C
C       ENTRY POINTS            CURVD
C
C       SPECIAL CONDITIONS      NONE
C
C       COMMON BLOCKS           NONE
C
C       I/O                     NONE
C
C       PRECISION               SINGLE
C
C       REQUIRED ULIB           NONE
C       ROUTINES
C
C       SPECIALIST              RUSSELL K. REW, NCAR, BOULDER, COLORADO  80302
C
C       LANGUAGE                FORTRAN
C
C       HISTORY                 ORIGINALLY WRITTEN BY A. K. CLINE, AUGUST 1973.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION CURVD(T,N,X,Y,YP,SIGMA,IT)
      DIMENSION X(N),Y(N),YP(N)
      S = X(N) - X(1)
      SIGMAP = (ABS(SIGMA)*FLOAT(N-1)) / S
      IF (IT .EQ. 1) I1=2
   20 DO 40 I=I1,N
        IF (X(I)-T) 40,40,60
   40 CONTINUE
      I = N
   60 IF (X(I-1).LE.T .OR. T.LE.X(1)) GO TO 80
      I1 = 2
      GO TO 20
   80 DEL1 = T - X(I-1)
      DEL2 = X(I) - T
      DELS = X(I) - X(I-1)
      EXPS1 = EXP(SIGMAP*DEL1)
      COSHD1 = 0.5*(EXPS1 + 1.0/EXPS1)
      EXPS = EXP(SIGMAP*DEL2)
      COSHD2 = 0.5*(EXPS + 1.0/EXPS)
      EXPS = EXPS1*EXPS
      SINHS = 0.5*(EXPS - 1.0/EXPS)/SIGMAP
      CURVD = (YP(I)*COSHD1 - YP(I-1)*COSHD2)/SINHS
     2      + ((Y(I)-YP(I)) - (Y(I-1)-YP(I-1)))/DELS
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FUNCTION CURVI (XL,XU,N,X,Y,YP,SIGMA)
C
C
C       DIMENSION OF            X(N),Y(N),YP(N)
C       ARGUMENTS
C
C       LATEST REVISION         FEBRUARY 5, 1974
C
C       PURPOSE                 CURVI INTEGRATES A CURVE BETWEEN TWO GIVEN
C                               LIMITS USING AN INTERPOLATORY SPLINE UNDER
C                               TENSION.  THE SUBROUTINE CURV1 SHOULD BE CALLED
C                               EARLIER TO DETERMINE CERTAIN NECESSARY
C                               PARAMETERS.
C
C       USAGE                   YTI = CURVI (XL,XU,N,X,Y,YP,SIGMA)
C
C       ARGUMENTS
C
C       ON INPUT                XL AND XU
C                                 CONTAIN THE LOWER AND UPPER LIMITS OF
C                                 INTEGRATION, RESPECTIVELY.  (XL NEED NOT BE
C                                 LESS THAN OR EQUAL TO XU,
C                                 CURVI (XL,XU,...) = -CURVI (XU,XL,...)).
C
C                               N
C                                 CONTAINS THE NUMBER OF POINTS WHICH WERE
C                                 INTERPOLATED TO DETERMINE THE CURVE.
C
C                               X AND Y
C                                 ARRAYS CONTAINING THE ABSCISSAE AND ORDINATES
C                                 OF THE INTERPOLATED POINTS.
C
C                               YP
C                                 AN ARRAY WITH VALUES PROPORTIONAL TO THE
C                                 SECOND DERIVATIVE OF THE CURVE AT THE NODES.
C
C                               SIGMA
C                                 CONTAINS THE TENSION FACTOR (ITS SIGN IS
C                                 IGNORED).
C
C                               THE PARAMETERS N, X, Y, YP AND SIGMA SHOULD BE
C                               INPUT UNALTERED FROM THE OUTPUT OF CURV1.
C
C       ON OUTPUT               CURVI
C                                 CONTAINS THE INTEGRAL VALUE.
C
C                               NONE OF THE INPUT PARAMETERS ARE ALTERED.
C
C       ENTRY POINTS            CURVI
C
C       SPECIAL CONDITIONS      NONE
C
C       COMMON BLOCKS           NONE
C
C       I/O                     NONE
C
C       PRECISION               SINGLE
C
C       REQUIRED ULIB           NONE
C       ROUTINES
C
C       SPECIALIST              RUSSELL K. REW, NCAR, BOULDER, COLORADO 80302
C
C       LANGUAGE                FORTRAN
C
C       HISTORY                 ORIGINALLY WRITTEN BY A. K. CLINE, AUGUST 1973.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION CURVI( XL, XU, N, X, Y, YP, SIGMA )
      DIMENSION   X(N),  Y(N),  YP(N)
      NN = N
      S  = X(NN) - X(1)
C
C DENORMALIZE SIGMA
C
      SIGMAP = ABS(SIGMA) * FLOAT(NN-1) / S
C
C DETERMINE ACTUAL UPPER AND LOWER BOUNDS
C
      XXL = XL
      XXU = XU
      SSIGN = 1.0
      IF ( XL .LT. XU )  GO TO 101
      XXL = XU
      XXU = XL
      SSIGN = -1.0
      IF ( XL .GT. XU )  GO TO 101
C
C RETURN ZERO IF UPPER AND LOWER LIMITS ARE EQUAL
C
      CURVI = 0.0
      RETURN
C
C SEARCH FOR PROPER INTERVALS
C
  101 DO 102  I=2,NN
         IF ( X(I) -XXL )  102,103,103
  102 CONTINUE
      I = NN
  103 IL = I
      NP1 = NN + 1
      DO 104  I=2,NN
         NP1MI = NP1 - I
         IF ( X(NP1MI) - XXU )  105,105,104
  104 CONTINUE
      I = NN
  105 IU = NN + 2 - I
      IF ( IL .EQ. IU )  GO TO 110
C
C INTEGRATE FROM XXL TO X(IL)
C
      SUM = 0.0
      IF ( XXL .EQ. X(IL) )  GO TO 106
      DEL1 = XXL - X(IL-1)
      DEL2 = X(IL) - XXL
      DELS = X(IL) - X(IL-1)
      EXPS1 = EXP( SIGMAP*DEL1 )
      COSHD1 = 0.5 * ( EXPS1 + 1.0/EXPS1 )
      EXPS = EXP( SIGMAP*DEL2 )
      COSHD2 = 0.5 * ( EXPS + 1.0/EXPS )
      EXPS = EXPS1 * EXPS
      SINHS = SIGMAP * 0.5 * ( EXPS - 1.0/EXPS )
      COSHS = 0.5 * (EXPS + 1.0/EXPS )
      SUM = ( YP(IL) * ( COSHS-COSHD1 )  - YP(IL-1) * (1.0-COSHD2 ) )
     2    / SINHS  +  0.5 * ( ( Y(IL) - YP(IL) ) * ( DELS*DELS -
     3      DEL1*DEL1 )  + ( Y(IL-1) - YP(IL-1) ) * DEL2*DEL2 ) / DELS
  106 IF ( IU-IL .EQ. 1 )  GO TO 108
C
C INTEGRATE OVER INTERIOR ANGLES
C
      ILP1 = IL + 1
      IUM1 = IU - 1
      DO 107  I=ILP1,IUM1
         DELS = X(I) - X(I-1)
         EXPS = EXP( SIGMAP*DELS )
         SINHS = SIGMAP * 0.5 * ( EXPS - 1.0/EXPS )
         COSHS = 0.5 * ( EXPS + 1.0/EXPS )
         SUM  =  SUM  +  ( YP(I) + YP(I-1) ) * ( COSHS - 1.0 ) / SINHS
     2         - 0.5 * ( ( YP(I) + YP(I-1) ) - ( Y(I) + Y(I-1) ) )*DELS
  107 CONTINUE
  108 IF ( XXU .EQ. X(IU-1) )  GO TO 109
C
C INTEGRATE FROM X(IU-1) TAO XXU
C
      DEL1 = XXU - X(IU-1)
      DEL2 = X(IU) - XXU
      DELS = X(IU) - X(IU-1)
      EXPS1 = EXP( SIGMAP*DEL1 )
      COSHD1 = 0.5 * ( EXPS1 + 1.0/EXPS1 )
      EXPS = EXP( SIGMAP*DEL2 )
      COSHD2 = 0.5 * ( EXPS + 1.0/EXPS )
      EXPS = EXPS1 * EXPS
      SINHS = SIGMAP * 0.5 * ( EXPS - 1.0/EXPS )
      COSHS = 0.5 * (EXPS + 1.0/EXPS )
      SUM  =  SUM  +  ( YP(IU) * ( COSHD1 - 1.0 )  -
     2                  YP(IU-1) * ( COSHD2 - COSHS ) ) / SINHS
     3             +  0.5 * ( ( Y(IU) - YP(IU) ) * DEL1 * DEL1 -
     4                        ( Y(IU-1) - YP(IU-1) ) * ( DEL2 * DEL2 -
     5                          DELS * DELS ) ) / DELS
  109 CURVI = SSIGN * SUM
      RETURN
C
C INTEGRATE FROM XXL TO XXU
C
  110 DELU1 = XXU - X(IU-1)
      DELU2 = X(IU) - XXU
      DELL1 = XXL - X(IU-1)
      DELL2 = X(IU) - XXL
      DELS  = X(IU) - X(IU-1)
      EXPS1 = EXP( SIGMAP * DELU1 )
      COSHU1 = 0.5 * ( EXPS1 + 1.0/EXPS1 )
      EXPS  = EXP( SIGMAP * DELU2 )
      COSHU2 = 0.5 * ( EXPS + 1.0/EXPS )
      EXPS  = EXPS1 * EXPS
      SINHS = 0.5 * SIGMAP * ( EXPS - 1.0/EXPS )
      EXPS1 = EXP( SIGMAP * DELL1 )
      COSHL1 = 0.5 * ( EXPS1 + 1.0/EXPS1 )
      EXPS  = EXP( SIGMAP * DELL2 )
      COSHL2 = 0.5 * ( EXPS + 1.0/EXPS )
      SUM  =  ( YP(IU) * ( COSHU1 - COSHL1 )  -
     2          YP(IU-1) * ( COSHU2 - COSHL2 ) ) / SINHS  +
     3        0.5 * ( ( Y(IU) - YP(IU) ) * ( DELU1*DELU1 - DELL1*DELL1 )
     4          - ( Y(IU-1) - YP(IU-1) ) * ( DELU2*DELU2 - DELL2*DELL2 )
     5                                         ) / DELS
      GO TO 109
      END

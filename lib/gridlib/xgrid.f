C+---------------------------------------------------------------------
C
      SUBROUTINE XGRID (NPTS, MODE, XMIN, XMAX, X)
C
C  ONE-LINER:  Simple 1-D grid point distributions
C
C  PURPOSE:  XGRID generates NPTS abscissas in the range [XMIN, XMAX],
C            distributed according to MODE.
C
C  ARGUMENTS:
C    ARG    TYPE  I/O/S   DIM     DESCRIPTION
C    NPTS     I     I      -      Desired number of abscissas.
C    MODE     I     I      -      MODE = 0 means equally-spaced points;
C                                      = 1 means sinusoidal bunching
C                                          towards XMIN;
C                                      = 2 means sinusoidal bunching
C                                          towards XMIN and XMAX;
C                                      = 3 means sinusoidal bunching
C                                          towards XMAX.
C  XMIN,XMAX  R     I      -      First and last values to include.
C    X        R     O     NPTS    Desired distribution of abscissas.
C
C  ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C  HISTORY:
C    04/30/82    PJT    Original design and coding.
C    01/20/84    DAS    Added MODE=2 option (for circular arc case).
C    02/03/86    DAS    Changed loops from 1:NPTS to 2:NPTS-1, in case
C                       X(1) and X(NPTS) are passed as XMIN, XMAX.
C    12/02/92    DAS    Added bunching at XMAX option.
C
C  AUTHOR:     Phillip J. Trosin, Informatics General, Palo Alto, CA.
C
C----------------------------------------------------------------------

      DIMENSION X (NPTS)

      PIBY2  = ASIN (1.0)
      DX     = XMAX - XMIN
      DTHETA = PIBY2 / (NPTS - 1)

      IF (MODE .EQ. 0) THEN   ! Uniform

         DX = DX / (NPTS - 1)

         DO 100 I = 2, NPTS - 1
            X (I) = XMIN + DX * (I - 1)
  100    CONTINUE

      ELSE IF (MODE .EQ. 1) THEN  ! Sinusoidal bunching towards XMIN

         DO 200 I = 2, NPTS - 1
            X (I) = XMIN + DX * (1.0 - COS ((I - 1) * DTHETA))
  200    CONTINUE

      ELSE IF (MODE .EQ. 2) THEN  ! Sinusoidal bunching at both ends:

         DTHETA = DTHETA + DTHETA
         DX = DX * 0.5

         DO 300 I = 2, NPTS - 1
            X (I) = XMIN + DX * (1.0 - COS ((I - 1) * DTHETA))
  300    CONTINUE

      ELSE IF (MODE. EQ. 3) THEN  ! Sinusoidal bunching towards XMAX

         DO 400 I = 2, NPTS - 1
            X (I) = XMIN + DX * SIN ((I - 1) * DTHETA)
  400    CONTINUE

      END IF

      X (1) = XMIN
      X (NPTS) = XMAX

      RETURN
      END

C+----------------------------------------------------------------------
C
      SUBROUTINE RFFT2N (N, F)
C
C
C     Description:
C
C           This routine performs efficient calculation of the first N
C        Fourier coefficients of a sequence of 2N values by fast Fourier
C        transform, for N a power of 2.  The coefficients are unnormalized
C        in that if the original sequence is reconstructed by Fourier
C        synthesis, the result will be multiplied by 2N.  The rest of the
C        2N coefficients may be obtained with a simple modification if
C        necessary using:
C
C              F (N) = CMPLX (REAL (F (0)) - AIMAG (F (0)), ZERO) and,
C
C              F (2 * N - I) = CONJG (F (I)), for I = 1, 2, ... N
C
C        (more storage will be required).
C
C           Note that at this level, we consider the input array F to be
C        COMPLEX, but the machine representation should be the same as for
C        a calling routine which considers F to be REAL, with Re and Im
C        parts stored in alternate array elements.
C
C           The algorithm was adapted from Brigham (see Bibliography), and
C        rewritten in terms of complex variables in the (so far unrealized)
C        hope that it would be clearer.  Source and object code are in any
C        event more compact than the original.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        N                   I    I      Length of the data sequence in
C                                        complex form (i.e., there are
C                                        actually 2N real values).  N
C                                        must be a power of 2.
C        F        0:N-1      Z    I/O    Input data sequence (either a
C                                        REAL array of length 2N, or a
C                                        COMPLEX array of length N with
C                                        the input data packed into the
C                                        Re and Im parts.  On output, F
C                                        contains the Fourier coefs.
C
C
C     Development environments:  Digital VAX-11/780  VMS/V4.1   FORTRAN
C                                Cray X-MP/48        COS/V1.14  CFT/V1.11
C
C
C     Standards violations:  IMPLICIT NONE is non-standard.
C
C
C     Bibliography:
C
C        (1)  Brigham, E. Oran.  The Fast Fourier Transform.
C                Englewood Cliffs, New Jersey:  Prentice-Hall, 1974.
C                Chap. 10.  NOTE: The algorithm from section 10-10 was
C                amended to handle the first elements of A, B correctly.
C
C
C     Author:  Robert Kennelly, NASA-Ames/Sterling Software.
C
C
C     History:
C
C         2 Apr. 1986    RAK    Modified form of FOUCF from FLO6 - uses
C                               COMPLEX notation consistently (as prompted
C                               by CFFT2).
C        23 June 1986    RAK    Rewritten in standard form for 2N-point
C                               REAL transform by N-point COMPLEX FFT.
C        16 July 1986    RAK    A, B redimensioned to (LC) from (LC+1) -
C                               the last element wasn't used anywhere.
C                               First elements of A, B are special case
C                               and are handled separately.
C        18 July 1986    RAK    Removed normalization since it is
C                               application-dependent.  Return Nth coefs.
C                               since a similar IMSL routine does.
C        25 July 1986    RAK    Drop Nth coefs. so that the result can
C                               be returned in place.  Use "rotated"
C                               form for FM.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   ZERO, HALF, PI, HALFPI
      PARAMETER
     >   (ZERO   = 0.0E+0,
     >    HALF   = 0.5E+0,
     >    PI     = 3.14159265358979E+0,
     >    HALFPI = HALF * PI)

C     Arguments.

      INTEGER
     >   N
      COMPLEX
     >   F (0:N - 1)

C     Local variables.

      INTEGER
     >   I
      REAL
     >   PHASE, PIN
      COMPLEX
     >   FMROT, FP

C     Procedures.

      EXTERNAL
     >   CFFT2


C     Execution.
C     ----------

      PIN = PI / REAL (N)

C     Perform unnormalized Fourier analysis on the real sequence packed in F.

      CALL CFFT2 (N, F, 'ANALYSIS')

C     Rearrange the result to recover the Re and Im parts of the transform.

      F (0) = CMPLX (REAL (F (0)) + AIMAG (F (0)), ZERO)

C     Work from both ends so that the results can overwrite F.  The second
C     calculation of F (N / 2) is redundant but harmless.

      DO 10, I = 1, N / 2
         FP = F (I) + CONJG (F (N - I))
         PHASE = REAL (I) * PIN
         FMROT = (F (I) - CONJG (F (N - I))) *
     >      CMPLX (SIN (PHASE), COS (PHASE))

         F (I)     = HALF *       (FP - FMROT)
         F (N - I) = HALF * CONJG (FP + FMROT)
   10 CONTINUE


C     Termination.
C     ------------

      RETURN
      END

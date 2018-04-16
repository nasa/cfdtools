C+----------------------------------------------------------------------
C
      SUBROUTINE CFFT2 (N, C, MODE)
C
C
C     Description and usage:
C
C           Computes the Fourier transform of a complex sequence "in place."
C        The length of the sequence must be a power of 2.  The notation
C        used for the (output) function F(x) defined at N equally spaced
C        points in the interval [0, 2 pi] and the associated coefficients
C        C(n) (input) is as follows for Fourier "synthesis" (default):
C
C                          ---
C                          \                   +2 pi i j n / N
C                 F(j) =   /    C(2 pi n / N) e
C                          ---
C                        n=0,N-1
C                                         for  j = 0, 1, 2, ... , N-1
C
C        The Fourier coefficients C(n) may be obtained by conjugating the
C        input sequence, transforming, and then conjugating the result
C        ("analysis" mode).  This routine is unnormalized  - an analysis
C        followed by synthesis yields the original sequence multiplied by N.
C
C           Note that the convention used here is the same as that used
C        by the IMSL library, but variations are possible - check both the
C        normalization and the sign of the exponential term!  The original
C        authors of this routine do not make clear in their paper that
C        this program must be run "backwards" (as described above) to
C        obtain Fourier coefficients.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        N                   I    I      Length of the data sequence.
C        C         N         Z    I/O    Input AND output sequences (result
C                                        overwrites the original data).
C        MODE               C*1   I      Character flag indicating transform
C                                        direction.  If MODE = 'A' or 'a'
C                                        an analysis is performed by pre-
C                                        and post-conjugating the data.
C                                        Default mode is synthesis (+ve
C                                        imaginary exponent).
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
C        (1)  Cooley, J. W., P. A. W. Lewis, and P. D. Welch.  "The Fast
C                Fourier Transform and its Applications" in IEEE Trans.
C                on Education, Vol. 12, No. 1.  (March 1969)  Pp. 27-34.
C
C        (2)  Dahlquist, G., and A. Bjork.  Numerical Methods (Englewood
C                Cliffs, NJ: Prentice-Hall, 1974).  Pp. 413-417.
C
C
C     Author:  Robert Kennelly, NASA-Ames/Sterling Software.
C
C
C     History:
C
C         4 May  1984    RAK    Translation into FORTRAN 77, from code
C                               quoted by Dahlquist and Bjork.  Order of
C                               input parameters is reversed.
C        26 Mar. 1986    RAK    Tidied up the code and header.
C         7 Apr. 1986    RAK    Use successive multiplies by 2 for LE,
C                               and lag LE to form LEV2 (instead of **).
C        19 June 1986    RAK    Added synthesis/analysis switch (MODE).
C        17 July 1986    RAK    Use N as input parameter, compute log2(N).
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   ZERO, ONE, TWO, PI
      PARAMETER
     >   (ZERO = 0.0E+0,
     >    ONE  = 1.0E+0,
     >    TWO  = 2.0E+0,
     >    PI   = 3.14159265358979E+0)

C     Arguments.

      INTEGER
     >   N
      COMPLEX
     >   C (N)
      CHARACTER
     >   MODE

C     Local variables.

      INTEGER
     >   I, IP, J, K, L, LE, LEV2, NV2
      REAL
     >   ANGLE
      COMPLEX
     >   TEMP, U, W


C     Execution.
C     ----------

      NV2 = N / 2

C     If "analysis" mode, we must conjugate the data and the result.

      IF (MODE .EQ. 'A' .OR. MODE .EQ. 'a') THEN
         DO 10, I = 1, N
            C (I) = CONJG (C (I))
   10    CONTINUE
      END IF

C     Reorder the data.

      J = 1
      DO 30, I = 1, N - 1
         IF (I .LT. J) THEN
            TEMP = C (J)
            C (J) = C (I)
            C (I) = TEMP
         END IF

         K = NV2
   20    CONTINUE
         IF (K .LT. J) THEN
            J = J - K
            K = K / 2
            GO TO 20

         END IF
         J = J + K
   30 CONTINUE

C     Transform.

      LE = 1
      DO 60, L = 1, NINT (LOG (REAL (N)) / LOG (TWO))
         LEV2 = LE
         LE = 2 * LE
         U = CMPLX (ONE, ZERO)
         ANGLE = PI / REAL (LEV2)
         W = CMPLX (COS (ANGLE), SIN (ANGLE))
         DO 50, J = 1, LEV2

C           Cray FORTRAN compiler directive - permits the loop to vectorize
C           by "ignoring vector dependencies."
CDIR$       IVDEP

            DO 40, I = J, N, LE
               IP = I + LEV2
               TEMP = C (IP) * U
               C (IP) = C (I) - TEMP
               C (I) = C (I) + TEMP
   40       CONTINUE
            U = U * W
   50    CONTINUE
   60 CONTINUE

      IF (MODE .EQ. 'A' .OR. MODE .EQ. 'a') THEN
         DO 70, I = 1, N
            C (I) = CONJG (C (I))
   70    CONTINUE
      END IF


C     Termination.
C     ------------

      RETURN
      END

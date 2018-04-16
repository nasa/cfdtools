C+----------------------------------------------------------------------
C
      SUBROUTINE FSCALC ( MODE, N, F, A, B, SINCOS, DX, FX, IER )
C
C ACRONYM: Discrete Fourier Series CALCulations
C                   -       -      ----
C PURPOSE: FSCALC combines in one module several of the computations
C          associated with Fourier series for periodic uniform data:
C
C     MODE='a':  analyze - determine all of the Fourier coefficients;
C          's':  synthesize - recover the data points from the coefs.;
C          '1':  calculate the 1st derivatives at the data points;
C          '2':   "    "    "  2nd  "    "    "    "    "    " .
C          
C          It performs just one of these functions per call, but takes
C          advantage of repeated calls with the same value of N by
C          avoiding repeated sine/cosine evaluations.
C
C          To clarify:  FSCALC calculates all of the coefficients a(j),
C          b(j) (j=0:M-1) in the discrete Fourier series represented by
C
C          f(x) = a(0)/2 + SUM  [ a(j)*cos(jx) + b(j)*sin(jx) ]
C                        j=1:M-1
C
C          where N = 2M and the real values F(1:N) represent periodic
C          function f(x) at uniformly spaced intervals, with F(N+1)=F(1)
C          understood and therefore suppressed.  Note that b(0) = 0 is
C          assigned here.
C
C          Requiring that N be even halves the number of sines and cosines
C          needed.  If N were required to be a multiple of 4, this could
C          be halved again, but the saving is not considered worthwhile
C          since the bulk of the computation is done in the recursions.
C
C          The units of x are indicated by the interval DX - needed only
C          for the derivative calculations.  (The abscissas are assumed
C          to be 0, h, 2h, ..., (N-1)h for the analysis and synthesis
C          modes, where h = 2*pi / N.)
C
C METHOD:  FSCALC packages the algorithms used by FSCOEF (for the jth
C          coefficients) and FSEVAL (for summing a Fourier series by
C          efficient recursion) in a form appropriate for the case where
C          ALL of the coefficients are desired (as opposed to the first
C          so many) and where evaluations at the data points (as opposed
C          to arbitrary abscissas) are all that are needed.  (Use FSEVAL
C          for other evaluations of a series given the coefficients.)
C
C          The same recursion applies to both a(j) and b(j) - just a
C          different combination of the last iterates distinguishes them.
C          (If FSCOEF/FSEVAL were used, this advantage would be lost.)
C
C          Recovering the function values involves two similar but
C          different recursions for the sums of the "even" (cosine) terms
C          and the "odd" (sine) terms.  The values for the range [pi,2*pi]
C          are related to those for [0,pi] as different combinations of
C          the odds and evens.  No other symmetries are available.
C
C          For derivatives at the data points, the same recursions are
C          applied to appropriately adjusted forms of the coefficients.
C          Doing this in-line here is considered preferable to storing
C          adjusted a(j)s and b(j)s then invoking FSEVAL.
C
C          Since N defines the trig. functions needed here, it is assumed
C          that if the current value of N is the same as the value of N
C          on the previous call, then the SINCOS(*) argument contains the
C          appropriate sine/cosine values from some earlier call.
C
C          Note that the notation here is somewhat different from that
C          of FSCOEF/FSEVAL in an attempt to clarify usage.  The details
C          of the basic recursion used are given in "Numerical Methods
C          That (Usually) Work" by F.S.Acton (Harper and Row, 1970).
C          Briefly,
C
C           SUM  a(j)*cos(jx)  =  alpha(1)*cos(x) - alpha(2)
C          j=1:N
C
C          where  alpha(j)  =  2*cos(x)*alpha(j+1) - alpha(j+2) + a(j)
C
C          for j = N, N-1, N-2, ..., 1  with  alpha(N+1) = alpha(N+2) = 0.
C
C          A similar sum of sines degenerates to just  alpha(1)*sin(x).
C
C ARGUMENTS:
C    ARG    DIM  TYPE I/O/S DESCRIPTION
C   MODE     -    C*1   I   See PURPOSE above.  Briefly:
C                           'A' or 'a': "analyze"
C                           'S' or 's': "synthesize"
C                           '1' or '2': 1st or 2nd derivatives resp.
C    N       -     I    I   Number of data points (x,f(x)) not counting
C                           the repeated point implied by periodicity -
C                           it would be point N+1.  N >= 2 should be EVEN.
C    F      1:N    R    I   Discretized function values at uniformly-spaced
C                           intervals 0, DX, 2*DX, ..., (N-1)*DX, where DX
C                           is needed only for derivatives;
C                           input for "analyze" mode;
C                           unused by other modes
C    A   0:(N/2)-1 R   I/O  Coefficients a(j), j=0:M-1 where M=N/2;
C    B   0:(N/2)-1 R   I/O  Coefficients b(j), j=0:M-1 with b(0)=0;
C                           output for "analyze" mode;
C                           input for other modes
C  SINCOS  0:N+1   R   I/O  Storage for cos(jh), j=0:(N/2) followed by
C                               "    "  sin(jh),   "    "   (h = 2*pi/N);
C                           output if N is not the same as NLAST (local
C                           variable here, DATAd to 0 and SAVEd);
C                           input if N is the same as NLAST
C   DX       -     R    I   Distance (in some units) between data abscissas
C   FX     1:N+1   R    O   FX(1:N) are the quantities computed at the data
C                           points according to MODE (not used for "analyze"
C                           mode); FX(N+1) is needed to match FX(1) when
C                           advantage is taken of certain symmetry
C  IER       -     I    O   Error return code:
C                           IER = 0: no error was detected;
C                                 1: N is not even. (Odd N may be handled
C                                    with minor modifications not done yet.)
C                                 2: N is too small (less than 4);
C                                 3: MODE is not a valid character.
C
C ENVIRONMENT: VAX/VMS FORTRAN 77; IMPLICIT NONE is the only extension used.
C
C DEVELOPMENT HISTORY:
C       DATE   INITIALS   DESCRIPTION
C    06/16/87    DAS      Adaptation of FSCOEF/FSEVAL for the all-coef-
C                         ficients/data-point-evaluations-only cases.
C    09/18/87    DAS      Description for F(*) was unclear.  (Not used
C                         for modes other than 'analyze' mode.)
C
C AUTHOR: David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   N, IER
      REAL
     >   A (0 : (N/2) - 1), B (0 : (N/2) - 1), DX, F (N), FX (N),
     >   SINCOS (0 : N + 1)
      CHARACTER
     >   MODE * 1
      
C     Local constants:

      REAL
     >   HALF, PI, TWO, ZERO
      PARAMETER
     >  (HALF = 0.5E+0,
     >   PI   = 3.1415927E+0,
     >   TWO  = 2.E+0,
     >   ZERO = 0.E+0)

C     Local variables:

      INTEGER
     >   I, J, M, NLAST
      REAL
     >   AI, AIP1, AIP2, BI, BIP1, BIP2, C, EVENS, H, HJ, ODDS, SCALE,
     >   TWOCOSX

C     Data (initialized and SAVEd):

      SAVE
     >   NLAST
      DATA
     >   NLAST /0/

C     Execution:

      M = N / 2
      IF (M + M .NE. N) THEN
         IER = 1
         GO TO 99
      END IF

      IF (N .LT. 4) THEN
         IER = 2
         GO TO 99
      END IF

      IER = 0
      C = TWO / N
      H = C * PI

      IF (N .NE. NLAST) THEN
         NLAST = N

C        Even N means the [pi,2*pi] range repeats values from [0,pi].
C        It is easier to duplicate them here than select from half as
C        many below...
        
         DO 10, J = 0, M / 2
            HJ = H * J
            SINCOS (J) = TWO * COS (HJ)
            SINCOS (M - J) = -SINCOS (J)
            SINCOS (M + 1 + J) = SIN (HJ)
            SINCOS (N + 1 - J) = SINCOS (M + 1 + J)
   10    CONTINUE
      END IF

      IF (MODE .EQ. 'A' .OR. MODE .EQ. 'a') THEN

C        "Analyze" mode: compute the coefficients by recursion.

         DO 30, J = 0, M - 1
            TWOCOSX = SINCOS (J)
            AIP2 = ZERO
            AIP1 = ZERO

            DO 20, I = N, 2, -1
               AI = F (I) + TWOCOSX * AIP1 - AIP2
               AIP2 = AIP1
               AIP1 = AI
   20       CONTINUE

            A (J) = C * (F (1) + AI * TWOCOSX * HALF - AIP2)
            B (J) = C * AI * SINCOS (M + 1 + J)
   30    CONTINUE

      ELSE IF (MODE .EQ. 'S' .OR. MODE .EQ. 's') THEN

C        "Synthesize" mode: recover the data values from the coefficients.
C        Since SIN (pi + x) = -sin (pi - x) etc., the second half of
C        the values are related to the first half (and the abscissas are
C        the same as needed for computing the coefficients above).

         DO 50, J = 0, M
            TWOCOSX = SINCOS (J)
            AIP2 = ZERO
            AIP1 = ZERO
            BIP2 = ZERO
            BIP1 = ZERO

            DO 40, I = M - 1, 1, -1
               AI = A (I) + TWOCOSX * AIP1 - AIP2
               BI = B (I) + TWOCOSX * BIP1 - BIP2
               BIP2 = BIP1
               BIP1 = BI
               AIP2 = AIP1
               AIP1 = AI
   40       CONTINUE

            EVENS = (A (0) + AI * TWOCOSX) * HALF - AIP2
            ODDS  = BI * SINCOS (M + 1 + J)
            FX (J + 1) = EVENS + ODDS
            FX (N + 1 - J) = EVENS - ODDS
   50    CONTINUE

      ELSE IF (MODE .EQ. '1') THEN

C        First derivatives mode: do it in-line as above, with suitable
C        adjustments to the Fourier coefficients:

         SCALE = H / DX
         DO 70, J = 0, M
            TWOCOSX = SINCOS (J)
            AIP2 = ZERO
            AIP1 = ZERO
            BIP2 = ZERO
            BIP1 = ZERO

            DO 60, I = M - 1, 1, -1
               AI =  I * B (I) + TWOCOSX * AIP1 - AIP2
               BI = -I * A (I) + TWOCOSX * BIP1 - BIP2
               BIP2 = BIP1
               BIP1 = BI
               AIP2 = AIP1
               AIP1 = AI
   60       CONTINUE

            EVENS = AI * TWOCOSX * HALF - AIP2
            ODDS  = BI * SINCOS (M + 1 + J)
            FX (J + 1) = SCALE * (EVENS + ODDS)
            FX (N + 1 - J) = SCALE * (EVENS - ODDS)
   70    CONTINUE

      ELSE IF (MODE .EQ. '2') THEN

C        Second derivatives mode:

         SCALE = (H / DX) ** 2
         DO 90, J = 0, M
            TWOCOSX = SINCOS (J)
            AIP2 = ZERO
            AIP1 = ZERO
            BIP2 = ZERO
            BIP1 = ZERO

            DO 80, I = M - 1, 1, -1
               AI = (-I * I) * A (I) + TWOCOSX * AIP1 - AIP2
               BI = (-I * I) * B (I) + TWOCOSX * BIP1 - BIP2
               BIP2 = BIP1
               BIP1 = BI
               AIP2 = AIP1
               AIP1 = AI
   80       CONTINUE

            EVENS = AI * TWOCOSX * HALF - AIP2
            ODDS  = BI * SINCOS (M + 1 + J)
            FX (J + 1) = SCALE * (EVENS + ODDS)
            FX (N + 1 - J) = SCALE * (EVENS - ODDS)
   90    CONTINUE

      ELSE

C        Invalid mode:

         IER = 3

      END IF

   99 RETURN
      END

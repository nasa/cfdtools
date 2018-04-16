C+------------------------------------------------------------------------------
C
      SUBROUTINE BSZEROMN (MODE, IV, C, U1, U2, UOPTIM, V, LUNERR, IER)
C
C     One-liner:
C
C        Find a B-spline curve zero, minimum, or maximum in some range
C
C     Description:
C
C           BSZEROMN finds a zero, minimum or maximum of the specifed dependent
C        variable (X, Y, ...) in the specified range of the parametric variable
C        U for a B-spline curve in n-space (rational or nonrational).  It
C        returns the optimal value of U and all corresponding coordinate values.
C
C           Reentrant forms of the univariate zero-finding and minimization
C        utilities by Richard Brent are used (function values only).  These
C        typically converge in a dozen or so iterations, so there is little
C        point in attempting to use the available analytic derivatives.  Full
C        machine precision is specified to the utilities via zero tolerances.
C
C     Arguments:
C
C        Name      Dim    Type   I/O/S   Description
C
C        MODE             C*(*)  I       Only the first two characters are
C                                        checked.  They should be upper case:
C                                        'ZE' means find a zero;
C                                        'MI' means find a minimum;
C                                        'MA' means find a maximum.
C
C        IV               I      I       IV = 1, 2, ... indicates the
C                                        coordinate being zeroed or optimized.
C
C        C         (*)    R      I       B-spline vector in DTNURBS form.
C
C        U1,              R      I       Lower and upper range of the
C        U2               R      I       parametric variable U to search.
C
C        UOPTIM           R        O     Optimal value of U found.
C
C        V       |C(2)|   R        O     V(IV) is the optimal value of the
C                                        specified coordinate; corresponding
C                                        values of the other coordinates
C                                        are also returned.
C
C        LUNERR           I      I       Logical unit for error messages.
C
C        IER              I        O     Success/error code:
C                                        0 = no error;
C                                        1 = 100-iteration limit hit;
C                                        2 = some other fatal error in ZERORC
C                                            or FMINRC, which give diagnostics.
C                                        3 = polynomial degree too high for
C                                            local storage;
C                                      < 0 = error code reported by DTSPVL, q.v.
C     Procedures:
C
C        DTSPVL    Evaluates a univariate B-spline at a given U.
C        FMINRC    Seeks a local minimum of a function in a specified interval.
C        ZERORC    Seeks a zero of a function in a specified interval.
C
C     Environment:  VAX/VMS; FORTRAN 77 with minor extensions.
C
C     History:
C
C     04/21/92   D.A.Saunders   Initial implementation, with normalizing the
C                               chord of a B-spline airfoil in mind.
C     02/01/94        "         Dispensed with DOUBLE PRECISION.  Compile with
C                               /REAL_LENGTH=64 or -r8 switches as necessary.
C
C     Author:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      CHARACTER
     >   MODE * 2
      INTEGER
     >   IV, LUNERR, IER
      REAL
     >   C (*), U1, U2, UOPTIM, V (*)

C     Local constants:

      INTEGER
     >   MAXFUN, MAXORDER, NWORK
      REAL
     >   ONE, TOL
      CHARACTER
     >   NAME * 8

      PARAMETER
     >  (MAXFUN   = 100,               ! Maximum number of function evaluations
     >   MAXORDER = 7,                 ! Local storage handles polynomial degree
     >   NAME     = 'BSZEROMN',        ! <= 6 (order <= 7).
     >   NWORK    = 5 * MAXORDER - 2,  ! This is why we limit the degree.
     >   ONE      = 1.E+0,
     >   TOL      = 0.E+0)

C     Local variables:

      INTEGER
     >   ICASE, ISTAT, NUMFUN
      REAL
     >   FSIGN, FU, HOLD (13), WORK (NWORK)  ! ZERORC can't use DTSPVL's WORK.
      CHARACTER
     >   UTILITY (2) * 6

C     Procedures:
 
      EXTERNAL
     >   DTSPVL, FMINRC, ZERORC

C     Storage:

      DATA
     >   UTILITY /'ZERORC', 'FMINRC'/


C     Execution:

C     Distinguish the cases here - may be needed for an error message too.

      FSIGN = ONE
      ICASE = 2
      IF (MODE .EQ. 'ZE') THEN
         ICASE = 1
      ELSE IF (MODE .EQ. 'MA') THEN
         FSIGN = -ONE
      END IF

      ISTAT = 2          ! 2 initializes the iteration
      NUMFUN = MAXFUN

C     Most of the error checking is left to the lower level utilities.

      IF (C (3) .GT. MAXORDER) THEN
         IER = 3
         WRITE (LUNERR, 1000) 'Local storage exceeded.'
         GO TO 900
      END IF

C     Zero-finding and minimum-finding are so similar they can be done
C     in the same loop.

  200 CONTINUE

         IF (ICASE .EQ. 1) THEN

            CALL ZERORC (U1, U2, UOPTIM, FU, TOL, NUMFUN, NAME,
     >                   -LUNERR, HOLD, ISTAT)
         ELSE

            CALL FMINRC (U1, U2, UOPTIM, FU, TOL, NUMFUN, NAME,
     >                   -LUNERR, ISTAT)
         END IF

         IF (ISTAT .EQ. -1) THEN  ! Probable application error

            IER = 1
            WRITE (LUNERR, 1000) 'Iteration limit reached.'
            GO TO 900

         ELSE IF (ISTAT .LT. 0) THEN  ! Some other fatal error

            IER = 2
            WRITE (LUNERR, 1000) 'Fatal return from a utility.'
            GO TO 900

         ELSE IF (ISTAT .GT. 0) THEN  ! Evaluate the function

            CALL DTSPVL (UOPTIM, C, WORK, NWORK, V, IER)

            IF (IER .NE. 0) THEN
               WRITE (LUNERR, 1000) 'Fatal return from DTSPVL.'
               GO TO 900
            END IF

            FU = V (IV) * FSIGN
            GO TO 200

         ELSE  ! ISTAT = 0 (success).  Ensure all coordinates are at U = UOPTIM:

            CALL DTSPVL (UOPTIM, C, WORK, NWORK, V, IER)

         END IF

      GO TO 999


C     Error handling:

  900 WRITE (LUNERR, 1010) UTILITY (ICASE), ISTAT, IER


C     Termination:

  999 RETURN


C     Formats:

 1000 FORMAT (/, ' *** Trouble: ', A)
 1010 FORMAT (' BSZEROMN was using ', A, '.  ISTAT:', I3, '  IER:', I4)

      END
